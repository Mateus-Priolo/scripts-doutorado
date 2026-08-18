#!/usr/bin/env Rscript
# =============================================================================
# 02_cnv_qc_prepare.R — QC de células scWGS + geração de inputs MEDICC2
# Entrada : data/raw/GSE173279_scWGS_{JK136,JK142,JK153}_coarse_cn.tsv.gz
#           data/raw/GSE173279_scWGS_{JK136,JK142,JK153}_metadata.csv.gz
# Saída   : data/processed/cnv_qc_object.rds
#           data/medicc2_input/{all_cells,JK136,JK142,JK153}_medicc2.tsv
#           results/qc/qc_metrics.csv
#           results/qc/qc_plots.pdf
# =============================================================================

set.seed(42)

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  library(future)
  library(future.apply)
})

WORKDIR <- "/home/renanomete/projetos/matdata/evo_clonal"
setwd(WORKDIR)

n_workers <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "4"))
plan(multicore, workers = n_workers)
options(future.globals.maxSize = 32 * 1024^3)

# Carregar apenas funções utilitárias básicas (ensure_dir, log_*)
source("scripts/functions.R")
ensure_dir("data/processed")
ensure_dir("data/medicc2_input")
ensure_dir("results/qc")

cat("=================================================================\n")
cat(" ETAPA 2 — QC e preparação de input para MEDICC2\n")
cat(sprintf(" Início   : %s\n", format(Sys.time())))
cat(sprintf(" Workers  : %d\n", n_workers))
cat(sprintf(" Workdir  : %s\n", WORKDIR))
cat("=================================================================\n")

# =============================================================================
# Parâmetros globais
# =============================================================================
PATIENTS     <- c("JK136", "JK142", "JK153")
QC_MAD_N     <- 3L
QC_MIN_CELLS <- 5L
DIPLOID_CN   <- 2L

# =============================================================================
# FUNÇÕES AUXILIARES DE QC (reimplementadas para garantir robustez)
# =============================================================================

#' Coeficiente de Gini normalizado (0 = todos iguais, 1 = máxima desigualdade)
gini_coefficient <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 2) return(NA_real_)
  x <- sort(x)
  n <- length(x)
  sum_x <- sum(x)
  if (sum_x == 0) return(NA_real_)
  2 * sum(seq_len(n) * x) / (n * sum_x) - (n + 1) / n
}

#' Ploidia modal (valor inteiro mais frequente)
estimate_ploidy_modal <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_integer_)
  tab <- table(round(x))
  as.integer(names(tab)[which.max(tab)])
}

#' Score de aberração CNV: soma dos desvios absolutos do diploide (2)
cnv_aberration_score <- function(mat, diploid = 2) {
  # mat: bins x células
  apply(mat, 2, function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) return(NA_real_)
    sum(abs(x - diploid), na.rm = TRUE)
  })
}

#' Detecta outliers superiores usando MAD (retorna TRUE se outlier, FALSE caso contrário)
is_outlier_high <- function(x, nmad = 3) {
  x_clean <- x[!is.na(x)]
  if (length(x_clean) == 0) return(rep(TRUE, length(x)))  # se não há dados, considera outlier
  med <- median(x_clean)
  mad_val <- mad(x_clean, constant = 1)
  threshold <- med + nmad * mad_val
  ifelse(is.na(x) | x > threshold, TRUE, FALSE)
}

#' Detecta outliers inferiores usando MAD
is_outlier_low <- function(x, nmad = 3) {
  x_clean <- x[!is.na(x)]
  if (length(x_clean) == 0) return(rep(TRUE, length(x)))
  med <- median(x_clean)
  mad_val <- mad(x_clean, constant = 1)
  threshold <- med - nmad * mad_val
  ifelse(is.na(x) | x < threshold, TRUE, FALSE)
}

# =============================================================================
# FUNÇÃO: Parse de coordenadas genômicas (formato: chr1_1040001_2080000)
# =============================================================================
parse_bin_coords <- function(bin_names) {
  parts <- strsplit(bin_names, "[_:\\-]")
  chrom <- sapply(parts, `[`, 1)
  start <- as.integer(sapply(parts, `[`, 2))
  end   <- as.integer(sapply(parts, `[`, 3))
  
  if (any(is.na(start)) || any(is.na(end))) return(NULL)
  data.frame(chrom = chrom, start = start, end = end, stringsAsFactors = FALSE)
}

# =============================================================================
# FUNÇÃO: detectar orientação e ler coarse_cn.tsv.gz + associar metadata
# =============================================================================
read_cn_file <- function(cn_path, meta_path = NULL) {
  if (!file.exists(cn_path)) stop("Arquivo não encontrado: ", cn_path)
  log_info(sprintf("Lendo %s ...", basename(cn_path)))

  hdr <- fread(cn_path, nrows = 3, data.table = FALSE,
               header = TRUE, check.names = FALSE)
  col_names <- colnames(hdr)

  is_coord <- function(x) grepl("^(chr)?[0-9XYMxy]+[:\\._\\-][0-9]+", x)

  frac_coord_cols  <- mean(is_coord(col_names))
  first_col_coords <- mean(is_coord(as.character(hdr[[1]])))

  log_info(sprintf("  Detecção: %.0f%% colunas com padrão genômico | %.0f%% valores col1 com padrão",
                   frac_coord_cols * 100, first_col_coords * 100))

  dt <- fread(cn_path, data.table = FALSE, header = TRUE, check.names = FALSE)

  if (frac_coord_cols > 0.3) {
    log_info("  Orientação: células × bins  →  transpondo para bins × células")
    id_col_idx <- which(!sapply(dt, is.numeric))[1]
    if (!is.na(id_col_idx) && !is_coord(dt[[id_col_idx]][1])) {
      cell_ids <- as.character(dt[[id_col_idx]])
      cn_vals  <- as.matrix(dt[, -id_col_idx, drop = FALSE])
    } else {
      cell_ids <- paste0("cell_", seq_len(nrow(dt)))
      cn_vals  <- as.matrix(dt)
    }
    bin_names <- colnames(cn_vals)
    cn_mat    <- t(cn_vals)
    colnames(cn_mat) <- cell_ids
    rownames(cn_mat) <- bin_names
  } else {
    log_info("  Orientação: bins × células")
    coord_idx <- which(!sapply(dt, is.numeric))[1]
    if (is.na(coord_idx)) {
      coord_idx <- 1L
      log_info("  AVISO: nenhuma coluna não-numérica detectada; usando col 1 como bin ID")
    }
    bin_names <- as.character(dt[[coord_idx]])
    cn_mat    <- as.matrix(dt[, -coord_idx, drop = FALSE])
    rownames(cn_mat) <- bin_names
  }

  mode(cn_mat) <- "integer"

  # Associação com metadata (se existir)
  meta <- NULL
  if (!is.null(meta_path) && file.exists(meta_path)) {
    meta <- fread(meta_path, data.table = FALSE)
    n_cn_cells <- ncol(cn_mat)
    n_meta_rows <- nrow(meta)

    if (n_cn_cells != n_meta_rows) {
      stop(sprintf(
        "Número de colunas na matriz (%d) difere do número de linhas no metadata (%d).\n",
        n_cn_cells, n_meta_rows),
        "  Verifique se a matriz e o metadata correspondem exatamente (ordem 1:1)."
      )
    }

    if ("barcode" %in% colnames(meta)) {
      cell_ids <- meta$barcode
    } else if ("cell_id" %in% colnames(meta)) {
      cell_ids <- as.character(meta$cell_id)
    } else {
      cell_ids <- paste0("cell_", seq_len(n_meta_rows))
    }

    colnames(cn_mat) <- cell_ids
    log_info(sprintf("  Células renomeadas com barcode. Exemplo: %s", cell_ids[1]))
  } else {
    log_info("  Metadata não fornecido; mantendo nomes originais das colunas.")
  }

  # Remover colunas 100% NA
  all_na_cols <- colSums(is.na(cn_mat)) == nrow(cn_mat)
  if (any(all_na_cols)) {
    log_info(sprintf("  Removendo %d células 100%% NA", sum(all_na_cols)))
    cn_mat <- cn_mat[, !all_na_cols, drop = FALSE]
    if (!is.null(meta)) meta <- meta[!all_na_cols, , drop = FALSE]
  }

  # Parse de coordenadas
  bin_coords <- parse_bin_coords(rownames(cn_mat))
  if (is.null(bin_coords)) {
    sample_bins <- head(rownames(cn_mat), 5)
    stop("Não foi possível parsear coordenadas dos bins.")
  }

  log_info(sprintf("  OK: %d bins × %d células", nrow(cn_mat), ncol(cn_mat)))
  log_info(sprintf("  Exemplo de bin: %s", rownames(cn_mat)[1]))
  log_info(sprintf("  Exemplo de célula: %s", colnames(cn_mat)[1]))

  list(cn_mat = cn_mat, bin_coords = bin_coords, meta = meta)
}

# =============================================================================
# FUNÇÃO: calcular métricas QC por célula
# =============================================================================
compute_qc_metrics <- function(cn_mat, patient_id) {
  cells <- colnames(cn_mat)
  log_info(sprintf("  Calculando QC: %d células (%s)...", length(cells), patient_id))

  aberr    <- cnv_aberration_score(cn_mat, DIPLOID_CN)
  n_meas   <- colSums(!is.na(cn_mat))
  gini_v   <- vapply(cells, function(c) gini_coefficient(cn_mat[, c]), numeric(1))
  ploidy_v <- vapply(cells, function(c) estimate_ploidy_modal(cn_mat[, c]), integer(1))

  metrics <- data.frame(
    cell             = cells,
    patient          = patient_id,
    n_bins           = nrow(cn_mat),
    n_measured       = n_meas,
    pct_measured     = n_meas / nrow(cn_mat),
    aberration_score = aberr,
    gini             = gini_v,
    mean_cn          = colMeans(cn_mat, na.rm = TRUE),
    ploidy           = ploidy_v,
    stringsAsFactors = FALSE,
    row.names        = NULL
  )

  # NA nas métricas = falha automática
  metrics$fail_low_coverage <- is_outlier_low( metrics$pct_measured,     QC_MAD_N)
  metrics$fail_high_gini    <- is_outlier_high(metrics$gini,             QC_MAD_N)
  metrics$fail_high_aberr   <- is_outlier_high(metrics$aberration_score, QC_MAD_N)
  
  # Se qualquer métrica for NA, considera falha no respectivo critério
  metrics$fail_low_coverage[is.na(metrics$pct_measured)]     <- TRUE
  metrics$fail_high_gini[is.na(metrics$gini)]                <- TRUE
  metrics$fail_high_aberr[is.na(metrics$aberration_score)]   <- TRUE

  metrics$pass_qc <- !(metrics$fail_low_coverage |
                       metrics$fail_high_gini   |
                       metrics$fail_high_aberr)

  n_pass <- sum(metrics$pass_qc, na.rm = TRUE)
  log_info(sprintf("  %s: %d/%d células passaram QC (%.0f%%)",
                   patient_id, n_pass, nrow(metrics),
                   100 * n_pass / nrow(metrics)))

  # Estatísticas com na.rm
  log_info(sprintf("  Aberration score: mediana=%.1f, min=%.1f, max=%.1f",
                   median(aberr, na.rm = TRUE),
                   min(aberr, na.rm = TRUE),
                   max(aberr, na.rm = TRUE)))
  log_info(sprintf("  Gini: mediana=%.3f | Ploidy modal mais comum: %d",
                   median(gini_v, na.rm = TRUE),
                   as.integer(names(sort(table(ploidy_v), decreasing=TRUE))[1])))

  metrics
}

# =============================================================================
# FUNÇÃO: gerar TSV input para MEDICC2
# =============================================================================
write_medicc2_tsv <- function(cn_mat, bin_coords, out_path,
                               cells_pass = colnames(cn_mat)) {
  cells_use <- intersect(cells_pass, colnames(cn_mat))

  log_info(sprintf("  [%s] cells_pass=%d | em cn_mat=%d | interseção=%d",
                   basename(out_path),
                   length(cells_pass), ncol(cn_mat), length(cells_use)))

  if (length(cells_use) == 0) {
    log_info(sprintf("  [SKIP] %s — nenhuma célula passou QC encontrada na matriz",
                     basename(out_path)))
    return(invisible(NULL))
  }
  if (length(cells_use) < QC_MIN_CELLS) {
    log_info(sprintf("  [SKIP] %s — %d células < mínimo de %d para MEDICC2",
                     basename(out_path), length(cells_use), QC_MIN_CELLS))
    return(invisible(NULL))
  }

  cn_use  <- cn_mat[, cells_use, drop = FALSE]
  n_bins  <- nrow(cn_use)
  n_cells <- ncol(cn_use)
  log_info(sprintf("  Escrevendo %s (%d células × %d bins)...",
                   basename(out_path), n_cells, n_bins))

  normal_dt <- data.table(
    sample_id   = "diploid_normal",
    chrom       = bin_coords$chrom,
    start       = bin_coords$start,
    end         = bin_coords$end,
    copy_number = DIPLOID_CN
  )

  cell_dts <- lapply(seq_len(n_cells), function(j) {
    cv <- cn_use[, j]
    cv[is.na(cv)] <- as.integer(round(median(cv, na.rm = TRUE)))
    cv[cv < 0L]  <- 0L
    data.table(
      sample_id   = colnames(cn_use)[j],
      chrom       = bin_coords$chrom,
      start       = bin_coords$start,
      end         = bin_coords$end,
      copy_number = as.integer(cv)
    )
  })

  out_dt <- rbindlist(c(list(normal_dt), cell_dts))
  fwrite(out_dt, out_path, sep = "\t", quote = FALSE, col.names = TRUE)
  log_info(sprintf("  OK %s (%.1f MB)", out_path, file.size(out_path) / 1024^2))
  invisible(out_path)
}

# =============================================================================
# STEP 1 — Carregar dados e metadados
# =============================================================================
log_step(1, "Carregando dados brutos e metadados...")

cn_data   <- vector("list", length(PATIENTS))
meta_list <- vector("list", length(PATIENTS))
names(cn_data) <- names(meta_list) <- PATIENTS

for (pat in PATIENTS) {
  cn_path   <- sprintf("data/raw/GSE173279_scWGS_%s_coarse_cn.tsv.gz", pat)
  meta_path <- sprintf("data/raw/GSE173279_scWGS_%s_metadata.csv.gz", pat)

  res <- read_cn_file(cn_path, meta_path)
  cn_data[[pat]]   <- list(cn_mat = res$cn_mat, bin_coords = res$bin_coords)
  meta_list[[pat]] <- res$meta

  if (!is.null(res$meta)) {
    meta_list[[pat]]$patient <- pat
    log_info(sprintf("  Metadata %s: %d linhas", pat, nrow(res$meta)))
  }
  gc()
}

for (pat in PATIENTS) {
  if (!is.null(meta_list[[pat]]) && !"cell" %in% colnames(meta_list[[pat]])) {
    meta_list[[pat]]$cell <- colnames(cn_data[[pat]]$cn_mat)
  }
}

cell_col    <- "cell"
patient_col <- "patient"

# =============================================================================
# STEP 2 — QC
# =============================================================================
log_step(2, "Calculando métricas QC por paciente...")

qc_list <- lapply(PATIENTS, function(pat)
  compute_qc_metrics(cn_data[[pat]]$cn_mat, patient_id = pat))
names(qc_list) <- PATIENTS

qc_metrics_all <- rbindlist(qc_list, fill = TRUE)
setDF(qc_metrics_all)

fwrite(qc_metrics_all, "results/qc/qc_metrics.csv")
log_info(sprintf("OK results/qc/qc_metrics.csv (%d células)", nrow(qc_metrics_all)))

qc_summary <- do.call(rbind, lapply(PATIENTS, function(pat) {
  m <- qc_list[[pat]]
  data.frame(
    patient      = pat,
    n_total      = nrow(m),
    n_pass       = sum(m$pass_qc, na.rm = TRUE),
    n_fail       = sum(!m$pass_qc, na.rm = TRUE),
    pct_pass     = round(mean(m$pass_qc, na.rm = TRUE) * 100, 1),
    median_aberr = round(median(m$aberration_score, na.rm = TRUE), 1),
    median_gini  = round(median(m$gini, na.rm = TRUE), 3)
  )
}))
cat("\nResumo QC por paciente:\n")
print(qc_summary)

# =============================================================================
# STEP 3 — Plots QC
# =============================================================================
log_step(3, "Gerando plots QC...")

make_qc_plots <- function(metrics_df, title_tag) {
  pal <- c("TRUE" = "#2196F3", "FALSE" = "#F44336")
  df_plot <- metrics_df[!is.na(metrics_df$pass_qc), ]
  if (nrow(df_plot) == 0) {
    return(ggplot() + annotate("text", x=1, y=1, label="Sem dados para plot") + theme_void())
  }

  p1 <- ggplot(df_plot, aes(aberration_score, fill = pass_qc)) +
    geom_histogram(bins = 50, alpha = 0.85, colour = "white", linewidth = 0.2) +
    scale_fill_manual(values = pal, name = "Pass QC") +
    labs(title = "Score de aberração CNV", x = "Σ |CN − 2|", y = "N células") +
    theme_bw(base_size = 9)

  p2 <- ggplot(df_plot, aes(gini, fill = pass_qc)) +
    geom_histogram(bins = 50, alpha = 0.85, colour = "white", linewidth = 0.2) +
    scale_fill_manual(values = pal, name = "Pass QC") +
    labs(title = "Índice de Gini CNV", x = "Gini", y = "N células") +
    theme_bw(base_size = 9)

  p3 <- ggplot(df_plot, aes(pct_measured, fill = pass_qc)) +
    geom_histogram(bins = 50, alpha = 0.85, colour = "white", linewidth = 0.2) +
    scale_fill_manual(values = pal, name = "Pass QC") +
    labs(title = "Cobertura de bins", x = "Fração bins medidos", y = "N células") +
    theme_bw(base_size = 9)

  p4 <- ggplot(df_plot, aes(aberration_score, gini, colour = pass_qc)) +
    geom_point(alpha = 0.55, size = 1.2) +
    scale_colour_manual(values = pal, name = "Pass QC") +
    labs(title = "Gini vs Aberração CNV") +
    theme_bw(base_size = 9)

  (p1 | p2) / (p3 | p4) +
    plot_annotation(
      title    = sprintf("QC scWGS — %s", title_tag),
      subtitle = sprintf("%d células | %d pass (%.0f%%) | %d fail",
                         nrow(df_plot), sum(df_plot$pass_qc),
                         mean(df_plot$pass_qc) * 100,
                         sum(!df_plot$pass_qc)),
      theme = theme(plot.title    = element_text(size = 13, face = "bold"),
                    plot.subtitle = element_text(size = 9, colour = "grey40"))
    )
}

pdf("results/qc/qc_plots.pdf", width = 13, height = 8)
print(make_qc_plots(qc_metrics_all, "Todas as células"))
for (pat in PATIENTS) print(make_qc_plots(qc_list[[pat]], pat))
invisible(dev.off())
log_info("OK results/qc/qc_plots.pdf")

# =============================================================================
# STEP 4 — Gerar inputs MEDICC2
# =============================================================================
log_step(4, "Gerando inputs MEDICC2...")

for (pat in PATIENTS) {
  cells_ok <- qc_list[[pat]]$cell[which(qc_list[[pat]]$pass_qc)]
  write_medicc2_tsv(
    cn_mat     = cn_data[[pat]]$cn_mat,
    bin_coords = cn_data[[pat]]$bin_coords,
    out_path   = sprintf("data/medicc2_input/%s_medicc2.tsv", pat),
    cells_pass = cells_ok
  )
  gc()
}

log_info("Identificando bins comuns entre pacientes...")
common_bins <- Reduce(intersect, lapply(cn_data, function(d) rownames(d$cn_mat)))
log_info(sprintf("  Bins comuns: %d", length(common_bins)))

if (length(common_bins) == 0)
  stop("Nenhum bin em comum entre os pacientes.")

all_cells_cn <- do.call(cbind, lapply(PATIENTS, function(pat)
  cn_data[[pat]]$cn_mat[common_bins, , drop = FALSE]))

all_coords <- cn_data[[PATIENTS[1]]]$bin_coords[
  match(common_bins, rownames(cn_data[[PATIENTS[1]]]$cn_mat)), ]

all_cells_pass <- qc_metrics_all$cell[which(qc_metrics_all$pass_qc)]

write_medicc2_tsv(
  cn_mat     = all_cells_cn,
  bin_coords = all_coords,
  out_path   = "data/medicc2_input/all_cells_medicc2.tsv",
  cells_pass = all_cells_pass
)
rm(all_cells_cn); gc()

# =============================================================================
# STEP 5 — Salvar objeto QC
# =============================================================================
log_step(5, "Salvando objeto QC...")

meta_all <- rbindlist(Filter(Negate(is.null), meta_list), fill = TRUE)
setDF(meta_all)

meta_final <- merge(
  meta_all,
  qc_metrics_all[, c("cell", "patient", "aberration_score",
                      "gini", "ploidy", "pass_qc")],
  by.x = c(cell_col, "patient"),
  by.y = c("cell",   "patient"),
  all.x = TRUE
)

qc_obj <- list(
  meta        = meta_final,
  qc_metrics  = qc_metrics_all,
  patient_col = patient_col,
  cell_col    = cell_col,
  qc_summary  = qc_summary,
  params      = list(mad_n      = QC_MAD_N,
                     min_cells  = QC_MIN_CELLS,
                     diploid_cn = DIPLOID_CN,
                     patients   = PATIENTS)
)

saveRDS(qc_obj, "data/processed/cnv_qc_object.rds")
log_info("OK data/processed/cnv_qc_object.rds")

# =============================================================================
# Relatório final
# =============================================================================
cat("\n=================================================================\n")
cat(sprintf(" Etapa 2 concluída: %s\n", format(Sys.time())))
cat("=================================================================\n\n")

cat("Arquivos gerados:\n")
expected <- c(
  "data/processed/cnv_qc_object.rds",
  "results/qc/qc_metrics.csv",
  "results/qc/qc_plots.pdf",
  sprintf("data/medicc2_input/%s_medicc2.tsv", c("all_cells", PATIENTS))
)
for (f in expected) {
  if (file.exists(f))
    cat(sprintf("  [OK] %-60s  %.1f MB\n", f, file.size(f) / 1024^2))
  else
    cat(sprintf("  [!!] %-60s  NAO GERADO\n", f))
}
cat("\n")
