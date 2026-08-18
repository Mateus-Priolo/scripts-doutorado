#!/usr/bin/env Rscript
# =============================================================================
# functions.R — Utilitários compartilhados | Pipeline scWGS + MEDICC2
# GSE173279 | Evolução Clonal Tumoral | GridUNESP
# Diretório de trabalho: /home/renanomete/projetos/matdata/evo_clonal
# =============================================================================

suppressPackageStartupMessages(library(Matrix))

WORKDIR <- "/home/renanomete/projetos/matdata/evo_clonal"

log_step <- function(step_num, msg) {
  ts <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
  cat(sprintf("\n%s [STEP %02d] %s\n", ts, step_num, msg))
  flush.console()
}

log_info <- function(msg) {
  cat(sprintf("  \u21b3 %s\n", msg))
  flush.console()
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  invisible(path)
}

# QC — Coeficiente de Gini
gini_coefficient <- function(x) {
  x <- sort(x[!is.na(x) & x > 0])
  n <- length(x)
  if (n < 2) return(NA_real_)
  sum((2L * seq_along(x) - n - 1L) * x) / (n * sum(x))
}

# QC — Score de aberração CNV: sum(|CN - 2|) por célula
cnv_aberration_score <- function(cn_mat, diploid_cn = 2L) {
  if (inherits(cn_mat, "dgCMatrix"))
    Matrix::colSums(abs(cn_mat - diploid_cn))
  else
    colSums(abs(cn_mat - diploid_cn))
}

# Filtragem QC via n x MAD
is_outlier_low  <- function(x, n = 3) x < median(x, na.rm=TRUE) - n * mad(x, na.rm=TRUE)
is_outlier_high <- function(x, n = 3) x > median(x, na.rm=TRUE) + n * mad(x, na.rm=TRUE)

# Estimativa de ploidy via moda KDE
estimate_ploidy_modal <- function(cn_vec) {
  cn_pos <- cn_vec[cn_vec > 0 & !is.na(cn_vec)]
  if (length(cn_pos) < 10) return(2L)
  d <- tryCatch(density(cn_pos, bw = 0.5, from = 1, to = max(cn_pos) + 1),
                error = function(e) NULL)
  if (is.null(d)) return(2L)
  max(1L, as.integer(round(d$x[which.max(d$y)])))
}

# Selecao de bins variaveis por MAD
select_variable_bins <- function(cn_mat, n_bins = 3000, exclude_sex_chr = TRUE) {
  if (exclude_sex_chr && !is.null(rownames(cn_mat))) {
    sex <- grepl("^chrX|^chrY|^X:|^Y:", rownames(cn_mat))
    log_info(sprintf("Excluindo %d bins em chr sexuais.", sum(sex)))
    cn_mat <- cn_mat[!sex, , drop = FALSE]
  }
  bin_mad  <- apply(cn_mat, 1, mad, na.rm = TRUE)
  variable <- bin_mad > 0
  log_info(sprintf("Bins com MAD > 0: %d / %d", sum(variable), length(variable)))
  cn_mat  <- cn_mat[variable, , drop = FALSE]
  bin_mad <- bin_mad[variable]
  n_sel   <- min(n_bins, nrow(cn_mat))
  top     <- order(bin_mad, decreasing = TRUE)[seq_len(n_sel)]
  log_info(sprintf("Bins selecionados: %d", n_sel))
  cn_mat[top, , drop = FALSE]
}

# Parsear coordenadas genomicas dos bins
parse_bin_coords <- function(bin_names) {
  if (all(grepl("^[^:]+:[0-9]+-[0-9]+$", head(bin_names, 20)))) {
    parts <- strsplit(bin_names, ":|-")
    return(do.call(rbind, lapply(parts, function(x)
      data.frame(chrom=x[1], start=as.integer(x[2]), end=as.integer(x[3]),
                 stringsAsFactors=FALSE))))
  }
  if (all(grepl("^[^_]+_[0-9]+_[0-9]+$", head(bin_names, 20)))) {
    parts <- strsplit(bin_names, "_")
    return(do.call(rbind, lapply(parts, function(x)
      data.frame(chrom=x[1], start=as.integer(x[2]), end=as.integer(x[3]),
                 stringsAsFactors=FALSE))))
  }
  warning("Formato de nome de bin nao reconhecido.")
  return(NULL)
}

cat("  [functions.R carregado | MEDICC2 pipeline]\n")
