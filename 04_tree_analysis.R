#!/usr/bin/env Rscript
# =============================================================================
# 04_tree_analysis.R — Análise da árvore filogenética e eventos de CNV
# Entrada : results/medicc2/  +  data/processed/cnv_qc_object.rds
# Saída   : data/processed/tree_analysis.rds
#           results/tree_analysis/{evo_metrics_per_cell,clone_summary,
#                                  patient_summary}.csv
# =============================================================================

set.seed(42)

suppressPackageStartupMessages({
  library(ape)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(viridis)
  library(RColorBrewer)
})

WORKDIR <- "/home/renanomete/projetos/matdata/evo_clonal"
setwd(WORKDIR)

source("scripts/functions.R")
ensure_dir("results/tree_analysis")

cat("=================================================================\n")
cat(" ETAPA 4 — Análise da Árvore Filogenética MEDICC2\n")
cat(sprintf(" Diretório : %s\n", WORKDIR))
cat(sprintf(" Início    : %s\n", format(Sys.time())))
cat("=================================================================\n")

# =============================================================================
# Função: carregar outputs de uma rodada MEDICC2
# =============================================================================
load_medicc2_outputs <- function(dir_path, run_name) {
  if (!dir.exists(dir_path)) {
    warning("Diretório não encontrado: ", dir_path)
    return(NULL)
  }

  files <- list.files(dir_path, full.names = TRUE)
  if (length(files) == 0) {
    warning("Diretório vazio: ", dir_path)
    return(NULL)
  }

  # CORREÇÃO: MEDICC2 gera árvores com extensão .new (Newick), não apenas
  # .tree / .nwk / .newick
  tree_file <- files[grepl("\\.(tree|nwk|newick|new)$", files,
                            ignore.case = TRUE)][1]
  dist_file <- files[grepl("distance", files, ignore.case = TRUE) &
                     grepl("\\.(tsv|csv|txt)$", files)][1]
  evt_file  <- files[grepl("event",    files, ignore.case = TRUE) &
                     grepl("\\.(tsv|csv|txt)$", files)][1]

  if (is.na(tree_file)) {
    warning("[", run_name, "] Arquivo de árvore não encontrado em ", dir_path,
            "\n  Arquivos presentes: ", paste(basename(files), collapse = ", "))
    return(NULL)
  }

  res <- list(run_name = run_name, dir = dir_path)
  res$tree <- ape::read.tree(tree_file)
  log_info(sprintf("[%s] Árvore: %d tips, %d nós internos",
                   run_name, ape::Ntip(res$tree), ape::Nnode(res$tree)))

  if (!is.na(dist_file)) {
    dm <- fread(dist_file, data.table = FALSE)
    rownames(dm) <- dm[[1]]
    res$dist_mat <- as.matrix(dm[, -1, drop = FALSE])
  }
  if (!is.na(evt_file)) res$events <- fread(evt_file, data.table = FALSE)

  return(res)
}

# =============================================================================
# 1. Carregar dados
# =============================================================================
log_step(1, "Carregando dados QC e resultados MEDICC2...")

qc_obj     <- readRDS("data/processed/cnv_qc_object.rds")
meta       <- qc_obj$meta
qc_metrics <- qc_obj$qc_metrics
patient_col <- qc_obj$patient_col

runs <- list(
  global = load_medicc2_outputs("results/medicc2/global", "global"),
  JK136  = load_medicc2_outputs("results/medicc2/JK136",  "JK136"),
  JK142  = load_medicc2_outputs("results/medicc2/JK142",  "JK142"),
  JK153  = load_medicc2_outputs("results/medicc2/JK153",  "JK153")
)
runs <- Filter(Negate(is.null), runs)
log_info(sprintf("Runs carregados: %s", paste(names(runs), collapse = ", ")))

if (length(runs) == 0)
  stop("Nenhum resultado MEDICC2 encontrado em results/medicc2/\n",
       "Verifique se o script 03 foi concluído com sucesso.")

# =============================================================================
# 2. Métricas por célula
# =============================================================================
log_step(2, "Calculando métricas de evolução clonal por célula...")

compute_tree_metrics <- function(run, qc_metrics, patient_col) {
  tree     <- run$tree
  run_name <- run$run_name

  cell_tips <- tree$tip.label[tree$tip.label != "diploid_normal"]
  if (length(cell_tips) == 0) {
    warning("[", run_name, "] Sem células tumorais nas tips.")
    return(NULL)
  }

  # Distância da raiz por célula
  root_dists <- ape::dist.nodes(tree)
  root_node  <- ape::Ntip(tree) + 1L
  tip_idx    <- match(cell_tips, tree$tip.label)
  evo_dist   <- root_dists[tip_idx, root_node]
  names(evo_dist) <- cell_tips

  # Clustering clonal via hclust na matriz de distâncias MEDICC2
  clone_clusters <- NULL
  if (!is.null(run$dist_mat)) {
    dm_cells <- intersect(cell_tips, rownames(run$dist_mat))
    if (length(dm_cells) >= 5) {
      dm_sub <- as.dist(run$dist_mat[dm_cells, dm_cells])
      hc     <- hclust(dm_sub, method = "ward.D2")
      k      <- max(2L, min(10L, as.integer(sqrt(length(dm_cells) / 2))))
      clone_clusters <- paste0("clone_", cutree(hc, k = k))
      names(clone_clusters) <- dm_cells
      log_info(sprintf("[%s] %d clones (k = %d)", run_name, k,
                       length(unique(clone_clusters))))
    }
  }

  res_df <- data.frame(
    cell             = cell_tips,
    run              = run_name,
    evo_dist_to_root = evo_dist[cell_tips],
    clone_cluster    = if (!is.null(clone_clusters))
                         clone_clusters[cell_tips] else NA_character_,
    patient          = qc_metrics[[patient_col]][
                         match(cell_tips, qc_metrics$cell)],
    aberration_score = qc_metrics$aberration_score[
                         match(cell_tips, qc_metrics$cell)],
    stringsAsFactors = FALSE,
    row.names        = NULL
  )

  # Validação biológica: correlação Spearman entre distância evolutiva e
  # aberração CNV (esperado > 0.2 para sinal real)
  corr <- tryCatch(
    cor(res_df$evo_dist_to_root, res_df$aberration_score,
        use = "complete.obs", method = "spearman"),
    error = function(e) NA_real_
  )
  log_info(sprintf("[%s] Spearman(evo_dist, aberração CNV) = %.3f%s",
                   run_name, corr,
                   if (!is.na(corr) && corr < 0.1)
                     "  ← ATENÇÃO: correlação baixa — verifique referência diploide"
                   else "  OK"))

  # CORREÇÃO: remover apenas o atributo extra, sem destruir o data.frame
  # (attributes(xa) <- NULL strips names/class e quebra o rbindlist)
  attr(res_df, "validation_corr") <- corr
  res_df
}

tree_metrics <- lapply(runs, compute_tree_metrics,
                       qc_metrics = qc_metrics, patient_col = patient_col)
tree_metrics <- Filter(Negate(is.null), tree_metrics)

# CORREÇÃO: remover apenas o atributo "validation_corr" antes do rbindlist,
# preservando names, row.names e class do data.frame
metrics_all <- rbindlist(
  lapply(tree_metrics, function(x) {
    attr(x, "validation_corr") <- NULL
    x
  }),
  fill = TRUE
)
setDF(metrics_all)

fwrite(metrics_all, "results/tree_analysis/evo_metrics_per_cell.csv")
log_info(sprintf("OK results/tree_analysis/evo_metrics_per_cell.csv (%d células)",
                 nrow(metrics_all)))

# =============================================================================
# 3. Resumos
# =============================================================================
log_step(3, "Resumos por clone e paciente...")

summary_clone <- metrics_all %>%
  filter(!is.na(clone_cluster)) %>%
  group_by(run, clone_cluster) %>%
  summarise(
    n_cells          = n(),
    evo_dist_median  = median(evo_dist_to_root, na.rm = TRUE),
    evo_dist_mean    = mean(evo_dist_to_root,   na.rm = TRUE),
    abr_median       = median(aberration_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(run, evo_dist_median)

fwrite(summary_clone, "results/tree_analysis/clone_summary.csv")
cat("\nResumo por clone:\n")
print(as.data.frame(summary_clone))

summary_patient <- metrics_all %>%
  group_by(run, patient) %>%
  summarise(
    n_cells          = n(),
    evo_dist_median  = median(evo_dist_to_root, na.rm = TRUE),
    evo_dist_max     = max(evo_dist_to_root,    na.rm = TRUE),
    .groups = "drop"
  )

fwrite(summary_patient, "results/tree_analysis/patient_summary.csv")
cat("\nResumo por paciente:\n")
print(as.data.frame(summary_patient))

# =============================================================================
# 4. Salvar
# =============================================================================
log_step(4, "Salvando objeto de análise...")

saveRDS(
  list(
    runs         = runs,
    tree_metrics = tree_metrics,
    metrics_all  = metrics_all,
    qc_obj       = qc_obj
  ),
  "data/processed/tree_analysis.rds"
)
log_info("OK data/processed/tree_analysis.rds")

cat("\n=================================================================\n")
cat(sprintf(" Etapa 4 concluída: %s\n", format(Sys.time())))
cat("=================================================================\n\n")

cat("Arquivos gerados:\n")
for (f in c(
  "data/processed/tree_analysis.rds",
  "results/tree_analysis/evo_metrics_per_cell.csv",
  "results/tree_analysis/clone_summary.csv",
  "results/tree_analysis/patient_summary.csv"
)) {
  if (file.exists(f))
    cat(sprintf("  [OK] %-60s  %.1f MB\n", f, file.size(f) / 1024^2))
  else
    cat(sprintf("  [!!] %-60s  NAO GERADO\n", f))
}
cat("\n")
