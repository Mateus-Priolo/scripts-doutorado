#!/bin/bash
# =============================================================================
# SCRIPT: 02_qc_and_normalize_v3.sh
# PROPOSITO:
#   - QC/SCT para GSE182109 e Synapse (serao integrados via Harmony)
#   - QC leve para GSE173278 (apenas para label transfer downstream)
# =============================================================================
#SBATCH --job-name=gbm_qc_v3
#SBATCH --output=/home/renanomete/projetos/matdata/met_mat_data/singlecell_novo/logs_v3/02_qc_v3_%j.out
#SBATCH --error=/home/renanomete/projetos/matdata/met_mat_data/singlecell_novo/logs_v3/02_qc_v3_%j.err
#SBATCH --chdir=/home/renanomete/projetos/matdata/met_mat_data/singlecell_novo
#SBATCH -t 24:00:00
#SBATCH -c 4
#SBATCH --mem=180G
 
set -euo pipefail
 
BASE=/home/renanomete/projetos/matdata/met_mat_data/singlecell_novo
RESULTS=${BASE}/results_v3
LOGS=${BASE}/logs_v3
 
cd "${BASE}"
mkdir -p "${RESULTS}" "${LOGS}"
 
export INPUT="${RESULTS}/gse_raw.rds ${RESULTS}/synapse_raw.rds ${RESULTS}/gse173278_primary_raw.rds"
export OUTPUT="${RESULTS}/gse_qc.rds ${RESULTS}/synapse_qc.rds ${RESULTS}/gse173278_primary_qc.rds"
 
module load miniconda/24.4.0-libmamba
source activate gbm_scrnaseq
 
echo "[$(date '+%F %T')] Job iniciado em $(hostname)"
echo "[$(date '+%F %T')] Diretorio de execucao: $(pwd)"
 
job-nanny Rscript - <<'REOF'
suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(scDblFinder)
  library(SingleCellExperiment)
  library(BiocParallel)
  library(Matrix)
  library(patchwork)
})
 
BASE    <- "/home/renanomete/projetos/matdata/met_mat_data/singlecell_novo"
RESULTS <- file.path(BASE, "results_v3")
 
N_WORKERS <- 4
options(future.globals.maxSize = 60 * 1024^3)
 
get_counts_matrix <- function(so, assay = "RNA") {
  m <- tryCatch(GetAssayData(so, assay = assay, layer = "counts"), error = function(e) NULL)
  if (!is.null(m)) return(m)
  m <- tryCatch(suppressWarnings(GetAssayData(so, assay = assay, slot = "counts")), error = function(e) NULL)
  if (!is.null(m)) return(m)
  m <- tryCatch(slot(so[[assay]], "counts"), error = function(e) NULL)
  if (!is.null(m) && length(m) > 0) return(m)
  NULL
}
 
prepare_rna_assay <- function(so, label) {
  if (!"RNA" %in% names(so@assays)) stop("Assay RNA nao encontrado em ", label)
  DefaultAssay(so) <- "RNA"
  if (inherits(so[["RNA"]], "Assay5")) {
    cat(sprintf("  [%s] Assay5 detectado; JoinLayers...\n", label))
    so <- JoinLayers(so, assay = "RNA")
  }
  so
}
 
compute_qc_metrics <- function(so) {
  counts_mat <- get_counts_matrix(so, assay = "RNA")
  if (is.null(counts_mat)) stop("Nao foi possivel acessar contagens no assay RNA")
  so$nCount_RNA   <- Matrix::colSums(counts_mat)
  so$nFeature_RNA <- Matrix::colSums(counts_mat > 0)
  so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
  so[["percent.rb"]] <- PercentageFeatureSet(so, pattern = "^RP[SL]")
  if ("mitoRatio" %in% colnames(so@meta.data)) {
    mito_num <- suppressWarnings(as.numeric(so$mitoRatio))
    idx <- !is.na(mito_num)
    if (any(idx)) so$percent.mt[idx] <- mito_num[idx] * 100
  }
  if ("pct_counts_mitochondrial" %in% colnames(so@meta.data)) {
    mito_num2 <- suppressWarnings(as.numeric(so$pct_counts_mitochondrial))
    idx2 <- !is.na(mito_num2)
    if (any(idx2)) so$percent.mt[idx2] <- mito_num2[idx2]
  }
  so
}
 
plot_qc_violin <- function(so, out_pdf, title_text) {
  plot_features <- c("nFeature_RNA", "nCount_RNA")
  if ("percent.mt" %in% colnames(so@meta.data)) {
    vals <- so$percent.mt
    if (length(unique(vals[!is.na(vals)])) > 1) {
      plot_features <- c(plot_features, "percent.mt")
    }
  }
  pdf(out_pdf, width = 4 * length(plot_features), height = 5)
  print(VlnPlot(so, features = plot_features, ncol = length(plot_features), pt.size = 0) +
        plot_annotation(title = title_text))
  dev.off()
}
 
run_doublet_detection <- function(so, n_workers = 4) {
  tryCatch({
    sce <- as.SingleCellExperiment(so)
    sample_col <- if ("orig.ident" %in% colnames(colData(sce))) "orig.ident" else NULL
    bp  <- MulticoreParam(workers = n_workers)
    sce <- if (!is.null(sample_col)) scDblFinder(sce, samples = sample_col, BPPARAM = bp) else
                                     scDblFinder(sce, BPPARAM = bp)
    so$dbl_class <- sce$scDblFinder.class
    so$dbl_score <- sce$scDblFinder.score
    n_dbl <- sum(so$dbl_class == "doublet", na.rm = TRUE)
    cat(sprintf("    Doublets: %d (%.2f%%)\n", n_dbl, 100 * n_dbl / ncol(so)))
    subset(so, subset = dbl_class == "singlet")
  }, error = function(e) {
    cat(sprintf("    [AVISO] scDblFinder falhou: %s\n", conditionMessage(e)))
    so
  })
}
 
run_sct <- function(so) {
  DefaultAssay(so) <- "RNA"
  so <- NormalizeData(so, verbose = FALSE)
  so <- SCTransform(so, assay = "RNA", vst.flavor = "v2",
                    vars.to.regress = "percent.mt", variable.features.n = 3000,
                    method = "glmGamPoi", conserve.memory = TRUE, verbose = FALSE)
  so
}
 
run_qc_single_object <- function(so, label, min_genes=300, max_genes=8000,
                                  max_mito=0.25, do_doublets=TRUE, n_workers=4,
                                  make_plot=TRUE) {
  cat(sprintf("\n==== QC: %s ====\n", label))
  cat(sprintf("  Entrada: %d celulas, %d genes\n", ncol(so), nrow(so)))
  so <- prepare_rna_assay(so, label)
  so <- compute_qc_metrics(so)
 
  keep_na <- !is.na(so$nFeature_RNA) & !is.na(so$nCount_RNA) & !is.na(so$percent.mt)
  n_before_na <- ncol(so)
  so <- subset(so, cells = colnames(so)[keep_na])
  cat(sprintf("  Apos remocao de NAs: %d celulas\n", ncol(so)))
  if (ncol(so) == 0) stop("Nenhuma celula restou apos NAs para ", label)
 
  if (make_plot) {
    plot_qc_violin(so,
      file.path(RESULTS, paste0(label, "_qc_antes.pdf")),
      paste("QC antes da filtragem -", label))
  }
 
  n_before <- ncol(so)
  so <- subset(so, subset = nFeature_RNA > min_genes &
                             nFeature_RNA < max_genes &
                             percent.mt < (max_mito * 100))
  cat(sprintf("  Apos filtros QC: %d celulas (removidas: %d)\n",
              ncol(so), n_before - ncol(so)))
  if (ncol(so) == 0) stop("Nenhuma celula restou apos filtros QC para ", label)
 
  if (do_doublets) {
    cat("  Detectando doublets...\n")
    so <- run_doublet_detection(so, n_workers = n_workers)
    cat(sprintf("  Celulas apos doublets: %d\n", ncol(so)))
  }
  if (ncol(so) == 0) stop("Nenhuma celula restou apos doublets para ", label)
 
  cat("  Rodando SCTransform v2...\n")
  so <- run_sct(so)
  cat(sprintf("  SCTransform concluido. HVGs: %d\n", length(VariableFeatures(so))))
  so
}
 
# -------------------------------------------------------------------
# GSE182109
# -------------------------------------------------------------------
cat("\n==============================\n")
cat("Processando GSE182109\n")
cat("==============================\n")
 
gse <- readRDS(file.path(RESULTS, "gse_raw.rds"))
gse <- prepare_rna_assay(gse, "GSE182109")
gse <- compute_qc_metrics(gse)
plot_qc_violin(gse, file.path(RESULTS, "GSE182109_qc_antes.pdf"),
               "QC antes da filtragem - GSE182109")
 
if (!"orig.ident" %in% colnames(gse@meta.data)) stop("orig.ident nao encontrado em gse_raw.rds")
 
gse_list    <- SplitObject(gse, split.by = "orig.ident")
cat(sprintf("  Amostras para QC/SCT: %d\n", length(gse_list)))
gse_qc_list <- vector("list", length(gse_list))
names(gse_qc_list) <- names(gse_list)
 
for (nm in names(gse_list)) {
  cat(sprintf("\n---- Amostra: %s ----\n", nm))
  gse_qc_list[[nm]] <- run_qc_single_object(
    gse_list[[nm]], label = nm, min_genes = 300, max_genes = 8000,
    max_mito = 0.25, do_doublets = TRUE, n_workers = N_WORKERS, make_plot = FALSE)
  gc()
}
 
cat("\nUnindo amostras do GSE...\n")
gse_qc <- merge(gse_qc_list[[1]], y = gse_qc_list[-1],
                project = "GSE182109_QC", merge.data = TRUE)
saveRDS(gse_qc, file.path(RESULTS, "gse_qc.rds"))
cat(sprintf("Salvo: results_v3/gse_qc.rds (%d celulas)\n", ncol(gse_qc)))
rm(gse, gse_list, gse_qc_list, gse_qc); gc()
 
# -------------------------------------------------------------------
# Synapse
# -------------------------------------------------------------------
cat("\n==============================\n")
cat("Processando Synapse\n")
cat("==============================\n")
 
syn <- readRDS(file.path(RESULTS, "synapse_raw.rds"))
if (!"Patient" %in% colnames(syn@meta.data)) stop("Patient nao encontrado em synapse_raw.rds")
cat("Patient antes do QC:\n")
print(table(syn$Patient, useNA = "ifany"))
 
syn_qc <- run_qc_single_object(
  syn, label = "Synapse", min_genes = 300, max_genes = 9000,
  max_mito = 0.20, do_doublets = TRUE, n_workers = N_WORKERS, make_plot = TRUE)
 
cat("Patient apos QC:\n")
print(table(syn_qc$Patient, useNA = "ifany"))
saveRDS(syn_qc, file.path(RESULTS, "synapse_qc.rds"))
cat(sprintf("Salvo: results_v3/synapse_qc.rds (%d celulas)\n", ncol(syn_qc)))
rm(syn, syn_qc); gc()
 
# -------------------------------------------------------------------
# GSE173278 primary — QC leve, sem SCTransform
# Sera usado APENAS para label transfer no job 04
# NAO entra na integracao Harmony
# -------------------------------------------------------------------
cat("\n==============================\n")
cat("Processando GSE173278 (label transfer only)\n")
cat("==============================\n")
 
g173 <- readRDS(file.path(RESULTS, "gse173278_primary_raw.rds"))
g173 <- prepare_rna_assay(g173, "GSE173278_primary")
 
# Metricas de QC usando o slot counts (matriz binaria placeholder)
g173$nCount_RNA   <- Matrix::colSums(GetAssayData(g173, assay = "RNA", layer = "counts"))
g173$nFeature_RNA <- Matrix::colSums(GetAssayData(g173, assay = "RNA", layer = "counts") > 0)
g173[["percent.mt"]] <- PercentageFeatureSet(g173, pattern = "^MT-")
 
plot_qc_violin(g173, file.path(RESULTS, "GSE173278_primary_qc_antes.pdf"),
               "QC before filtering - GSE173278 primary (label transfer only)")
 
# Filtro leve — mantemos mais celulas pois e referencia para label transfer
n_before <- ncol(g173)
keep_cells <- !is.na(g173$nFeature_RNA) & g173$nFeature_RNA > 200
g173 <- subset(g173, cells = colnames(g173)[keep_cells])
cat(sprintf("  Celulas mantidas: %d de %d\n", ncol(g173), n_before))
 
# FindVariableFeatures no slot data (normalizacao original do provedor)
DefaultAssay(g173) <- "RNA"
g173 <- FindVariableFeatures(g173, assay = "RNA", selection.method = "vst",
                              nfeatures = 3000, verbose = FALSE)
cat(sprintf("  HVGs: %d\n", length(VariableFeatures(g173))))
 
saveRDS(g173, file.path(RESULTS, "gse173278_primary_qc.rds"))
cat(sprintf("Salvo: results_v3/gse173278_primary_qc.rds (%d celulas)\n", ncol(g173)))
 
cat("\n[DONE] QC V3 concluido.\n")
REOF
