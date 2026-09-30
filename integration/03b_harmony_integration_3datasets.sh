#!/bin/bash
# =============================================================================
# SCRIPT: 03_harmony_integration_v3.sh
# PROPOSITO:
#   - integrar GSE182109 + Synapse (ambos com counts brutos reais)
#   - PCA + Harmony + UMAP + clustering
#   - GSE173278 NAO entra aqui (sera usado no job 04 para label transfer)
# =============================================================================
#SBATCH --job-name=gbm_harmony_v3
#SBATCH --output=/home/renanomete/projetos/matdata/met_mat_data/singlecell_novo/logs_v3/03_harmony_v3_%j.out
#SBATCH --error=/home/renanomete/projetos/matdata/met_mat_data/singlecell_novo/logs_v3/03_harmony_v3_%j.err
#SBATCH --chdir=/home/renanomete/projetos/matdata/met_mat_data/singlecell_novo
#SBATCH -t 24:00:00
#SBATCH -c 4
#SBATCH --mem=120G
 
set -euo pipefail
 
BASE=/home/renanomete/projetos/matdata/met_mat_data/singlecell_novo
RESULTS=${BASE}/results_v3
LOGS=${BASE}/logs_v3
 
cd "${BASE}"
mkdir -p "${RESULTS}" "${LOGS}"
 
export INPUT="${RESULTS}/gse_qc.rds ${RESULTS}/synapse_qc.rds"
export OUTPUT="${RESULTS}/integrated_harmony_v3.rds ${RESULTS}/pca_elbow_v3.pdf"
 
module load miniconda/24.4.0-libmamba
source activate gbm_scrnaseq
 
echo "[$(date '+%F %T')] Job iniciado em $(hostname)"
echo "[$(date '+%F %T')] Diretorio de execucao: $(pwd)"
 
job-nanny Rscript - <<'REOF'
suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
})
 
BASE    <- "/home/renanomete/projetos/matdata/met_mat_data/singlecell_novo"
RESULTS <- file.path(BASE, "results_v3")
 
# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------
ensure_metadata_columns <- function(meta, dataset_label = NULL) {
  meta <- as.data.frame(meta)
  if (!"dataset" %in% colnames(meta))
    meta$dataset <- if (!is.null(dataset_label)) dataset_label else "unknown_dataset"
  if (!"Patient" %in% colnames(meta)) {
    meta$Patient <- if ("orig.ident" %in% colnames(meta)) meta$orig.ident else "unknown_patient"
  } else {
    idx <- is.na(meta$Patient) | meta$Patient == ""
    if (any(idx)) meta$Patient[idx] <- if ("orig.ident" %in% colnames(meta)) meta$orig.ident[idx] else "unknown_patient"
  }
  if (!"Type" %in% colnames(meta))              meta$Type              <- NA_character_
  if (!"idh_codel_subtype" %in% colnames(meta)) meta$idh_codel_subtype <- NA_character_
  meta
}
 
get_counts_matrix <- function(so, assay = "RNA", label = "objeto") {
  if (!assay %in% names(so@assays)) stop("Assay ", assay, " nao encontrado em ", label)
  DefaultAssay(so) <- assay
  if (inherits(so[[assay]], "Assay5")) so <- tryCatch(JoinLayers(so, assay=assay), error=function(e) so)
  m <- tryCatch(GetAssayData(so, assay=assay, layer="counts"), error=function(e) NULL)
  if (!is.null(m)) return(m)
  m <- tryCatch(suppressWarnings(GetAssayData(so, assay=assay, slot="counts")), error=function(e) NULL)
  if (!is.null(m)) return(m)
  m <- tryCatch(slot(so[[assay]], "counts"), error=function(e) NULL)
  if (!is.null(m)) return(m)
  stop("Nao foi possivel acessar counts no assay ", assay, " para ", label)
}
 
standardize_gene_names <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("\\.[0-9]+$", "", x)
  x <- gsub("_", "-", x, fixed = TRUE)
  toupper(x)
}
 
rebuild_minimal_rna_object <- function(so, dataset_label, label_print) {
  cat(sprintf("[LOAD] Reconstruindo %s...\n", label_print))
  meta   <- ensure_metadata_columns(so@meta.data, dataset_label = dataset_label)
  counts <- get_counts_matrix(so, assay = "RNA", label = label_print)
  if (!inherits(counts, "dgCMatrix")) counts <- as(counts, "CsparseMatrix")
 
  new_genes <- standardize_gene_names(rownames(counts))
  keep <- !duplicated(new_genes) & !is.na(new_genes) & new_genes != ""
  if (sum(!keep) > 0)
    cat(sprintf("  [%s] Removendo %d genes duplicados/invalidos.\n", label_print, sum(!keep)))
  counts <- counts[keep, , drop = FALSE]
  rownames(counts) <- new_genes[keep]
  meta <- meta[colnames(counts), , drop = FALSE]
 
  so_min <- CreateSeuratObject(counts=counts, meta.data=meta,
                               project=dataset_label, min.cells=1, min.features=1)
  DefaultAssay(so_min) <- "RNA"
  cat(sprintf("  %s: %d celulas, %d genes\n", label_print, ncol(so_min), nrow(so_min)))
  so_min
}
 
prepare_object <- function(path, dataset_label, label_print) {
  cat(sprintf("[READ] Lendo %s...\n", basename(path)))
  so_raw <- readRDS(path)
  cat(sprintf("  %s original: %d celulas, %d genes\n", label_print, ncol(so_raw), nrow(so_raw)))
 
  so_min <- rebuild_minimal_rna_object(so_raw, dataset_label, label_print)
  rm(so_raw); gc()
 
  DefaultAssay(so_min) <- "RNA"
  so_min <- NormalizeData(so_min, verbose = FALSE)
  so_min <- FindVariableFeatures(so_min, assay="RNA", selection.method="vst",
                                  nfeatures=3000, verbose=FALSE)
 
  hvgs <- VariableFeatures(so_min)
  hvgs <- unique(hvgs[!is.na(hvgs) & hvgs != ""])
  if (length(hvgs) == 0) stop("Nenhuma VariableFeature para ", label_print)
  cat(sprintf("  %s preparado: %d celulas, %d HVGs\n", label_print, ncol(so_min), length(hvgs)))
  list(obj = so_min, hvgs = hvgs)
}
 
# ------------------------------------------------------------------
# [1/7] Preparando objetos — apenas GSE182109 e Synapse
# ------------------------------------------------------------------
cat("[1/7] Preparando objetos (GSE182109 + Synapse)...\n")
 
gse_prep <- prepare_object(file.path(RESULTS, "gse_qc.rds"),    "GSE182109", "GSE182109")
syn_prep <- prepare_object(file.path(RESULTS, "synapse_qc.rds"), "Synapse",   "Synapse")
 
gse_obj  <- gse_prep$obj; gse_hvgs <- gse_prep$hvgs; rm(gse_prep); gc()
syn_obj  <- syn_prep$obj; syn_hvgs <- syn_prep$hvgs; rm(syn_prep); gc()
 
cat(sprintf("  GSE182109: %d celulas\n", ncol(gse_obj)))
cat(sprintf("  Synapse:   %d celulas\n", ncol(syn_obj)))
 
# ------------------------------------------------------------------
# [2/7] Selecionando features
# ------------------------------------------------------------------
cat("[2/7] Selecionando features...\n")
 
obj_list <- list(gse_obj, syn_obj)
hvgs <- tryCatch(SelectIntegrationFeatures(object.list=obj_list, nfeatures=3000),
                 error=function(e) character(0))
cat(sprintf("  HVGs por SelectIntegrationFeatures: %d\n", length(hvgs)))
 
if (length(hvgs) < 500) {
  cat("  [AVISO] Usando uniao dos HVGs individuais.\n")
  hvgs <- unique(c(gse_hvgs, syn_hvgs))
}
if (length(hvgs) < 500) {
  cat("  [AVISO] Usando genes presentes em ambos os datasets.\n")
  hvgs <- intersect(rownames(gse_obj), rownames(syn_obj))
}
 
hvgs <- unique(hvgs[!is.na(hvgs) & hvgs != ""])
hvgs <- hvgs[hvgs %in% c(rownames(gse_obj), rownames(syn_obj))]
hvgs <- head(hvgs, 3000)
cat(sprintf("  Features finais: %d\n", length(hvgs)))
if (length(hvgs) < 500) stop("Poucas features apos fallback.")
 
# ------------------------------------------------------------------
# [3/7] Merge
# ------------------------------------------------------------------
cat("[3/7] Merge dos objetos...\n")
 
combined <- merge(gse_obj, syn_obj)
rm(obj_list, gse_obj, syn_obj, gse_hvgs, syn_hvgs); gc()
 
cat(sprintf("  Combinado: %d celulas, %d genes\n", ncol(combined), nrow(combined)))
DefaultAssay(combined) <- "RNA"
 
if (inherits(combined[["RNA"]], "Assay5")) {
  combined <- tryCatch(JoinLayers(combined, assay="RNA"),
    error = function(e) { cat(sprintf("  [AVISO] JoinLayers: %s\n", conditionMessage(e))); combined })
}
 
# ------------------------------------------------------------------
# [4/7] Batch + IDH status
# ------------------------------------------------------------------
cat("[4/7] Definindo batch e IDH status...\n")
 
synapse_idh_map <- c(
  "SM001"="IDH-mut","SM002"="IDH-mut","SM004"="IDH-mut","SM006"="IDH-wt",
  "SM008"="IDH-mut","SM011"="IDH-wt","SM012"="IDH-wt","SM015"="IDH-mut",
  "SM017"="IDH-wt","SM018"="IDH-wt","SM019"="IDH-mut"
)
 
combined$IDH_status <- dplyr::case_when(
  combined$dataset == "GSE182109" & grepl("^LGG$",  combined$Type, ignore.case=TRUE) ~ "IDH-mut",
  combined$dataset == "GSE182109" & grepl("GBM",    combined$Type, ignore.case=TRUE) ~ "IDH-wt",
  combined$dataset == "Synapse"   & combined$Patient %in% names(synapse_idh_map) ~
    unname(synapse_idh_map[combined$Patient]),
  TRUE ~ "Unknown"
)
combined$IDH_status <- factor(combined$IDH_status, levels=c("IDH-mut","IDH-wt","Unknown"))
cat("  IDH_status:\n"); print(table(combined$IDH_status, useNA="ifany"))
 
combined$batch <- paste(combined$dataset, combined$Patient, sep="_")
cat("  Batches: "); cat(length(unique(combined$batch)), "\n")
 
# ------------------------------------------------------------------
# [5/7] ScaleData + PCA
# ------------------------------------------------------------------
cat("[5/7] ScaleData + PCA...\n")
 
features_use <- intersect(hvgs, rownames(combined))
cat(sprintf("  Features usadas: %d\n", length(features_use)))
if (length(features_use) < 500) stop("Poucas features para ScaleData/PCA.")
 
VariableFeatures(combined) <- features_use
 
candidate_regress <- intersect(c("percent.mt","nCount_RNA"), colnames(combined@meta.data))
vars_regress <- candidate_regress[sapply(candidate_regress, function(v) {
  vals <- combined@meta.data[[v]][!is.na(combined@meta.data[[v]])]
  length(unique(vals)) > 1
})]
cat("  Variaveis a regredir:\n"); print(vars_regress)
 
combined <- ScaleData(combined, features=features_use, vars.to.regress=vars_regress,
                      do.scale=TRUE, do.center=TRUE, verbose=FALSE)
combined <- RunPCA(combined, features=features_use, npcs=50, verbose=FALSE)
 
pdf(file.path(RESULTS, "pca_elbow_v3.pdf"), width=8, height=5)
print(ElbowPlot(combined, ndims=50) + ggtitle("Elbow Plot - GSE182109 + Synapse V3"))
dev.off()
 
# ------------------------------------------------------------------
# [6/7] Harmony
# ------------------------------------------------------------------
cat("[6/7] Rodando Harmony...\n")
 
combined <- RunHarmony(
  object=combined, group.by.vars="batch",
  reduction.use="pca", dims.use=1:30,
  reduction.save="harmony", plot_convergence=FALSE, verbose=TRUE
)
 
# ------------------------------------------------------------------
# [7/7] Vizinhanca, clustering, UMAP
# ------------------------------------------------------------------
cat("[7/7] Vizinhanca, clustering e UMAP...\n")
 
N_DIMS <- 30
combined <- FindNeighbors(combined, reduction="harmony", dims=1:N_DIMS, verbose=FALSE)
 
for (res in c(0.2, 0.4, 0.6, 0.8, 1.0)) {
  combined <- FindClusters(combined, resolution=res, algorithm=1, verbose=FALSE)
  cc <- paste0("RNA_snn_res.", res)
  cat(sprintf("  Resolucao %.1f -> %d clusters\n", res, length(unique(combined@meta.data[[cc]]))))
}
Idents(combined) <- "RNA_snn_res.0.6"
 
combined <- RunUMAP(combined, reduction="harmony", dims=1:N_DIMS,
                    min.dist=0.3, spread=1, verbose=FALSE)
 
saveRDS(combined, file.path(RESULTS, "integrated_harmony_v3.rds"))
cat(sprintf("Salvo: results_v3/integrated_harmony_v3.rds (%d celulas)\n", ncol(combined)))
cat("\n[DONE] Integracao Harmony V3 concluida.\n")
REOF
