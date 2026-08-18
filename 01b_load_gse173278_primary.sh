#!/bin/bash
# =============================================================================
# SCRIPT: 01_prepare_raw_v3.sh
# PROPOSITO:
#   - copiar gse_raw.rds ja validado do pipeline anterior
#   - reconstruir synapse_raw.rds com Patient correto (SM001...SM019)
#   - carregar GSE173278 primary GBM tissue (APENAS para label transfer,
#     NAO sera integrado via Harmony — dados ja normalizados pelo provedor)
# =============================================================================
#SBATCH --job-name=gbm_raw_v3
#SBATCH --output=/home/renanomete/projetos/matdata/met_mat_data/singlecell_novo/logs_v2/01_raw_v3_%j.out
#SBATCH --error=/home/renanomete/projetos/matdata/met_mat_data/singlecell_novo/logs_v2/01_raw_v3_%j.err
#SBATCH --chdir=/home/renanomete/projetos/matdata/met_mat_data/singlecell_novo
#SBATCH -t 12:00:00
#SBATCH -c 4
#SBATCH --mem=96G
 
set -euo pipefail
 
BASE=/home/renanomete/projetos/matdata/met_mat_data/singlecell_novo
RESULTS_OLD=${BASE}/results
RESULTS=${BASE}/results_v3
LOGS=${BASE}/logs_v3
SYNAPSE_DIR=${BASE}/synapse_data
G173_DIR=${BASE}/gse173278
 
cd "${BASE}"
mkdir -p "${RESULTS}" "${LOGS}" "${G173_DIR}" "${SYNAPSE_DIR}"
 
export INPUT="${RESULTS_OLD}/gse_raw.rds \
${SYNAPSE_DIR}/analysis_scRNAseq_tumor_counts.h5ad \
${SYNAPSE_DIR}/analysis_scRNAseq_tumor_gene_expression.tsv.gz \
${SYNAPSE_DIR}/41588_2021_926_MOESM2_ESM.xlsx \
${G173_DIR}/GSE173278_scRNAseq_filtered_cells_barcodes.tsv.gz \
${G173_DIR}/GSE173278_scRNAseq_filtered_cells_genes.tsv.gz \
${G173_DIR}/GSE173278_scRNAseq_filtered_cells_metadata.csv.gz \
${G173_DIR}/GSE173278_scRNAseq_filtered_cells_norm_counts_matrix.mtx.gz"
export OUTPUT="${RESULTS}/gse_raw.rds \
${RESULTS}/synapse_raw.rds \
${RESULTS}/gse173278_primary_raw.rds"
 
module load miniconda/24.4.0-libmamba
source activate gbm_scrnaseq
 
echo "[$(date '+%F %T')] Job iniciado em $(hostname)"
echo "[$(date '+%F %T')] Diretorio de execucao: $(pwd)"
 
python - <<'PY'
import os
import anndata as ad
 
base = "/home/renanomete/projetos/matdata/met_mat_data/singlecell_novo"
synapse_dir = os.path.join(base, "synapse_data")
h5ad_path = os.path.join(synapse_dir, "analysis_scRNAseq_tumor_counts.h5ad")
obs_csv = os.path.join(synapse_dir, "analysis_scRNAseq_tumor_counts_obs.csv")
 
if not os.path.exists(h5ad_path):
    raise SystemExit(f"Arquivo nao encontrado: {h5ad_path}")
 
adata = ad.read_h5ad(h5ad_path, backed="r")
obs = adata.obs.copy()
obs["cell_name"] = obs.index.astype(str)
obs.to_csv(obs_csv, index=False)
 
print(f"OBS exportado para: {obs_csv}")
print("Colunas do obs:")
print(list(obs.columns))
PY
 
job-nanny Rscript - <<'REOF'
suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(Matrix)
  library(dplyr)
  library(readxl)
  library(stringr)
})
 
BASE        <- "/home/renanomete/projetos/matdata/met_mat_data/singlecell_novo"
RESULTS_OLD <- file.path(BASE, "results")
RESULTS     <- file.path(BASE, "results_v3")
SYNAPSE_DIR <- file.path(BASE, "synapse_data")
G173_DIR    <- file.path(BASE, "gse173278")
 
dir.create(RESULTS, showWarnings = FALSE, recursive = TRUE)
 
find_col <- function(df, candidates) {
  cn <- colnames(df)
  idx <- match(tolower(candidates), tolower(cn), nomatch = 0)
  if (any(idx > 0)) return(cn[idx[idx > 0][1]])
  for (cand in candidates) {
    hit <- cn[grepl(tolower(cand), tolower(cn), fixed = TRUE)]
    if (length(hit) > 0) return(hit[1])
  }
  NULL
}
 
normalize_cell_key <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- sub("-[0-9]+$", "", x)
  x
}
 
# -------------------------------------------------------------------
# 1. Reaproveita GSE182109 raw ja validado
# -------------------------------------------------------------------
cat("[1/3] Copiando gse_raw.rds validado para results_v3...\n")
src_gse <- file.path(RESULTS_OLD, "gse_raw.rds")
dst_gse <- file.path(RESULTS, "gse_raw.rds")
if (!file.exists(src_gse)) stop("Arquivo nao encontrado: ", src_gse)
ok <- file.copy(src_gse, dst_gse, overwrite = TRUE)
if (!ok) stop("Falha ao copiar gse_raw.rds para results_v3/")
cat("  Copiado: results_v3/gse_raw.rds\n")
 
# -------------------------------------------------------------------
# 2. Reconstroi Synapse raw com Patient correto
# -------------------------------------------------------------------
cat("\n[2/3] Reconstruindo Synapse raw com Patient correto...\n")
 
counts_h5ad <- file.path(SYNAPSE_DIR, "analysis_scRNAseq_tumor_counts.h5ad")
counts_tsv  <- file.path(SYNAPSE_DIR, "analysis_scRNAseq_tumor_gene_expression.tsv.gz")
obs_csv     <- file.path(SYNAPSE_DIR, "analysis_scRNAseq_tumor_counts_obs.csv")
 
if (!file.exists(obs_csv)) stop("OBS CSV do Synapse nao encontrado: ", obs_csv)
 
if (file.exists(counts_h5ad) && requireNamespace("SeuratDisk", quietly = TRUE)) {
  suppressPackageStartupMessages(library(SeuratDisk))
  h5seurat_path <- file.path(SYNAPSE_DIR, "analysis_scRNAseq_tumor_counts.h5seurat")
  if (!file.exists(h5seurat_path)) {
    cat("  Convertendo h5ad -> h5seurat...\n")
    Convert(counts_h5ad, dest = "h5seurat", overwrite = TRUE)
  }
  syn_so <- LoadH5Seurat(h5seurat_path, meta.data = TRUE, reductions = FALSE)
  cat("  Synapse lido via SeuratDisk (h5ad)\n")
} else if (file.exists(counts_tsv)) {
  cat("  Lendo TSV de contagens do Synapse...\n")
  mat_raw <- fread(counts_tsv, sep = "\t", header = TRUE, data.table = FALSE)
  genes   <- mat_raw[[1]]
  mat_raw <- mat_raw[, -1, drop = FALSE]
  rownames(mat_raw) <- make.unique(genes)
  syn_so <- CreateSeuratObject(
    counts = as(as.matrix(mat_raw), "sparseMatrix"),
    project = "Synapse_GBM", min.cells = 3, min.features = 200
  )
  cat("  Synapse lido via TSV\n")
} else {
  stop("Nenhum arquivo de contagens do Synapse encontrado em synapse_data/")
}
 
syn_so$dataset <- "Synapse"
syn_so$Type    <- "GBM_Tumor"
cat(sprintf("    %d celulas, %d genes\n", ncol(syn_so), nrow(syn_so)))
 
obs_df <- fread(obs_csv, data.table = FALSE)
obs_df$cell_name <- as.character(obs_df$cell_name)
if (!"sampleid" %in% colnames(obs_df)) stop("Coluna 'sampleid' nao encontrada.")
 
cell_df  <- data.frame(cell_name = colnames(syn_so), stringsAsFactors = FALSE)
meta_obs <- left_join(cell_df, obs_df, by = "cell_name")
n_exact  <- sum(!is.na(meta_obs$sampleid))
cat(sprintf("  Match exato: %d / %d\n", n_exact, nrow(meta_obs)))
 
if (n_exact == 0) {
  cell_df$cell_key <- normalize_cell_key(cell_df$cell_name)
  obs_df$cell_key  <- normalize_cell_key(obs_df$cell_name)
  obs_df2  <- obs_df %>% group_by(cell_key) %>% slice(1) %>% ungroup()
  meta_obs <- left_join(cell_df, obs_df2, by = "cell_key", suffix = c("", ".obs"))
  cat(sprintf("  Match normalizado: %d / %d\n", sum(!is.na(meta_obs$sampleid)), nrow(meta_obs)))
}
 
if (sum(!is.na(meta_obs$sampleid)) == 0 && nrow(obs_df) == ncol(syn_so)) {
  cat("  [AVISO] Fallback por ordem.\n")
  meta_obs <- obs_df
  meta_obs$cell_name <- colnames(syn_so)
}
 
rownames(meta_obs) <- meta_obs$cell_name
meta_obs <- meta_obs[colnames(syn_so), , drop = FALSE]
if (all(is.na(meta_obs$sampleid))) stop("sampleid todo NA apos juncao.")
 
syn_meta <- syn_so@meta.data
syn_meta$cell_name <- rownames(syn_meta)
meta_obs_keep <- meta_obs
if ("cell_name" %in% colnames(meta_obs_keep)) meta_obs_keep$cell_name <- NULL
meta_obs_keep <- meta_obs_keep %>%
  tibble::rownames_to_column("cell_name") %>%
  dplyr::select(-any_of(c("dataset", "Type")))
syn_meta2 <- left_join(syn_meta, meta_obs_keep, by = "cell_name", suffix = c("", ".obs"))
rownames(syn_meta2) <- syn_meta2$cell_name
syn_meta2$cell_name <- NULL
syn_so@meta.data <- syn_meta2
 
cat("  sampleid final:\n")
print(table(syn_so$sampleid, useNA = "ifany"))
 
sample_counts <- table(as.character(syn_so$sampleid))
expected_counts <- c(
  "SM001"=5774L,"SM002"=4819L,"SM004"=4986L,"SM006"=5548L,"SM008"=4906L,
  "SM011"=5105L,"SM012"=6865L,"SM015"=4648L,"SM017"=7617L,"SM018"=2995L,"SM019"=2021L
)
 
if (length(sample_counts) != length(expected_counts)) {
  stop("Numero de grupos em sampleid (", length(sample_counts),
       ") difere do esperado (", length(expected_counts), ").")
}
 
sample_to_sm <- sapply(sample_counts, function(n) {
  hit <- names(expected_counts)[expected_counts == as.integer(n)]
  if (length(hit) == 1) return(hit)
  NA_character_
}, USE.NAMES = TRUE)
 
if (any(is.na(sample_to_sm))) stop("Falha no mapeamento sampleid -> SM.")
 
patient_vec <- unname(sample_to_sm[as.character(syn_so$sampleid)])
syn_so$Patient    <- patient_vec
syn_so$orig.ident <- patient_vec
cat("  Patient no Synapse:\n")
print(table(syn_so$Patient))
 
moesm_candidates <- c(
  file.path(SYNAPSE_DIR, "41588_2021_926_MOESM2_ESM.xlsx"),
  file.path(BASE, "41588_2021_926_MOESM2_ESM.xlsx")
)
moesm_candidates <- moesm_candidates[file.exists(moesm_candidates)]
if (length(moesm_candidates) > 0) {
  moesm_raw <- read_excel(moesm_candidates[1], sheet = "STable1", skip = 1)
  colnames(moesm_raw) <- as.character(moesm_raw[1, ])
  moesm_df <- moesm_raw[-1, , drop = FALSE]
  if ("case_barcode" %in% colnames(moesm_df)) {
    meta_joined <- left_join(syn_so@meta.data, moesm_df, by = c("Patient" = "case_barcode"))
    rownames(meta_joined) <- colnames(syn_so)
    syn_so@meta.data <- meta_joined
    cat("  Metadata clinica MOESM anexada.\n")
  } else {
    cat("  [AVISO] case_barcode nao encontrado no MOESM.\n")
  }
} else {
  cat("  [AVISO] MOESM nao encontrado.\n")
}
 
saveRDS(syn_so, file.path(RESULTS, "synapse_raw.rds"))
cat(sprintf("  Synapse salvo: %d celulas, %d genes\n", ncol(syn_so), nrow(syn_so)))
 
# -------------------------------------------------------------------
# 3. Carrega GSE173278 primary — APENAS para label transfer downstream
#    NAO sera integrado via Harmony (dados ja normalizados pelo provedor)
# -------------------------------------------------------------------
cat("\n[3/3] Carregando GSE173278 primary (somente para label transfer)...\n")
 
bc_file   <- file.path(G173_DIR, "GSE173278_scRNAseq_filtered_cells_barcodes.tsv.gz")
gene_file <- file.path(G173_DIR, "GSE173278_scRNAseq_filtered_cells_genes.tsv.gz")
meta_file <- file.path(G173_DIR, "GSE173278_scRNAseq_filtered_cells_metadata.csv.gz")
mtx_file  <- file.path(G173_DIR, "GSE173278_scRNAseq_filtered_cells_norm_counts_matrix.mtx.gz")
 
if (!all(file.exists(c(bc_file, gene_file, meta_file, mtx_file)))) {
  stop("Arquivos do GSE173278 nao encontrados em ", G173_DIR)
}
 
barcodes <- fread(cmd = paste("zcat", shQuote(bc_file)), header = FALSE)$V1
genes    <- fread(cmd = paste("zcat", shQuote(gene_file)), header = FALSE)
meta     <- fread(cmd = paste("zcat", shQuote(meta_file)), data.table = FALSE)
mat      <- as(readMM(mtx_file), "CsparseMatrix")
storage.mode(mat@x) <- "double"
rownames(mat) <- make.unique(as.character(genes[[1]]))
colnames(mat) <- barcodes
 
barcode_col <- find_col(meta, c("barcode","barcodes","cell","cell_id","cell_name","cell_barcode"))
patient_col <- find_col(meta, c("patient","patient_id","sample","sample_id","case","case_id"))
if (is.null(barcode_col)) stop("Coluna de barcode nao detectada.")
if (is.null(patient_col)) stop("Coluna de patient nao detectada.")
 
meta <- meta[meta[[barcode_col]] %in% barcodes, , drop = FALSE]
rownames(meta) <- meta[[barcode_col]]
meta <- meta[barcodes, , drop = FALSE]
 
cell_name <- as.character(meta[[barcode_col]])
patient   <- as.character(meta[[patient_col]])
 
x <- tolower(cell_name)
is_gs        <- grepl("gs|gliomasphere|sphere|line", x)
is_pde       <- grepl("pde|explant", x)
is_recurrent <- grepl("recurrent|recurrence|\\brec\\b|\\-r\\b|_r\\b", x) |
                grepl("recurrent|recurrence|\\brec\\b|\\-r\\b|_r\\b", tolower(patient))
is_tissue    <- !is_pde & !is_gs
keep         <- is_tissue & !is_recurrent
 
cat("  Resumo GSE173278:\n")
print(c(tissue=sum(is_tissue), pde=sum(is_pde), gs=sum(is_gs),
        recurrent=sum(is_recurrent), keep_primary_tissue=sum(keep)))
if (sum(keep) == 0) stop("Nenhuma celula primary tissue no GSE173278.")
 
mat_sub  <- mat[, keep, drop = FALSE]
meta_sub <- meta[keep, , drop = FALSE]
rm(mat); gc()
 
# Downsampling para 30k antes de qualquer operacao Seurat
MAX_CELLS <- 30000
n_total   <- ncol(mat_sub)
cat(sprintf("  Celulas apos filtro: %d\n", n_total))
if (n_total > MAX_CELLS) {
  set.seed(42)
  idx <- sample(n_total, MAX_CELLS)
  mat_sub  <- mat_sub[, idx, drop = FALSE]
  meta_sub <- meta_sub[idx, , drop = FALSE]
  cat(sprintf("  Downsampling: %d -> %d\n", n_total, MAX_CELLS))
}
 
patient_vec <- as.character(meta_sub[[patient_col]])
patient_vec[is.na(patient_vec) | patient_vec == ""] <- "G173278"
orig_vec  <- gsub("[[:space:]/]+", "_", patient_vec)
new_cells <- paste0("G173_", orig_vec, "__", colnames(mat_sub))
colnames(mat_sub) <- new_cells
rownames(meta_sub) <- new_cells
 
# Filtra genes/celulas antes de criar o objeto
gene_ncells <- Matrix::rowSums(mat_sub > 0)
cell_ngenes <- Matrix::colSums(mat_sub > 0)
genes_pass  <- names(gene_ncells)[gene_ncells >= 3]
cells_pass  <- names(cell_ngenes)[cell_ngenes >= 200]
mat_norm    <- mat_sub[genes_pass, cells_pass, drop = FALSE]
meta_sub_f  <- meta_sub[cells_pass, , drop = FALSE]
rm(mat_sub); gc()
 
# Cria objeto Seurat com a matriz normalizada diretamente no slot data
# O slot counts fica com uma matriz binaria (placeholder)
mat_binary <- mat_norm
mat_binary@x <- rep(1.0, length(mat_binary@x))
assay_obj <- CreateAssayObject(counts = mat_binary, min.cells = 0, min.features = 0)
assay_obj <- SetAssayData(assay_obj, layer = "data", new.data = mat_norm)
rm(mat_binary, mat_norm); gc()
 
g173_so <- CreateSeuratObject(
  counts    = assay_obj,
  meta.data = meta_sub_f,
  project   = "GSE173278_primary"
)
rm(assay_obj, meta_sub_f); gc()
 
cells_keep <- colnames(g173_so)
g173_so$orig.ident         <- orig_vec[match(cells_keep, new_cells)]
g173_so$Patient            <- patient_vec[match(cells_keep, new_cells)]
g173_so$dataset            <- "GSE173278_primary"
g173_so$Type               <- "Primary GBM"
g173_so$ModelSystem        <- "Tissue"
g173_so$TimePoint          <- "Primary"
g173_so$idh_codel_subtype  <- "IDHwt"
g173_so$IDH_status         <- "IDH-wt"
g173_so$already_normalized <- TRUE
 
cat(sprintf("  GSE173278: %d celulas, %d genes\n", ncol(g173_so), nrow(g173_so)))
saveRDS(g173_so, file.path(RESULTS, "gse173278_primary_raw.rds"))
cat(sprintf("  Salvo: results_v3/gse173278_primary_raw.rds\n"))
 
cat("\n[DONE] Raw V3 concluido.\n")
REOF
