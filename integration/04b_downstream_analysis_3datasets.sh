#!/bin/bash
# =============================================================================
# SCRIPT: 04_downstream_v3.sh
# PROPOSITO:
#   - analises downstream do objeto integrado GSE182109 + Synapse
#   - encontra marcadores por cluster
#   - anotacao automatica por modulo de genes
#   - separa IDH-mut vs IDH-wt
#   - label transfer do GSE173278 para o objeto integrado
#     (valida anotacoes de tipo celular usando os dados de referencia)
# =============================================================================
#SBATCH --job-name=gbm_downstream_v3
#SBATCH --output=/home/renanomete/projetos/matdata/met_mat_data/singlecell_novo/logs_v3/04_downstream_v3_%j.out
#SBATCH --error=/home/renanomete/projetos/matdata/met_mat_data/singlecell_novo/logs_v3/04_downstream_v3_%j.err
#SBATCH --chdir=/home/renanomete/projetos/matdata/met_mat_data/singlecell_novo
#SBATCH -t 24:00:00
#SBATCH -c 4
#SBATCH --mem=128G

set -euo pipefail

BASE=/home/renanomete/projetos/matdata/met_mat_data/singlecell_novo
RESULTS=${BASE}/results_v3
LOGS=${BASE}/logs_v3

cd "${BASE}"
mkdir -p "${RESULTS}" "${LOGS}"

export INPUT="${RESULTS}/integrated_harmony_v3.rds \
${RESULTS}/gse173278_primary_qc.rds"
export OUTPUT="${RESULTS}/annotated_integrated_v3.rds \
${RESULTS}/markers_all_clusters_v3.csv \
${RESULTS}/umap_clusters_v3.pdf \
${RESULTS}/umap_dataset_v3.pdf \
${RESULTS}/umap_patient_v3.pdf \
${RESULTS}/umap_idh_v3.pdf \
${RESULTS}/idh_composition_by_cluster_v3.pdf \
${RESULTS}/cluster_annotation_auto_v3.csv \
${RESULTS}/label_transfer_g173_v3.csv \
${RESULTS}/umap_label_transfer_v3.pdf \
${RESULTS}/idh_mut_integrated_v3.rds \
${RESULTS}/idh_wt_integrated_v3.rds"

module load miniconda/24.4.0-libmamba
source activate gbm_scrnaseq

echo "[$(date '+%F %T')] Job iniciado em $(hostname)"
echo "[$(date '+%F %T')] Diretorio de execucao: $(pwd)"

job-nanny Rscript - <<'REOF'
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(data.table)
  library(patchwork)
})

BASE    <- "/home/renanomete/projetos/matdata/met_mat_data/singlecell_novo"
RESULTS <- file.path(BASE, "results_v3")

join_rna_layers_safely <- function(so) {
  if (!"RNA" %in% names(so@assays)) return(so)
  DefaultAssay(so) <- "RNA"
  so_try <- tryCatch(JoinLayers(so, assay="RNA"), error=function(e) NULL)
  if (!is.null(so_try)) return(so_try)
  so
}

safe_write_empty_markers <- function(path) {
  fwrite(data.frame(cluster=character(), gene=character(), p_val=numeric(),
                    avg_log2FC=numeric(), pct.1=numeric(), pct.2=numeric(),
                    p_val_adj=numeric()), path)
}

# ------------------------------------------------------------------
# [1/9] Carrega objeto integrado
# ------------------------------------------------------------------
cat("[1/9] Carregando objeto integrado...\n")
so <- readRDS(file.path(RESULTS, "integrated_harmony_v3.rds"))
if (!"RNA" %in% names(so@assays)) stop("Assay RNA nao encontrado.")
DefaultAssay(so) <- "RNA"
cat(sprintf("  %d celulas, %d genes\n", ncol(so), nrow(so)))

# ------------------------------------------------------------------
# [2/9] Padroniza metadata
# ------------------------------------------------------------------
cat("[2/9] Padronizando metadata...\n")
for (col in c("dataset","Patient","Type","idh_codel_subtype")) {
  if (!col %in% colnames(so@meta.data))
    so@meta.data[[col]] <- NA_character_
}

# ------------------------------------------------------------------
# [3/9] IDH status
# ------------------------------------------------------------------
cat("[3/9] Definindo IDH_status...\n")
synapse_idh_map <- c(
  "SM001"="IDH-mut","SM002"="IDH-mut","SM004"="IDH-mut","SM006"="IDH-wt",
  "SM008"="IDH-mut","SM011"="IDH-wt","SM012"="IDH-wt","SM015"="IDH-mut",
  "SM017"="IDH-wt","SM018"="IDH-wt","SM019"="IDH-mut"
)

so$IDH_status <- dplyr::case_when(
  so$dataset == "Synapse" & so$Patient %in% names(synapse_idh_map) ~
    unname(synapse_idh_map[so$Patient]),
  so$dataset == "Synapse" & grepl("^IDHmut", so$idh_codel_subtype, ignore.case=TRUE) ~ "IDH-mut",
  so$dataset == "Synapse" & grepl("^IDHwt",  so$idh_codel_subtype, ignore.case=TRUE) ~ "IDH-wt",
  grepl("^LGG$", so$Type, ignore.case=TRUE) ~ "IDH-mut",
  grepl("GBM",   so$Type, ignore.case=TRUE) ~ "IDH-wt",
  TRUE ~ "Unknown"
)
so$IDH_status <- factor(so$IDH_status, levels=c("IDH-mut","IDH-wt","Unknown"))
cat("IDH_status:\n"); print(table(so$IDH_status, useNA="ifany"))
cat("Dataset x IDH_status:\n"); print(table(so$dataset, so$IDH_status, useNA="ifany"))

# ------------------------------------------------------------------
# [4/9] UMAPs
# ------------------------------------------------------------------
cat("[4/9] Gerando UMAPs...\n")

cluster_col <- if ("RNA_snn_res.0.6" %in% colnames(so@meta.data)) "RNA_snn_res.0.6" else {
  cand <- grep("^RNA_snn_res\\.", colnames(so@meta.data), value=TRUE)
  if (length(cand) == 0) stop("Coluna de cluster nao encontrada.")
  cand[1]
}
if (!"umap" %in% names(so@reductions)) stop("UMAP nao encontrado.")
Idents(so) <- cluster_col

pdf(file.path(RESULTS, "umap_clusters_v3.pdf"), width=10, height=8)
print(DimPlot(so, reduction="umap", group.by=cluster_col, label=TRUE, repel=TRUE, raster=TRUE))
dev.off()

pdf(file.path(RESULTS, "umap_dataset_v3.pdf"), width=10, height=8)
print(DimPlot(so, reduction="umap", group.by="dataset", raster=TRUE))
dev.off()

pdf(file.path(RESULTS, "umap_patient_v3.pdf"), width=10, height=8)
print(DimPlot(so, reduction="umap", group.by="Patient", raster=TRUE))
dev.off()

pdf(file.path(RESULTS, "umap_idh_v3.pdf"), width=10, height=8)
print(DimPlot(so, reduction="umap", group.by="IDH_status", raster=TRUE))
dev.off()

# ------------------------------------------------------------------
# [5/9] Marcadores por cluster
# ------------------------------------------------------------------
cat("[5/9] Encontrando marcadores...\n")
DefaultAssay(so) <- "RNA"
Idents(so) <- cluster_col
so <- join_rna_layers_safely(so)

markers <- tryCatch(
  FindAllMarkers(so, assay="RNA", only.pos=TRUE, test.use="wilcox",
                 min.pct=0.25, logfc.threshold=0.25, max.cells.per.ident=3000, verbose=FALSE),
  error = function(e) { cat(sprintf("  [AVISO] FindAllMarkers falhou: %s\n", conditionMessage(e))); data.frame() }
)

if (!is.data.frame(markers) || nrow(markers) == 0) {
  safe_write_empty_markers(file.path(RESULTS, "markers_all_clusters_v3.csv"))
} else {
  if (!"p_val_adj" %in% colnames(markers)) markers$p_val_adj <- NA_real_
  if (!"avg_log2FC" %in% colnames(markers) && "avg_logFC" %in% colnames(markers))
    markers$avg_log2FC <- markers$avg_logFC
  if (!"avg_log2FC" %in% colnames(markers)) markers$avg_log2FC <- NA_real_
  markers_sig <- markers
  if (any(!is.na(markers_sig$p_val_adj)))
    markers_sig <- markers_sig %>% dplyr::filter(!is.na(p_val_adj) & p_val_adj < 0.05)
  if (nrow(markers_sig) > 0)
    markers_sig <- markers_sig %>% dplyr::arrange(cluster, dplyr::desc(avg_log2FC))
  fwrite(markers_sig, file.path(RESULTS, "markers_all_clusters_v3.csv"))
  cat(sprintf("  Marcadores significativos: %d\n", nrow(markers_sig)))
}

# ------------------------------------------------------------------
# [6/9] Anotacao automatica por modulo de genes
# ------------------------------------------------------------------
cat("[6/9] Anotacao automatica...\n")

marker_panels <- list(
  Glioma         = c("SOX2","OLIG1","OLIG2","GFAP","S100B","CHI3L1","PDGFRA"),
  Myeloid        = c("PTPRC","ITGAM","CD68","P2RY12","TMEM119","CD14","FCER1G"),
  Tcells         = c("CD3D","CD3E","CD4","CD8A","IL7R"),
  Bcells         = c("CD79A","MS4A1","CD19"),
  Endothelial    = c("PECAM1","VWF","KDR"),
  Pericytes      = c("PDGFRB","RGS5","ACTA2"),
  Oligodendrocytes = c("MBP","MOG","PLP1")
)
marker_panels <- lapply(marker_panels, function(g) intersect(g, rownames(so)))

for (nm in names(marker_panels)) {
  genes <- marker_panels[[nm]]
  if (length(genes) >= 2) {
    so <- AddModuleScore(so, features=list(genes), name=paste0(nm,"_score"),
                         assay="RNA", search=FALSE)
  }
}

score_cols  <- grep("_score1$", colnames(so@meta.data), value=TRUE)
score_names <- sub("_score1$", "", score_cols)
if (length(score_cols) > 0) {
  score_mat <- so@meta.data[, score_cols, drop=FALSE]
  so$predicted_celltype <- score_names[apply(score_mat, 1, which.max)]
} else {
  so$predicted_celltype <- "Unassigned"
}

cluster_annot <- so@meta.data %>%
  dplyr::group_by(.data[[cluster_col]], predicted_celltype) %>%
  dplyr::summarise(n=dplyr::n(), .groups="drop") %>%
  dplyr::group_by(.data[[cluster_col]]) %>%
  dplyr::slice_max(n, n=1, with_ties=FALSE) %>%
  dplyr::rename(cluster=!!cluster_col, annotation=predicted_celltype)
fwrite(cluster_annot, file.path(RESULTS, "cluster_annotation_auto_v3.csv"))

# ------------------------------------------------------------------
# [7/9] IDH composition + subset
# ------------------------------------------------------------------
cat("[7/9] Composicao IDH por cluster e subset...\n")

idh_comp <- so@meta.data %>%
  dplyr::filter(IDH_status %in% c("IDH-mut","IDH-wt")) %>%
  dplyr::group_by(.data[[cluster_col]], IDH_status) %>%
  dplyr::summarise(n=dplyr::n(), .groups="drop") %>%
  dplyr::group_by(IDH_status) %>%
  dplyr::mutate(pct=100*n/sum(n)) %>%
  dplyr::ungroup()

pdf(file.path(RESULTS, "idh_composition_by_cluster_v3.pdf"), width=12, height=6)
print(ggplot(idh_comp, aes(x=.data[[cluster_col]], y=pct, fill=IDH_status)) +
  geom_bar(stat="identity", position="dodge") +
  theme_bw(base_size=12) +
  labs(x="Cluster", y="% de celulas", fill="IDH status"))
dev.off()

saveRDS(subset(so, subset=IDH_status=="IDH-mut"), file.path(RESULTS, "idh_mut_integrated_v3.rds"))
saveRDS(subset(so, subset=IDH_status=="IDH-wt"),  file.path(RESULTS, "idh_wt_integrated_v3.rds"))

# ------------------------------------------------------------------
# [8/9] Label transfer do GSE173278 -> objeto integrado
# ------------------------------------------------------------------
cat("[8/9] Label transfer do GSE173278...\n")

g173_path <- file.path(RESULTS, "gse173278_primary_qc.rds")
if (!file.exists(g173_path)) {
  cat("  [AVISO] gse173278_primary_qc.rds nao encontrado. Pulando label transfer.\n")
} else {
  g173 <- readRDS(g173_path)
  cat(sprintf("  GSE173278: %d celulas, %d genes\n", ncol(g173), nrow(g173)))

  # Padroniza nomes de genes no GSE173278 para bater com o objeto integrado
  standardize_gene_names <- function(x) {
    x <- trimws(as.character(x))
    x <- sub("\\.[0-9]+$", "", x)
    x <- gsub("_", "-", x, fixed=TRUE)
    toupper(x)
  }

  DefaultAssay(g173) <- "RNA"

  # Usa o slot data (normalizacao original) para FindVariableFeatures e ScaleData
  g173 <- ScaleData(g173, features=VariableFeatures(g173),
                    do.scale=TRUE, do.center=TRUE, verbose=FALSE)
  g173 <- RunPCA(g173, features=VariableFeatures(g173), npcs=30, verbose=FALSE)

  # Genes em comum entre referencia e query
  genes_common <- intersect(rownames(so), rownames(g173))
  cat(sprintf("  Genes em comum para label transfer: %d\n", length(genes_common)))

  if (length(genes_common) < 500) {
    cat("  [AVISO] Poucos genes em comum. Pulando label transfer.\n")
  } else {
    # Encontra anchors usando o objeto integrado como referencia
    # e o GSE173278 como query
    transfer_anchors <- tryCatch(
      FindTransferAnchors(
        reference   = so,
        query       = g173,
        dims        = 1:30,
        reference.reduction = "pca",
        verbose     = FALSE
      ),
      error = function(e) {
        cat(sprintf("  [AVISO] FindTransferAnchors falhou: %s\n", conditionMessage(e)))
        NULL
      }
    )

    if (!is.null(transfer_anchors)) {
      # Transfere as anotacoes de tipo celular e cluster
      predictions <- TransferData(
        anchorset   = transfer_anchors,
        refdata     = list(
          celltype = so$predicted_celltype,
          cluster  = as.character(so@meta.data[[cluster_col]])
        ),
        dims        = 1:30,
        verbose     = FALSE
      )

      g173$predicted_celltype_transfer <- predictions$predicted.celltype
      g173$predicted_cluster_transfer  <- predictions$predicted.cluster
      g173$prediction_score            <- predictions$prediction.score.max

      # Salva tabela de resultados
      lt_df <- data.frame(
        cell                     = colnames(g173),
        Patient                  = g173$Patient,
        predicted_celltype       = g173$predicted_celltype_transfer,
        predicted_cluster        = g173$predicted_cluster_transfer,
        prediction_score         = g173$prediction_score,
        stringsAsFactors         = FALSE
      )
      fwrite(lt_df, file.path(RESULTS, "label_transfer_g173_v3.csv"))
      cat("  Label transfer concluido.\n")
      cat("  Distribuicao de tipos celulares transferidos:\n")
      print(table(g173$predicted_celltype_transfer))

      # UMAP do label transfer
      g173 <- tryCatch({
        MapQuery(
          anchorset       = transfer_anchors,
          query           = g173,
          reference       = so,
          refdata         = list(celltype = "predicted_celltype"),
          reference.reduction = "pca",
          reduction.model = "umap",
          verbose         = FALSE
        )
      }, error = function(e) {
        cat(sprintf("  [AVISO] MapQuery falhou: %s\n", conditionMessage(e)))
        g173
      })

      if ("ref.umap" %in% names(g173@reductions)) {
        pdf(file.path(RESULTS, "umap_label_transfer_v3.pdf"), width=14, height=6)
        p1 <- DimPlot(so, reduction="umap", group.by="predicted_celltype",
                      label=TRUE, repel=TRUE, raster=TRUE) + ggtitle("Referencia (GSE182109+Synapse)")
        p2 <- DimPlot(g173, reduction="ref.umap", group.by="predicted_celltype_transfer",
                      label=TRUE, repel=TRUE, raster=TRUE) + ggtitle("Query projetado (GSE173278)")
        print(p1 + p2)
        dev.off()
      } else {
        # Fallback: UMAP proprio do GSE173278 colorido por label transfer
        g173 <- RunUMAP(g173, dims=1:30, verbose=FALSE)
        pdf(file.path(RESULTS, "umap_label_transfer_v3.pdf"), width=10, height=8)
        print(DimPlot(g173, reduction="umap", group.by="predicted_celltype_transfer",
                      label=TRUE, repel=TRUE, raster=TRUE) +
              ggtitle("GSE173278 - tipos celulares transferidos"))
        dev.off()
      }

    } else {
      cat("  Label transfer nao realizado (FindTransferAnchors falhou).\n")
      fwrite(data.frame(), file.path(RESULTS, "label_transfer_g173_v3.csv"))
    }
  }
  rm(g173); gc()
}

# ------------------------------------------------------------------
# [9/9] Salva objeto final
# ------------------------------------------------------------------
cat("[9/9] Salvando objeto final anotado...\n")
saveRDS(so, file.path(RESULTS, "annotated_integrated_v3.rds"))
cat("\n[DONE] Downstream V3 concluido.\n")
REOF
