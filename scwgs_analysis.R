#!/usr/bin/env Rscript
# =============================================================================
# scwgs_analysis.R - Análise de scWGS pós-MEDICC2 (GridUNESP)
# CORREÇÃO: filtrar diploid_normal de pw_mat, cn_profiles e bl_df
#           para evitar grupo SCNA espúrio no heatmap
# =============================================================================

suppressPackageStartupMessages({
  library(ape)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(pheatmap)
  library(RColorBrewer)
  library(viridisLite)
  library(scales)
  library(umap)
})

args <- commandArgs(trailingOnly = TRUE)
PATIENT <- if (length(args) >= 1) args[1] else "JK136"

BASE_DIR <- "/home/renanomete/projetos/matdata/new_evo_anal/scWGS/GSE173279/results/evo_clonal"
MEDICC_DIR <- file.path(BASE_DIR, "results", "medicc2", PATIENT)
OUT_DIR    <- file.path(BASE_DIR, "results", "scwgs_analysis", PATIENT)
LOG_DIR    <- file.path(BASE_DIR, "logs")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

LOG_FILE <- file.path(LOG_DIR, paste0("scwgs_", PATIENT, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
con_log <- file(LOG_FILE, open = "wt")

log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", paste0(...))
  message(msg)
  writeLines(msg, con = con_log)
}

log_msg("========================================")
log_msg(" scWGS Analysis - Patient: ", PATIENT)
log_msg(" Base dir: ", BASE_DIR)
log_msg(" Output dir: ", OUT_DIR)
log_msg("========================================")

# --- Localizar arquivos ---
find_file <- function(pattern) {
  files <- list.files(MEDICC_DIR, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) return(NA)
  return(files[1])
}

TREE_FILE    <- find_file("_final_tree\\.new$")
BL_FILE      <- find_file("_branch_lengths\\.tsv$")
CN_FILE      <- find_file("_final_cn_profiles\\.tsv$")
PW_FILE      <- find_file("_pairwise_distances\\.tsv$")
SUMMARY_FILE <- find_file("_summary\\.tsv$")

log_msg("Arquivos encontrados:")
log_msg("  Tree: ", ifelse(is.na(TREE_FILE), "FALTANDO", basename(TREE_FILE)))
log_msg("  Branch lengths: ", ifelse(is.na(BL_FILE), "FALTANDO", basename(BL_FILE)))
log_msg("  CN profiles: ", ifelse(is.na(CN_FILE), "FALTANDO", basename(CN_FILE)))
log_msg("  Pairwise: ", ifelse(is.na(PW_FILE), "FALTANDO", basename(PW_FILE)))

if (any(is.na(c(TREE_FILE, BL_FILE, CN_FILE, PW_FILE)))) {
  log_msg("[ERRO] Arquivos essenciais faltando. MEDICC2 não rodou completamente.")
  close(con_log)
  quit(save = "no", status = 1)
}

# --- Carregar dados ---
log_msg("Carregando dados...")
if (!is.na(SUMMARY_FILE)) {
  sum_df <- read.table(SUMMARY_FILE, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  if ("nsamples" %in% colnames(sum_df)) {
    log_msg("  Células: ", sum_df$nsamples[1])
    log_msg("  Tree length: ", sum_df$tree_length[1])
  }
}

# --- Branch lengths: excluir nós internos e diploid_normal ---
bl_df <- read.table(BL_FILE, header = FALSE, sep = "\t",
                    col.names = c("node_id", "branch_length"), stringsAsFactors = FALSE)
bl_cells <- bl_df %>%
  filter(!grepl("^internal|diploid_normal|diploid", node_id, ignore.case = TRUE))
log_msg("  Células com branch length (após filtro): ", nrow(bl_cells))

# --- CN profiles: excluir diploid_normal ---
cn_profiles <- read.table(CN_FILE, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cn_profiles <- cn_profiles %>%
  filter(!grepl("diploid_normal|diploid", sample_id, ignore.case = TRUE))
cn_col <- if ("cn_a" %in% colnames(cn_profiles)) "cn_a" else "copy_number"
log_msg("  CN profiles (após filtro): ", nrow(cn_profiles), " linhas, coluna CN = ", cn_col)

# --- Pairwise distances: excluir diploid_normal ---
pw_mat_raw <- read.table(PW_FILE, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
pw_mat_raw <- as.matrix(pw_mat_raw)
diploid_rows <- grepl("diploid_normal|diploid", rownames(pw_mat_raw), ignore.case = TRUE)
diploid_cols <- grepl("diploid_normal|diploid", colnames(pw_mat_raw), ignore.case = TRUE)
pw_mat <- pw_mat_raw[!diploid_rows, !diploid_cols]
pw_dist <- as.dist(pw_mat)
log_msg("  Pairwise distances (após filtro): ", nrow(pw_mat), "x", ncol(pw_mat))

tree <- read.tree(TREE_FILE)
log_msg("  Tree: ", length(tree$tip.label), " folhas")

# --- Genes de interesse ---
GENE_LOCI <- data.frame(
  gene   = c("TYMS","PCLAF","BIRC5","PBK","TPX2","EZH2","MYBL2","NEK2",
             "PLP1","MBP","COL1A2","COL1A1","CTHRC1","PCOLCE"),
  chrom  = c("chr18","chr1","chr17","chr8","chr20","chr7","chr20","chr1",
             "chrX","chr18","chr7","chr17","chr8","chr7"),
  pos_mb = c(0.7,35.1,76.0,19.0,32.0,148.2,42.0,212.4,
             103.3,76.4,94.0,48.2,114.3,94.8),
  axis   = c(rep("Proliferativo",8), rep("Diferenciado",6)),
  stringsAsFactors = FALSE
)
BIN_SIZE <- 1040000
GENE_LOCI$bin_start <- floor(GENE_LOCI$pos_mb * 1e6 / BIN_SIZE) * BIN_SIZE + 1
GENE_LOCI$bin_end   <- GENE_LOCI$bin_start + BIN_SIZE - 1

# --- Subclonagem ---
log_msg("Definindo subclones...")
hc <- hclust(pw_dist, method = "ward.D2")
K_VALS <- c(3,5,8,10,15)
sub_assign <- data.frame(cell_id = hc$labels)
for (k in K_VALS) {
  sub_assign[[paste0("k", k)]] <- cutree(hc, k = k)
}
K_DEFAULT <- 5
sub_assign$subclone <- sub_assign[[paste0("k", K_DEFAULT)]]

write.table(sub_assign, file.path(OUT_DIR, paste0(PATIENT, "_subclones.tsv")),
            sep = "\t", row.names = FALSE, quote = FALSE)
log_msg("  Subclones (k=5): ", length(unique(sub_assign$subclone)), " grupos")

calc_div <- function(x) {
  p <- table(x)/length(x)
  shannon <- -sum(p * log(p))
  simpson <- 1 - sum(p^2)
  c(shannon = shannon, simpson = simpson, n_clones = length(p), max_prop = max(p))
}
div_stats <- calc_div(sub_assign$subclone)
log_msg("  Diversidade: Shannon = ", round(div_stats["shannon"], 3),
        ", Simpson = ", round(div_stats["simpson"], 3))

# --- Extrair CN dos genes ---
log_msg("Extraindo CN dos genes...")
gene_cn_list <- list()
for (i in 1:nrow(GENE_LOCI)) {
  g <- GENE_LOCI$gene[i]
  chr <- GENE_LOCI$chrom[i]
  bs <- GENE_LOCI$bin_start[i]
  be <- GENE_LOCI$bin_end[i]
  bins <- cn_profiles %>%
    filter(chrom == chr, start <= be + BIN_SIZE, end >= bs - BIN_SIZE)
  if (nrow(bins) > 0) {
    bins$dist <- abs((bins$start + bins$end)/2 - (bs+be)/2)
    best <- bins %>% arrange(dist) %>% slice(1)
    gene_cn <- cn_profiles %>%
      filter(chrom == best$chrom, start == best$start, end == best$end) %>%
      select(sample_id, !!cn_col) %>%
      rename(!!g := !!cn_col)
    gene_cn_list[[g]] <- gene_cn
  }
}
if (length(gene_cn_list) > 0) {
  gene_cn_cells <- Reduce(function(x,y) full_join(x,y, by="sample_id"), gene_cn_list)
  gene_cn_cells <- left_join(gene_cn_cells, sub_assign, by=c("sample_id"="cell_id"))
  # Remover células que não foram atribuídas a nenhum subclone (segurança extra)
  gene_cn_cells <- gene_cn_cells %>% filter(!is.na(subclone))
  write.table(gene_cn_cells, file.path(OUT_DIR, paste0(PATIENT, "_gene_cn_per_cell.tsv")),
              sep="\t", row.names=FALSE, quote=FALSE)
  gene_cn_sub <- gene_cn_cells %>%
    group_by(subclone) %>%
    summarise(across(any_of(GENE_LOCI$gene), \(x) mean(x, na.rm = TRUE)))
  write.table(gene_cn_sub, file.path(OUT_DIR, paste0(PATIENT, "_gene_cn_per_subclone.tsv")),
              sep="\t", row.names=FALSE, quote=FALSE)
  log_msg("  CN genes extraído para ", ncol(gene_cn_sub)-1, " genes")
} else {
  log_msg("  Nenhum gene encontrado nos bins.")
}

# --- Perfil CNA médio por subclone ---
log_msg("Calculando perfil CNA médio por subclone...")
cn_sub <- cn_profiles %>%
  left_join(sub_assign, by=c("sample_id"="cell_id")) %>%
  filter(!is.na(subclone)) %>%
  group_by(subclone, chrom, start, end) %>%
  summarise(cn_mean = mean(!!sym(cn_col), na.rm=TRUE), .groups="drop")

cn_sub <- cn_sub %>%
  mutate(bin_id = paste0(chrom, ":", start, "-", end)) %>%
  select(subclone, bin_id, cn_mean) %>%
  pivot_wider(names_from=bin_id, values_from=cn_mean, values_fill=2)
cn_mat <- as.matrix(cn_sub[,-1])
rownames(cn_mat) <- paste0("SC", cn_sub$subclone)

write.table(cn_mat, file.path(OUT_DIR, paste0(PATIENT, "_subclone_cn_profiles.tsv")),
            sep="\t", quote=FALSE, row.names=TRUE)

# --- UMAP ---
if (nrow(pw_mat) > 50) {
  log_msg("Executando UMAP...")
  set.seed(42)
  umap_cfg <- umap.defaults
  umap_cfg$n_neighbors <- min(30, nrow(pw_mat)-1)
  umap_cfg$min_dist <- 0.1
  umap_res <- umap(pw_mat, config = umap_cfg)
  umap_df <- data.frame(UMAP1 = umap_res$layout[,1],
                        UMAP2 = umap_res$layout[,2],
                        cell_id = rownames(pw_mat)) %>%
    left_join(sub_assign, by="cell_id")

  p <- ggplot(umap_df, aes(UMAP1, UMAP2, color=factor(subclone))) +
    geom_point(size=0.8, alpha=0.7) +
    scale_color_brewer(palette="Set1") +
    labs(title=paste0(PATIENT, " - UMAP CNA (k=5 subclones)")) +
    theme_minimal()
  ggsave(file.path(OUT_DIR, paste0(PATIENT, "_umap_cna.pdf")), p, width=8, height=6)
  write.table(umap_df, file.path(OUT_DIR, paste0(PATIENT, "_umap_coords.tsv")),
              sep="\t", row.names=FALSE, quote=FALSE)
  log_msg("  UMAP gerado com sucesso.")
}

# --- Árvore de subclones ---
if (exists("cn_mat") && nrow(cn_mat) >= 3) {
  log_msg("Construindo árvore de subclones...")
  sub_dist <- dist(cn_mat)
  hc_sub <- hclust(sub_dist, method="ward.D2")
  tree_sub <- as.phylo(hc_sub)
  pdf(file.path(OUT_DIR, paste0(PATIENT, "_subclone_tree.pdf")), width=10, height=8)
  plot(tree_sub, main=paste0(PATIENT, " - Árvore de subclones"))
  add.scale.bar()
  dev.off()
  write.tree(tree_sub, file=file.path(OUT_DIR, paste0(PATIENT, "_subclone_tree.new")))
  log_msg("  Árvore de subclones salva.")
}

# --- Heatmap genes ---
if (exists("gene_cn_sub") && nrow(gene_cn_sub) >= 2) {
  mat <- as.matrix(gene_cn_sub[,-1])
  rownames(mat) <- paste0("SC", gene_cn_sub$subclone)
  pdf(file.path(OUT_DIR, paste0(PATIENT, "_gene_cn_heatmap.pdf")), width=12, height=8)
  pheatmap(mat, main=paste0(PATIENT, " - CN genes por subclone"),
           display_numbers=TRUE, number_format="%.1f",
           color=colorRampPalette(c("blue","white","red"))(50))
  dev.off()
  log_msg("  Heatmap de genes salvo.")
}

log_msg("========================================")
log_msg(" Análise concluída!")
log_msg(" Resultados em: ", OUT_DIR)
log_msg("========================================")
close(con_log)
