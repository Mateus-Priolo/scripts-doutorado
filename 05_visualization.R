#!/usr/bin/env Rscript
# =============================================================================
# 05_visualization.R — Todas as visualizacoes do pipeline MEDICC2
# Arvores filogeneticas: ape::plot.phylo (sem dependencia de ggtree)
# Demais figuras: ggplot2
#
# Entrada : data/processed/tree_analysis.rds
# Saida   : results/figures/
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
  library(scales)
  library(irlba)
  library(uwot)
  library(harmony)
  library(future)
  library(Matrix)
})

# Carregar ggtree se disponivel — usa ape como fallback transparente
HAS_GGTREE <- requireNamespace("ggtree", quietly = TRUE) &&
              requireNamespace("treeio", quietly = TRUE)
if (HAS_GGTREE) {
  suppressPackageStartupMessages({
    library(ggtree)
    library(treeio)
  })
  cat("  [ggtree disponivel — usando para arvores]\n")
} else {
  cat("  [ggtree indisponivel — usando ape::plot.phylo como fallback]\n")
}

WORKDIR <- "/home/renanomete/projetos/matdata/evo_clonal"
setwd(WORKDIR)

n_workers <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "4"))
plan(multicore, workers = n_workers)
options(future.globals.maxSize = 8 * 1024^3)

source("scripts/functions.R")
ensure_dir("results/figures")

cat("=================================================================\n")
cat(" ETAPA 5 — Visualizacoes\n")
cat(sprintf(" Diretorio: %s\n", WORKDIR))
cat(sprintf(" Inicio: %s\n", format(Sys.time())))
cat("=================================================================\n")

# =============================================================================
# 1. Carregar dados
# =============================================================================
log_step(1, "Carregando dados...")

ta          <- readRDS("data/processed/tree_analysis.rds")
runs        <- ta$runs
tree_metrics<- ta$tree_metrics
metrics_all <- ta$metrics_all
qc_obj      <- ta$qc_obj
meta        <- qc_obj$meta
patient_col <- qc_obj$patient_col
cn_int      <- qc_obj$cn_int
coords      <- qc_obj$coords

patients <- sort(unique(metrics_all$patient))

# Paletas
patient_pal <- setNames(
  brewer.pal(max(3, length(patients)), "Set1")[seq_along(patients)],
  patients
)

clone_levels <- sort(unique(na.omit(
  metrics_all$clone_cluster[metrics_all$run == "global"])))
n_clones <- max(1, length(clone_levels))
clone_pal <- setNames(
  colorRampPalette(c("#1F77B4","#FF7F0E","#2CA02C","#D62728","#9467BD",
                     "#8C564B","#E377C2","#7F7F7F","#BCBD22","#17BECF",
                     "#AEC7E8","#FFBB78"))(n_clones),
  clone_levels
)

# =============================================================================
# 2. Arvores filogeneticas
# =============================================================================
log_step(2, "Gerando arvores filogeneticas...")

# Funcao de plot usando ape (funciona com e sem ggtree)
plot_tree_ape <- function(tree, tip_colors, title = "", cex = 0.3) {
  old_par <- par(mar = c(1, 1, 2, 1))
  ape::plot.phylo(
    tree,
    type        = "phylogram",
    tip.color   = tip_colors,
    cex         = cex,
    main        = title,
    edge.width  = 0.8,
    label.offset = max(tree$edge.length, na.rm=TRUE) * 0.01
  )
  ape::add.scale.bar(cex = 0.6, col = "black")
  par(old_par)
}

# Funcao de plot usando ggtree (se disponivel)
plot_tree_ggtree <- function(tree, tip_data, colour_by, pal, title = "") {
  p <- ggtree(tree, layout = "rectangular", linewidth = 0.3) %<+% tip_data
  p <- p +
    geom_tippoint(aes(colour = .data[[colour_by]], size = is_tumor), alpha = 0.8) +
    scale_colour_manual(values = pal, na.value = "#CCCCCC") +
    scale_size_manual(values = c("TRUE"=1.5,"FALSE"=3), guide="none") +
    geom_treescale(x=0, y=0, offset=1, fontsize=3) +
    labs(title=title, colour=colour_by) +
    theme_tree2() +
    theme(legend.position="right",
          plot.title=element_text(size=10, face="bold"))
  return(p)
}

pdf("results/figures/01_phylo_tree_global.pdf", width=14, height=10)

for (colour_var in c("patient", "clone")) {
  run <- runs[["global"]]
  if (is.null(run)) next
  tree <- run$tree

  tip_labels <- tree$tip.label
  is_tumor   <- tip_labels != "diploid_normal"

  if (colour_var == "patient") {
    patient_tips <- meta[[patient_col]][match(tip_labels, rownames(meta))]
    patient_tips[!is_tumor] <- "Reference"
    pal_use <- c(patient_pal, "Reference"="#888888")
    tip_colors <- pal_use[ifelse(is.na(patient_tips), "Reference", patient_tips)]
    title_str <- "Arvore filogenetica global — por paciente"
    legend_var <- "patient"
    tip_data <- data.frame(label=tip_labels, is_tumor=is_tumor,
                           patient=patient_tips, stringsAsFactors=FALSE)
  } else {
    gm <- metrics_all[metrics_all$run=="global",]
    clone_tips <- gm$clone_cluster[match(tip_labels, gm$cell)]
    clone_tips[!is_tumor] <- "Reference"
    pal_use <- c(clone_pal, "Reference"="#888888")
    tip_colors <- pal_use[ifelse(is.na(clone_tips), "Reference", clone_tips)]
    title_str <- "Arvore filogenetica global — por clone"
    legend_var <- "clone_cluster"
    tip_data <- data.frame(label=tip_labels, is_tumor=is_tumor,
                           clone_cluster=clone_tips, stringsAsFactors=FALSE)
  }

  if (HAS_GGTREE) {
    print(plot_tree_ggtree(tree, tip_data, legend_var, pal_use, title_str))
  } else {
    plot_tree_ape(tree, tip_colors, title_str)
    # Legenda manual
    legend("topright",
           legend = names(pal_use),
           col    = pal_use,
           pch    = 16, cex = 0.6, bty = "n")
  }
}
dev.off()
log_info("OK results/figures/01_phylo_tree_global.pdf")

# Arvores por paciente
pdf("results/figures/02_phylo_trees_by_patient.pdf", width=18, height=8)
patient_tree_list <- list()
for (pat in patients) {
  run <- runs[[pat]]
  if (is.null(run)) next
  tree <- run$tree
  tip_labels <- tree$tip.label
  is_tumor   <- tip_labels != "diploid_normal"
  patient_tips <- meta[[patient_col]][match(tip_labels, rownames(meta))]
  patient_tips[!is_tumor] <- "Reference"
  pal_use    <- c(patient_pal, "Reference"="#888888")
  tip_colors <- pal_use[ifelse(is.na(patient_tips), "Reference", patient_tips)]

  if (HAS_GGTREE) {
    tip_data <- data.frame(label=tip_labels, is_tumor=is_tumor,
                           patient=patient_tips, stringsAsFactors=FALSE)
    patient_tree_list[[pat]] <- plot_tree_ggtree(
      tree, tip_data, "patient", pal_use, sprintf("Filogenia — %s", pat))
  } else {
    plot_tree_ape(tree, tip_colors, sprintf("Filogenia — %s", pat), cex=0.4)
  }
}
if (HAS_GGTREE && length(patient_tree_list) > 0) {
  print(patchwork::wrap_plots(patient_tree_list, nrow=1))
}
dev.off()
log_info("OK results/figures/02_phylo_trees_by_patient.pdf")

# =============================================================================
# 3. UMAP sobre CN (Harmony por paciente)
# =============================================================================
log_step(3, "Calculando UMAP...")

cn_var <- select_variable_bins(as.matrix(cn_int), n_bins=2000)

bin_med <- apply(cn_var, 1, median)
bin_mad <- apply(cn_var, 1, mad); bin_mad[bin_mad==0] <- 1
cn_sc   <- sweep(cn_var, 1, bin_med, "-")
cn_sc   <- sweep(cn_sc,  1, bin_mad, "/")
cn_sc[cn_sc >  5] <-  5
cn_sc[cn_sc < -5] <- -5

n_pcs <- min(30L, nrow(cn_sc)-1L, ncol(cn_sc)-1L)
pca   <- irlba(t(cn_sc), nv=n_pcs)
pca_emb <- pca$u %*% diag(pca$d)
rownames(pca_emb) <- colnames(cn_sc)

patient_vec <- meta[[patient_col]][match(rownames(pca_emb), rownames(meta))]
harm_emb <- HarmonyMatrix(
  data_mat  = pca_emb,
  meta_data = data.frame(patient=patient_vec, row.names=rownames(pca_emb)),
  vars_use  = "patient",
  do_pca    = FALSE,
  verbose   = FALSE
)
rownames(harm_emb) <- rownames(pca_emb)

umap_coords <- umap(harm_emb, n_neighbors=30, min_dist=0.3,
                    n_threads=n_workers, seed=42, verbose=FALSE)
rownames(umap_coords) <- rownames(harm_emb)
colnames(umap_coords) <- c("UMAP_1","UMAP_2")

global_m <- metrics_all %>%
  filter(run=="global") %>%
  distinct(cell, .keep_all=TRUE) %>%
  select(cell, evo_dist_to_root, clone_cluster)

umap_df <- data.frame(
  umap_coords,
  cell    = rownames(umap_coords),
  patient = meta[[patient_col]][match(rownames(umap_coords), rownames(meta))],
  stringsAsFactors = FALSE
) %>% left_join(global_m, by="cell")

p_umap_pat <- ggplot(umap_df, aes(UMAP_1, UMAP_2, colour=patient)) +
  geom_point(size=.6, alpha=.7) +
  scale_colour_manual(values=patient_pal, name="Paciente") +
  labs(title="UMAP — Paciente (pos-Harmony)") +
  theme_bw(base_size=10) +
  guides(colour=guide_legend(override.aes=list(size=3)))

p_umap_dist <- ggplot(umap_df, aes(UMAP_1, UMAP_2, colour=evo_dist_to_root)) +
  geom_point(size=.6, alpha=.8) +
  scale_colour_viridis_c(option="plasma", name="Dist.\nevolutiva",
                          na.value="#CCCCCC") +
  labs(title="UMAP — Distancia evolutiva (MEDICC2)",
       subtitle="Maior distancia = mais alteracoes desde o clone ancestral") +
  theme_bw(base_size=10)

p_umap_clone <- ggplot(umap_df, aes(UMAP_1, UMAP_2, colour=clone_cluster)) +
  geom_point(size=.6, alpha=.8) +
  scale_colour_manual(values=clone_pal, na.value="#CCCCCC", name="Clone") +
  labs(title="UMAP — Clone cluster (MEDICC2 + hclust)") +
  theme_bw(base_size=10) +
  guides(colour=guide_legend(override.aes=list(size=3)))

pdf("results/figures/03_umap_panels.pdf", width=16, height=5)
print(p_umap_pat | p_umap_dist | p_umap_clone)
dev.off()
log_info("OK results/figures/03_umap_panels.pdf")

# =============================================================================
# 4. Heatmap CN por clone x cromossomo
# =============================================================================
log_step(4, "Heatmap CN por clone x cromossomo...")

gm_full <- metrics_all %>% filter(run=="global", !is.na(clone_cluster))
if (nrow(gm_full) > 0) {
  cn_cells <- intersect(colnames(cn_int), gm_full$cell)
  cn_sub   <- as.matrix(cn_int[, cn_cells, drop=FALSE])
  clone_vec<- gm_full$clone_cluster[match(cn_cells, gm_full$cell)]
  chrom_f  <- coords$chrom[match(rownames(cn_sub), rownames(coords))]
  valid    <- !is.na(chrom_f)
  cn_sub   <- cn_sub[valid,, drop=FALSE]
  chrom_f  <- chrom_f[valid]

  hm_long <- data.frame(chrom=chrom_f, t(cn_sub), check.names=FALSE) %>%
    pivot_longer(-chrom, names_to="cell", values_to="cn") %>%
    left_join(data.frame(cell=cn_cells, clone=clone_vec), by="cell") %>%
    filter(!is.na(clone)) %>%
    group_by(clone, chrom) %>%
    summarise(median_cn=median(cn, na.rm=TRUE), .groups="drop")

  chr_order <- c(paste0("chr",1:22),"chrX","chrY",as.character(1:22),"X","Y")
  hm_long$chrom <- factor(hm_long$chrom,
                            levels=chr_order[chr_order %in% hm_long$chrom])

  p_hm <- ggplot(hm_long, aes(chrom, clone, fill=median_cn)) +
    geom_tile() +
    scale_fill_gradient2(low="#2166AC", mid="white", high="#D73027",
                          midpoint=2, limits=c(0,6), oob=squish,
                          name="CN\n(mediana)") +
    labs(title="Perfil de CN por clone x cromossomo",
         subtitle="Azul=delecao, Vermelho=ganho, Branco=diploide",
         x="Cromossomo", y="Clone") +
    theme_bw(base_size=10) +
    theme(axis.text.x=element_text(angle=90,hjust=1,size=8),
          panel.grid=element_blank())

  pdf("results/figures/04_cn_heatmap_by_clone.pdf", width=14, height=5)
  print(p_hm)
  dev.off()
  log_info("OK results/figures/04_cn_heatmap_by_clone.pdf")
}

# =============================================================================
# 5. Boxplot distancia evolutiva
# =============================================================================
log_step(5, "Boxplot de distancia evolutiva...")

pt <- metrics_all %>% filter(!is.na(evo_dist_to_root))

p_box_clone <- ggplot(pt %>% filter(!is.na(clone_cluster)),
                       aes(reorder(clone_cluster, evo_dist_to_root, median),
                           evo_dist_to_root, fill=clone_cluster)) +
  geom_violin(alpha=.7, draw_quantiles=c(.25,.5,.75), trim=TRUE) +
  geom_jitter(width=.1, alpha=.15, size=.4) +
  scale_fill_manual(values=clone_pal, guide="none") +
  facet_wrap(~run, scales="free_x") +
  labs(title="Distancia Evolutiva por Clone (MEDICC2)",
       subtitle="Ordenados pela mediana | Menor distancia = clone ancestral",
       x="Clone", y="Distancia evolutiva ate a raiz") +
  theme_bw(base_size=10) +
  theme(axis.text.x=element_text(angle=45,hjust=1))

p_box_pat <- ggplot(pt, aes(patient, evo_dist_to_root, fill=patient)) +
  geom_violin(alpha=.7, draw_quantiles=.5, trim=TRUE) +
  geom_jitter(width=.1, alpha=.15, size=.4) +
  scale_fill_manual(values=patient_pal, guide="none") +
  facet_wrap(~run, scales="free_x") +
  labs(title="Distancia Evolutiva por Paciente",
       x="Paciente", y="Distancia evolutiva") +
  theme_bw(base_size=10)

pdf("results/figures/05_evo_dist_boxplots.pdf", width=14, height=10)
print(p_box_clone / p_box_pat)
dev.off()
log_info("OK results/figures/05_evo_dist_boxplots.pdf")

# =============================================================================
# 6. Composicao clonal por paciente
# =============================================================================
log_step(6, "Composicao clonal por paciente...")

comp_df <- metrics_all %>%
  filter(run=="global", !is.na(clone_cluster), !is.na(patient)) %>%
  count(patient, clone_cluster) %>%
  group_by(patient) %>%
  mutate(frac=n/sum(n)) %>%
  ungroup()

p_comp <- ggplot(comp_df, aes(patient, frac, fill=clone_cluster)) +
  geom_col(position="stack", width=.7) +
  scale_fill_manual(values=clone_pal, name="Clone") +
  scale_y_continuous(labels=percent) +
  labs(title="Composicao Clonal por Paciente",
       subtitle="Clones definidos por MEDICC2 + hclust (analise global)",
       x="Paciente", y="Proporcao de celulas") +
  theme_bw(base_size=11)

pdf("results/figures/06_clonal_composition.pdf", width=8, height=6)
print(p_comp)
dev.off()
log_info("OK results/figures/06_clonal_composition.pdf")

# =============================================================================
# 7. Validacao: distancia evolutiva x aberracao CNV
# =============================================================================
log_step(7, "Plot de validacao...")

for (run_name in unique(metrics_all$run)) {
  m <- metrics_all %>% filter(run==run_name, !is.na(evo_dist_to_root))
  if (nrow(m) < 10) next
  corr <- cor(m$evo_dist_to_root, m$aberration_score,
              use="complete.obs", method="spearman")
  p_v <- ggplot(m %>% sample_n(min(3000,nrow(m))),
                 aes(evo_dist_to_root, aberration_score, colour=patient)) +
    geom_point(alpha=.5, size=.7) +
    geom_smooth(method="loess", colour="black", se=TRUE) +
    scale_colour_manual(values=patient_pal) +
    labs(title=sprintf("Validacao: Dist. Evolutiva vs Aberracao CNV — %s", run_name),
         subtitle=sprintf("Spearman rho = %.3f (esperado: positivo)", corr),
         x="Distancia evolutiva MEDICC2", y="Score de aberracao CNV") +
    theme_bw(base_size=10)
  pdf(sprintf("results/figures/07_validation_%s.pdf", run_name), width=9, height=6)
  print(p_v)
  dev.off()
}
log_info("OK results/figures/07_validation_*.pdf")

# =============================================================================
# 8. Relatorio final
# =============================================================================
log_step(8, "Relatorio final...")

n_cells_tot <- nrow(metrics_all %>% filter(run=="global") %>% distinct(cell))
n_cl_global <- length(unique(na.omit(
  metrics_all$clone_cluster[metrics_all$run=="global"])))
corr_global <- tryCatch({
  m <- metrics_all %>% filter(run=="global", !is.na(evo_dist_to_root))
  cor(m$evo_dist_to_root, m$aberration_score, use="complete.obs", method="spearman")
}, error=function(e) NA)

report <- c(
  "=================================================================",
  "  RELATORIO FINAL — Filogenia Clonal CNV com MEDICC2",
  sprintf("  Dataset: GSE173279 | %s", format(Sys.time(), "%Y-%m-%d")),
  sprintf("  Diretorio: %s", WORKDIR),
  "=================================================================",
  "",
  sprintf("Pacientes:           %s", paste(patients, collapse=", ")),
  sprintf("Celulas pos-QC:      %d", n_cells_tot),
  sprintf("Clones (global):     %d", n_cl_global),
  sprintf("Validacao Spearman:  %.3f", corr_global),
  "",
  "Figuras:",
  "  results/figures/01_phylo_tree_global.pdf",
  "  results/figures/02_phylo_trees_by_patient.pdf",
  "  results/figures/03_umap_panels.pdf",
  "  results/figures/04_cn_heatmap_by_clone.pdf",
  "  results/figures/05_evo_dist_boxplots.pdf",
  "  results/figures/06_clonal_composition.pdf",
  "  results/figures/07_validation_*.pdf",
  "",
  "Tabelas:",
  "  results/tree_analysis/evo_metrics_per_cell.csv",
  "  results/tree_analysis/clone_summary.csv",
  "  results/tree_analysis/patient_summary.csv",
  "  results/qc/qc_metrics.csv",
  "================================================================="
)

writeLines(report, "results/SUMMARY_REPORT.txt")
cat(paste(report, collapse="\n"), "\n")

cat("\n=================================================================\n")
cat(sprintf(" Etapa 5 concluida: %s\n", format(Sys.time())))
cat(" Pipeline MEDICC2 completo!\n")
cat("=================================================================\n")
