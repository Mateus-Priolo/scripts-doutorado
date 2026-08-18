#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(R.utils)
  library(ggplot2)
  library(irlba)
  library(uwot)
  library(igraph)
})

parse_args <- function() {
  x <- commandArgs(trailingOnly = TRUE)
  out <- list()
  i <- 1
  while (i <= length(x)) {
    key <- x[i]
    if (!startsWith(key, "--")) stop("Unexpected argument: ", key)
    key <- sub("^--", "", key)
    if (i == length(x) || startsWith(x[i + 1], "--")) {
      out[[key]] <- TRUE
      i <- i + 1
    } else {
      out[[key]] <- x[i + 1]
      i <- i + 2
    }
  }
  out
}

args <- parse_args()

`%||%` <- function(a, b) if (!is.null(a)) a else b

mode <- args$mode %||% "all"
results_dir <- args$results %||% stop("--results is required")
qc_file <- args$qc %||% stop("--qc is required")
clinical_file <- args$clinical %||% stop("--clinical is required")
vmr_dir <- file.path(results_dir, "VMR_matrix")
promoter_dir <- file.path(results_dir, "promoter_matrix")
gene_dir <- file.path(results_dir, "gene_matrix")
outdir <- file.path(results_dir, "downstream")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "cell_groups"), showWarnings = FALSE, recursive = TRUE)

npcs_use <- as.integer(args$npcs %||% 15)
n_neighbors <- as.integer(args$n_neighbors %||% 30)
min_dist <- as.numeric(args$min_dist %||% 0.10)
leiden_resolution <- as.numeric(args$leiden_resolution %||% 0.60)
set.seed(2)

message("=== MethSCAn downstream clustering ===")
message("Mode: ", mode)
message("Results dir: ", results_dir)

clean_barcode <- function(x) {
  x <- as.character(x)
  x <- gsub("\\.", "-", x)
  trimws(x)
}

safe_read_matrix <- function(path) {
  if (!file.exists(path)) stop("Matrix file not found: ", path)
  dt <- fread(path, sep = ",")
  rn <- dt[[1]]
  dt[[1]] <- NULL
  mat <- as.matrix(dt)
  rownames(mat) <- rn
  storage.mode(mat) <- "numeric"
  mat
}

prcomp_iterative <- function(x, n = 10, n_iter = 50, min_gain = 0.001) {
  mse <- rep(NA_real_, n_iter)
  na_loc <- is.na(x)
  x[na_loc] <- 0
  for (i in seq_len(n_iter)) {
    prev_imp <- x[na_loc]
    pr <- irlba::prcomp_irlba(x, center = FALSE, scale. = FALSE, n = n)
    new_imp <- (pr$x %*% t(pr$rotation))[na_loc]
    x[na_loc] <- new_imp
    mse[i] <- mean((prev_imp - new_imp) ^ 2)
    gain <- mse[i] / max(mse, na.rm = TRUE)
    if (is.finite(gain) && gain < min_gain) break
  }
  pr$mse_iter <- mse[seq_len(i)]
  pr
}

plot_umap <- function(df, color_col, outfile, title_text) {
  p <- ggplot(df, aes(x = UMAP1, y = UMAP2, color = .data[[color_col]])) +
    geom_point(size = 0.6, alpha = 0.9) +
    coord_fixed() +
    theme_bw() +
    labs(title = title_text, color = color_col)
  ggsave(outfile, p, width = 8, height = 6)
}

write_group_file_vs_rest <- function(meta, cluster_value, out_csv) {
  grp <- meta %>%
    transmute(cell, group = ifelse(leiden_cluster == cluster_value, "group_A", "group_B"))
  fwrite(grp, out_csv, sep = ",", col.names = FALSE)
}

write_group_file_pairwise <- function(meta, cluster_a, cluster_b, out_csv) {
  grp <- meta %>%
    transmute(cell,
              group = case_when(
                leiden_cluster == cluster_a ~ "group_A",
                leiden_cluster == cluster_b ~ "group_B",
                TRUE ~ "-"
              ))
  fwrite(grp, out_csv, sep = ",", col.names = FALSE)
}

vmr_matrix_file <- file.path(vmr_dir, "mean_shrunken_residuals.csv.gz")
if (!file.exists(vmr_matrix_file)) stop("Missing VMR mean_shrunken_residuals matrix: ", vmr_matrix_file)

meth_mtx <- safe_read_matrix(vmr_matrix_file)

qc <- fread(qc_file) %>% as.data.frame()
qc$cell_barcode <- clean_barcode(qc$cell_barcode)
qc$case_barcode <- clean_barcode(qc$case_barcode)

clin <- fread(clinical_file) %>% as.data.frame()
clin$case_barcode <- clean_barcode(clin$case_barcode)
if (!"idh_status" %in% colnames(clin)) {
  clin$idh_status <- ifelse(grepl("IDHwt", clin$idh_codel_subtype), "IDHwt", "IDHmut")
}

for (nm in c("idh_status", "time_point", "who_grade", "histological_classification")) {
  if (!nm %in% colnames(qc)) qc[[nm]] <- NA
}

clin_sub <- clin %>%
  transmute(
    case_barcode,
    idh_status_clin = idh_status,
    time_point_clin = if ("time_point" %in% colnames(clin)) time_point else NA,
    who_grade_clin = if ("who_grade" %in% colnames(clin)) who_grade else NA,
    histological_classification_clin = if ("histological_classification" %in% colnames(clin)) histological_classification else NA
  )

meta_raw <- qc %>%
  left_join(clin_sub, by = "case_barcode") %>%
  distinct(cell_barcode, .keep_all = TRUE) %>%
  filter(cell_barcode %in% rownames(meth_mtx)) %>%
  mutate(
    idh_status = dplyr::coalesce(as.character(idh_status), as.character(idh_status_clin)),
    time_point = dplyr::coalesce(as.character(time_point), as.character(time_point_clin)),
    who_grade = dplyr::coalesce(as.character(who_grade), as.character(who_grade_clin)),
    histological_classification = dplyr::coalesce(
      as.character(histological_classification),
      as.character(histological_classification_clin)
    )
  )

meta <- meta_raw %>%
  transmute(
    cell = cell_barcode,
    case_barcode,
    idh_status,
    time_point,
    who_grade,
    histological_classification
  )
meth_mtx <- meth_mtx[meta$cell, , drop = FALSE]

npc <- min(npcs_use, ncol(meth_mtx) - 1, nrow(meth_mtx) - 1)
if (npc < 2) stop("Too few cells/features after filtering to run PCA.")

message("Running iterative PCA on VMR matrix...")
pca <- meth_mtx %>%
  scale(center = TRUE, scale = FALSE) %>%
  prcomp_iterative(n = npc)

pca_tbl <- as.data.frame(pca$x) %>%
  rownames_to_column("cell")

message("Running UMAP...")
umap_obj <- uwot::umap(
  pca$x[, seq_len(npc), drop = FALSE],
  min_dist = min_dist,
  n_neighbors = min(n_neighbors, nrow(pca$x) - 1),
  seed = 2,
  ret_nn = TRUE,
  verbose = TRUE
)

umap_tbl <- as.data.frame(umap_obj$embedding)
colnames(umap_tbl) <- c("UMAP1", "UMAP2")
umap_tbl$cell <- rownames(meth_mtx)

neighbor_graph_edges <- tibble(
  from = rep(seq_len(nrow(umap_obj$nn$euclidean$idx)), times = ncol(umap_obj$nn$euclidean$idx)),
  to = as.vector(umap_obj$nn$euclidean$idx),
  weight = as.vector(umap_obj$nn$euclidean$dist)
) %>%
  filter(from != to) %>%
  mutate(
    from = rownames(meth_mtx)[from],
    to = rownames(meth_mtx)[to]
  )

message("Running Leiden clustering...")
clust_obj <- neighbor_graph_edges %>%
  igraph::graph_from_data_frame(directed = FALSE) %>%
  igraph::cluster_leiden(resolution_parameter = leiden_resolution)

clust_tbl <- tibble(
  leiden_cluster = as.character(clust_obj$membership),
  cell = clust_obj$names
) %>%
  full_join(umap_tbl, by = "cell") %>%
  full_join(meta, by = "cell")

write.table(clust_tbl, file.path(outdir, "cell_cluster_assignments.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(as.data.frame(pca$x) %>% rownames_to_column("cell"), file.path(outdir, "pca_embeddings.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

cluster_sizes <- clust_tbl %>% count(leiden_cluster, name = "n_cells")
write.table(cluster_sizes, file.path(outdir, "cluster_sizes.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

cluster_case <- clust_tbl %>% count(case_barcode, leiden_cluster, name = "n_cells")
write.table(cluster_case, file.path(outdir, "cluster_by_case.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

cluster_case_frac <- clust_tbl %>%
  count(case_barcode, leiden_cluster, name = "n_cells") %>%
  group_by(case_barcode) %>%
  mutate(frac = n_cells / sum(n_cells)) %>%
  ungroup()
write.table(cluster_case_frac, file.path(outdir, "cluster_by_case_fraction.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

plot_umap(clust_tbl, "leiden_cluster", file.path(outdir, "UMAP_clusters.pdf"), paste0("MethSCAn UMAP - ", mode))
if ("case_barcode" %in% colnames(clust_tbl)) plot_umap(clust_tbl, "case_barcode", file.path(outdir, "UMAP_cases.pdf"), paste0("MethSCAn UMAP by case - ", mode))
if ("idh_status" %in% colnames(clust_tbl)) plot_umap(clust_tbl, "idh_status", file.path(outdir, "UMAP_IDH.pdf"), paste0("MethSCAn UMAP by IDH - ", mode))

p_pca <- pca_tbl %>%
  full_join(meta, by = "cell") %>%
  ggplot(aes(x = PC1, y = PC2, color = idh_status)) +
  geom_point(size = 0.6) +
  coord_fixed() +
  theme_bw() +
  labs(title = paste0("PCA on VMR matrix - ", mode))
ggsave(file.path(outdir, "PCA_PC1_PC2.pdf"), p_pca, width = 8, height = 6)

clusters <- sort(unique(clust_tbl$leiden_cluster))
for (cl in clusters) {
  write_group_file_vs_rest(clust_tbl, cl, file.path(outdir, "cell_groups", paste0("cluster_", cl, "_vs_rest.csv")))
}
if (length(clusters) >= 2) {
  comb <- t(combn(clusters, 2))
  for (i in seq_len(nrow(comb))) {
    a <- comb[i, 1]
    b <- comb[i, 2]
    write_group_file_pairwise(clust_tbl, a, b, file.path(outdir, "cell_groups", paste0("cluster_", a, "_vs_", b, ".csv")))
  }
}

summarize_optional_matrix <- function(dir_path, label) {
  f <- file.path(dir_path, "mean_shrunken_residuals.csv.gz")
  if (!file.exists(f)) return(invisible(NULL))
  mat <- safe_read_matrix(f)
  common <- intersect(rownames(mat), clust_tbl$cell)
  if (length(common) == 0) return(invisible(NULL))
  mat <- mat[common, , drop = FALSE]
  meta_sub <- clust_tbl %>% filter(cell %in% common)
  cluster_means <- lapply(sort(unique(meta_sub$leiden_cluster)), function(cl) {
    cells <- meta_sub$cell[meta_sub$leiden_cluster == cl]
    tibble(feature = colnames(mat), leiden_cluster = cl, mean_value = colMeans(mat[cells, , drop = FALSE], na.rm = TRUE))
  }) %>% bind_rows()
  write.table(cluster_means, file.path(outdir, paste0(label, "_cluster_means.tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
}

summarize_optional_matrix(promoter_dir, "promoter")
summarize_optional_matrix(gene_dir, "gene")

sink(file.path(outdir, "sessionInfo.txt"))
print(sessionInfo())
sink()

message("Done.")
