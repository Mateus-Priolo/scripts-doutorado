#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(parallel)
  library(pheatmap)
  library(PDclust)
})

plot_theme <- theme_bw(base_size = 12) +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    axis.line = element_blank(),
    strip.background = element_rect(fill = "white")
  )

option_list <- list(
  make_option("--cov_dir", type = "character", help = "Directory with per-cell .cov.gz files"),
  make_option("--qc_file", type = "character", help = "analysis_scRRBS_sequencing_qc.tsv"),
  make_option("--clinical_file", type = "character", help = "clinical_metadata.tsv"),
  make_option("--outdir", type = "character", help = "Output directory"),
  make_option("--n_clusters", type = "integer", default = 6, help = "Number of PDclust clusters [default %default]"),
  make_option("--cores_pairwise", type = "integer", default = 4, help = "Cores for create_pairwise_master [default %default]"),
  make_option("--cores_read", type = "integer", default = 8, help = "Cores for reading coverage files [default %default]"),
  make_option("--include_sex_scaffold", action = "store_true", default = FALSE,
              help = "Keep GL/X/Y/MT chromosomes; by default they are removed to match Fig1b."),
  make_option("--pattern", type = "character", default = "\\.cov\\.gz$", help = "Filename regex [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

stopifnot(!is.null(opt$cov_dir), !is.null(opt$qc_file), !is.null(opt$clinical_file), !is.null(opt$outdir))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

auto_read <- function(path) {
  if (grepl("\\.tsv$", path, ignore.case = TRUE)) {
    fread(path, sep = "\t", header = TRUE, data.table = FALSE)
  } else if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    fread(path, sep = ",", header = TRUE, data.table = FALSE)
  } else {
    fread(path, header = TRUE, data.table = FALSE)
  }
}

message("[1/8] Reading QC and clinical metadata")
qc <- auto_read(opt$qc_file)
clin <- auto_read(opt$clinical_file)

needed_qc <- c("cell_barcode", "case_barcode", "cpg_unique", "bisulfite_conversion_rate", "tumor_status")
missing_qc <- setdiff(needed_qc, colnames(qc))
if (length(missing_qc) > 0) stop("Missing QC columns: ", paste(missing_qc, collapse = ", "))
if (!"case_barcode" %in% colnames(clin)) stop("clinical_file must contain case_barcode")
if (!"idh_codel_subtype" %in% colnames(clin)) stop("clinical_file must contain idh_codel_subtype")

clin <- clin %>%
  mutate(
    case_barcode = as.character(case_barcode),
    idh_status = ifelse(grepl("IDHwt", idh_codel_subtype), "IDHwt", "IDHmut")
  )

qc_pass <- qc %>%
  mutate(
    cell_barcode = as.character(cell_barcode),
    case_barcode = as.character(case_barcode)
  ) %>%
  filter(cpg_unique > 40000, bisulfite_conversion_rate > 95, tumor_status == 1) %>%
  left_join(clin, by = "case_barcode")

message("QC-pass tumor cells: ", nrow(qc_pass))

message("[2/8] Locating .cov.gz files")
files <- list.files(opt$cov_dir, pattern = opt$pattern, full.names = TRUE)
if (length(files) == 0) stop("No coverage files found in ", opt$cov_dir)
file_tbl <- data.frame(
  cell_barcode = sub("\\.cov\\.gz$", "", basename(files)),
  file_path = files,
  stringsAsFactors = FALSE
)

qc_pass <- qc_pass %>% inner_join(file_tbl, by = "cell_barcode")
if (nrow(qc_pass) == 0) stop("No QC-pass cells matched files in cov_dir")
message("QC-pass cells with coverage files: ", nrow(qc_pass))

files_use <- qc_pass$file_path
names(files_use) <- qc_pass$cell_barcode

message("[3/8] Reading coverage files and converting to PDclust input")
read_cov_to_pd <- function(f) {
  dat <- tryCatch(
    fread(cmd = sprintf("zcat < %s", shQuote(f)), showProgress = FALSE, data.table = FALSE),
    error = function(e) e
  )
  if (inherits(dat, "error")) {
    message("Read failed for ", f, ": ", dat$message)
    return(NULL)
  }
  if (ncol(dat) < 4) {
    message("Unexpected coverage format in ", f)
    return(NULL)
  }

  # Support both standard bismark columns and the prepared MethSCAn .cov.gz
  cn <- colnames(dat)
  if (all(c("chr", "start", "end", "methylation_percentage") %in% cn)) {
    dat <- dat[, c("chr", "start", "end", "methylation_percentage")]
    colnames(dat) <- c("chr", "start", "end", "meth")
  } else {
    dat <- dat[, 1:4]
    colnames(dat) <- c("chr", "start", "end", "meth")
  }

  dat$chr <- as.character(dat$chr)
  dat$start <- as.integer(dat$start)
  dat$end <- as.integer(dat$end)
  dat$meth <- as.numeric(dat$meth)
  dat <- dat %>% filter(!is.na(chr), !is.na(start), !is.na(end), !is.na(meth))

  if (!opt$include_sex_scaffold) {
    dat <- dat %>% filter(!grepl("GL|X|Y|MT", chr))
  }

  dat$chr <- ifelse(grepl("^chr", dat$chr), dat$chr, paste0("chr", dat$chr))
  dat
}

scgp_files <- mclapply(files_use, read_cov_to_pd, mc.cores = opt$cores_read)
keep_idx <- !vapply(scgp_files, is.null, logical(1))
scgp_files <- scgp_files[keep_idx]
qc_pass <- qc_pass[match(names(scgp_files), qc_pass$cell_barcode), , drop = FALSE]
message("Cells successfully loaded for PDclust: ", length(scgp_files))
if (length(scgp_files) < 10) stop("Too few cells loaded for PDclust")

message("[4/8] Pairwise dissimilarity with PDclust")
pairwise <- create_pairwise_master(scgp_files, cores_to_use = opt$cores_pairwise, digital = FALSE)
saveRDS(pairwise, file.path(opt$outdir, "pdclust_pairwise.rds"))

message("[5/8] Converting to dissimilarity matrix")
pd_mat <- convert_to_dissimilarity_matrix(pairwise)
saveRDS(pd_mat, file.path(opt$outdir, "pdclust_dissimilarity_matrix.rds"))
write.table(pd_mat, file.path(opt$outdir, "pdclust_dissimilarity_matrix.tsv"), sep = "\t", quote = FALSE)

message("[6/8] Clustering")
cluster_res <- cluster_dissimilarity(pd_mat, num_clusters = opt$n_clusters)
saveRDS(cluster_res, file.path(opt$outdir, "pdclust_cluster_results.rds"))

cluster_assignments <- as.data.frame(cluster_res$cluster_assignments)
cluster_assignments$cell_barcode <- rownames(cluster_assignments)
colnames(cluster_assignments)[1] <- "pd_cluster"
cluster_assignments <- cluster_assignments %>%
  left_join(qc_pass, by = "cell_barcode")
write.table(cluster_assignments, file.path(opt$outdir, "pdclust_cluster_assignments.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

message("[7/8] Heatmap and MDS visualization")
heatmap_pal <- colorRampPalette(RColorBrewer::brewer.pal(8, name = "YlOrRd"))(21)
pdf(file.path(opt$outdir, "PDclust_heatmap.pdf"), width = 8, height = 7)
pheatmap(
  pd_mat,
  cluster_rows = cluster_res$hclust_obj,
  cluster_cols = cluster_res$hclust_obj,
  treeheight_row = 0,
  border_color = NA,
  color = heatmap_pal,
  show_colnames = FALSE,
  show_rownames = FALSE,
  annotation_col = cluster_res$cluster_assignments
)
dev.off()

viz_df <- visualize_clusters(pd_mat, cluster_labels = cluster_res$cluster_assignments)
viz_df <- as.data.frame(viz_df)
mds_df <- viz_df %>%
  tibble::rownames_to_column("cell_barcode") %>%
  left_join(qc_pass, by = "cell_barcode") %>%
  left_join(cluster_assignments %>% select(cell_barcode, pd_cluster), by = "cell_barcode")
write.table(mds_df, file.path(opt$outdir, "pdclust_mds_embeddings.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

pdf(file.path(opt$outdir, "PDclust_MDS_by_case.pdf"), width = 7, height = 5)
print(
  ggplot(mds_df, aes(V1, V2, color = case_barcode)) +
    geom_point(alpha = 0.8, size = 1) +
    labs(x = "MDS Dimension 1", y = "MDS Dimension 2", color = "Subject") +
    plot_theme
)
dev.off()

pdf(file.path(opt$outdir, "PDclust_MDS_by_IDH.pdf"), width = 6, height = 5)
print(
  ggplot(mds_df, aes(V1, V2, color = idh_status)) +
    geom_point(alpha = 0.8, size = 1) +
    labs(x = "MDS Dimension 1", y = "MDS Dimension 2", color = "IDH") +
    plot_theme
)
dev.off()

pdf(file.path(opt$outdir, "PDclust_MDS_by_cluster.pdf"), width = 6, height = 5)
print(
  ggplot(mds_df, aes(V1, V2, color = factor(pd_cluster))) +
    geom_point(alpha = 0.8, size = 1) +
    labs(x = "MDS Dimension 1", y = "MDS Dimension 2", color = "PDclust cluster") +
    plot_theme
)
dev.off()

message("[8/8] Cluster composition summaries")
cluster_case <- mds_df %>% count(pd_cluster, case_barcode)
cluster_case_frac <- cluster_case %>% group_by(pd_cluster) %>% mutate(frac = n / sum(n)) %>% ungroup()
write.table(cluster_case, file.path(opt$outdir, "pdclust_cluster_by_case.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(cluster_case_frac, file.path(opt$outdir, "pdclust_cluster_by_case_fraction.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

message("Done. Outputs in: ", opt$outdir)
