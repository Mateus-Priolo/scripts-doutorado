#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
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
cluster_file <- args$clusters %||% stop("--clusters is required")
outdir <- args$outdir %||% stop("--outdir is required")
context_file <- args$context %||% NA_character_
promoter_file <- args$promoter %||% NA_character_
tfbs_file <- args$tfbs %||% NA_character_
replication_file <- args$replication %||% NA_character_
epiallele_context_file <- args$epiallele_context %||% NA_character_
epiallele_summary_file <- args$epiallele_summary %||% NA_character_
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

clean_barcode <- function(x) {
  x <- as.character(x)
  x <- gsub("\\.", "-", x)
  trimws(x)
}

pick_first <- function(cols, candidates) {
  idx <- match(tolower(candidates), tolower(cols))
  idx <- idx[!is.na(idx)]
  if (length(idx) == 0) return(NA_character_)
  cols[idx[1]]
}

safe_fread <- function(path) {
  if (is.na(path) || !file.exists(path)) return(NULL)
  fread(path) %>% as.data.frame()
}

cluster_meta <- fread(cluster_file) %>% as.data.frame()
cluster_meta$cell <- clean_barcode(cluster_meta$cell)

context_dt <- safe_fread(context_file)
if (!is.null(context_dt)) {
  cell_col <- pick_first(colnames(context_dt), c("cell_barcode", "cell", "barcode", "cell_id"))
  if (is.na(cell_col)) stop("Could not identify cell column in context disorder table.")
  context_dt[[cell_col]] <- clean_barcode(context_dt[[cell_col]])
  context_join <- context_dt %>%
    rename(cell = !!cell_col) %>%
    inner_join(cluster_meta, by = "cell")

  if (all(c("PDR", "promoter_PDR", "enhancer_PDR", "cgi_PDR") %in% colnames(context_join))) {
    long_ctx <- context_join %>%
      select(cell, leiden_cluster, case_barcode, idh_status, PDR, promoter_PDR, enhancer_PDR, cgi_PDR) %>%
      pivot_longer(cols = c(PDR, promoter_PDR, enhancer_PDR, cgi_PDR), names_to = "context", values_to = "value")
  } else {
    context_col <- pick_first(colnames(context_join), c("context", "annotation", "feature_class", "feature"))
    value_col <- pick_first(colnames(context_join), c("PDR", "DNAme_disorder", "disorder", "value"))
    if (is.na(context_col) || is.na(value_col)) {
      warning("Context disorder table found, but no compatible columns were identified. Skipping context summary.")
      long_ctx <- NULL
    } else {
      long_ctx <- context_join %>%
        transmute(cell, leiden_cluster, case_barcode, idh_status, context = .data[[context_col]], value = as.numeric(.data[[value_col]]))
    }
  }

  if (!is.null(long_ctx)) {
    cluster_summary <- long_ctx %>%
      group_by(leiden_cluster, context) %>%
      summarise(n = sum(is.finite(value)), mean_value = mean(value, na.rm = TRUE), median_value = median(value, na.rm = TRUE), .groups = "drop")
    fwrite(cluster_summary, file.path(outdir, "context_disorder_by_cluster.tsv"), sep = "\t")

    case_summary <- long_ctx %>%
      group_by(case_barcode, context) %>%
      summarise(n = sum(is.finite(value)), mean_value = mean(value, na.rm = TRUE), median_value = median(value, na.rm = TRUE), .groups = "drop")
    fwrite(case_summary, file.path(outdir, "context_disorder_by_case.tsv"), sep = "\t")

    p <- ggplot(long_ctx, aes(x = factor(leiden_cluster), y = value)) +
      geom_boxplot(outlier.size = 0.2) +
      facet_wrap(~ context, scales = "free_y") +
      theme_bw() +
      xlab("MethSCAn cluster") +
      ylab("DNAme disorder / PDR") +
      ggtitle(paste0("DNAme disorder by cluster - ", mode))
    ggsave(file.path(outdir, "Boxplot_context_disorder_by_cluster.pdf"), p, width = 11, height = 7)
  }
}

summarize_feature_disorder <- function(path, feature_label, out_prefix) {
  dt <- safe_fread(path)
  if (is.null(dt)) return(invisible(NULL))
  cell_col <- pick_first(colnames(dt), c("cell_barcode", "cell", "barcode", "cell_id"))
  feature_col <- pick_first(colnames(dt), c("promoter_id", "gene_name", "gene_symbol", "tfbs_motif", "motif", "feature", "context"))
  value_col <- pick_first(colnames(dt), c("PDR", "DNAme_disorder", "disorder", "value"))
  if (is.na(cell_col) || is.na(feature_col) || is.na(value_col)) {
    warning("Skipping ", feature_label, ": unable to infer columns.")
    return(invisible(NULL))
  }
  dt[[cell_col]] <- clean_barcode(dt[[cell_col]])
  joined <- dt %>%
    rename(cell = !!cell_col) %>%
    inner_join(cluster_meta, by = "cell") %>%
    transmute(cell, leiden_cluster, case_barcode, feature = .data[[feature_col]], value = as.numeric(.data[[value_col]]))

  feat_cluster <- joined %>%
    group_by(leiden_cluster, feature) %>%
    summarise(n = sum(is.finite(value)), mean_value = mean(value, na.rm = TRUE), .groups = "drop")
  fwrite(feat_cluster, file.path(outdir, paste0(out_prefix, "_by_cluster.tsv")), sep = "\t")

  top_var <- feat_cluster %>%
    group_by(feature) %>%
    summarise(sd_cluster = sd(mean_value, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(sd_cluster)) %>%
    slice_head(n = 100)
  fwrite(top_var, file.path(outdir, paste0("top100_variable_", out_prefix, ".tsv")), sep = "\t")
}

summarize_feature_disorder(promoter_file, "promoter", "promoter_disorder")
summarize_feature_disorder(tfbs_file, "tfbs", "tfbs_disorder")
summarize_feature_disorder(replication_file, "replication", "replication_disorder")

if (!is.na(epiallele_context_file) && file.exists(epiallele_context_file)) {
  epi <- fread(epiallele_context_file) %>% as.data.frame()
  cell_col <- pick_first(colnames(epi), c("cell_barcode", "cell", "barcode", "cell_id"))
  if (!is.na(cell_col)) {
    epi[[cell_col]] <- clean_barcode(epi[[cell_col]])
    epi_join <- epi %>% rename(cell = !!cell_col) %>% inner_join(cluster_meta, by = "cell")
    fwrite(epi_join, file.path(outdir, "epiallele_context_joined.tsv"), sep = "\t")
  }
}

if (!is.na(epiallele_summary_file) && file.exists(epiallele_summary_file)) {
  epi_sum <- fread(epiallele_summary_file) %>% as.data.frame()
  cell_col <- pick_first(colnames(epi_sum), c("cell_barcode", "cell", "barcode", "cell_id"))
  if (!is.na(cell_col)) {
    epi_sum[[cell_col]] <- clean_barcode(epi_sum[[cell_col]])
    epi_sum_join <- epi_sum %>% rename(cell = !!cell_col) %>% inner_join(cluster_meta, by = "cell")
    fwrite(epi_sum_join, file.path(outdir, "epiallele_density_summary_joined.tsv"), sep = "\t")
  }
}

sink(file.path(outdir, "sessionInfo.txt"))
print(sessionInfo())
sink()

message("Done disorder downstream.")
