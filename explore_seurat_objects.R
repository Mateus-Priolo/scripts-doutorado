#!/usr/bin/env Rscript

# Explorar objetos Seurat integrados: colunas do objeto e meta.data, com sumários
# Uso: Rscript explore_seurat_objects.R --object1 <path> --object2 <path> --outdir <dir>

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
})

# Função para gerar sumário de uma coluna (vetor)
summarize_column <- function(x, colname) {
  cat("\n### Column:", colname, "\n")
  cat("Class:", class(x), "\n")
  if (is.numeric(x)) {
    cat("Summary (min, 1st Qu., median, mean, 3rd Qu., max, NA count):\n")
    print(summary(x))
    cat("Number of NAs:", sum(is.na(x)), "\n")
  } else if (is.factor(x) || is.character(x) || is.logical(x)) {
    tab <- table(x, useNA = "ifany")
    cat("Table (top 20 levels if many):\n")
    if (length(tab) > 20) {
      print(head(tab, 20))
      cat("... and", length(tab) - 20, "more levels\n")
    } else {
      print(tab)
    }
  } else {
    cat("Non-standard type, showing structure:\n")
    str(x, vec.len = 5)
  }
  cat("\n")
}

# Parsing simples de argumentos
args <- commandArgs(trailingOnly = TRUE)
object1_path <- NULL
object2_path <- NULL
outdir <- "."

i <- 1
while (i <= length(args)) {
  if (args[i] == "--object1") {
    object1_path <- args[i+1]
    i <- i + 2
  } else if (args[i] == "--object2") {
    object2_path <- args[i+1]
    i <- i + 2
  } else if (args[i] == "--outdir") {
    outdir <- args[i+1]
    i <- i + 2
  } else {
    i <- i + 1
  }
}

if (is.null(object1_path) || is.null(object2_path)) {
  stop("Please provide both --object1 and --object2 paths")
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Redirecionar saída para um arquivo de log
log_file <- file.path(outdir, "exploration_summary.txt")
sink(log_file, split = TRUE)  # envia para o arquivo e também para o console

cat("========================================\n")
cat("Exploração de objetos Seurat integrados\n")
cat("Data:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================\n\n")

# Função para explorar um objeto
explore_object <- function(obj_path, obj_name) {
  cat("\n####################################################################\n")
  cat("## OBJETO:", obj_name, "\n")
  cat("## Caminho:", obj_path, "\n")
  cat("####################################################################\n")
  
  if (!file.exists(obj_path)) {
    cat("ERRO: Arquivo não encontrado!\n")
    return()
  }
  
  cat("\nCarregando objeto...\n")
  obj <- readRDS(obj_path)
  
  cat("\n--- INFORMAÇÕES GERAIS ---\n")
  cat("Classe do objeto:", class(obj), "\n")
  cat("Número de células (ncol):", ncol(obj), "\n")
  cat("Número de genes (nrow):", nrow(obj), "\n")
  cat("Assays disponíveis:", names(obj@assays), "\n")
  cat("Reduções disponíveis:", names(obj@reductions), "\n")
  
  # Explorar meta.data
  cat("\n--- META.DATA ---\n")
  meta <- obj@meta.data
  cat("Número de colunas em meta.data:", ncol(meta), "\n")
  cat("Nomes das colunas:\n")
  print(colnames(meta))
  
  # Para cada coluna, gerar sumário
  cat("\n--- SUMÁRIO DETALHADO POR COLUNA ---\n")
  for (col in colnames(meta)) {
    summarize_column(meta[[col]], col)
  }
  
  # Verificar se existe a coluna 'seurat_clusters' ou 'integrated_snn_res.0.8' etc.
  cat("\n--- COLUNAS ADICIONAIS DE INTERESSE (se existirem) ---\n")
  possible_cluster_cols <- grep("res\\.|cluster", colnames(meta), value = TRUE, ignore.case = TRUE)
  if (length(possible_cluster_cols) > 0) {
    cat("Possíveis colunas de cluster:", paste(possible_cluster_cols, collapse = ", "), "\n")
  } else {
    cat("Nenhuma coluna de cluster óbvia encontrada.\n")
  }
  
  # Informações sobre identidades ativas
  cat("\n--- Identidades ativas (Idents) ---\n")
  print(table(Idents(obj), useNA = "ifany"))
  
  cat("\n")
}

# Explorar os dois objetos
explore_object(object1_path, "annotated_integrated_v3")
explore_object(object2_path, "integrated_harmony_v3")

cat("\n========================================\n")
cat("Fim da exploração.\n")
cat("Arquivo de saída:", log_file, "\n")
cat("========================================\n")

sink()  # fechar redirecionamento
