#!/usr/bin/env Rscript
# =============================================================================
# 02_diagnose.R — Inspeciona estrutura dos arquivos coarse_cn e metadata
# Rodar ANTES de submeter o script 02 para confirmar orientação da matriz
# Uso: Rscript scripts/02_diagnose.R
# =============================================================================

suppressPackageStartupMessages(library(data.table))

WORKDIR <- "/home/renanomete/projetos/matdata/evo_clonal"
setwd(WORKDIR)

cat("=================================================================\n")
cat(" DIAGNÓSTICO — Estrutura dos arquivos brutos\n")
cat("=================================================================\n\n")

PATIENTS <- c("JK136", "JK142", "JK153")

for (pat in PATIENTS) {
  cn_path   <- sprintf("data/raw/GSE173279_scWGS_%s_coarse_cn.tsv.gz", pat)
  meta_path <- sprintf("data/raw/GSE173279_scWGS_%s_metadata.csv.gz", pat)

  cat(sprintf("━━━ %s ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n", pat))

  # Ler apenas as 4 primeiras linhas e 6 primeiras colunas
  hdr <- fread(cn_path, nrows = 4, data.table = FALSE,
               header = TRUE, check.names = FALSE)

  cat(sprintf("  Dimensão (4 linhas lidas): %d linhas × %d colunas\n",
              nrow(hdr), ncol(hdr)))
  cat(sprintf("  Nomes das 6 primeiras colunas:\n"))
  cat(sprintf("    %s\n", paste(head(colnames(hdr), 6), collapse = " | ")))
  cat(sprintf("  Primeiros valores da coluna 1 (linhas 1-4):\n"))
  cat(sprintf("    %s\n", paste(hdr[[1]], collapse = " | ")))
  cat(sprintf("  Primeiros valores da coluna 2 (linhas 1-4):\n"))
  cat(sprintf("    %s\n", paste(hdr[[2]], collapse = " | ")))
  cat(sprintf("  Classe das colunas 1-6:\n"))
  cat(sprintf("    %s\n\n",
              paste(sapply(hdr[, 1:min(6,ncol(hdr))], class), collapse = " | ")))

  # Contar total de linhas e colunas sem carregar tudo
  n_total_lines <- as.integer(system(
    sprintf("zcat %s | wc -l", cn_path), intern = TRUE))
  n_total_cols  <- ncol(hdr)
  cat(sprintf("  Total de linhas (incluindo header): %d\n", n_total_lines))
  cat(sprintf("  Total de colunas: %d\n", n_total_cols))

  # Heurística de orientação
  frac_coord_cols <- mean(grepl("^(chr)?[0-9XYMxy]+[:\\-_][0-9]+",
                                colnames(hdr)))
  first_col_char  <- is.character(hdr[[1]])

  cat(sprintf("  Fração de colunas com padrão genômico: %.2f\n", frac_coord_cols))
  cat(sprintf("  Coluna 1 é character: %s\n", first_col_char))

  if (frac_coord_cols > 0.5) {
    cat("  ➜ ORIENTAÇÃO DETECTADA: células × bins (linhas=células, colunas=bins)\n")
    cat(sprintf("  ➜ ~%d células, ~%d bins\n",
                n_total_lines - 1, n_total_cols - as.integer(first_col_char)))
  } else {
    cat("  ➜ ORIENTAÇÃO DETECTADA: bins × células (linhas=bins, colunas=células)\n")
    cat(sprintf("  ➜ ~%d bins, ~%d células\n",
                n_total_lines - 1, n_total_cols - as.integer(first_col_char)))
  }
  cat("\n")

  # Metadata
  if (file.exists(meta_path)) {
    meta_hdr <- fread(meta_path, nrows = 3, data.table = FALSE)
    cat(sprintf("  Metadata colunas: %s\n",
                paste(colnames(meta_hdr), collapse = " | ")))
    cat(sprintf("  Metadata linha 1: %s\n\n",
                paste(as.character(meta_hdr[1, ]), collapse = " | ")))
  }
}

cat("=================================================================\n")
cat(" Copie a saída acima e ajuste o script 02 conforme necessário.\n")
cat("=================================================================\n")
