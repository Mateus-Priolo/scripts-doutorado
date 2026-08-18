#!/bin/bash

# ==============================================================
# Configuração central do pipeline scDNAme multimodal
# Ajuste estes caminhos/recursos antes de rodar no servidor.
# ==============================================================

CONFIG_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Projeto
export BASE="${BASE:-/home/matpg/scDNAme}"
export TAB="${TAB:-$BASE/tables}"
export PIPELINE_TAG="${PIPELINE_TAG:-v2_nobam_safe}"
export PROJECT_OUT_BASE="${PROJECT_OUT_BASE:-$BASE/results_scDNAme_multimodal_${PIPELINE_TAG}}"
export LOG_DIR="${LOG_DIR:-$PROJECT_OUT_BASE/logs}"

# Ambientes
export CONDA_SH="${CONDA_SH:-$HOME/miniconda3/etc/profile.d/conda.sh}"
export CONDA_ENV="${CONDA_ENV:-methscan}"
# Compatibilidade retroativa: todos apontam para o mesmo ambiente
export PY_ENV="${PY_ENV:-$CONDA_ENV}"
export METHSCAN_ENV="${METHSCAN_ENV:-$CONDA_ENV}"
export R_ENV="${R_ENV:-$CONDA_ENV}"

# Recursos
export THREADS="${THREADS:-8}"
export MEM_GB="${MEM_GB:-64}"
export UMAP_NEIGHBORS="${UMAP_NEIGHBORS:-30}"
export UMAP_MIN_DIST="${UMAP_MIN_DIST:-0.10}"
export LEIDEN_RESOLUTION="${LEIDEN_RESOLUTION:-0.60}"
export PCA_NPCS="${PCA_NPCS:-15}"

# MethSCAn
export METHSCAN_CHUNKSIZE_BP="${METHSCAN_CHUNKSIZE_BP:-5000000}"
export METHSCAN_SMOOTH_BW="${METHSCAN_SMOOTH_BW:-1000}"
export METHSCAN_SCAN_BW="${METHSCAN_SCAN_BW:-2000}"
export METHSCAN_SCAN_STEP="${METHSCAN_SCAN_STEP:-100}"
export METHSCAN_VAR_THRESHOLD="${METHSCAN_VAR_THRESHOLD:-0.02}"
export METHSCAN_MIN_CELLS="${METHSCAN_MIN_CELLS:-6}"
export METHSCAN_DIFF_BW="${METHSCAN_DIFF_BW:-2000}"
export METHSCAN_DIFF_STEP="${METHSCAN_DIFF_STEP:-1000}"
export METHSCAN_DIFF_THRESHOLD="${METHSCAN_DIFF_THRESHOLD:-0.02}"
export METHSCAN_DIFF_MIN_CELLS="${METHSCAN_DIFF_MIN_CELLS:-6}"

# Filtros QC do paper/pipeline atual
export MIN_UNIQUE_CPG="${MIN_UNIQUE_CPG:-40000}"
export MIN_BS_CONVERSION="${MIN_BS_CONVERSION:-95}"
export REQUIRE_TUMOR_STATUS="${REQUIRE_TUMOR_STATUS:-1}"

# Entradas conhecidas do projeto
export QC_FILE="${QC_FILE:-$TAB/analysis_scRRBS_sequencing_qc.tsv}"
export CLINICAL_FILE="${CLINICAL_FILE:-$TAB/clinical_metadata.tsv}"
export BISMARK_COVERAGE_AGG="/home/matpg/scDNAme/analysis_scRRBS_bismark_coverage.txt.gz"
export REF_PROMOTERS="${REF_PROMOTERS:-$TAB/ref_promoters.tsv}"
export REF_GENES="${REF_GENES:-$TAB/ref_genes.tsv}"

# Tabelas de disorder / downstream do paper
export CONTEXT_DISORDER_FILE="${CONTEXT_DISORDER_FILE:-$TAB/analysis_scRRBS_context_specific_DNAme_disorder.tsv}"
export PROMOTER_DISORDER_FILE="${PROMOTER_DISORDER_FILE:-$TAB/analysis_scRRBS_individual_promoter_DNAme_disorder.tsv}"
export TFBS_DISORDER_FILE="${TFBS_DISORDER_FILE:-$TAB/analysis_scRRBS_individual_TFBS_motif_DNAme_disorder.tsv}"
export REPLICATION_DISORDER_FILE="${REPLICATION_DISORDER_FILE:-$TAB/analysis_scRRBS_replication_timing_DNAme_disorder.tsv}"
export EPIALLELE_CONTEXT_FILE="${EPIALLELE_CONTEXT_FILE:-$TAB/analysis_RRBS_context_specific_epiallele_methylation.tsv}"
export EPIALLELE_SUMMARY_FILE="${EPIALLELE_SUMMARY_FILE:-$TAB/analysis_scRRBS_epiallele_CpG_density_summary.tsv}"

# SCNA / Ginkgo
export BAM_MANIFEST="${BAM_MANIFEST:-$BASE/bam_manifest_scRRBS.tsv}"
export GINKGO_DIR="${GINKGO_DIR:-$HOME/tools/ginkgo}"
export GINKGO_RUN_CLI="${GINKGO_RUN_CLI:-0}"
export ENABLE_SCNA="${ENABLE_SCNA:-0}"
export GINKGO_GENOME="${GINKGO_GENOME:-hg19}"
export GINKGO_BINNING="${GINKGO_BINNING:-variable_500000_101_bowtie}"

mkdir -p "$PROJECT_OUT_BASE" "$LOG_DIR"
