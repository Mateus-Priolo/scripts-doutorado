#!/bin/bash
#SBATCH -t 06:00:00
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH --mail-type=END,FAIL

set -euo pipefail

MODE="${1:-all}"
if [[ "$MODE" != "all" && "$MODE" != "IDHwt" && "$MODE" != "IDHmut" ]]; then
  echo "Uso: sbatch run_disorder_syn22257780_v1.sh [all|IDHwt|IDHmut]"
  exit 1
fi

if [[ -n "${SLURM_SUBMIT_DIR:-}" && -f "${SLURM_SUBMIT_DIR}/config_scDNAme_multimodal_v1.sh" ]]; then
  SCRIPT_DIR="$SLURM_SUBMIT_DIR"
else
  SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi

source "$SCRIPT_DIR/config_scDNAme_multimodal_v1.sh"

RESULTS_DIR="$PROJECT_OUT_BASE/methscan_${MODE}"
CLUSTERS="$RESULTS_DIR/downstream/cell_cluster_assignments.tsv"
OUTDIR="$RESULTS_DIR/disorder_downstream"
mkdir -p "$OUTDIR"

if [[ ! -f "$CLUSTERS" ]]; then
  echo "Arquivo de clusters não encontrado: $CLUSTERS"
  echo "Rode primeiro o módulo MethSCAn para este modo."
  exit 1
fi

if [[ -f "$CONDA_SH" ]]; then
  source "$CONDA_SH"
else
  echo "Conda profile não encontrado em: $CONDA_SH"
  exit 1
fi
echo "Ativando ambiente conda único: $CONDA_ENV"
conda activate "$CONDA_ENV"

if command -v job-nanny >/dev/null 2>&1; then
  job-nanny Rscript "$SCRIPT_DIR/script_disorder_downstream_syn22257780_v1.R" \
    --mode "$MODE" \
    --clusters "$CLUSTERS" \
    --context "$CONTEXT_DISORDER_FILE" \
    --promoter "$PROMOTER_DISORDER_FILE" \
    --tfbs "$TFBS_DISORDER_FILE" \
    --replication "$REPLICATION_DISORDER_FILE" \
    --epiallele_context "$EPIALLELE_CONTEXT_FILE" \
    --epiallele_summary "$EPIALLELE_SUMMARY_FILE" \
    --outdir "$OUTDIR"
else
  Rscript "$SCRIPT_DIR/script_disorder_downstream_syn22257780_v1.R" \
    --mode "$MODE" \
    --clusters "$CLUSTERS" \
    --context "$CONTEXT_DISORDER_FILE" \
    --promoter "$PROMOTER_DISORDER_FILE" \
    --tfbs "$TFBS_DISORDER_FILE" \
    --replication "$REPLICATION_DISORDER_FILE" \
    --epiallele_context "$EPIALLELE_CONTEXT_FILE" \
    --epiallele_summary "$EPIALLELE_SUMMARY_FILE" \
    --outdir "$OUTDIR"
fi
