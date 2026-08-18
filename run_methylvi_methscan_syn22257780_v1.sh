#!/bin/bash
#SBATCH -t 48:00:00
#SBATCH -c 8
#SBATCH --mem=96G
#SBATCH --mail-type=END,FAIL
#SBATCH -J methylvi_scRRBS

set -euo pipefail

source ~/miniconda3/etc/profile.d/conda.sh
conda activate ${SCVI_ENV:-methylvi}

SCRIPT_DIR="/home/matpg/scDNAme/scriptsdname/scDNAme_multimodal_v2"
BASE_OUT="/home/matpg/scDNAme/results_scDNAme_multimodal_v2_nobam_safe"
CLINICAL="/home/matpg/scDNAme/tables/clinical_metadata.tsv"

MODE_ARG="${1:-all}"
MATRIX_DIR="${2:-VMR_matrix}"
RUN_DIFF="${RUN_DIFF:-1}"
N_LATENT="${N_LATENT:-10}"
MAX_EPOCHS="${MAX_EPOCHS:-250}"
BATCH_SIZE="${BATCH_SIZE:-128}"
MIN_CELLS_PER_FEATURE="${MIN_CELLS_PER_FEATURE:-10}"
MAX_FEATURES="${MAX_FEATURES:-5000}"
NNEI="${NNEI:-30}"
MIN_DIST="${MIN_DIST:-0.10}"
RES="${RES:-0.12}"
OUT_SUBDIR="${OUT_SUBDIR:-methylvi_casebatch_v1}"

if [[ "$MODE_ARG" == "all_modes" ]]; then
  MODES=(all IDHwt IDHmut)
else
  MODES=($MODE_ARG)
fi

for MODE in "${MODES[@]}"; do
  OUTDIR="${BASE_OUT}/methscan_${MODE}"
  FILTERED_DIR="${OUTDIR}/methscan_filtered"

  python "${SCRIPT_DIR}/script_methylvi_methscan_syn22257780_v1.py" \
    --mode "${MODE}" \
    --results "${OUTDIR}" \
    --matrix_dir "${MATRIX_DIR}" \
    --qc "${OUTDIR}/prep/qc_keep.tsv" \
    --clinical "${CLINICAL}" \
    --out_subdir "${OUT_SUBDIR}" \
    --batch_key case_barcode \
    --n_latent "${N_LATENT}" \
    --max_epochs "${MAX_EPOCHS}" \
    --batch_size "${BATCH_SIZE}" \
    --min_cells_per_feature "${MIN_CELLS_PER_FEATURE}" \
    --max_features "${MAX_FEATURES}" \
    --n_neighbors "${NNEI}" \
    --min_dist "${MIN_DIST}" \
    --resolution "${RES}"

  if [[ "$RUN_DIFF" == "1" ]]; then
    mkdir -p "${OUTDIR}/${OUT_SUBDIR}/dmrs"
    while IFS= read -r group_file; do
      stem="$(basename "$group_file" .csv)"
      methscan diff \
        -bw 2000 \
        --stepsize 1000 \
        --threshold 0.02 \
        --min-cells 6 \
        --threads 8 \
        --write-header \
        "${FILTERED_DIR}" \
        "$group_file" \
        "${OUTDIR}/${OUT_SUBDIR}/dmrs/${stem}.bed"
    done < <(find "${OUTDIR}/${OUT_SUBDIR}/cell_groups" -type f -name '*.csv' | sort)
  fi
done

echo "MethylVI package run finished."
