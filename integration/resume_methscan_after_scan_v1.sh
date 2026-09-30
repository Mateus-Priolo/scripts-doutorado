#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH --mail-type=END,FAIL

set -euo pipefail

source ~/miniconda3/etc/profile.d/conda.sh
conda activate methscan

SCRIPT_DIR="/home/matpg/scDNAme/scriptsdname/scDNAme_multimodal_v2"
BASE_OUT="/home/matpg/scDNAme/results_scDNAme_multimodal_v2_nobam_safe"
CLINICAL="/home/matpg/scDNAme/tables/clinical_metadata.tsv"
PROM_BED="/home/matpg/scDNAme/annotations_methscan/promoters.noheader.bed"
GENE_BED="/home/matpg/scDNAme/annotations_methscan/genes.noheader.bed"

for MODE in all IDHwt IDHmut; do
  OUTDIR="${BASE_OUT}/methscan_${MODE}"
  PREP_DIR="${OUTDIR}/prep"
  FILTERED_DIR="${OUTDIR}/methscan_filtered"
  DMR_DIR="${OUTDIR}/dmrs"

  mkdir -p "$DMR_DIR"

  echo "======================================"
  echo "MODE: ${MODE}"
  echo "OUTDIR: ${OUTDIR}"
  echo "======================================"

  echo "[6a] matrix em VMRs"
  methscan matrix --threads 8 \
    "${OUTDIR}/VMRs.noheader.bed" \
    "${FILTERED_DIR}" \
    "${OUTDIR}/VMR_matrix"

  echo "[6b] matrix em promotores"
  methscan matrix --threads 8 \
    "${PROM_BED}" \
    "${FILTERED_DIR}" \
    "${OUTDIR}/promoter_matrix"

  echo "[6c] matrix em genes"
  methscan matrix --threads 8 \
    "${GENE_BED}" \
    "${FILTERED_DIR}" \
    "${OUTDIR}/gene_matrix"

  echo "[7] downstream em R"
  Rscript "${SCRIPT_DIR}/script_methscan_downstream_syn22257780_v1.R" \
    --mode "${MODE}" \
    --results "${OUTDIR}" \
    --qc "${PREP_DIR}/qc_keep.tsv" \
    --clinical "${CLINICAL}" \
    --npcs 15 \
    --n_neighbors 30 \
    --min_dist 0.10 \
    --leiden_resolution 0.60

  echo "[8] methscan diff"
  if compgen -G "${OUTDIR}/downstream/cell_groups/*.csv" > /dev/null; then
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
        "${DMR_DIR}/${stem}.bed"
    done < <(find "${OUTDIR}/downstream/cell_groups" -type f -name '*.csv' | sort)
  fi

done

echo "Resume MethSCAn finished."
