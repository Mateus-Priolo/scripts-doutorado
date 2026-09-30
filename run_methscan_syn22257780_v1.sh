#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH --mail-type=END,FAIL

set -euo pipefail

MODE="${1:-all}"
if [[ "$MODE" != "all" && "$MODE" != "IDHwt" && "$MODE" != "IDHmut" ]]; then
  echo "Uso: sbatch run_methscan_syn22257780_v1.sh [all|IDHwt|IDHmut]"
  exit 1
fi

if [[ -n "${SLURM_SUBMIT_DIR:-}" && -f "${SLURM_SUBMIT_DIR}/config_scDNAme_multimodal_v1.sh" ]]; then
  SCRIPT_DIR="$SLURM_SUBMIT_DIR"
else
  SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi

source "$SCRIPT_DIR/config_scDNAme_multimodal_v1.sh"

OUTDIR="$PROJECT_OUT_BASE/methscan_${MODE}"
PREP_DIR="$OUTDIR/prep"
RAW_DIR="$OUTDIR/methscan_raw"
FILTERED_DIR="$OUTDIR/methscan_filtered"
VMR_BED="$OUTDIR/VMRs.bed"
DMR_DIR="$OUTDIR/dmrs"
mkdir -p "$OUTDIR" "$DMR_DIR"

run_rscript() {
  if command -v job-nanny >/dev/null 2>&1; then
    job-nanny Rscript "$@"
  else
    Rscript "$@"
  fi
}

if [[ -f "$CONDA_SH" ]]; then
  source "$CONDA_SH"
else
  echo "Conda profile não encontrado em: $CONDA_SH"
  exit 1
fi

echo "=== MethSCAn pipeline ==="
echo "Modo: $MODE"
echo "Saída: $OUTDIR"

echo "Ativando ambiente conda único: $CONDA_ENV"
conda activate "$CONDA_ENV"
echo "python: $(command -v python)"
echo "methscan: $(command -v methscan || true)"
echo "Rscript: $(command -v Rscript || true)"

echo "[1/7] Preparando inputs por célula..."
python "$SCRIPT_DIR/script_prepare_methscan_inputs_syn22257780_v1.py" \
  --coverage "$BISMARK_COVERAGE_AGG" \
  --qc "$QC_FILE" \
  --clinical "$CLINICAL_FILE" \
  --ref-promoters "$REF_PROMOTERS" \
  --ref-genes "$REF_GENES" \
  --mode "$MODE" \
  --outdir "$PREP_DIR" \
  --min-unique-cpg "$MIN_UNIQUE_CPG" \
  --min-bs "$MIN_BS_CONVERSION" \
  --require-tumor "$REQUIRE_TUMOR_STATUS" \
  --chunksize 500000

mapfile -t COV_FILES < <(find "$PREP_DIR/cov_by_cell" -type f -name '*.cov.gz' | sort)
if [[ ${#COV_FILES[@]} -eq 0 ]]; then
  echo "Nenhum arquivo .cov.gz por célula foi gerado em $PREP_DIR/cov_by_cell"
  exit 1
fi

echo "[2/7] methscan prepare..."
ulimit -n 65535 || true
methscan prepare \
  --input-format bismark \
  --round-sites \
  --chunksize "$METHSCAN_CHUNKSIZE_BP" \
  "${COV_FILES[@]}" \
  "$RAW_DIR"

echo "[3/7] methscan filter..."
methscan filter \
  --cell-names "$PREP_DIR/qc_keep_cell_names.txt" \
  --keep \
  "$RAW_DIR" \
  "$FILTERED_DIR"

echo "[4/7] methscan smooth..."
methscan smooth \
  -bw "$METHSCAN_SMOOTH_BW" \
  "$FILTERED_DIR"

echo "[5/7] methscan scan..."
methscan scan \
  -bw "$METHSCAN_SCAN_BW" \
  --stepsize "$METHSCAN_SCAN_STEP" \
  --var-threshold "$METHSCAN_VAR_THRESHOLD" \
  --min-cells "$METHSCAN_MIN_CELLS" \
  --threads "$THREADS" \
  --write-header \
  "$FILTERED_DIR" \
  "$VMR_BED"

echo "[6/7] methscan matrix (VMRs + promotores + genes)..."
methscan matrix --threads "$THREADS" "$VMR_BED" "$FILTERED_DIR" "$OUTDIR/VMR_matrix"
if [[ -s "$PREP_DIR/bed/promoters.bed" ]]; then
  methscan matrix --threads "$THREADS" "$PREP_DIR/bed/promoters.bed" "$FILTERED_DIR" "$OUTDIR/promoter_matrix"
fi
if [[ -s "$PREP_DIR/bed/genes.bed" ]]; then
  methscan matrix --threads "$THREADS" "$PREP_DIR/bed/genes.bed" "$FILTERED_DIR" "$OUTDIR/gene_matrix"
fi

echo "[7/7] Downstream em R + grupos para DMR..."
run_rscript "$SCRIPT_DIR/script_methscan_downstream_syn22257780_v1.R" \
  --mode "$MODE" \
  --results "$OUTDIR" \
  --qc "$PREP_DIR/qc_keep.tsv" \
  --clinical "$CLINICAL_FILE" \
  --npcs "$PCA_NPCS" \
  --n_neighbors "$UMAP_NEIGHBORS" \
  --min_dist "$UMAP_MIN_DIST" \
  --leiden_resolution "$LEIDEN_RESOLUTION"

if compgen -G "$OUTDIR/downstream/cell_groups/*.csv" > /dev/null; then
  echo "[extra] Running methscan diff for cluster contrasts..."
  while IFS= read -r group_file; do
    stem="$(basename "$group_file" .csv)"
    methscan diff \
      -bw "$METHSCAN_DIFF_BW" \
      --stepsize "$METHSCAN_DIFF_STEP" \
      --threshold "$METHSCAN_DIFF_THRESHOLD" \
      --min-cells "$METHSCAN_DIFF_MIN_CELLS" \
      --threads "$THREADS" \
      --write-header \
      "$FILTERED_DIR" \
      "$group_file" \
      "$DMR_DIR/${stem}.bed"
  done < <(find "$OUTDIR/downstream/cell_groups" -type f -name '*.csv' | sort)
fi

echo "Finalizado: $MODE"
