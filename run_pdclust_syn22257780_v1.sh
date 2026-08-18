#!/bin/bash
#SBATCH -t 48:00:00
#SBATCH -c 8
#SBATCH --mem=128G
#SBATCH --mail-type=END,FAIL
#SBATCH -J pdclust_fig1b

set -euo pipefail

source ~/miniconda3/etc/profile.d/conda.sh
conda activate methscan

WORKDIR="/home/matpg/scDNAme/scriptsdname/scDNAme_multimodal_v2"
OUTDIR="/home/matpg/scDNAme/results_scDNAme_multimodal_v2_nobam_safe/pdclust_all"
COV_DIR="/home/matpg/scDNAme/results_scDNAme_multimodal_v2_nobam_safe/methscan_all/prep/cov_by_cell"
QC_FILE="/home/matpg/scDNAme/tables/analysis_scRRBS_sequencing_qc.tsv"
CLINICAL_FILE="/home/matpg/scDNAme/tables/clinical_metadata.tsv"

mkdir -p "$OUTDIR"
cd "$WORKDIR"

Rscript "$WORKDIR/script_pdclust_syn22257780_v1.R" \
  --cov_dir "$COV_DIR" \
  --qc_file "$QC_FILE" \
  --clinical_file "$CLINICAL_FILE" \
  --outdir "$OUTDIR" \
  --n_clusters 6 \
  --cores_pairwise 4 \
  --cores_read 8
