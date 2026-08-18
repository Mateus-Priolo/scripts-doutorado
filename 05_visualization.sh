#!/bin/bash
#SBATCH --job-name=viz
#SBATCH --chdir=/home/renanomete/projetos/matdata/evo_clonal
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=logs/05_viz_%j.out
#SBATCH --error=logs/05_viz_%j.err

export INPUT="data/processed/tree_analysis.rds \
              data/processed/cnv_qc_object.rds \
              scripts/05_visualization.R \
              scripts/functions.R"
export OUTPUT="results/figures/ \
               results/SUMMARY_REPORT.txt"

set -euo pipefail

WORKDIR="/home/renanomete/projetos/matdata/evo_clonal"
cd "${WORKDIR}"

module load miniconda/24.4.0-libmamba
# Inicializar Conda para o shell
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate evo_clonal_medicc2

mkdir -p results/figures logs

echo "[$(date)] Iniciando etapa 5: visualizacoes"

job-nanny Rscript scripts/05_visualization.R

echo "[$(date)] Pipeline MEDICC2 completo!"
