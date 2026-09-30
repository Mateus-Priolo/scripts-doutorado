#!/bin/bash
#SBATCH --job-name=scwgs_analysis
#SBATCH --output=/home/renanomete/projetos/matdata/new_evo_anal/scWGS/GSE173279/results/evo_clonal/logs/scwgs_%a_%j.out
#SBATCH --error=/home/renanomete/projetos/matdata/new_evo_anal/scWGS/GSE173279/results/evo_clonal/logs/scwgs_%a_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --array=0-2   # JK136, JK153, global

PATIENTS=("JK136" "JK153" "global")
PATIENT="${PATIENTS[$SLURM_ARRAY_TASK_ID]}"

BASE_DIR="/home/renanomete/projetos/matdata/new_evo_anal/scWGS/GSE173279/results/evo_clonal"
SCRIPT="${BASE_DIR}/scripts/scwgs_analysis.R"

echo "=========================================="
echo " Running scWGS analysis for $PATIENT"
echo " Job ID: $SLURM_JOB_ID"
echo " Host: $(hostname)"
echo " Date: $(date)"
echo "=========================================="

# Carregar módulo Conda
module load miniconda/24.4.0-libmamba

# Ativar ambiente
conda activate evo_clonal_medicc2

# Executar script R
Rscript "$SCRIPT" "$PATIENT"

echo "Done."
