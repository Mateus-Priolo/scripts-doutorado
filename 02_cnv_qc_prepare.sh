#!/bin/bash
# =============================================================================
# 02_cnv_qc_prepare.sh — QC scWGS + geração de inputs MEDICC2
# =============================================================================
#SBATCH --job-name=cnv_qc_prep
#SBATCH --chdir=/home/renanomete/projetos/matdata/evo_clonal
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=logs/02_qc_%j.out
#SBATCH --error=logs/02_qc_%j.err

set -euo pipefail

WORKDIR="/home/renanomete/projetos/matdata/evo_clonal"
cd "${WORKDIR}"

module load miniconda/24.4.0-libmamba
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate evo_clonal_medicc2

mkdir -p data/processed data/medicc2_input results/qc logs

echo "======================================================"
echo " ETAPA 2 — QC scWGS + preparação MEDICC2"
echo " $(date) | CPUs: ${SLURM_CPUS_PER_TASK} | MEM: 64G"
echo "======================================================"

for PAT in JK136 JK142 JK153; do
  CN_FILE="data/raw/GSE173279_scWGS_${PAT}_coarse_cn.tsv.gz"
  META_FILE="data/raw/GSE173279_scWGS_${PAT}_metadata.csv.gz"
  if [[ ! -f "${CN_FILE}" ]]; then
    echo "  [ERRO] Arquivo não encontrado: ${CN_FILE}" >&2
    exit 1
  fi
  if [[ ! -f "${META_FILE}" ]]; then
    echo "  [AVISO] Metadata não encontrado: ${META_FILE} (continuando sem ele)"
  fi
done

echo "[$(date)] Iniciando Rscript 02..."
Rscript scripts/02_cnv_qc_prepare.R

echo ""
echo "======================================================"
echo " Etapa 2 concluída: $(date)"
echo "======================================================"

echo ""
echo "Verificando outputs para etapa 3:"
ALL_OK=1
for F in \
  data/processed/cnv_qc_object.rds \
  data/medicc2_input/all_cells_medicc2.tsv \
  data/medicc2_input/JK136_medicc2.tsv \
  data/medicc2_input/JK142_medicc2.tsv \
  data/medicc2_input/JK153_medicc2.tsv; do
  if [[ -f "${F}" ]]; then
    SIZE=$(du -sh "${F}" | cut -f1)
    echo "  [OK] ${F}  (${SIZE})"
  else
    echo "  [!!] AUSENTE: ${F}" >&2
    ALL_OK=0
  fi
done

[[ "${ALL_OK}" -eq 0 ]] && { echo "[ERRO] Outputs incompletos."; exit 1; }

echo ""
echo "Pronto para submeter: sbatch scripts/03_medicc2.sh"
