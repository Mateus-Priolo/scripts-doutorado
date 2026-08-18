#!/bin/bash
#SBATCH --job-name=dl_GSE173279
#SBATCH --chdir=/home/renanomete/projetos/matdata/evo_clonal
#SBATCH --time=01:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --output=logs/01_download_%j.out
#SBATCH --error=logs/01_download_%j.err

export INPUT=""
export OUTPUT="data/raw/GSE173279_scWGS_all_cells_cn.tsv.gz \
               data/raw/GSE173279_scWGS_all_cells_metadata.csv.gz \
               data/raw/GSE173279_scWGS_JK136_coarse_cn.tsv.gz \
               data/raw/GSE173279_scWGS_JK136_metadata.csv.gz \
               data/raw/GSE173279_scWGS_JK142_coarse_cn.tsv.gz \
               data/raw/GSE173279_scWGS_JK142_metadata.csv.gz \
               data/raw/GSE173279_scWGS_JK153_coarse_cn.tsv.gz \
               data/raw/GSE173279_scWGS_JK153_metadata.csv.gz"

set -euo pipefail

WORKDIR="/home/renanomete/projetos/matdata/evo_clonal"
cd "${WORKDIR}"
mkdir -p data/raw logs

echo "======================================================"
echo " Download GSE173279 | $(date)"
echo "======================================================"

BASEURL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE173nnn/GSE173279/suppl"

FILES=(
  "GSE173279_scWGS_all_cells_cn.tsv.gz"
  "GSE173279_scWGS_all_cells_metadata.csv.gz"
  "GSE173279_scWGS_JK136_coarse_cn.tsv.gz"
  "GSE173279_scWGS_JK136_metadata.csv.gz"
  "GSE173279_scWGS_JK142_coarse_cn.tsv.gz"
  "GSE173279_scWGS_JK142_metadata.csv.gz"
  "GSE173279_scWGS_JK153_coarse_cn.tsv.gz"
  "GSE173279_scWGS_JK153_metadata.csv.gz"
)

for f in "${FILES[@]}"; do
  DEST="data/raw/${f}"
  if [[ -f "${DEST}" && $(stat -c%s "${DEST}") -gt 1000 ]]; then
    echo "  [SKIP] ${f} ja existe"
    continue
  fi
  echo "  Baixando ${f}..."
  wget --quiet --tries=5 --timeout=120 --retry-connrefused \
       -O "${DEST}" "${BASEURL}/${f}" \
    || { echo "  [ERRO] Falha: ${f}"; exit 1; }
  echo "  OK ${f} ($(du -sh "${DEST}" | cut -f1))"
done

echo ""
echo "Verificacao:"
for f in "${FILES[@]}"; do
  SIZE=$(stat -c%s "data/raw/${f}" 2>/dev/null || echo 0)
  [[ "${SIZE}" -lt 1000 ]] && { echo "[ERRO] Arquivo vazio: ${f}"; exit 1; }
  echo "  OK data/raw/${f}"
done

echo ""
echo " Download concluido: $(date)"
echo "======================================================"
