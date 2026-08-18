nalysis.sh — Análise da árvore filogenética MEDICC2
# Depende de: 03_medicc2.sh  +  data/processed/cnv_qc_object.rds
# =============================================================================
#SBATCH --job-name=tree_analysis
#SBATCH --chdir=/home/renanomete/projetos/matdata/evo_clonal
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=logs/04_tree_%j.out
#SBATCH --error=logs/04_tree_%j.err

set -euo pipefail

WORKDIR="/home/renanomete/projetos/matdata/evo_clonal"
cd "${WORKDIR}"

module load miniconda/24.4.0-libmamba
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate evo_clonal_medicc2

mkdir -p results/tree_analysis logs

echo "======================================================"
echo " ETAPA 4 — Análise da Árvore MEDICC2"
echo " $(date) | CPUs: ${SLURM_CPUS_PER_TASK}"
echo "======================================================"

# Verificar pré-requisitos
for F in \
  data/processed/cnv_qc_object.rds \
  scripts/04_tree_analysis.R \
  scripts/functions.R; do
  if [[ ! -f "${F}" ]]; then
    echo "[ERRO] Arquivo não encontrado: ${F}" >&2
    exit 1
  fi
done

# Verificar se ao menos um resultado MEDICC2 existe
N_MEDICC2=$(find results/medicc2/ -name "*.new" -o -name "*.nwk" \
            2>/dev/null | wc -l)
if [[ "${N_MEDICC2}" -eq 0 ]]; then
  echo "[ERRO] Nenhuma árvore MEDICC2 encontrada em results/medicc2/" >&2
  echo "       Execute o script 03 antes de submeter o 04." >&2
  exit 1
fi
echo "  Árvores MEDICC2 encontradas: ${N_MEDICC2}"

echo "[$(date)] Iniciando Rscript 04..."
job-nanny Rscript scripts/04_tree_analysis.R

echo ""
echo "======================================================"
echo " Etapa 4 concluída: $(date)"
echo "======================================================"
