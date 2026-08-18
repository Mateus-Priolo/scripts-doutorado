#!/bin/bash
#SBATCH -t 00:10:00
#SBATCH -c 2
#SBATCH --mem=64G
#SBATCH --mail-type=END,FAIL
#SBATCH -J explore_seurat

set -eo pipefail

source ~/miniconda3/etc/profile.d/conda.sh
set +u
conda activate R
set -u

# Diretório onde estão os objetos
OBJECT_DIR="/home/matpg/scRNA/singlecell/integrated/results_v3"
OUTDIR="${OBJECT_DIR}/exploration_output"
mkdir -p "$OUTDIR"

# CAMINHO ABSOLUTO PARA O SCRIPT R
RSCRIPT_PATH="/home/matpg/scRNA/singlecell/integrated/scripts/explore_seurat_objects.R"

if [ ! -f "$RSCRIPT_PATH" ]; then
    echo "ERRO: Script R não encontrado em $RSCRIPT_PATH"
    exit 1
fi

Rscript "$RSCRIPT_PATH" \
  --object1 "${OBJECT_DIR}/annotated_integrated_v3.rds" \
  --object2 "${OBJECT_DIR}/integrated_harmony_v3.rds" \
  --outdir "$OUTDIR"

echo "Exploração concluída. Resultados em: $OUTDIR/exploration_summary.txt"
