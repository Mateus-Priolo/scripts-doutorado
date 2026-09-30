#!/bin/bash
# =============================================================================
# 00_setup_conda.sh — Ambiente Conda para pipeline MEDICC2 (scWGS / CNV)
# GSE173279 | GridUNESP
#
# Execute UMA VEZ em sessão interativa ANTES de submeter qualquer job:
#   srun --pty --mem=8G --cpus-per-task=2 --time=02:00:00 bash
#   cd /home/renanomete/projetos/matdata/evo_clonal
#   bash 00_setup_conda.sh
# =============================================================================

set -euo pipefail

WORKDIR="/home/renanomete/projetos/matdata/evo_clonal"
ENV_NAME="evo_clonal_medicc2"

cd "${WORKDIR}"

echo "======================================================"
echo " Setup: ${ENV_NAME}"
echo " Diretório: ${WORKDIR}"
echo " $(date)"
echo "======================================================"

module load miniconda/24.4.0-libmamba

# Remover ambiente anterior se existir
if conda env list | grep -q "^${ENV_NAME}"; then
  echo "[!] Removendo ambiente anterior '${ENV_NAME}'..."
  conda env remove -n "${ENV_NAME}" -y
fi

# ── 1. Criar ambiente Python + R base ─────────────────────────────────────────
echo ""
echo "[1] Criando ambiente base (Python 3.10 + R 4.3)..."
conda create -n "${ENV_NAME}" \
  -c conda-forge -c bioconda \
  python=3.10 \
  r-base=4.3 \
  r-essentials \
  -y

source activate "${ENV_NAME}"

# ── 2. MEDICC2 — via conda-forge (instalação oficial) ─────────────────────────
# MEDICC2: Turajlic et al. 2021, Genome Biology
# Computa distâncias de edição ponderadas em CN inteiro e infere árvore
# filogenética por Maximum Parsimony / Neighbor-Joining
echo ""
echo "[2] Instalando MEDICC2..."
conda install -n "${ENV_NAME}" \
  -c conda-forge -c bioconda \
  medicc2 \
  -y

# Verificar instalação
medicc2 --version || { echo "ERRO: medicc2 não instalado corretamente"; exit 1; }

# ── 3. Dependências Python para pré/pós-processamento ─────────────────────────
echo ""
echo "[3] Instalando dependências Python..."
conda install -n "${ENV_NAME}" \
  -c conda-forge \
  pandas \
  numpy \
  scipy \
  matplotlib \
  seaborn \
  ete3 \
  -y

# ── 4. Pacotes R para QC, visualização e análise da árvore ───────────────────
echo ""
echo "[4] Instalando pacotes R..."
conda install -n "${ENV_NAME}" \
  -c conda-forge \
  r-data.table \
  r-matrix \
  r-ggplot2 \
  r-patchwork \
  r-dplyr \
  r-tidyr \
  r-viridis \
  r-rcolorbrewer \
  r-scales \
  r-irlba \
  r-uwot \
  r-rann \
  r-igraph \
  r-harmony \
  r-future \
  r-future.apply \
  -y

# ── 5. Pacotes R para árvore filogenética e genomics ──────────────────────────
echo ""
echo "[5] Instalando pacotes R para filogenia e genômica..."
conda install -n "${ENV_NAME}" \
  -c conda-forge -c bioconda \
  r-ape \
  r-phangorn \
  r-ggtree \
  bioconductor-genomicranges \
  bioconductor-iranges \
  -y

# ── 6. Pacotes R via CRAN (não disponíveis no conda) ─────────────────────────
echo ""
echo "[6] Instalando pacotes R adicionais via CRAN..."
Rscript - <<'REOF'
options(repos = c(CRAN = "https://cran.r-project.org"))
pkgs <- c("ggtreeExtra", "aplot", "tidytree", "treeio")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    tryCatch(
      install.packages(p, quiet = TRUE),
      error = function(e) {
        # Tentar via BiocManager
        if (!requireNamespace("BiocManager", quietly = TRUE))
          install.packages("BiocManager")
        BiocManager::install(p, ask = FALSE, quiet = TRUE)
      }
    )
  }
}
REOF

# ── 7. Verificação final ───────────────────────────────────────────────────────
echo ""
echo "[7] Verificação final..."

echo "--- Python / MEDICC2 ---"
python -c "import medicc2; print('medicc2:', medicc2.__version__)" 2>/dev/null || \
  medicc2 --version

python -c "import pandas, numpy, scipy, matplotlib, ete3; print('Deps Python: OK')"

echo "--- R packages ---"
Rscript - <<'REOF'
pkgs <- c("data.table", "Matrix", "ggplot2", "patchwork", "dplyr", "tidyr",
          "harmony", "irlba", "uwot", "ape", "igraph",
          "GenomicRanges", "future")
for (p in pkgs) {
  ok <- requireNamespace(p, quietly = TRUE)
  cat(sprintf("  %-20s %s\n", p, if (ok) "OK" else "FALHOU"))
}
REOF

echo ""
echo "======================================================"
echo " Ambiente '${ENV_NAME}' pronto!"
echo " Ative com: source activate ${ENV_NAME}"
echo " $(date)"
echo "======================================================"
