#!/bin/bash
# =============================================================================
# submit_all.sh — Submissão encadeada do pipeline MEDICC2
# GSE173279 | GridUNESP
#
# Uso:
#   bash submit_all.sh           # submete tudo do início
#   bash submit_all.sh --from 3  # retoma a partir da etapa 3
# =============================================================================

set -euo pipefail

WORKDIR="/home/renanomete/projetos/matdata/evo_clonal"
cd "${WORKDIR}"

mkdir -p logs data/raw data/processed data/medicc2_input \
          results/qc results/medicc2/{global,JK136,JK142,JK153} \
          results/tree_analysis results/figures

echo "======================================================"
echo " Pipeline MEDICC2 — Filogenia Clonal CNV"
echo " Dataset: GSE173279 | GridUNESP"
echo " Diretório: ${WORKDIR}"
echo " $(date)"
echo "======================================================"

# Parsear --from
START_FROM=1
while [[ $# -gt 0 ]]; do
  case $1 in
    --from) START_FROM="$2"; shift 2 ;;
    *) echo "Argumento desconhecido: $1"; exit 1 ;;
  esac
done
echo "Iniciando a partir da etapa: ${START_FROM}"
echo ""

submit_job() {
  local script="$1"
  local dep="$2"
  if [[ "${dep}" == "none" ]]; then
    sbatch --parsable "${script}"
  else
    sbatch --parsable --dependency=afterok:"${dep}" "${script}"
  fi
}

check_file() {
  local f="$1"
  if [[ ! -f "${f}" ]]; then
    echo "  ERRO: arquivo necessário não encontrado: ${f}"
    exit 1
  fi
}

# ── Etapa 1: Download ──────────────────────────────────────────────────────
JOB_DL="none"
if [[ "${START_FROM}" -le 1 ]]; then
  echo "[1] Submetendo download..."
  JOB_DL=$(submit_job "scripts/01_download.sh" "none")
  echo "    Job ID: ${JOB_DL}"
else
  echo "[1] Pulando download."
  check_file "data/raw/GSE173279_scWGS_all_cells_cn.tsv.gz"
fi

# ── Etapa 2: QC e preparação MEDICC2 ──────────────────────────────────────
if [[ "${START_FROM}" -le 2 ]]; then
  echo "[2] Submetendo QC e preparação de input..."
  JOB_QC=$(submit_job "scripts/02_cnv_qc_prepare.sh" "${JOB_DL}")
  echo "    Job ID: ${JOB_QC} (dep: ${JOB_DL})"
else
  echo "[2] Pulando QC."
  check_file "data/processed/cnv_qc_object.rds"
  check_file "data/medicc2_input/all_cells_medicc2.tsv"
  JOB_QC="none"
fi

# ── Etapa 3: MEDICC2 ───────────────────────────────────────────────────────
if [[ "${START_FROM}" -le 3 ]]; then
  echo "[3] Submetendo MEDICC2..."
  JOB_M2=$(submit_job "scripts/03_medicc2_run.sh" "${JOB_QC}")
  echo "    Job ID: ${JOB_M2} (dep: ${JOB_QC})"
  echo "    ⚠  MEDICC2 pode levar até 24h para datasets grandes."
else
  echo "[3] Pulando MEDICC2."
  JOB_M2="none"
fi

# ── Etapa 4: Análise da árvore ─────────────────────────────────────────────
if [[ "${START_FROM}" -le 4 ]]; then
  echo "[4] Submetendo análise da árvore..."
  JOB_TREE=$(submit_job "scripts/04_tree_analysis.sh" "${JOB_M2}")
  echo "    Job ID: ${JOB_TREE} (dep: ${JOB_M2})"
else
  echo "[4] Pulando análise da árvore."
  check_file "data/processed/tree_analysis.rds"
  JOB_TREE="none"
fi

# ── Etapa 5: Visualizações ─────────────────────────────────────────────────
if [[ "${START_FROM}" -le 5 ]]; then
  echo "[5] Submetendo visualizações..."
  JOB_VIZ=$(submit_job "scripts/05_visualization.sh" "${JOB_TREE}")
  echo "    Job ID: ${JOB_VIZ} (dep: ${JOB_TREE})"
else
  echo "[5] Nada a submeter."
fi

echo ""
echo "======================================================"
echo " Encadeamento completo:"
[[ "${START_FROM}" -le 1 ]] && echo "  Etapa 1 — Download:              ${JOB_DL}"
[[ "${START_FROM}" -le 2 ]] && echo "  Etapa 2 — QC + input MEDICC2:   ${JOB_QC}"
[[ "${START_FROM}" -le 3 ]] && echo "  Etapa 3 — MEDICC2 (filo.):      ${JOB_M2}"
[[ "${START_FROM}" -le 4 ]] && echo "  Etapa 4 — Análise da árvore:    ${JOB_TREE}"
[[ "${START_FROM}" -le 5 ]] && echo "  Etapa 5 — Visualizações:        ${JOB_VIZ}"
echo ""
echo " Monitorar:"
echo "   squeue -u renanomete"
echo "   tail -f logs/03_medicc2_<JOB_ID>.out"
echo ""
echo " Retomar de uma etapa após falha:"
echo "   bash submit_all.sh --from <ETAPA>"
echo "======================================================"
