#!/bin/bash
#SBATCH -t 04:20:00
#SBATCH -c 4
#SBATCH --mem=64G
#SBATCH --mail-type=FAIL

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config_scDNAme_multimodal_v1.sh"

MODULE="${1:-all}"

submit_job() {
  local script="$1"
  shift
  sbatch --parsable "$script" "$@"
}

submit_dependent_job() {
  local dep="$1"
  local script="$2"
  shift 2
  sbatch --parsable --dependency="afterok:${dep}" "$script" "$@"
}

case "$MODULE" in
  methscan)
    submit_job "$SCRIPT_DIR/run_methscan_syn22257780_v1.sh" all
    submit_job "$SCRIPT_DIR/run_methscan_syn22257780_v1.sh" IDHwt
    submit_job "$SCRIPT_DIR/run_methscan_syn22257780_v1.sh" IDHmut
    ;;
  disorder)
    submit_job "$SCRIPT_DIR/run_disorder_syn22257780_v1.sh" all
    submit_job "$SCRIPT_DIR/run_disorder_syn22257780_v1.sh" IDHwt
    submit_job "$SCRIPT_DIR/run_disorder_syn22257780_v1.sh" IDHmut
    ;;
  scna)
    if [[ "$ENABLE_SCNA" != "1" ]]; then
      echo "SCNA está desabilitado na config (ENABLE_SCNA=$ENABLE_SCNA)."
      echo "Sem BAMs, este módulo não deve ser rodado. Ajuste ENABLE_SCNA=1 apenas se você obtiver BAM_MANIFEST funcional."
      exit 1
    fi
    submit_job "$SCRIPT_DIR/run_scna_ginkgo_prep_syn22257780_v1.sh" all
    ;;
  all)
    jid_all=$(submit_job "$SCRIPT_DIR/run_methscan_syn22257780_v1.sh" all)
    jid_wt=$(submit_job "$SCRIPT_DIR/run_methscan_syn22257780_v1.sh" IDHwt)
    jid_mut=$(submit_job "$SCRIPT_DIR/run_methscan_syn22257780_v1.sh" IDHmut)

    d_all=$(submit_dependent_job "$jid_all" "$SCRIPT_DIR/run_disorder_syn22257780_v1.sh" all)
    d_wt=$(submit_dependent_job "$jid_wt" "$SCRIPT_DIR/run_disorder_syn22257780_v1.sh" IDHwt)
    d_mut=$(submit_dependent_job "$jid_mut" "$SCRIPT_DIR/run_disorder_syn22257780_v1.sh" IDHmut)

    echo "MethSCAn jobs: all=$jid_all IDHwt=$jid_wt IDHmut=$jid_mut"
    echo "Disorder jobs: all=$d_all IDHwt=$d_wt IDHmut=$d_mut"
    echo "SCNA não foi submetido neste modo seguro sem BAM."
    ;;
  all_with_scna)
    if [[ "$ENABLE_SCNA" != "1" ]]; then
      echo "SCNA está desabilitado na config (ENABLE_SCNA=$ENABLE_SCNA)."
      echo "Ative ENABLE_SCNA=1 somente se você tiver BAM_MANIFEST e BAMs válidos."
      exit 1
    fi

    jid_all=$(submit_job "$SCRIPT_DIR/run_methscan_syn22257780_v1.sh" all)
    jid_wt=$(submit_job "$SCRIPT_DIR/run_methscan_syn22257780_v1.sh" IDHwt)
    jid_mut=$(submit_job "$SCRIPT_DIR/run_methscan_syn22257780_v1.sh" IDHmut)

    d_all=$(submit_dependent_job "$jid_all" "$SCRIPT_DIR/run_disorder_syn22257780_v1.sh" all)
    d_wt=$(submit_dependent_job "$jid_wt" "$SCRIPT_DIR/run_disorder_syn22257780_v1.sh" IDHwt)
    d_mut=$(submit_dependent_job "$jid_mut" "$SCRIPT_DIR/run_disorder_syn22257780_v1.sh" IDHmut)

    s_all=$(submit_job "$SCRIPT_DIR/run_scna_ginkgo_prep_syn22257780_v1.sh" all)

    echo "MethSCAn jobs: all=$jid_all IDHwt=$jid_wt IDHmut=$jid_mut"
    echo "Disorder jobs: all=$d_all IDHwt=$d_wt IDHmut=$d_mut"
    echo "SCNA job: all=$s_all"
    ;;
  *)
    echo "Uso: bash run_master_scDNAme_multimodal_v1.sh [methscan|disorder|scna|all|all_with_scna]"
    exit 1
    ;;
esac
