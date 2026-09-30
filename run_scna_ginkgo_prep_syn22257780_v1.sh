#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH --mail-type=END,FAIL

set -euo pipefail

MODE="${1:-all}"
if [[ "$MODE" != "all" && "$MODE" != "IDHwt" && "$MODE" != "IDHmut" ]]; then
  echo "Uso: sbatch run_scna_ginkgo_prep_syn22257780_v1.sh [all|IDHwt|IDHmut]"
  exit 1
fi

if [[ -n "${SLURM_SUBMIT_DIR:-}" && -f "${SLURM_SUBMIT_DIR}/config_scDNAme_multimodal_v1.sh" ]]; then
  SCRIPT_DIR="$SLURM_SUBMIT_DIR"
else
  SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi

source "$SCRIPT_DIR/config_scDNAme_multimodal_v1.sh"

if [[ "$ENABLE_SCNA" != "1" ]]; then
  echo "SCNA/Ginkgo prep está desabilitado na config (ENABLE_SCNA=$ENABLE_SCNA)."
  echo "Este pacote foi ajustado para o cenário sem BAM. Ative ENABLE_SCNA=1 apenas se você tiver BAM_MANIFEST válido."
  exit 1
fi

OUTDIR="$PROJECT_OUT_BASE/scna_${MODE}"
PREP_DIR="$OUTDIR/prep"
BED_DIR="$OUTDIR/ginkgo_bed"
mkdir -p "$PREP_DIR" "$BED_DIR"

if [[ ! -f "$BAM_MANIFEST" ]]; then
  echo "BAM manifest não encontrado: $BAM_MANIFEST"
  echo "Preencha o arquivo template_bam_manifest_scRRBS.tsv e ajuste BAM_MANIFEST na config."
  exit 1
fi

if [[ -f "$CONDA_SH" ]]; then
  source "$CONDA_SH"
else
  echo "Conda profile não encontrado em: $CONDA_SH"
  exit 1
fi
echo "Ativando ambiente conda único: $CONDA_ENV"
conda activate "$CONDA_ENV"

python "$SCRIPT_DIR/script_prepare_ginkgo_inputs_syn22257780_v1.py" \
  --bam-manifest "$BAM_MANIFEST" \
  --qc "$QC_FILE" \
  --clinical "$CLINICAL_FILE" \
  --mode "$MODE" \
  --outdir "$PREP_DIR" \
  --min-unique-cpg "$MIN_UNIQUE_CPG" \
  --min-bs "$MIN_BS_CONVERSION" \
  --require-tumor "$REQUIRE_TUMOR_STATUS"

if ! command -v bedtools >/dev/null 2>&1; then
  echo "bedtools não encontrado no PATH."
  exit 1
fi

if ! command -v samtools >/dev/null 2>&1; then
  echo "samtools não encontrado no PATH."
  exit 1
fi

echo "Gerando arquivos .bed por célula para Ginkgo..."
{
  read -r header
  while IFS=$'\t' read -r cell bam case idh exists; do
    [[ -z "$cell" ]] && continue
    if [[ "$exists" != "True" && "$exists" != "TRUE" && "$exists" != "1" ]]; then
      echo "Pulando $cell: BAM ausente em $bam"
      continue
    fi
    out_bed="$BED_DIR/${cell}.bed"
    samtools view -b -F 260 "$bam" | bedtools bamtobed -i stdin | awk 'BEGIN{OFS="\t"} {print $1,$2,$3}' > "$out_bed"
  done
} < "$PREP_DIR/bam_manifest_filtered.tsv"

cat > "$OUTDIR/README_SCNA_NEXT_STEPS.txt" <<TXT
SCNA/Ginkgo prep concluído.

Arquivos BED por célula:
  $BED_DIR

Manifest filtrado:
  $PREP_DIR/bam_manifest_filtered.tsv

Observações:
1) O Ginkgo fornece CLI e instalação local, mas o repositório/documentação pública lista explicitamente dados de binning prontos para hg19 e outros genomas antigos; hg38 não aparece nessa lista pública.
2) Se seus BAMs estão em hg38, valide primeiro se você tem bins customizados compatíveis no Ginkgo local; caso contrário, pare aqui e não rode a chamada automática.
3) Se você tiver uma instalação local funcional do Ginkgo e binning compatível, rode algo no estilo:
   $GINKGO_DIR/cli/ginkgo.sh --input $BED_DIR --genome $GINKGO_GENOME --binning $GINKGO_BINNING
TXT

if [[ "$GINKGO_RUN_CLI" == "1" ]]; then
  if [[ ! -x "$GINKGO_DIR/cli/ginkgo.sh" ]]; then
    echo "Ginkgo CLI não encontrado em $GINKGO_DIR/cli/ginkgo.sh"
    exit 1
  fi
  if [[ "$GINKGO_GENOME" != "hg19" ]]; then
    echo "GINKGO_RUN_CLI=1, mas GINKGO_GENOME=$GINKGO_GENOME. A documentação pública do Ginkgo lista bins prontos para hg19 e outros genomas antigos, não hg38."
    echo "Abortando execução automática; mantendo apenas a preparação dos BEDs."
    exit 0
  fi

  echo "Rodando Ginkgo CLI..."
  "$GINKGO_DIR/cli/ginkgo.sh" \
    --input "$BED_DIR" \
    --genome "$GINKGO_GENOME" \
    --binning "$GINKGO_BINNING"
fi
