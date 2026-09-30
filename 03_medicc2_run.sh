#!/bin/bash
#SBATCH --job-name=medicc2_final
#SBATCH --output=/home/renanomete/projetos/matdata/evo_clonal/logs/medicc2_final_%j.out
#SBATCH --error=/home/renanomete/projetos/matdata/evo_clonal/logs/medicc2_final_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00

set -o pipefail
set -e

source ~/.bashrc
conda activate evo_clonal_medicc2 || { echo "[ERRO] Falha ao ativar ambiente"; exit 1; }
python -c "import joblib" 2>/dev/null || { echo "[ERRO] Instale joblib: pip install joblib"; exit 1; }

BASE_DIR="/home/renanomete/projetos/matdata/evo_clonal"
INPUT_DIR="${BASE_DIR}/data/medicc2_input"
RESULTS_DIR="${BASE_DIR}/results/medicc2"
LOG_DIR="${BASE_DIR}/logs"
DIAG_DIR="${RESULTS_DIR}/diagnostics"
mkdir -p "${RESULTS_DIR}" "${LOG_DIR}" "${DIAG_DIR}"
cd "${BASE_DIR}" || exit 1

echo "======================================================"
echo " MEDICC2 iniciado: $(date)"
echo "======================================================"
echo "[INFO] Host: $(hostname)"
echo "[INFO] Python: $(which python)"
echo "[INFO] medicc2: $(which medicc2)"
echo "[INFO] Input dir: ${INPUT_DIR}"
ls -lh "${INPUT_DIR}"
echo ""

# =============================================================================
# Função de diagnóstico detalhado
# =============================================================================
diagnose_input() {
    local INPUT_FILE="$1"
    local LABEL="$2"
    local DIAG_FILE="${DIAG_DIR}/${LABEL}_diagnosis.txt"

    {
        echo "=============================================="
        echo " DIAGNÓSTICO DO ARQUIVO: ${LABEL}"
        echo " Arquivo: ${INPUT_FILE}"
        echo " Data: $(date)"
        echo "=============================================="
        echo ""

        if [[ ! -f "${INPUT_FILE}" ]]; then
            echo "[ERRO] Arquivo não encontrado."
            return 1
        fi

        local total_lines=$(wc -l < "${INPUT_FILE}")
        local total_cols=$(head -1 "${INPUT_FILE}" | awk -F'\t' '{print NF}')
        echo "Total de linhas: ${total_lines}"
        echo "Total de colunas: ${total_cols}"
        echo ""

        echo "--- Cabeçalho completo (com índices) ---"
        head -1 "${INPUT_FILE}" | tr '\t' '\n' | cat -n
        echo ""

        echo "--- Primeiras 5 linhas ---"
        head -5 "${INPUT_FILE}" | cat -A
        echo ""

        echo "--- Últimas 5 linhas ---"
        tail -5 "${INPUT_FILE}" | cat -A
        echo ""

        # Análise da 5ª coluna
        echo "--- Análise da 5ª coluna (valores de cópia/alelo) ---"
        local col5_values=$(mktemp)
        tail -n +2 "${INPUT_FILE}" | cut -f5 | sort | uniq -c | sort -rn > "${col5_values}"
        local unique_count=$(wc -l < "${col5_values}")

        echo "Número de valores únicos: ${unique_count}"
        echo "Top 20 valores mais frequentes:"
        head -20 "${col5_values}"
        echo ""

        echo "--- Verificação de integridade dos valores ---"
        local empty_count=$(tail -n +2 "${INPUT_FILE}" | awk -F'\t' '$5 == "" {print}' | wc -l)
        local nonint_count=$(tail -n +2 "${INPUT_FILE}" | awk -F'\t' '$5 !~ /^[0-9]+$/ {print}' | wc -l)
        local pipe_count=$(tail -n +2 "${INPUT_FILE}" | awk -F'\t' '$5 ~ /\|/ {print}' | wc -l)
        local slash_count=$(tail -n +2 "${INPUT_FILE}" | awk -F'\t' '$5 ~ /\// {print}' | wc -l)
        local comma_count=$(tail -n +2 "${INPUT_FILE}" | awk -F'\t' '$5 ~ /,/ {print}' | wc -l)

        echo "Linhas com valor vazio: ${empty_count}"
        echo "Linhas com valor não inteiro (incluindo vazios): ${nonint_count}"
        echo "Linhas contendo '|': ${pipe_count}"
        echo "Linhas contendo '/': ${slash_count}"
        echo "Linhas contendo ',': ${comma_count}"
        echo ""

        if [[ ${nonint_count} -gt 0 ]]; then
            echo "Exemplos de linhas com valores problemáticos (até 5):"
            tail -n +2 "${INPUT_FILE}" | awk -F'\t' '$5 !~ /^[0-9]+$/' | head -5 | cat -A
            echo ""
        fi

        rm -f "${col5_values}"
    } > "${DIAG_FILE}"

    echo "[DIAG] Relatório salvo em: ${DIAG_FILE}"
    grep -E "Total de linhas|Total de colunas" "${DIAG_FILE}"
}

# =============================================================================
# Função para preparar dados no formato de alelos (cn_a, cn_b)
# =============================================================================
prepare_allele_file() {
    local INPUT_FILE="$1"
    local OUT_FILE="$2"
    local LABEL="$3"

    local header=$(head -1 "${INPUT_FILE}")

    if echo "${header}" | grep -q "cn_b"; then
        # Já tem cn_a e cn_b
        cp "${INPUT_FILE}" "${OUT_FILE}"
        echo "cn_a cn_b"
    elif echo "${header}" | grep -q "cn_a"; then
        # Apenas cn_a
        cp "${INPUT_FILE}" "${OUT_FILE}"
        echo "cn_a"
    else
        # Formato copy_number: precisa converter
        # Detecta separador
        local separator="|"
        if tail -n +2 "${INPUT_FILE}" | cut -f5 | grep -q '/'; then
            separator="/"
        elif tail -n +2 "${INPUT_FILE}" | cut -f5 | grep -q ','; then
            separator=","
        fi
        echo "[INFO] Preparando ${LABEL}: convertendo copy_number -> cn_a,cn_b (separador '${separator}')"

        # Verifica se há linhas com separador
        local has_sep=$(tail -n +2 "${INPUT_FILE}" | cut -f5 | grep -E '[|/,]' | head -1)
        if [[ -n "${has_sep}" ]]; then
            awk -F'\t' -v sep="${separator}" 'BEGIN {OFS=FS}
                NR==1 {print "sample_id","chrom","start","end","cn_a","cn_b"; next}
                {
                    split($5, a, sep)
                    if (length(a) == 2) {
                        print $1, $2, $3, $4, a[1], a[2]
                    } else {
                        print $1, $2, $3, $4, $5, $5
                    }
                }' "${INPUT_FILE}" > "${OUT_FILE}"
        else
            # Valores inteiros: duplica coluna
            echo "[INFO] Valores inteiros detectados: duplicando coluna para cn_a e cn_b."
            awk -F'\t' 'BEGIN {OFS=FS}
                NR==1 {print "sample_id","chrom","start","end","cn_a","cn_b"; next}
                {print $1, $2, $3, $4, $5, $5}' "${INPUT_FILE}" > "${OUT_FILE}"
        fi
        echo "cn_a cn_b"
    fi
}

# =============================================================================
# Função principal de execução
# =============================================================================
run_medicc() {
    local INPUT_FILE="$1"
    local OUT_DIR="$2"
    local LABEL="$3"

    mkdir -p "${OUT_DIR}"
    echo "------------------------------------------------------"
    echo "[INFO] Iniciando run: ${LABEL}"
    echo "[INFO] Input original: ${INPUT_FILE}"
    echo "[INFO] Output: ${OUT_DIR}"
    echo "[INFO] Data/hora: $(date)"
    echo "------------------------------------------------------"

    if [[ ! -f "${INPUT_FILE}" ]]; then
        echo "[ERRO] Arquivo não encontrado: ${INPUT_FILE}"
        return 1
    fi

    # Diagnóstico do arquivo original
    diagnose_input "${INPUT_FILE}" "${LABEL}"

    # Preparar arquivo de alelos
    local ALLELE_FILE="${OUT_DIR}/${LABEL}_alleles.tsv"
    local ALLELE_COLS=$(prepare_allele_file "${INPUT_FILE}" "${ALLELE_FILE}" "${LABEL}")

    echo "[INFO] Colunas de alelo usadas: ${ALLELE_COLS}"
    echo "[INFO] Primeiras linhas do arquivo preparado:"
    head -5 "${ALLELE_FILE}" | cat -A
    echo ""

    # Executar MEDICC2
    local RUN_LOG="${OUT_DIR}/${LABEL}_medicc2.log"
    local ALLELE_COLS_COMMA="${ALLELE_COLS// /,}"

    medicc2 \
        "${ALLELE_FILE}" \
        "${OUT_DIR}" \
        --input-type tsv \
        --input-allele-columns "${ALLELE_COLS_COMMA}" \
        --normal-name diploid_normal \
        -j "${SLURM_CPUS_PER_TASK}" \
        --verbose \
        > "${RUN_LOG}" 2>&1

    local STATUS=$?
    if [[ ${STATUS} -ne 0 ]]; then
        echo "[ERRO] MEDICC2 falhou para ${LABEL} (exit code ${STATUS})"
        echo "[INFO] Últimas linhas do log:"
        tail -n 40 "${RUN_LOG}"
        return 1
    fi

    echo "[INFO] MEDICC2 finalizado para ${LABEL}"
    echo "[INFO] Log: ${RUN_LOG}"

    local TREE_FILE
    TREE_FILE=$(find "${OUT_DIR}" -type f \( -name "*.new" -o -name "*.nwk" -o -name "*.tree" \) 2>/dev/null | head -1)
    if [[ -n "${TREE_FILE}" ]]; then
        echo "[OK] Árvore encontrada: ${TREE_FILE}"
    else
        echo "[AVISO] Run terminou, mas não encontrei arquivo de árvore em ${OUT_DIR}"
        find "${OUT_DIR}" -maxdepth 2 -type f | sort
    fi
    echo ""
    return 0
}

# =============================================================================
# Execução principal
# =============================================================================
OVERALL_OK=1

# 1. GLOBAL
echo "[INFO] ========== Processando dataset global =========="
GLOBAL_INPUT="${INPUT_DIR}/all_cells_medicc2_clean.tsv"
if [[ -f "${GLOBAL_INPUT}" ]]; then
    run_medicc "${GLOBAL_INPUT}" "${RESULTS_DIR}/global" "global" || OVERALL_OK=0
else
    echo "[AVISO] Arquivo global não encontrado: ${GLOBAL_INPUT}"
    OVERALL_OK=0
fi

# 2. JK142
echo "[INFO] ========== Processando dataset JK142 =========="
JK142_INPUT="${INPUT_DIR}/JK142_medicc2.tsv"
if [[ -f "${JK142_INPUT}" ]]; then
    run_medicc "${JK142_INPUT}" "${RESULTS_DIR}/JK142" "JK142" || OVERALL_OK=0
else
    echo "[AVISO] Arquivo JK142 não encontrado"
    OVERALL_OK=0
fi

# 3. JK136 e JK153
for PAT in JK136 JK153; do
    echo "[INFO] ========== Processando dataset ${PAT} =========="
    INPUT_PAT="${INPUT_DIR}/${PAT}_medicc2.tsv"
    if [[ -f "${INPUT_PAT}" ]]; then
        run_medicc "${INPUT_PAT}" "${RESULTS_DIR}/${PAT}" "${PAT}" || OVERALL_OK=0
    else
        echo "[AVISO] Arquivo ${PAT} não encontrado"
        OVERALL_OK=0
    fi
done

# =============================================================================
# Relatório final
# =============================================================================
echo ""
echo "======================================================"
echo " MEDICC2 finalizado: $(date)"
echo "======================================================"
echo ""
echo "Árvores geradas:"
for DIR in "${RESULTS_DIR}/global" "${RESULTS_DIR}/JK136" "${RESULTS_DIR}/JK142" "${RESULTS_DIR}/JK153"; do
    if [[ -d "${DIR}" ]]; then
        TREE_FILE=$(find "${DIR}" -type f \( -name "*.new" -o -name "*.nwk" -o -name "*.tree" \) 2>/dev/null | head -1)
        if [[ -n "${TREE_FILE}" ]]; then
            echo "  [OK] $(basename "${DIR}"): ${TREE_FILE} ($(du -sh "${TREE_FILE}" | cut -f1))"
        else
            echo "  [!!] $(basename "${DIR}"): árvore não encontrada"
        fi
    else
        echo "  [!!] $(basename "${DIR}"): diretório não existe"
    fi
done

echo ""
echo "Relatórios de diagnóstico salvos em: ${DIAG_DIR}"
if [[ "${OVERALL_OK}" -eq 1 ]]; then
    echo "Todos os runs concluíram com outputs válidos."
else
    echo "[AVISO] Um ou mais runs falharam — verifique os logs."
    exit 1
fi
