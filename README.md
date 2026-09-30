# scripts-doutorado

Repositório de scripts e pipelines do doutorado, organizado em três branches
órfãs independentes — cada uma representa um projeto distinto.

## 🌿 Estrutura das branches

| Branch | Projeto | Conteúdo |
|---|---|---|
| **`scRNA`** | Análise de single-cell RNA-seq | `gse182109/` (pipeline GSE182109) + `integration/` (integração GSE182109 + syn22257780) |
| **`scDNAme`** | Análise de single-cell DNA methylation | Pipelines de methscan, methylVI, pdclust e ginkgo para syn22257780 |
| **`scWGS`** | Análise de single-cell whole-genome sequencing | Pipelines de CNV, MEDICC2, árvores filogenéticas e visualização |

## 📂 Estrutura interna

### `scRNA`
gse182109/ # Pipeline GSE182109
integration/ # Pipeline de integração com syn22257780

### `scDNAme`
Arquivos na raiz cobrindo:
- `config_*`, `run_*`, `script_*` → pipelines de methscan, methylVI, pdclust
- `README_*` → documentação de cada pipeline
- `DEPENDENCIAS_*` → dependências de ambiente

### `scWGS`
Arquivos na raiz cobrindo:
- `00_setup_conda.sh` a `05_visualization.*` → pipeline sequencial
- `scwgs_analysis.R/.sh` → análise principal
- `functions.R` → funções auxiliares
- `submit_*.sh` → scripts de submissão

## Uso rápido
```bash
# Clonar o repositório
git clone https://github.com/<usuario>/scripts-doutorado.git
cd scripts-doutorado
# Escolher o projeto desejado
git checkout scRNA      # ou scDNAme ou scWGS
