# Pipeline scWGS + MEDICC2 — Filogenia clonal por CNV (GSE173279)

Este README documenta os scripts do pipeline **scWGS + MEDICC2** aplicado ao dataset **GSE173279** (evolução clonal tumoral em glioma) executado no GridUNESP.

O pipeline cobre todo o percurso:

- Download dos dados públicos do GEO.
- QC das células scWGS e preparação dos inputs do MEDICC2.
- Inferência da filogenia clonal via MEDICC2.
- Análise da árvore (métricas por célula, clones, validação).
- Visualizações finais (árvores, UMAP, heatmaps, composição clonal).
- Análise pós-MEDICC2 (`scwgs_analysis.R`) por paciente, com subclonagem e CN de genes de interesse.

> ⚠️ Caminhos absolutos estão hardcoded nos scripts. Ajuste conforme o ambiente.
>
> ⚠️ Existem **dois `WORKDIR` diferentes** nos scripts: o pipeline principal usa `evo_clonal`, enquanto `scwgs_analysis.R`/`submit_scwgs_analysis.sh` usam `new_evo_anal/scWGS/GSE173279/results/evo_clonal`. Unifique antes de rodar em produção.

---

## 1. Visão geral do fluxo

```
00_setup_conda.sh              (executar UMA VEZ em sessão interativa)

submit_all.sh                  (orquestrador SLURM)
      │
      ├─► 01_download.sh                    baixa arquivos brutos do GEO
      │
      ├─► 02_cnv_qc_prepare.sh              QC scWGS + inputs MEDICC2
      │        └── 02_cnv_qc_prepare.R
      │
      ├─► 03_medicc2_run.sh                 roda MEDICC2 (global + JK136 + JK142 + JK153)
      │
      ├─► 04_tree_analysis.sh               análise da árvore e clones
      │        └── 04_tree_analysis.R
      │
      └─► 05_visualization.sh               figuras finais
               └── 05_visualization.R

00_setup_conda.sh (opcional, diagnóstico)
      └── 02_diagnose.R                     inspeciona estrutura dos brutos

submit_scwgs_analysis.sh / scwgs_analysis.sh   (job array alternativo, por paciente)
      └── scwgs_analysis.R                  análise pós-MEDICC2 (JK136, JK153, global)
```

---

## 2. Estrutura de diretórios

| Diretório | Conteúdo |
|---|---|
| `scripts/` | Scripts R e wrappers SLURM |
| `data/raw/` | Arquivos brutos baixados do GEO |
| `data/processed/` | Objetos R intermediários (`cnv_qc_object.rds`, `tree_analysis.rds`) |
| `data/medicc2_input/` | TSVs de entrada do MEDICC2 |
| `results/qc/` | Métricas e plots de QC |
| `results/medicc2/` | Árvores e distâncias MEDICC2 (`global`, `JK136`, `JK142`, `JK153`) |
| `results/medicc2/diagnostics/` | Relatórios de diagnóstico dos inputs |
| `results/tree_analysis/` | Métricas por célula, clones, pacientes |
| `results/figures/` | Figuras finais |
| `logs/` | Logs SLURM |

Base alternativa usada pelo `scwgs_analysis.R`:

| Diretório | Conteúdo |
|---|---|
| `${BASE_DIR}/scripts/` | Scripts da análise pós-MEDICC2 |
| `${BASE_DIR}/results/medicc2/<PATIENT>/` | Saídas do MEDICC2 por paciente |
| `${BASE_DIR}/results/scwgs_analysis/<PATIENT>/` | Saídas da análise pós-MEDICC2 |
| `${BASE_DIR}/logs/` | Logs por paciente |

---

## 3. Dependências

### Sistema

| Recurso | Observação |
|---|---|
| SLURM | Gerenciador de jobs (`sbatch`, `squeue`) |
| `module load miniconda/24.4.0-libmamba` | Módulo do cluster |
| `conda` | Ativação do ambiente |
| `wget` | Download dos dados |
| `Rscript` | Scripts R |
| `python` | Verificação do MEDICC2 |
| `job-nanny` | Opcional; usado nos wrappers 04 e 05 |

### Ambiente conda `evo_clonal_medicc2`

Criado por `00_setup_conda.sh`.

| Grupo | Pacotes |
|---|---|
| Base | `python=3.10`, `r-base=4.3`, `r-essentials` |
| MEDICC2 | `medicc2` (conda-forge) |
| Python | `pandas`, `numpy`, `scipy`, `matplotlib`, `seaborn`, `ete3` |
| R (core) | `data.table`, `Matrix`, `ggplot2`, `patchwork`, `dplyr`, `tidyr`, `viridis`, `RColorBrewer`, `scales`, `irlba`, `uwot`, `rann`, `igraph`, `harmony`, `future`, `future.apply` |
| R (filogenia/genômica) | `ape`, `phangorn`, `ggtree`, `GenomicRanges`, `IRanges` |
| R (CRAN extra) | `ggtreeExtra`, `aplot`, `tidytree`, `treeio` |
| R (análise pós) | `pheatmap`, `viridisLite`, `umap` |

---

## 4. Entradas esperadas

### Dados brutos (baixados por `01_download.sh`)

| Arquivo | Descrição |
|---|---|
| `data/raw/GSE173279_scWGS_all_cells_cn.tsv.gz` | CN de todas as células |
| `data/raw/GSE173279_scWGS_all_cells_metadata.csv.gz` | Metadados de todas as células |
| `data/raw/GSE173279_scWGS_JK136_coarse_cn.tsv.gz` | CN coarse JK136 |
| `data/raw/GSE173279_scWGS_JK136_metadata.csv.gz` | Metadados JK136 |
| `data/raw/GSE173279_scWGS_JK142_coarse_cn.tsv.gz` | CN coarse JK142 |
| `data/raw/GSE173279_scWGS_JK142_metadata.csv.gz` | Metadados JK142 |
| `data/raw/GSE173279_scWGS_JK153_coarse_cn.tsv.gz` | CN coarse JK153 |
| `data/raw/GSE173279_scWGS_JK153_metadata.csv.gz` | Metadados JK153 |

Base URL: `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE173nnn/GSE173279/suppl`

### Arquivos intermediários consumidos por etapas posteriores

| Arquivo | Consumido por |
|---|---|
| `data/processed/cnv_qc_object.rds` | `04_tree_analysis.R`, `05_visualization.R` |
| `data/medicc2_input/all_cells_medicc2.tsv` | `03_medicc2_run.sh` (global) |
| `data/medicc2_input/JK136_medicc2.tsv` | `03_medicc2_run.sh` |
| `data/medicc2_input/JK142_medicc2.tsv` | `03_medicc2_run.sh` |
| `data/medicc2_input/JK153_medicc2.tsv` | `03_medicc2_run.sh` |
| `data/processed/tree_analysis.rds` | `05_visualization.R` |

### Arquivos MEDICC2 por paciente (consumidos por `scwgs_analysis.R`)

| Arquivo | Padrão |
|---|---|
| Árvore final | `*_final_tree.new` |
| Branch lengths | `*_branch_lengths.tsv` |
| CN profiles | `*_final_cn_profiles.tsv` |
| Distâncias pairwise | `*_pairwise_distances.tsv` |
| Sumário | `*_summary.tsv` |

---

## 5. Saídas principais

### QC (etapa 2)

| Arquivo | Conteúdo |
|---|---|
| `data/processed/cnv_qc_object.rds` | Objeto R com metadados + QC |
| `results/qc/qc_metrics.csv` | Métricas por célula |
| `results/qc/qc_plots.pdf` | Histogramas e scatter de QC |
| `data/medicc2_input/all_cells_medicc2.tsv` | Input global |
| `data/medicc2_input/JK136_medicc2.tsv` | Input JK136 |
| `data/medicc2_input/JK142_medicc2.tsv` | Input JK142 |
| `data/medicc2_input/JK153_medicc2.tsv` | Input JK153 |

### MEDICC2 (etapa 3)

| Arquivo | Conteúdo |
|---|---|
| `results/medicc2/<RUN>/*.new` | Árvore filogenética |
| `results/medicc2/<RUN>/*_alleles.tsv` | Arquivo de alelos preparado |
| `results/medicc2/<RUN>/*_medicc2.log` | Log da execução |
| `results/medicc2/diagnostics/<RUN>_diagnosis.txt` | Diagnóstico do input |

### Análise da árvore (etapa 4)

| Arquivo | Conteúdo |
|---|---|
| `data/processed/tree_analysis.rds` | Objeto consolidado |
| `results/tree_analysis/evo_metrics_per_cell.csv` | Distância à raiz + clone por célula |
| `results/tree_analysis/clone_summary.csv` | Resumo por clone |
| `results/tree_analysis/patient_summary.csv` | Resumo por paciente |

### Visualizações (etapa 5)

| Arquivo | Conteúdo |
|---|---|
| `results/figures/01_phylo_tree_global.pdf` | Árvore global colorida por paciente/clone |
| `results/figures/02_phylo_trees_by_patient.pdf` | Árvores por paciente |
| `results/figures/03_umap_panels.pdf` | UMAP por paciente/distância/clone |
| `results/figures/04_cn_heatmap_by_clone.pdf` | Heatmap CN × clone × cromossomo |
| `results/figures/05_evo_dist_boxplots.pdf` | Boxplot de distância evolutiva |
| `results/figures/06_clonal_composition.pdf` | Composição clonal por paciente |
| `results/figures/07_validation_<run>.pdf` | Distância evolutiva × aberração CNV |
| `results/SUMMARY_REPORT.txt` | Relatório final consolidado |

### Análise pós-MEDICC2 (`scwgs_analysis.R`)

| Arquivo | Conteúdo |
|---|---|
| `<OUT>/<PATIENT>_subclones.tsv` | Atribuição de subclones por célula |
| `<OUT>/<PATIENT>_gene_cn_per_cell.tsv` | CN dos genes de interesse por célula |
| `<OUT>/<PATIENT>_gene_cn_per_subclone.tsv` | CN médio dos genes por subclone |
| `<OUT>/<PATIENT>_subclone_cn_profiles.tsv` | Matriz subclone × bin |
| `<OUT>/<PATIENT>_umap_cna.pdf` | UMAP sobre a matriz de distâncias MEDICC2 |
| `<OUT>/<PATIENT>_umap_coords.tsv` | Coordenadas UMAP |
| `<OUT>/<PATIENT>_subclone_tree.pdf` / `.new` | Árvore dos subclones |
| `<OUT>/<PATIENT>_gene_cn_heatmap.pdf` | Heatmap dos genes por subclone |
| `logs/scwgs_<PATIENT>_<timestamp>.log` | Log da execução |

---

## 6. Como executar

### Passo 0 — Setup do ambiente (uma vez)

Em sessão interativa:

```bash
bash 00_setup_conda.sh
```

### Passo 1 — Pipeline completo

```bash
bash submit_all.sh
```

### Passo 2 — Retomar de uma etapa

```bash
bash submit_all.sh --from 3    # retoma a partir da etapa 3
```

### Passo 3 — Diagnóstico opcional antes do QC

Se estiver em dúvida sobre a orientação dos arquivos brutos:

```bash
Rscript scripts/02_diagnose.R
```

### Alternativa — análise pós-MEDICC2 por paciente

```bash
sbatch submit_scwgs_analysis.sh         # job array 0-2: JK136, JK153, global
```

### Alternativa — rodar scripts individuais

```bash
sbatch scripts/01_download.sh
sbatch scripts/02_cnv_qc_prepare.sh
sbatch scripts/03_medicc2_run.sh
sbatch scripts/04_tree_analysis.sh
sbatch scripts/05_visualization.sh
```

---

## 7. Descrição dos scripts

### 7.1 `00_setup_conda.sh`

Script de setup do ambiente `evo_clonal_medicc2`. Executado uma única vez em sessão interativa.

| Etapa | Descrição |
|---|---|
| 1 | Cria ambiente com Python 3.10 + R 4.3 |
| 2 | Instala `medicc2` via conda-forge/bioconda |
| 3 | Instala dependências Python (pandas, numpy, scipy, matplotlib, seaborn, ete3) |
| 4 | Instala pacotes R base |
| 5 | Instala pacotes R de filogenia (`ape`, `phangorn`, `ggtree`, `GenomicRanges`, `IRanges`) |
| 6 | Instala pacotes R via CRAN (`ggtreeExtra`, `aplot`, `tidytree`, `treeio`) |
| 7 | Verifica instalações (Python, MEDICC2, R) |

---

### 7.2 `submit_all.sh`

Orquestrador SLURM do pipeline MEDICC2.

| Etapa | Script | Dependência |
|---|---|---|
| 1 | `scripts/01_download.sh` | — |
| 2 | `scripts/02_cnv_qc_prepare.sh` | `afterok` da etapa 1 |
| 3 | `scripts/03_medicc2_run.sh` | `afterok` da etapa 2 |
| 4 | `scripts/04_tree_analysis.sh` | `afterok` da etapa 3 |
| 5 | `scripts/05_visualization.sh` | `afterok` da etapa 4 |

| Parâmetro | Efeito |
|---|---|
| `--from N` | Começa a partir da etapa N (1 a 5) |
| `--from 1` (default) | Roda pipeline inteiro |

---

### 7.3 `01_download.sh`

Download dos arquivos públicos do GSE173279.

| Característica | Descrição |
|---|---|
| Fonte | `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE173nnn/GSE173279/suppl` |
| Recursos | 2 CPUs, 4 GB, 1h30 |
| Idempotente | Sim — pula arquivos já presentes com `> 1000 bytes` |
| Validação | Falha se qualquer arquivo ficar vazio ao final |

---

### 7.4 `02_diagnose.R`

Script opcional de diagnóstico. Lê as primeiras linhas dos arquivos `_coarse_cn.tsv.gz` e `_metadata.csv.gz` para inferir a orientação da matriz (células × bins vs bins × células) e a estrutura de colunas.

| Saída | Descrição |
|---|---|
| stdout | Número de linhas/colunas, exemplos de valores, heurística de orientação |

---

### 7.5 `02_cnv_qc_prepare.R`

QC das células e preparo dos inputs do MEDICC2.

| Etapa | Descrição |
|---|---|
| 1 | Lê CN coarse + metadados para JK136, JK142, JK153 |
| 2 | Detecta orientação da matriz e associa barcodes |
| 3 | Calcula métricas de QC: aberração CNV, Gini, ploidia modal, cobertura |
| 4 | Aplica filtros MAD (n=3) para `pct_measured`, `gini`, `aberration_score` |
| 5 | Gera `qc_metrics.csv` e `qc_plots.pdf` |
| 6 | Escreve TSVs de input MEDICC2 (normal diploide + células QC-pass) |
| 7 | Identifica bins comuns entre pacientes e gera input global |
| 8 | Salva `cnv_qc_object.rds` |

Funções auxiliares:

| Função | Descrição |
|---|---|
| `gini_coefficient` | Coeficiente de Gini |
| `estimate_ploidy_modal` | Ploidia modal |
| `cnv_aberration_score` | Soma de `|CN − 2|` |
| `is_outlier_high` / `is_outlier_low` | Outliers via MAD |
| `parse_bin_coords` | Parse de `chr1_1040001_2080000` |
| `read_cn_file` | Detecção de orientação + associação com metadados |
| `write_medicc2_tsv` | Escreve TSV no formato esperado pelo MEDICC2 |

---

### 7.6 `02_cnv_qc_prepare.sh`

Wrapper SLURM da etapa 2. Verifica existência dos inputs e chama `02_cnv_qc_prepare.R`.

| Parâmetro | Valor |
|---|---|
| Tempo | 4 h |
| CPUs | 8 |
| Memória | 64 GB |
| Conda env | `evo_clonal_medicc2` |

---

### 7.7 `03_medicc2_run.sh`

Executa o MEDICC2 para os quatro conjuntos.

| Run | Input |
|---|---|
| global | `data/medicc2_input/all_cells_medicc2_clean.tsv` |
| JK142 | `data/medicc2_input/JK142_medicc2.tsv` |
| JK136 | `data/medicc2_input/JK136_medicc2.tsv` |
| JK153 | `data/medicc2_input/JK153_medicc2.tsv` |

Etapas por run:

| Etapa | Descrição |
|---|---|
| Diagnóstico | `diagnose_input` analisa coluna 5 (CN ou alelo) |
| Preparo | `prepare_allele_file` converte `copy_number` em `cn_a`, `cn_b` se necessário |
| Execução | `medicc2 <input> <outdir> --input-type tsv --input-allele-columns cn_a,cn_b --normal-name diploid_normal -j N` |
| Relatório | Busca árvore (`*.new`, `*.nwk`, `*.tree`) |

| Recurso | Valor |
|---|---|
| Tempo | 24 h |
| CPUs | 8 |
| Memória | 32 GB |
| Conda env | `evo_clonal_medicc2` |
| Dependência Python | `joblib` |

---

### 7.8 `04_tree_analysis.R`

Análise da árvore MEDICC2.

| Etapa | Descrição |
|---|---|
| 1 | Carrega `cnv_qc_object.rds` + saídas MEDICC2 |
| 2 | Para cada run: calcula distância à raiz por célula |
| 3 | Clustering clonal via `hclust` na matriz de distâncias (`ward.D2`) |
| 4 | Validação: Spearman(dist. evolutiva, aberração CNV) |
| 5 | Gera `evo_metrics_per_cell.csv` |
| 6 | Resumos por clone e por paciente |
| 7 | Salva `tree_analysis.rds` |

| Função | Descrição |
|---|---|
| `load_medicc2_outputs` | Carrega árvore, distâncias e eventos |
| `compute_tree_metrics` | Distância à raiz + cluster + validação |

---

### 7.9 `04_tree_analysis.sh`

Wrapper SLURM da etapa 4.

| Parâmetro | Valor |
|---|---|
| Tempo | 4 h |
| CPUs | 4 |
| Memória | 32 GB |
| Conda env | `evo_clonal_medicc2` |
| Verificação | Confirma que pelo menos uma árvore MEDICC2 existe |

> ⚠️ O arquivo anexado está **truncado no topo** — falta o shebang e as primeiras linhas. Reconstitua o cabeçalho antes de submeter.

---

### 7.10 `05_visualization.R`

Todas as figuras do pipeline.

| Etapa | Descrição |
|---|---|
| 1 | Carrega `tree_analysis.rds` |
| 2 | Árvores filogenéticas (`ggtree` se disponível, senão `ape`) |
| 3 | UMAP sobre CN (PCA + Harmony por paciente + UMAP) |
| 4 | Heatmap CN por clone × cromossomo |
| 5 | Boxplot de distância evolutiva |
| 6 | Composição clonal por paciente |
| 7 | Validação (dist. evolutiva × aberração CNV) |
| 8 | Relatório final (`SUMMARY_REPORT.txt`) |

| Característica | Descrição |
|---|---|
| Dependência opcional | `ggtree` + `treeio` (fallback: `ape::plot.phylo`) |
| Harmony | Corrige por paciente no UMAP |
| Amostragem | Até 3000 pontos por run no plot de validação |

---

### 7.11 `05_visualization.sh`

Wrapper SLURM da etapa 5.

| Parâmetro | Valor |
|---|---|
| Tempo | 4 h |
| CPUs | 4 |
| Memória | 32 GB |
| Conda env | `evo_clonal_medicc2` |

---

### 7.12 `functions.R`

Utilitários compartilhados.

| Função | Uso |
|---|---|
| `log_step` | Log de etapa numerada |
| `log_info` | Log informativo |
| `ensure_dir` | Cria diretório se não existir |
| `gini_coefficient` | Gini |
| `cnv_aberration_score` | `sum(|CN − 2|)` |
| `is_outlier_low` / `is_outlier_high` | Outliers via MAD |
| `estimate_ploidy_modal` | Ploidia modal via KDE |
| `select_variable_bins` | Bins variáveis por MAD (exclui cromossomos sexuais) |
| `parse_bin_coords` | Parse de coordenadas de bins (2 formatos) |

---

### 7.13 `submit_scwgs_analysis.sh` e `scwgs_analysis.sh`

Wrappers SLURM do job array pós-MEDICC2.

| Parâmetro | Valor |
|---|---|
| Array | 0-2 → `JK136`, `JK153`, `global` |
| Tempo | 4 h |
| CPUs | 4 |
| Memória | 32 GB |
| Conda env | `evo_clonal_medicc2` |

> ⚠️ Os dois arquivos são idênticos. Mantenha apenas um.

---

### 7.14 `scwgs_analysis.R`

Análise pós-MEDICC2 por paciente.

| Etapa | Descrição |
|---|---|
| 1 | Localiza arquivos MEDICC2 do paciente |
| 2 | Filtra `diploid_normal` das distâncias, CN profiles e branch lengths |
| 3 | Constrói subclones via `hclust` (`ward.D2`) com `k = 3, 5, 8, 10, 15` (default `k = 5`) |
| 4 | Estatísticas de diversidade clonal (Shannon, Simpson) |
| 5 | Extrai CN de genes de interesse (14 genes: proliferação e diferenciação) |
| 6 | Calcula CN médio por subclone |
| 7 | Gera UMAP sobre a matriz de distâncias |
| 8 | Constrói árvore de subclones (`hclust` + `as.phylo`) |
| 9 | Heatmap dos genes por subclone |

Genes analisados:

| Gene | Cromossomo | Eixo |
|---|---|---|
| TYMS | chr18 | Proliferativo |
| PCLAF | chr1 | Proliferativo |
| BIRC5 | chr17 | Proliferativo |
| PBK | chr8 | Proliferativo |
| TPX2 | chr20 | Proliferativo |
| EZH2 | chr7 | Proliferativo |
| MYBL2 | chr20 | Proliferativo |
| NEK2 | chr1 | Proliferativo |
| PLP1 | chrX | Diferenciado |
| MBP | chr18 | Diferenciado |
| COL1A2 | chr7 | Diferenciado |
| COL1A1 | chr17 | Diferenciado |
| CTHRC1 | chr8 | Diferenciado |
| PCOLCE | chr7 | Diferenciado |

---

## 8. Observações e cuidados

| Item | Observação |
|---|---|
| Caminhos hardcoded | Ajustar `WORKDIR` e `BASE_DIR` em todos os scripts |
| Dois WORKDIRs | `evo_clonal` (pipeline principal) e `new_evo_anal/scWGS/GSE173279/results/evo_clonal` (análise pós). Unificar |
| `04_tree_analysis.sh` truncado | Falta o shebang `#!/bin/bash` no topo — reconstituir |
| Duplicação | `submit_scwgs_analysis.sh` e `scwgs_analysis.sh` são idênticos |
| `submit_all.sh` e etapa 3 | O orquestrador aponta para `scripts/03_medicc2_run.sh` — nome correto |
| Global input | Esperado como `all_cells_medicc2_clean.tsv` em `03_medicc2_run.sh`, mas gerado como `all_cells_medicc2.tsv` por `02_cnv_qc_prepare.R` — verificar reconciliação |
| MEDICC2 — árvore | Extensão real é `.new` (Newick), não apenas `.tree`/`.nwk` |
| MEDICC2 — alelos | `prepare_allele_file` aceita `copy_number` inteiro (duplica), `cn_a`/`cn_b` (mantém) ou separadores `|`, `/`, `,` |
| `diploid_normal` | Filtrado das distâncias, CN profiles e branch lengths em `scwgs_analysis.R` |
| `hclust` | Usa `ward.D2` consistentemente |
| Clones (etapa 4) | `k = max(2, min(10, sqrt(n_cells/2)))` |
| Subclones (pós) | `k = 5` por padrão, com alternativas 3, 8, 10, 15 |
| Validação biológica | Spearman(dist. evolutiva, aberração CNV) esperado > 0.2 |
| `job-nanny` | Opcional nas etapas 4 e 5 |
| `ggtree` | Opcional; fallback transparente para `ape::plot.phylo` |
| Conda em wrapper | Etapa 3 usa `~/.bashrc` + `conda activate`; demais usam `module load` + `source conda.sh` — padronizar |

---

## 9. Apêndice — tabela rápida de scripts

| Script | Tipo | Função |
|---|---|---|
| `00_setup_conda.sh` | Setup | Cria ambiente `evo_clonal_medicc2` |
| `submit_all.sh` | Orquestrador | Encadeia todas as etapas |
| `01_download.sh` | SLURM | Download dos dados |
| `02_diagnose.R` | R (opcional) | Diagnostica estrutura dos brutos |
| `02_cnv_qc_prepare.R` | R | QC + inputs MEDICC2 |
| `02_cnv_qc_prepare.sh` | SLURM | Wrapper da etapa 2 |
| `03_medicc2_run.sh` | SLURM | Roda MEDICC2 (global + 3 pacientes) |
| `04_tree_analysis.R` | R | Análise da árvore e clones |
| `04_tree_analysis.sh` | SLURM | Wrapper da etapa 4 |
| `05_visualization.R` | R | Figuras finais |
| `05_visualization.sh` | SLURM | Wrapper da etapa 5 |
| `functions.R` | R | Utilitários compartilhados |
| `submit_scwgs_analysis.sh` | SLURM (array) | Wrapper pós-MEDICC2 |
| `scwgs_analysis.sh` | SLURM (array) | Duplicata do anterior |
| `scwgs_analysis.R` | R | Análise pós-MEDICC2 por paciente |

---

## 10. Saída final esperada

Após `bash submit_all.sh`:

| Caminho | Conteúdo |
|---|---|
| `results/medicc2/global/` | Árvore global |
| `results/medicc2/JK136/` | Árvore JK136 |
| `results/medicc2/JK142/` | Árvore JK142 |
| `results/medicc2/JK153/` | Árvore JK153 |
| `results/tree_analysis/evo_metrics_per_cell.csv` | Distância e clone por célula |
| `results/tree_analysis/clone_summary.csv` | Resumo por clone |
| `results/tree_analysis/patient_summary.csv` | Resumo por paciente |
| `results/figures/*.pdf` | Figuras finais |
| `results/SUMMARY_REPORT.txt` | Relatório consolidado |

Após `sbatch submit_scwgs_analysis.sh`:

| Caminho | Conteúdo |
|---|---|
| `results/scwgs_analysis/JK136/` | Subclones, UMAP, árvore e heatmap de genes |
| `results/scwgs_analysis/JK153/` | Idem |
| `results/scwgs_analysis/global/` | Idem |
