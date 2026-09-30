# Pipeline scRNA-seq de gliomas — GSE182109 + Synapse

- Reconstrói e prepara três fontes de dados: **GSE182109**, **Synapse** e **GSE173278**.
- Integra via Harmony apenas **GSE182109 + Synapse**.
- Gera um objeto final anotado: `annotated_integrated_v3.rds`.

> ⚠️ Os caminhos são absolutos e estão hardcoded. Ajuste conforme o seu ambiente antes de executar. O GSE173278 foi removido da integração pois não consegui corrigir os batch effects dele com os demais datasets.

---

## 1. Visão geral do fluxo

```
01b_load_gse173278_primary.sh
        │
        ▼
02b_prepare_gse173278_primary.sh
        │
        ▼
03b_harmony_integration_3datasets.sh
        │
        ▼
04b_downstream_analysis_3datasets.sh
```

Fluxo de dados:

```
GSE182109 (gse_raw.rds) ─┐
                         ├─► QC/SCT ─► Harmony ─► clusters/UMAP ─► downstream
Synapse (synapse_raw.rds)┘

GSE173278 (raw) ─► QC leve ─► label transfer ─► validação de anotações
```

O pipeline é orquestrado por `submit_pipeline.sh`, que submete os jobs no SLURM com dependências `afterok`.

---

## 2. Estrutura de diretórios esperada

Base assumida pelos scripts:

```bash
BASE=/home/renanomete/projetos/matdata/met_mat_data/singlecell_novo
```

| Diretório | Conteúdo |
|---|---|
| `results/` | Resultados antigos, incluindo `gse_raw.rds` validado |
| `results_v3/` | Saídas do pipeline V3 |
| `logs_v3/` | Logs SLURM |
| `synapse_data/` | Dados brutos do Synapse |
| `gse173278/` | Dados brutos do GSE173278 |
| `scripts/` | Scripts do pipeline |

---

## 3. Dependências

### Sistema

| Recurso | Observação |
|---|---|
| SLURM | Gerenciador de jobs |
| `module load miniconda/24.4.0-libmamba` | Módulo do cluster |
| Conda env `gbm_scrnaseq` | Ambiente R/Python do pipeline |
| `job-nanny` | Wrapper para `Rscript`; se não existir, substituir por `Rscript` |

### R

| Pacote | Uso |
|---|---|
| Seurat, SeuratObject, SeuratDisk | Manipulação de objetos e conversão h5ad |
| harmony | Integração de batches |
| dplyr, tidyverse, tibble, stringr | Manipulação de dados |
| ggplot2, patchwork | Visualizações |
| data.table | Leitura/escrita rápida |
| Matrix | Matrizes esparsas |
| scDblFinder, SingleCellExperiment, BiocParallel | Detecção de doublets |
| readxl | Leitura de metadados clínicos (.xlsx) |
| future | Paralelização |

### Python

| Pacote | Uso |
|---|---|
| anndata | Leitura de h5ad do Synapse |

---

## 4. Entradas esperadas

| Arquivo | Descrição |
|---|---|
| `results/gse_raw.rds` | Objeto Seurat bruto do GSE182109 já validado |
| `synapse_data/analysis_scRNAseq_tumor_counts.h5ad` | Contagens Synapse em h5ad |
| `synapse_data/analysis_scRNAseq_tumor_gene_expression.tsv.gz` | Alternativa TSV de contagens Synapse |
| `synapse_data/41588_2021_926_MOESM2_ESM.xlsx` | Metadados clínicos do Synapse |
| `gse173278/GSE173278_scRNAseq_filtered_cells_barcodes.tsv.gz` | Barcodes GSE173278 |
| `gse173278/GSE173278_scRNAseq_filtered_cells_genes.tsv.gz` | Genes GSE173278 |
| `gse173278/GSE173278_scRNAseq_filtered_cells_metadata.csv.gz` | Metadados GSE173278 |
| `gse173278/GSE173278_scRNAseq_filtered_cells_norm_counts_matrix.mtx.gz` | Matriz normalizada GSE173278 |

O script `01b` também gera `analysis_scRNAseq_tumor_counts_obs.csv` a partir do h5ad do Synapse.

---

## 5. Saídas principais

| Arquivo | Conteúdo |
|---|---|
| `results_v3/gse_raw.rds` | Cópia validada do GSE182109 bruto |
| `results_v3/synapse_raw.rds` | Synapse bruto reconstruído com `Patient` correto |
| `results_v3/gse173278_primary_raw.rds` | GSE173278 bruto para label transfer |
| `results_v3/gse_qc.rds` | GSE182109 pós-QC/SCT |
| `results_v3/synapse_qc.rds` | Synapse pós-QC/SCT |
| `results_v3/gse173278_primary_qc.rds` | GSE173278 pós-QC leve |
| `results_v3/integrated_harmony_v3.rds` | Objeto integrado GSE182109 + Synapse |
| `results_v3/annotated_integrated_v3.rds` | Objeto final anotado |
| `results_v3/markers_all_clusters_v3.csv` | Marcadores por cluster |
| `results_v3/cluster_annotation_auto_v3.csv` | Anotação automática por cluster |
| `results_v3/label_transfer_g173_v3.csv` | Rótulos transferidos do GSE173278 |
| `results_v3/idh_mut_integrated_v3.rds` | Subset IDH-mut |
| `results_v3/idh_wt_integrated_v3.rds` | Subset IDH-wt |
| `results_v3/umap_*.pdf` | UMAPs por cluster, dataset, paciente, IDH |
| `results_v3/pca_elbow_v3.pdf` | Elbow plot da PCA |
| `exploration_output/exploration_summary.txt` | Sumário de exploração dos objetos |

---

## 6. Como executar

### Pipeline completo

```bash
cd /home/renanomete/projetos/matdata/met_mat_data/singlecell_novo/scripts
bash submit_pipeline.sh
```

Sequência de jobs:

| Ordem | Script | Dependência |
|---|---|---|
| 1 | `01b_load_gse173278_primary.sh` | — |
| 2 | `02b_prepare_gse173278_primary.sh` | `afterok` do job 1 |
| 3 | `03b_harmony_integration_3datasets.sh` | `afterok` do job 2 |
| 4 | `04b_downstream_analysis_3datasets.sh` | `afterok` do job 3 |

Execução individual:

```bash
sbatch 01b_load_gse173278_primary.sh
sbatch 02b_prepare_gse173278_primary.sh
# etc.
```

### Exploração dos objetos

```bash
sbatch run_explore_objects.sh
```

Gera `exploration_summary.txt` com informações de células, genes, colunas de metadados, clusters e identidades ativas.

---

## 7. Descrição dos scripts

### 7.1 `submit_pipeline.sh`

Orquestrador SLURM. Submete os quatro jobs principais em cadeia, usando `--dependency=afterok`. Exibe os IDs dos jobs.

---

### 7.2 `01b_load_gse173278_primary.sh`

**Propósito:** preparar os objetos brutos das três fontes.

| Etapa | O que faz |
|---|---|
| 1 | Copia `results/gse_raw.rds` para `results_v3/gse_raw.rds` |
| 2 | Reconstrói o Synapse (ver abaixo) e salva `synapse_raw.rds` |
| 3 | Carrega GSE173278 (ver abaixo) e salva `gse173278_primary_raw.rds` |

**Reconstrução do Synapse:**

| Passo | Detalhe |
|---|---|
| Leitura | `analysis_scRNAseq_tumor_counts.h5ad` via `SeuratDisk` ou TSV de contagens |
| OBS | Exporta `obs` do h5ad para CSV |
| Junção | Casa células com metadados do Synapse |
| Correção de Patient | Mapeia `sampleid` → `SM001...SM019` por contagem esperada de células |
| Metadados clínicos | Anexa `41588_2021_926_MOESM2_ESM.xlsx` |
| Saída | `synapse_raw.rds` |

**Carregamento do GSE173278:**

| Passo | Detalhe |
|---|---|
| Leitura | barcodes, genes, metadados e matriz `.mtx.gz` |
| Filtro | Apenas tecido primário; remove gliomasferas, PDEs e recorrentes |
| Downsampling | Máximo de 30.000 células (seed 42) |
| Objeto Seurat | Matriz normalizada no slot `data`; binária placeholder em `counts` |
| Metadados | `dataset = "GSE173278_primary"`, `Type = "Primary GBM"`, `ModelSystem = "Tissue"`, `IDH_status = "IDH-wt"` |
| Saída | `gse173278_primary_raw.rds` |

**Saídas:** `gse_raw.rds`, `synapse_raw.rds`, `gse173278_primary_raw.rds`.

---

### 7.3 `02b_prepare_gse173278_primary.sh`

**Propósito:** QC e normalização.

| Dataset | Estratégia de QC | Normalização |
|---|---|---|
| GSE182109 | Por amostra (`orig.ident`); filtros: `nFeature_RNA > 300`, `< 8000`, `percent.mt < 25%`; doublets com `scDblFinder` | `SCTransform` v2, 3000 HVGs, regride `percent.mt` |
| Synapse | Objeto único; filtros: `nFeature_RNA > 300`, `< 9000`, `percent.mt < 20%`; doublets com `scDblFinder` | `SCTransform` v2, 3000 HVGs, regride `percent.mt` |
| GSE173278 | QC leve: mantém células com `nFeature_RNA > 200` | Não roda `SCTransform` (dados já normalizados); apenas `FindVariableFeatures` (3000 genes) |

GSE182109 é processado por amostra e depois mesclado.

**Saídas:** `gse_qc.rds`, `synapse_qc.rds`, `gse173278_primary_qc.rds`, PDFs de QC.

---

### 7.4 `03b_harmony_integration_3datasets.sh`

**Propósito:** integrar GSE182109 + Synapse via Harmony.

> GSE173278 não entra aqui.

| Etapa | Descrição |
|---|---|
| 1 | Lê `gse_qc.rds` e `synapse_qc.rds` |
| 2 | Reconstrói objetos mínimos com assay RNA (padroniza nomes de genes: remove `.N`, troca `_`→`-`, `toupper`) |
| 3 | Normaliza e seleciona 3000 HVGs por dataset |
| 4 | `SelectIntegrationFeatures` |
| 5 | `merge` dos objetos |
| 6 | Define `IDH_status` (GSE182109: LGG→IDH-mut, GBM→IDH-wt; Synapse: mapa fixo `SM001` etc.) |
| 7 | Cria `batch = dataset_Patient` |
| 8 | `ScaleData` + `RunPCA` (50 PCs) |
| 9 | Harmony com `group.by.vars = "batch"`, `dims.use = 1:30` |
| 10 | `FindNeighbors` + `FindClusters` em resoluções 0.2, 0.4, 0.6, 0.8, 1.0 |
| 11 | Identidade ativa = `RNA_snn_res.0.6` |
| 12 | UMAP com redução `harmony` |

**Saídas:** `integrated_harmony_v3.rds`, `pca_elbow_v3.pdf`.

---

### 7.5 `04b_downstream_analysis_3datasets.sh`

**Propósito:** análises downstream, anotação e label transfer.

| Etapa | Descrição |
|---|---|
| 1 | Carrega `integrated_harmony_v3.rds` |
| 2 | Padroniza metadados e recalcula `IDH_status` |
| 3 | Gera UMAPs (clusters, dataset, paciente, IDH) |
| 4 | `FindAllMarkers` por cluster (Wilcoxon, `only.pos`, `min.pct=0.25`, `logfc.threshold=0.25`, máx 3000 células/identidade) |
| 5 | Anotação automática via `AddModuleScore` com painéis de marcadores |
| 6 | Composição IDH por cluster + subsets `idh_mut_integrated_v3.rds` e `idh_wt_integrated_v3.rds` |
| 7 | Label transfer do GSE173278 (referência = objeto integrado; query = GSE173278) |
| 8 | Salva `annotated_integrated_v3.rds` |

**Painéis de marcadores usados na anotação automática:**

| Tipo celular | Genes |
|---|---|
| Glioma | SOX2, OLIG1, OLIG2, GFAP, S100B, CHI3L1, PDGFRA |
| Myeloid | PTPRC, ITGAM, CD68, P2RY12, TMEM119, CD14, FCER1G |
| Tcells | CD3D, CD3E, CD4, CD8A, IL7R |
| Bcells | CD79A, MS4A1, CD19 |
| Endothelial | PECAM1, VWF, KDR |
| Pericytes | PDGFRB, RGS5, ACTA2 |
| Oligodendrocytes | MBP, MOG, PLP1 |

**Label transfer:**

| Passo | Detalhe |
|---|---|
| Âncoras | `FindTransferAnchors(dims=1:30, reference.reduction="pca")` |
| Transferência | `TransferData` de `predicted_celltype` e cluster |
| Projeção | `MapQuery` no UMAP da referência |
| Fallback | Se `MapQuery` falhar, roda UMAP próprio do GSE173278 |
| Saídas | `label_transfer_g173_v3.csv`, `umap_label_transfer_v3.pdf` |

**Saídas:** `annotated_integrated_v3.rds`, marcadores, UMAPs, composição IDH, subsets IDH, resultados de label transfer.

---

## 8. Scripts auxiliares

### 8.1 `explore_seurat_objects.R`

Script R para inspecionar objetos Seurat.

```bash
Rscript explore_seurat_objects.R \
  --object1 annotated_integrated_v3.rds \
  --object2 integrated_harmony_v3.rds \
  --outdir exploration_output
```

Conteúdo do `exploration_summary.txt`:

| Seção | Conteúdo |
|---|---|
| Informações gerais | Classe, número de células e genes |
| Assays e reduções | Disponíveis no objeto |
| `meta.data` | Colunas e sumários |
| Colunas de cluster | Detectadas automaticamente |
| Identidades ativas | Tabela de `Idents` |

### 8.2 `run_explore_objects.sh`

Wrapper SLURM para o script acima.

| Parâmetro | Valor |
|---|---|
| Tempo | 10 min |
| CPUs | 2 |
| Memória | 64 GB |
| Conda env | `R` |
| Caminhos | Apontam para `/home/matpg/scRNA/singlecell/integrated/results_v3` (ajustar) |

### 8.3 `scrna_integrated_objects.R` (legado)

Script antigo de integração focado no GSE182109. Assume que `merged_objects` já existe no ambiente.

| Etapa | Descrição |
|---|---|
| Normalização | `NormalizeData`, `FindVariableFeatures` (2000 genes), `ScaleData` |
| PCA | Combina genes variáveis + lista imune + lista estromal |
| Diagnóstico | `JackStraw`, `ElbowPlot` |
| Integração | Harmony corrigindo por `Patient` |
| Clusterização | `RunUMAP`, `FindNeighbors`, `FindClusters` (resolução 0.8021) |
| UMAPs | Por fragmento, paciente, tipo tumoral, tipo celular, subtipo e fase do ciclo |

> Este script não faz parte do pipeline V3. Está aqui como referência/legado.

---

## 9. Observações e cuidados

| Item | Observação |
|---|---|
| Caminhos hardcoded | Ajustar `BASE`, `RESULTS`, `SYNAPSE_DIR`, `G173_DIR` e caminhos em `run_explore_objects.sh` |
| Ambiente | `module load miniconda/24.4.0-libmamba` + `source activate gbm_scrnaseq` |
| `job-nanny` | Se indisponível, substituir `job-nanny Rscript` por `Rscript` |
| Seurat v5 | Scripts usam `JoinLayers` e checam `Assay5`; podem precisar de adaptação em v3/v4 |
| GSE173278 | Não é integrado; usado apenas como query no label transfer |
| Downsampling | GSE173278 reduzido para no máximo 30.000 células antes do QC |
| Anotação automática | `predicted_celltype` é baseado em módulos de genes, não em curadoria manual |
| IDH status | Derivado de regras e mapa fixo para Synapse; conferir antes de interpretar |
| Recursos | Jobs principais pedem de 96 GB a 180 GB e até 24 h; ajustar ao cluster |

---

## 10. Apêndice — relação com o `README_gse182109.md` original

O README original documentava um pipeline focado apenas no GSE182109.

| Script original | Função |
|---|---|
| `gse182109.R` | Carregamento e montagem do objeto Seurat |
| `gse182109_harmony.R` | Integração, clusterização e anotação |
| `denovoclustering.R` | Re-clusterização focada em glioma + oligodendrócitos |
| `cellstates_charts.R` | Visualização dos estados celulares |
| `repairgenes_findmarkers_clusterprofiler.Rmd` | Genes de reparo de DNA + enriquecimento |
| `DESeq2pipeline_Mateus.R` | Análise bulk RNA-seq com DESeq2 |
| `seuratutorial.R` | Tutorial Seurat (aprendizado) |

O pipeline V3 documentado aqui **não substitui** essas análises biológicas específicas; ele fornece uma etapa anterior de integração multi-dataset e anotação. Para reproduzir estados de glioma, enriquecimento funcional ou DESeq2, consulte o README original e os scripts correspondentes.

---

## 11. Saída final esperada

Objeto principal:

```r
so <- readRDS("results_v3/annotated_integrated_v3.rds")
```

Conteúdo do objeto:

| Item | Descrição |
|---|---|
| Células | GSE182109 + Synapse |
| Clusters | `RNA_snn_res.0.2`, `0.4`, `0.6`, `0.8`, `1.0` |
| Reduções | `pca`, `harmony`, `umap` |
| Metadados | `IDH_status`, `batch`, `dataset`, `Patient`, `Type` |
| Anotação | `predicted_celltype` por módulo de genes |
| Scores | Por tipo celular |
| Subsets | IDH-mut e IDH-wt salvos separadamente |

Label transfer do GSE173278:

```r
lt <- fread("results_v3/label_transfer_g173_v3.csv")
```

Colunas: `predicted_celltype`, `predicted_cluster`, `prediction_score` para cada célula do GSE173278.
