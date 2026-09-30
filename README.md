# Pipeline scDNAme multimodal — syn22257780 (V2, cenário sem BAM)

Este README documenta o pacote de scripts do pipeline **scDNAme multimodal** aplicado ao dataset Synapse `syn22257780`. Ele cobre desde o preparo dos inputs de metilação (coverage agregado por célula/CpG) até clustering, VMRs, DMRs, disorder downstream, MethylVI, PDclust e preparo opcional de SCNA/Ginkgo.

O pipeline foi **ajustado para o cenário real do projeto**:

- ✅ Existe `analysis_scRRBS_bismark_coverage.txt.gz` agregado por célula/CpG.
- ✅ Existem `analysis_scRRBS_sequencing_qc.tsv` e `clinical_metadata.tsv`.
- ❌ Não existem BAMs nem FASTQs. (não conseguimos acesso)

Por isso, o modo padrão é **seguro sem BAM**:

- **MethSCAn** principal (clusterização, VMR, DMR) — habilitado.
- **Disorder downstream** usando tabelas já processadas do paper — habilitado.
- **SCNA/Ginkgo** a partir de BAM — desabilitado por padrão (`ENABLE_SCNA=0`).
- **MethylVI** e **PDclust** — módulos opcionais, rodados separadamente.

> ⚠️ Caminhos absolutos (`/home/matpg/...`) estão hardcoded em vários scripts. Ajuste conforme o ambiente.

---

## 1. Visão geral do fluxo

```
run_master_scDNAme_multimodal_v1.sh
        │
        ├──► run_methscan_syn22257780_v1.sh <MODE>
        │           │
        │           ├── script_prepare_methscan_inputs_syn22257780_v1.py
        │           ├── methscan prepare / filter / smooth / scan / matrix
        │           └── script_methscan_downstream_syn22257780_v1.R
        │
        ├──► run_disorder_syn22257780_v1.sh <MODE>   (depende do MethSCAn)
        │           └── script_disorder_downstream_syn22257780_v1.R
        │
        └──► run_scna_ginkgo_prep_syn22257780_v1.sh <MODE>   (opcional, ENABLE_SCNA=1)
                    └── script_prepare_ginkgo_inputs_syn22257780_v1.py
```

Módulos **não** orquestrados pelo master (rodados separadamente):

```
run_methylvi_methscan_syn22257780_v1.sh
        └── script_methylvi_methscan_syn22257780_v1.py

run_pdclust_syn22257780_v1.sh
        └── script_pdclust_syn22257780_v1.R

resume_methscan_after_scan_v1.sh   (original)
resume_methscan_after_scan_v2.sh   (idempotente)
```

---

## 2. Estrutura de diretórios

Base assumida pelos scripts:

```bash
BASE=/home/matpg/scDNAme
TAB=$BASE/tables
PIPELINE_TAG=v2_nobam_safe
PROJECT_OUT_BASE=$BASE/results_scDNAme_multimodal_${PIPELINE_TAG}
```

| Diretório | Conteúdo |
|---|---|
| `$BASE/scriptsdname/scDNAme_multimodal_v2/` | Scripts do pacote |
| `$BASE/tables/` | Tabelas de QC, clínicas e disorder do paper |
| `$BASE/annotations_methscan/` | BEDs de promotores e genes |
| `$PROJECT_OUT_BASE/methscan_<MODE>/` | Saídas do MethSCAn por modo |
| `$PROJECT_OUT_BASE/logs/` | Logs SLURM |
| `$PROJECT_OUT_BASE/pdclust_all/` | Saídas do PDclust |
| `$BASE/results_scDNAme_multimodal_v2_nobam_safe/methscan_<MODE>/methylvi_casebatch_v1/` | Saídas do MethylVI |

---

## 3. Dependências

### Sistema

| Recurso | Observação |
|---|---|
| SLURM | Gerenciador de jobs |
| `conda` | Ativação de ambiente |
| `python` | Scripts `.py` |
| `Rscript` | Scripts `.R` |
| `find`, `sort` | Usados nos wrappers |
| `job-nanny` | Opcional; se ausente, os wrappers caem para `Rscript` puro |

### Ambiente conda

O pacote assume **um único ambiente conda**, por padrão chamado `methscan`:

```bash
CONDA_ENV=methscan
```

### Python

| Pacote | Uso |
|---|---|
| `pandas` | Leitura e manipulação de tabelas |
| `methscan` | CLI para prepare/filter/smooth/scan/matrix/diff |
| `numpy`, `scipy` | Álgebra numérica (MethylVI) |
| `anndata`, `mudata`, `scanpy`, `scvi-tools`, `torch`, `matplotlib` | MethylVI |

### R

| Pacote | Uso |
|---|---|
| `data.table` | Leitura rápida |
| `dplyr`, `tidyr`, `tibble` | Manipulação |
| `ggplot2` | Visualizações |
| `irlba` | PCA iterativa |
| `uwot` | UMAP |
| `igraph` | Leiden |
| `pheatmap` | Heatmaps (PDclust) |
| `PDclust` | Dissimilaridade pareada |
| `optparse`, `parallel` | CLI e paralelismo (PDclust) |

### Somente para SCNA futuro

| Ferramenta | Uso |
|---|---|
| `samtools` | Conversão BAM → BED |
| `bedtools` | `bamtobed` |

### Checagem rápida

```bash
bash check_methscan_env_v1.sh
```

---

## 4. Entradas esperadas

| Arquivo | Descrição |
|---|---|
| `$TAB/analysis_scRRBS_bismark_coverage.txt.gz` | Coverage agregado por célula/CpG (entrada do MethSCAn) |
| `$TAB/analysis_scRRBS_sequencing_qc.tsv` | QC das células scRRBS |
| `$TAB/clinical_metadata.tsv` | Metadados clínicos (IDH, time point, grau) |
| `$TAB/ref_promoters.tsv` | Referência de promotores (opcional) |
| `$TAB/ref_genes.tsv` | Referência de genes (opcional) |
| `$BASE/annotations_methscan/promoters.noheader.bed` | BED de promotores (usado no resume) |
| `$BASE/annotations_methscan/genes.noheader.bed` | BED de genes (usado no resume) |

### Tabelas de disorder do paper (usadas pelo `run_disorder_*`)

| Arquivo | Conteúdo |
|---|---|
| `analysis_scRRBS_context_specific_DNAme_disorder.tsv` | PDR por contexto |
| `analysis_scRRBS_individual_promoter_DNAme_disorder.tsv` | PDR por promotor |
| `analysis_scRRBS_individual_TFBS_motif_DNAme_disorder.tsv` | PDR por TFBS |
| `analysis_scRRBS_replication_timing_DNAme_disorder.tsv` | PDR por timing de replicação |
| `analysis_RRBS_context_specific_epiallele_methylation.tsv` | Epialelos por contexto |
| `analysis_scRRBS_epiallele_CpG_density_summary.tsv` | Densidade de epialelos |

### Manifesto de BAM (apenas para SCNA futuro)

| Arquivo | Conteúdo |
|---|---|
| `$BASE/bam_manifest_scRRBS.tsv` | Modelo preenchível: `cell_barcode`, `bam_path`, `case_barcode`, `idh_status` |

---

## 5. Saídas principais

### MethSCAn (por modo `all` / `IDHwt` / `IDHmut`)

| Arquivo | Conteúdo |
|---|---|
| `methscan_<MODE>/VMRs.bed` | VMRs detectados |
| `methscan_<MODE>/VMR_matrix/` | Matriz de médias por VMR |
| `methscan_<MODE>/promoter_matrix/` | Matriz por promotor |
| `methscan_<MODE>/gene_matrix/` | Matriz por gene |
| `methscan_<MODE>/prep/qc_keep.tsv` | Células que passaram no QC |
| `methscan_<MODE>/prep/qc_keep_cell_names.txt` | Lista de células mantidas |
| `methscan_<MODE>/prep/cov_by_cell/*.cov.gz` | Coverage por célula |
| `methscan_<MODE>/prep/bed/promoters.bed` | BED de promotores |
| `methscan_<MODE>/prep/bed/genes.bed` | BED de genes |
| `methscan_<MODE>/downstream/cell_cluster_assignments.tsv` | Clusters Leiden por célula |
| `methscan_<MODE>/downstream/pca_embeddings.tsv` | Embeddings da PCA |
| `methscan_<MODE>/downstream/cluster_sizes.tsv` | Tamanho dos clusters |
| `methscan_<MODE>/downstream/cluster_by_case.tsv` | Composição por caso |
| `methscan_<MODE>/downstream/cluster_by_case_fraction.tsv` | Fração por caso |
| `methscan_<MODE>/downstream/UMAP_clusters.pdf` | UMAP colorido por cluster |
| `methscan_<MODE>/downstream/UMAP_cases.pdf` | UMAP colorido por caso |
| `methscan_<MODE>/downstream/UMAP_IDH.pdf` | UMAP colorido por IDH |
| `methscan_<MODE>/downstream/PCA_PC1_PC2.pdf` | PCA PC1 × PC2 |
| `methscan_<MODE>/downstream/promoter_cluster_means.tsv` | Médias por cluster em promotores |
| `methscan_<MODE>/downstream/gene_cluster_means.tsv` | Médias por cluster em genes |
| `methscan_<MODE>/downstream/cell_groups/*.csv` | Grupos cluster vs rest e pairwise |
| `methscan_<MODE>/dmrs/*.bed` | DMRs por contraste |
| `methscan_<MODE>/downstream/sessionInfo.txt` | Info da sessão R |

### Disorder downstream (por modo)

| Arquivo | Conteúdo |
|---|---|
| `disorder_downstream/context_disorder_by_cluster.tsv` | PDR por cluster e contexto |
| `disorder_downstream/context_disorder_by_case.tsv` | PDR por caso e contexto |
| `disorder_downstream/promoter_disorder_by_cluster.tsv` | PDR por promotor e cluster |
| `disorder_downstream/top100_variable_promoter_disorder.tsv` | Top 100 promotores variáveis |
| `disorder_downstream/tfbs_disorder_by_cluster.tsv` | PDR por TFBS e cluster |
| `disorder_downstream/top100_variable_tfbs_disorder.tsv` | Top 100 TFBS variáveis |
| `disorder_downstream/replication_disorder_by_cluster.tsv` | PDR por timing e cluster |
| `disorder_downstream/top100_variable_replication_disorder.tsv` | Top 100 timings variáveis |
| `disorder_downstream/epiallele_context_joined.tsv` | Epialelos por contexto (cruzado) |
| `disorder_downstream/epiallele_density_summary_joined.tsv` | Sumário de densidade (cruzado) |
| `disorder_downstream/Boxplot_context_disorder_by_cluster.pdf` | Boxplot por cluster |
| `disorder_downstream/sessionInfo.txt` | Info da sessão R |

### MethylVI

| Arquivo | Conteúdo |
|---|---|
| `methylvi_casebatch_v1/methylvi_latent.tsv` | Latent space por célula |
| `methylvi_casebatch_v1/methylvi_umap.tsv` | Coordenadas UMAP |
| `methylvi_casebatch_v1/cell_cluster_assignments.tsv` | Clusters Leiden por célula |
| `methylvi_casebatch_v1/cluster_sizes.tsv` | Tamanho dos clusters |
| `methylvi_casebatch_v1/cluster_by_case.tsv` | Composição por caso |
| `methylvi_casebatch_v1/cluster_by_case_fraction.tsv` | Fração por caso |
| `methylvi_casebatch_v1/cell_groups/*.csv` | Grupos cluster vs rest |
| `methylvi_casebatch_v1/UMAP_clusters.pdf` | UMAP por cluster |
| `methylvi_casebatch_v1/UMAP_cases.pdf` | UMAP por caso |
| `methylvi_casebatch_v1/UMAP_IDH.pdf` | UMAP por IDH |
| `methylvi_casebatch_v1/methylvi_input_and_latent.h5mu` | MuData com inputs e latent |
| `methylvi_casebatch_v1/model/` | Modelo salvo |
| `methylvi_casebatch_v1/run_metadata.json` | Hiperparâmetros usados |
| `methylvi_casebatch_v1/dmrs/*.bed` | DMRs por contraste (opcional) |

### PDclust

| Arquivo | Conteúdo |
|---|---|
| `pdclust_pairwise.rds` | Objeto de dissimilaridade pareada |
| `pdclust_dissimilarity_matrix.rds` / `.tsv` | Matriz de dissimilaridade |
| `pdclust_cluster_results.rds` | Resultados do clustering |
| `pdclust_cluster_assignments.tsv` | Atribuições por célula |
| `PDclust_heatmap.pdf` | Heatmap ordenado |
| `PDclust_MDS_by_case.pdf` | MDS por caso |
| `PDclust_MDS_by_IDH.pdf` | MDS por IDH |
| `PDclust_MDS_by_cluster.pdf` | MDS por cluster |
| `pdclust_mds_embeddings.tsv` | Coordenadas MDS |
| `pdclust_cluster_by_case.tsv` | Composição por caso |
| `pdclust_cluster_by_case_fraction.tsv` | Fração por caso |

### SCNA/Ginkgo (futuro)

| Arquivo | Conteúdo |
|---|---|
| `scna_<MODE>/prep/bam_manifest_filtered.tsv` | Manifest filtrado por QC |
| `scna_<MODE>/prep/cells_to_keep.txt` | Lista de células |
| `scna_<MODE>/ginkgo_bed/*.bed` | BEDs por célula |
| `scna_<MODE>/README_SCNA_NEXT_STEPS.txt` | Instruções de follow-up |

---

## 6. Como executar

### Modo seguro (recomendado)

```bash
cd /home/matpg/scDNAme/scriptsdname/scDNAme_multimodal_v2
bash run_master_scDNAme_multimodal_v1.sh all
```

Submete:

| Job | Script | Dependência |
|---|---|---|
| MethSCAn `all` | `run_methscan_syn22257780_v1.sh all` | — |
| MethSCAn `IDHwt` | `run_methscan_syn22257780_v1.sh IDHwt` | — |
| MethSCAn `IDHmut` | `run_methscan_syn22257780_v1.sh IDHmut` | — |
| Disorder `all` | `run_disorder_syn22257780_v1.sh all` | `afterok` do MethSCAn `all` |
| Disorder `IDHwt` | `run_disorder_syn22257780_v1.sh IDHwt` | `afterok` do MethSCAn `IDHwt` |
| Disorder `IDHmut` | `run_disorder_syn22257780_v1.sh IDHmut` | `afterok` do MethSCAn `IDHmut` |

Não submete SCNA.

### Modos do master

| Comando | O que submete |
|---|---|
| `bash run_master_scDNAme_multimodal_v1.sh methscan` | Só MethSCAn (all, IDHwt, IDHmut) |
| `bash run_master_scDNAme_multimodal_v1.sh disorder` | Só disorder (all, IDHwt, IDHmut) |
| `bash run_master_scDNAme_multimodal_v1.sh scna` | Só SCNA (requer `ENABLE_SCNA=1`) |
| `bash run_master_scDNAme_multimodal_v1.sh all` | MethSCAn + disorder (seguro) |
| `bash run_master_scDNAme_multimodal_v1.sh all_with_scna` | MethSCAn + disorder + SCNA (requer `ENABLE_SCNA=1`) |

### Rodar módulos individuais

```bash
sbatch run_methscan_syn22257780_v1.sh all
sbatch run_methscan_syn22257780_v1.sh IDHwt
sbatch run_methscan_syn22257780_v1.sh IDHmut

sbatch run_disorder_syn22257780_v1.sh all
sbatch run_disorder_syn22257780_v1.sh IDHwt
sbatch run_disorder_syn22257780_v1.sh IDHmut

sbatch run_methylvi_methscan_syn22257780_v1.sh all VMR_matrix
sbatch run_pdclust_syn22257780_v1.sh
```

### Resume após o scan

Se o pipeline parou depois do `methscan scan`:

```bash
sbatch resume_methscan_after_scan_v2.sh   # idempotente
# ou
sbatch resume_methscan_after_scan_v1.sh   # original
```

---

## 7. Descrição dos scripts

### 7.1 `config_scDNAme_multimodal_v1.sh`

Configuração central. Exporta variáveis de caminho, ambiente e parâmetros.

| Grupo | Variáveis principais |
|---|---|
| Projeto | `BASE`, `TAB`, `PIPELINE_TAG`, `PROJECT_OUT_BASE`, `LOG_DIR` |
| Ambientes | `CONDA_SH`, `CONDA_ENV`, `PY_ENV`, `METHSCAN_ENV`, `R_ENV` |
| Recursos | `THREADS`, `MEM_GB`, `UMAP_NEIGHBORS`, `UMAP_MIN_DIST`, `LEIDEN_RESOLUTION`, `PCA_NPCS` |
| MethSCAn | `METHSCAN_CHUNKSIZE_BP`, `METHSCAN_SMOOTH_BW`, `METHSCAN_SCAN_BW`, `METHSCAN_SCAN_STEP`, `METHSCAN_VAR_THRESHOLD`, `METHSCAN_MIN_CELLS`, `METHSCAN_DIFF_*` |
| QC | `MIN_UNIQUE_CPG`, `MIN_BS_CONVERSION`, `REQUIRE_TUMOR_STATUS` |
| Entradas | `QC_FILE`, `CLINICAL_FILE`, `BISMARK_COVERAGE_AGG`, `REF_PROMOTERS`, `REF_GENES` |
| Disorder | `CONTEXT_DISORDER_FILE`, `PROMOTER_DISORDER_FILE`, `TFBS_DISORDER_FILE`, `REPLICATION_DISORDER_FILE`, `EPIALLELE_CONTEXT_FILE`, `EPIALLELE_SUMMARY_FILE` |
| SCNA | `BAM_MANIFEST`, `GINKGO_DIR`, `GINKGO_RUN_CLI`, `ENABLE_SCNA`, `GINKGO_GENOME`, `GINKGO_BINNING` |

---

### 7.2 `run_master_scDNAme_multimodal_v1.sh`

Orquestrador SLURM. Recebe o módulo como primeiro argumento (`methscan`, `disorder`, `scna`, `all`, `all_with_scna`). Submete jobs com `sbatch --parsable` e usa `--dependency=afterok` para encadear disorder após MethSCAn.

| Módulo | Comportamento |
|---|---|
| `methscan` | Submete os três modos do MethSCAn |
| `disorder` | Submete os três modos do disorder |
| `scna` | Só roda se `ENABLE_SCNA=1`; senão aborta com mensagem |
| `all` | MethSCAn + disorder (modo seguro) |
| `all_with_scna` | MethSCAn + disorder + SCNA (requer `ENABLE_SCNA=1`) |

---

### 7.3 `run_methscan_syn22257780_v1.sh`

Pipeline MethSCAn completo por modo.

| Etapa | Comando | Descrição |
|---|---|---|
| 1/7 | `script_prepare_methscan_inputs_*.py` | Filtra QC, separa coverage por célula, gera BEDs |
| 2/7 | `methscan prepare` | Converte `.cov.gz` em formato interno |
| 3/7 | `methscan filter` | Aplica whitelist de células |
| 4/7 | `methscan smooth` | Suaviza sinal |
| 5/7 | `methscan scan` | Detecta VMRs |
| 6/7 | `methscan matrix` | Gera matrizes VMR, promotor e gene |
| 7/7 | `script_methscan_downstream_*.R` | PCA, UMAP, Leiden, composição |
| extra | `methscan diff` | DMRs por contraste cluster vs rest e pairwise |

---

### 7.4 `script_prepare_methscan_inputs_syn22257780_v1.py`

Prepara os inputs do MethSCAn a partir do coverage agregado.

| Etapa | Descrição |
|---|---|
| 1 | Lê QC e tabela clínica, infere colunas por heurística |
| 2 | Filtra células: `cpg_unique > min_unique_cpg`, `bs > min_bs`, `tumor_status == require_tumor` |
| 3 | Filtra por modo (`all` / `IDHwt` / `IDHmut`) |
| 4 | Salva `qc_keep.tsv` e `qc_keep_cell_names.txt` |
| 5 | Lê o coverage em chunks e separa em `.cov.gz` por célula |
| 6 | Converte referências de promotores/genes em BEDs |

---

### 7.5 `script_methscan_downstream_syn22257780_v1.R`

Downstream do MethSCAn.

| Etapa | Descrição |
|---|---|
| 1 | Lê `mean_shrunken_residuals.csv.gz` da matriz VMR |
| 2 | Junta QC + clínica |
| 3 | PCA iterativa com `irlba::prcomp_irlba` (imputação de NAs) |
| 4 | UMAP via `uwot::umap` com grafo de vizinhos |
| 5 | Leiden via `igraph::cluster_leiden` |
| 6 | Escreve `cell_cluster_assignments.tsv`, `pca_embeddings.tsv`, `cluster_sizes.tsv`, `cluster_by_case*.tsv` |
| 7 | Gera UMAPs (cluster, caso, IDH) e PCA PC1×PC2 |
| 8 | Escreve grupos cluster vs rest e pairwise em `cell_groups/` |
| 9 | Sumariza matrizes de promotores e genes por cluster |
| 10 | Salva `sessionInfo.txt` |

---

### 7.6 `run_disorder_syn22257780_v1.sh`

Wrapper SLURM para o disorder downstream.

| Etapa | Descrição |
|---|---|
| 1 | Valida modo (`all`/`IDHwt`/`IDHmut`) |
| 2 | Carrega config e ativa ambiente conda |
| 3 | Verifica existência de `cell_cluster_assignments.tsv` |
| 4 | Chama `script_disorder_downstream_*.R` com tabelas de disorder |

---

### 7.7 `script_disorder_downstream_syn22257780_v1.R`

Cruzamento das tabelas de disorder com clusters MethSCAn.

| Etapa | Descrição |
|---|---|
| 1 | Lê atribuições de cluster |
| 2 | Cruza com tabela de contexto (PDR, promoter_PDR, enhancer_PDR, cgi_PDR) |
| 3 | Gera `context_disorder_by_cluster.tsv` e `_by_case.tsv` |
| 4 | Gera boxplot por cluster e contexto |
| 5 | Sumariza promotor, TFBS e replicação (top 100 variáveis) |
| 6 | Cruza epialelos (contexto e densidade) |
| 7 | Salva `sessionInfo.txt` |

---

### 7.8 `run_methylvi_methscan_syn22257780_v1.sh`

Wrapper SLURM para o MethylVI.

| Etapa | Descrição |
|---|---|
| 1 | Ativa ambiente `methylvi` (configurável via `SCVI_ENV`) |
| 2 | Chama `script_methylvi_methscan_*.py` por modo |
| 3 | Se `RUN_DIFF=1`, roda `methscan diff` para cada grupo |

Variáveis de ambiente aceitas:

| Variável | Default | Uso |
|---|---|---|
| `SCVI_ENV` | `methylvi` | Nome do ambiente conda |
| `RUN_DIFF` | `1` | Rodar DMRs pós-MethylVI |
| `N_LATENT` | `10` | Dimensão latente |
| `MAX_EPOCHS` | `250` | Épocas de treino |
| `BATCH_SIZE` | `128` | Tamanho do batch |
| `MIN_CELLS_PER_FEATURE` | `10` | Filtro mínimo por feature |
| `MAX_FEATURES` | `5000` | Máx. de features |
| `NNEI` | `30` | Vizinhos UMAP |
| `MIN_DIST` | `0.10` | `min_dist` UMAP |
| `RES` | `0.12` | Resolução Leiden |
| `OUT_SUBDIR` | `methylvi_casebatch_v1` | Subdiretório de saída |

---

### 7.9 `script_methylvi_methscan_syn22257780_v1.py`

Adaptação do METHYLVI para saídas do MethSCAn.

| Etapa | Descrição |
|---|---|
| 1 | Lê `methylated_sites.csv.gz` e `total_sites.csv.gz` da matriz escolhida |
| 2 | Detecta orientação (células × features) |
| 3 | Filtra features por cobertura mínima e mantém até `max_features` |
| 4 | Calcula fração de metilação = `mc / cov` |
| 5 | Constrói `AnnData` com `layers['mc']`, `layers['cov']` e `.X` como fração |
| 6 | Empacota em `MuData` com uma única modalidade `mCG` |
| 7 | `METHYLVI.setup_mudata()` com `batch_key=case_barcode` |
| 8 | Treina o modelo |
| 9 | Extrai latent, roda vizinhos, UMAP e Leiden |
| 10 | Salva tabelas, UMAPs e o modelo |

---

### 7.10 `run_pdclust_syn22257780_v1.sh`

Wrapper SLURM para o PDclust (reprodução da lógica da Fig. 1b).

| Parâmetro | Valor default |
|---|---|
| Tempo | 48 h |
| CPUs | 8 |
| Memória | 128 GB |
| `--n_clusters` | 6 |
| `--cores_pairwise` | 4 |
| `--cores_read` | 8 |

---

### 7.11 `script_pdclust_syn22257780_v1.R`

Análise PDclust/MDS.

| Etapa | Descrição |
|---|---|
| 1 | Lê QC e clínica; filtra `cpg_unique > 40000`, `bs > 95%`, `tumor_status == 1` |
| 2 | Localiza `.cov.gz` por célula em `cov_by_cell/` |
| 3 | Lê coverage e converte para formato PDclust (com remoção de GL/X/Y/MT por padrão) |
| 4 | `create_pairwise_master` para dissimilaridade pareada |
| 5 | `convert_to_dissimilarity_matrix` |
| 6 | `cluster_dissimilarity` com `num_clusters` configurável |
| 7 | Gera heatmap ordenado por `pheatmap` |
| 8 | Gera MDS (por caso, IDH e cluster) |
| 9 | Salva composições por cluster e por caso |

---

### 7.12 `run_scna_ginkgo_prep_syn22257780_v1.sh`

Wrapper SCNA/Ginkgo (desabilitado por padrão).

| Etapa | Descrição |
|---|---|
| 1 | Aborta se `ENABLE_SCNA != 1` |
| 2 | Valida `BAM_MANIFEST` |
| 3 | Chama `script_prepare_ginkgo_inputs_*.py` |
| 4 | Gera `.bed` por célula com `samtools view` + `bedtools bamtobed` |
| 5 | Escreve `README_SCNA_NEXT_STEPS.txt` com instruções |
| 6 | Se `GINKGO_RUN_CLI=1` e `GINKGO_GENOME=hg19`, roda Ginkgo CLI |

> ⚠️ A documentação pública do Ginkgo lista bins prontos para hg19 e outros genomas antigos — hg38 não aparece. Se seus BAMs forem hg38, valide antes de rodar.

---

### 7.13 `script_prepare_ginkgo_inputs_syn22257780_v1.py`

Filtra o manifesto de BAMs pelo QC e pelo modo.

| Etapa | Descrição |
|---|---|
| 1 | Lê manifesto, QC e clínica |
| 2 | Infere colunas por heurística |
| 3 | Normaliza barcodes (`.` → `-`) |
| 4 | Aplica filtros de QC |
| 5 | Filtra por modo (`all`/`IDHwt`/`IDHmut`) |
| 6 | Verifica existência física de cada BAM |
| 7 | Salva `bam_manifest_filtered.tsv` e `cells_to_keep.txt` |

---

### 7.14 `resume_methscan_after_scan_v1.sh` e `resume_methscan_after_scan_v2.sh`

Retomam o MethSCAn após a etapa `scan`.

| Script | Comportamento |
|---|---|
| `v1` | Sempre roda `matrix` (VMR, promotor, gene) + downstream + diff |
| `v2` | Idempotente: pula `matrix` se o `mean_shrunken_residuals.csv.gz` já existir |

Ambos iteram sobre os três modos (`all`, `IDHwt`, `IDHmut`).

| Parâmetro | `v1` | `v2` |
|---|---|---|
| Tempo | 24 h | 10 h |
| CPUs | 8 | 8 |
| Memória | 64 GB | 64 GB |
| Idempotente | ❌ | ✅ |

---

## 8. Observações e cuidados

| Item | Observação |
|---|---|
| Caminhos hardcoded | Ajustar `BASE`, `TAB`, `CONDA_ENV`, `BISMARK_COVERAGE_AGG` e caminhos nos wrappers |
| `job-nanny` | Opcional; wrappers caem para `Rscript` puro |
| Ambiente único | `methscan` é o ambiente único; `PY_ENV`, `METHSCAN_ENV` e `R_ENV` são aliases |
| Disorder | Não recalcula PDR a partir de reads; usa tabelas do paper |
| MethylVI | Assume CpG-only (RRBS); por isso, uma única modalidade `mCG` |
| SCNA | Bloqueado por padrão (`ENABLE_SCNA=0`); só habilitar com BAM_MANIFEST válido |
| Ginkgo | Bins públicos para hg19; validar antes de usar hg38 |
| `ulimit -n` | Ajustar para valores altos antes do `methscan prepare` (muitos `.cov.gz`) |
| Primeira rodada | Recomendado começar por `IDHwt` para calibrar `LEIDEN_RESOLUTION` |
| PDclust | Por padrão remove GL/X/Y/MT para casar com a Fig. 1b |
| Resume v2 | Idempotente: evita sobrescrever `matrix` já gerado |
| MuData | MethylVI salva `.h5mu` com inputs e latent |

---

## 9. Apêndice — tabela rápida de scripts

| Script | Tipo | Função |
|---|---|---|
| `config_scDNAme_multimodal_v1.sh` | Config | Configuração central |
| `run_master_scDNAme_multimodal_v1.sh` | Orquestrador | Dispara módulos via SLURM |
| `run_methscan_syn22257780_v1.sh` | Wrapper | Pipeline MethSCAn completo |
| `script_prepare_methscan_inputs_syn22257780_v1.py` | Python | Prepara inputs do MethSCAn |
| `script_methscan_downstream_syn22257780_v1.R` | R | PCA/UMAP/Leiden + grupos |
| `run_disorder_syn22257780_v1.sh` | Wrapper | Disorder downstream |
| `script_disorder_downstream_syn22257780_v1.R` | R | Cruza disorder com clusters |
| `run_methylvi_methscan_syn22257780_v1.sh` | Wrapper | MethylVI |
| `script_methylvi_methscan_syn22257780_v1.py` | Python | METHYLVI adaptado |
| `run_pdclust_syn22257780_v1.sh` | Wrapper | PDclust |
| `script_pdclust_syn22257780_v1.R` | R | PDclust/MDS |
| `run_scna_ginkgo_prep_syn22257780_v1.sh` | Wrapper | SCNA/Ginkgo (desabilitado) |
| `script_prepare_ginkgo_inputs_syn22257780_v1.py` | Python | Prepara BEDs para Ginkgo |
| `resume_methscan_after_scan_v1.sh` | Wrapper | Retoma pós-scan (não idempotente) |
| `resume_methscan_after_scan_v2.sh` | Wrapper | Retoma pós-scan (idempotente) |
| `check_methscan_env_v1.sh` | Check | Valida ambiente conda |
| `template_bam_manifest_scRRBS.tsv` | Template | Modelo de manifesto de BAMs |

---

## 10. Saída final esperada

Após rodar `bash run_master_scDNAme_multimodal_v1.sh all`:

| Diretório | Conteúdo |
|---|---|
| `$PROJECT_OUT_BASE/methscan_all/` | Clusters, VMRs, DMRs e matrizes (todas as células) |
| `$PROJECT_OUT_BASE/methscan_IDHwt/` | Mesmas saídas apenas para IDH-wt |
| `$PROJECT_OUT_BASE/methscan_IDHmut/` | Mesmas saídas apenas para IDH-mut |
| `$PROJECT_OUT_BASE/methscan_<MODE>/disorder_downstream/` | Sumários de disorder por cluster/caso |
| `$PROJECT_OUT_BASE/methscan_<MODE>/dmrs/` | DMRs por contraste |

Arquivos-chave para inspeção rápida:

```r
clusters <- read.delim(
  "$PROJECT_OUT_BASE/methscan_all/downstream/cell_cluster_assignments.tsv"
)

qc_keep <- read.delim(
  "$PROJECT_OUT_BASE/methscan_all/prep/qc_keep.tsv"
)

context_disorder <- read.delim(
  "$PROJECT_OUT_BASE/methscan_all/disorder_downstream/context_disorder_by_cluster.tsv"
)
```
