# Pacote operacional — scDNAme multimodal (syn22257780) — versão corrigida para cenário sem BAM

Este pacote foi ajustado para o cenário real do seu projeto:

- você **tem** `analysis_scRRBS_bismark_coverage.txt.gz` agregado por célula/CpG
- você **tem** `analysis_scRRBS_sequencing_qc.tsv` e `clinical_metadata.tsv`
- você **não tem** BAMs nem FASTQs

Por isso, o pacote agora assume um modo **seguro sem BAM**:

1. **MethSCAn principal** para clusterização/VMR/DMR.
2. **DNAme disorder downstream** usando as tabelas já processadas do repositório do artigo.
3. **SCNA/Ginkgo desabilitado por padrão**.

## O que roda neste cenário

### Roda

- `MethSCAn`
- `disorder downstream`

### Não roda neste cenário

- `SCNA/Ginkgo` a partir de BAM
- recálculo de `DNAme disorder` diretamente de epialleles em reads

## Principais correções desta versão

1. `BISMARK_COVERAGE_AGG` agora aponta por padrão para:

```bash
$TAB/analysis_scRRBS_bismark_coverage.txt.gz
```

2. `run_master_scDNAme_multimodal_v1.sh all` agora é um modo **seguro**, que submete apenas:
   - MethSCAn
   - disorder downstream dependente do MethSCAn

3. O módulo `SCNA` ficou bloqueado por padrão com:

```bash
ENABLE_SCNA=0
```

4. O módulo `all_with_scna` só deve ser usado se no futuro você obtiver BAMs e preencher `BAM_MANIFEST`.

## Estrutura dos arquivos

- `config_scDNAme_multimodal_v1.sh`
  - configuração central de caminhos, ambientes e parâmetros.
- `run_master_scDNAme_multimodal_v1.sh`
  - dispara os módulos em lote.
- `run_methscan_syn22257780_v1.sh`
  - pipeline MethSCAn por modo (`all`, `IDHwt`, `IDHmut`).
- `script_prepare_methscan_inputs_syn22257780_v1.py`
  - filtra QC, separa o coverage agregado em arquivos `.cov.gz` por célula e gera BEDs de promotores/genes.
- `script_methscan_downstream_syn22257780_v1.R`
  - PCA iterativa, UMAP, Leiden, composições e geração de grupos para `methscan diff`.
- `run_disorder_syn22257780_v1.sh`
  - usa os clusters do MethSCAn e resume disorder por cluster/caso/anotação.
- `script_disorder_downstream_syn22257780_v1.R`
  - lê as tabelas de disorder já processadas e cruza com os clusters MethSCAn.
- `run_scna_ginkgo_prep_syn22257780_v1.sh`
  - mantido no pacote, mas desabilitado por padrão neste cenário sem BAM.
- `script_prepare_ginkgo_inputs_syn22257780_v1.py`
  - mantido para uso futuro, caso existam BAMs.
- `template_bam_manifest_scRRBS.tsv`
  - modelo de manifesto de BAMs.

## Dependências do ambiente `methscan`

Este pacote agora assume **um único ambiente conda**, por padrão:

```bash
CONDA_ENV=methscan
```

### Python externo

Necessários pelos scripts do pacote:

- `pandas`
- `methscan` (CLI disponível no PATH do ambiente)

Os demais imports Python são da biblioteca padrão: `argparse`, `csv`, `gzip`, `os`, `re`, `sys`, `pathlib`, `typing`.

### R externo

Necessários pelos scripts downstream:

- `data.table`
- `dplyr`
- `tidyr`
- `tibble`
- `ggplot2`
- `irlba`
- `uwot`
- `igraph`

### Utilitários de sistema

Necessários neste cenário sem BAM:

- `bash`
- `sbatch`
- `conda`
- `python`
- `Rscript`
- `find`, `sort`

O `job-nanny` é **opcional**.

### Somente para SCNA futuro

Se um dia você ativar a trilha SCNA com BAM, o mesmo ambiente ou o PATH do job também precisará ter:

- `samtools`
- `bedtools`

### Checagem rápida do ambiente

Incluí o script:

```bash
check_methscan_env_v1.sh
```

Use antes da primeira submissão:

```bash
bash check_methscan_env_v1.sh
```

## Uso sugerido

### 1) Copiar os arquivos para o servidor

Copie a pasta inteira para algo como:

```bash
/home/matpg/scDNAme/scriptsdname/scDNAme_multimodal_v2/
```

### 2) Ajustar a configuração

Edite:

```bash
config_scDNAme_multimodal_v1.sh
```

Verifique principalmente:

- `BASE`
- `TAB`
- `CONDA_ENV` (agora ambiente único; `PY_ENV`, `METHSCAN_ENV` e `R_ENV` ficam só como aliases compatíveis)
- `BISMARK_COVERAGE_AGG`
- `THREADS`, `MEM_GB`
- `LEIDEN_RESOLUTION`, `UMAP_NEIGHBORS`, `PCA_NPCS`

Nesta versão, o padrão já é:

```bash
ENABLE_SCNA=0
```

## Modos de execução

### Opção mais segura: disparador principal

```bash
bash run_master_scDNAme_multimodal_v1.sh all
```

Esse modo submete:

- MethSCAn `all`
- MethSCAn `IDHwt`
- MethSCAn `IDHmut`
- disorder downstream para cada modo, com dependência `afterok`

E **não** tenta rodar SCNA.

### Rodar só o MethSCAn

```bash
bash run_master_scDNAme_multimodal_v1.sh methscan
```

ou

```bash
sbatch run_methscan_syn22257780_v1.sh all
sbatch run_methscan_syn22257780_v1.sh IDHwt
sbatch run_methscan_syn22257780_v1.sh IDHmut
```

### Rodar só disorder downstream

Depois que os clusters do MethSCAn existirem:

```bash
bash run_master_scDNAme_multimodal_v1.sh disorder
```

ou

```bash
sbatch run_disorder_syn22257780_v1.sh all
sbatch run_disorder_syn22257780_v1.sh IDHwt
sbatch run_disorder_syn22257780_v1.sh IDHmut
```

### SCNA no futuro

Se no futuro você obtiver BAMs:

1. preencha `template_bam_manifest_scRRBS.tsv`
2. aponte `BAM_MANIFEST` na config
3. mude na config:

```bash
ENABLE_SCNA=1
```

4. use:

```bash
bash run_master_scDNAme_multimodal_v1.sh all_with_scna
```

ou

```bash
bash run_master_scDNAme_multimodal_v1.sh scna
```

## Saídas principais

### MethSCAn

Por modo, em:

```bash
$PROJECT_OUT_BASE/methscan_<MODE>/
```

- `VMRs.bed`
- `VMR_matrix/`
- `promoter_matrix/`
- `gene_matrix/`
- `downstream/cell_cluster_assignments.tsv`
- `downstream/UMAP_clusters.pdf`
- `downstream/cluster_by_case_fraction.tsv`
- `dmrs/*.bed`

### Disorder

Por modo, em:

```bash
$PROJECT_OUT_BASE/methscan_<MODE>/disorder_downstream/
```

- `context_disorder_by_cluster.tsv`
- `context_disorder_by_case.tsv`
- `promoter_disorder_by_cluster.tsv`
- `tfbs_disorder_by_cluster.tsv`
- `Boxplot_context_disorder_by_cluster.pdf`

## Observações importantes

1. O parser do MethSCAn já está compatível com o layout mostrado para:

```bash
analysis_scRRBS_bismark_coverage.txt.gz
```

com colunas:

```text
cell_barcode chr start end methylation_percentage count_methylated count_unmethylated
```

2. O módulo de disorder **não recalcula** PDR/discordância a partir de reads; ele usa as tabelas já processadas presentes em `tables/`.

3. Antes de rodar o MethSCAn em muitas células, confira o limite de arquivos abertos (`ulimit -n`), porque o `prepare` pode precisar de um valor mais alto.

4. A primeira rodada recomendada continua sendo `IDHwt`, para verificar mistura por caso e calibrar `LEIDEN_RESOLUTION`.
