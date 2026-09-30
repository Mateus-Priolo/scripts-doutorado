Integração de Datasets scRNA-seq — GBM/LGG
Pipeline de integração de múltiplos datasets de single-cell RNA-seq de gliomas (GSE182109, Synapse syn22257780 e GSE173278), com QC, integração por Harmony, anotação automática, label transfer e análise downstream por subtipo IDH.

Visão geral
Este diretório reúne os scripts SLURM que orquestram a integração de datasets heterogêneos de glioma, unificando amostras de GBM primário, GBM recorrente e LGG provenientes de três fontes distintas. O pipeline foi desenhado para:

Carregar e padronizar metadados de múltiplas coortes

Aplicar QC rigoroso e SCTransform por amostra

Integrar por Harmony corrigindo efeito de lote por paciente/dataset

Anotar tipos celulares automaticamente por módulos de genes

Validar anotações via label transfer a partir de um dataset externo de referência

Separar objetos finais por status IDH (mutante vs wild-type)

Fluxo do pipeline
text
submit_pipeline.sh
    |
    +--> 01b_load_gse173278_primary.sh
    |        carrega raw: GSE182109 + Synapse + GSE173278
    |
    +--> 02b_prepare_gse173278_primary.sh
    |        QC + doublets + SCTransform por amostra
    |
    +--> 03b_harmony_integration_3datasets.sh
    |        Harmony (apenas GSE182109 + Synapse)
    |
    +--> 04b_downstream_analysis_3datasets.sh
             marcadores + anotacao + label transfer + IDH split

explore_seurat_objects.R + run_explore_objects.sh   (utilitario de inspecao)
scrna_integrated_objects.R                          (integracao antiga, baseada em Patient)
Descrição dos scripts
1. submit_pipeline.sh — Orquestrador SLURM
O que faz:

Submete os quatro jobs principais em sequência usando sbatch --dependency=afterok

Cada job só inicia após o anterior concluir com sucesso

Imprime os IDs dos jobs para rastreamento

Como usar:

text
bash submit_pipeline.sh
Saída: encadeamento de jobs 01b -> 02b -> 03b -> 04b no SLURM.

2. 01b_load_gse173278_primary.sh — Carregamento do raw
O que faz:

Copia o gse_raw.rds já validado do pipeline anterior para results_v3/

Reconstrói o objeto Synapse a partir do arquivo .h5ad exportado, corrigindo a coluna Patient (mapeia sampleid -> SM001...SM019 por contagem de células esperada)

Anexa metadados clínicos do artigo (MOESM2 xlsx)

Carrega o GSE173278 primary GBM (matriz normalizada pelo provedor) e filtra apenas células de tecido primário (exclui gliomasferas, explantes e recorrentes)

Aplica downsampling para 30.000 células antes de criar o objeto Seurat

Marca GSE173278 com already_normalized = TRUE (não será normalizado de novo)

Saída: gse_raw.rds, synapse_raw.rds, gse173278_primary_raw.rds

Observação: o GSE173278 NÃO será integrado via Harmony — serve apenas como referência para label transfer no passo 04.

3. 02b_prepare_gse173278_primary.sh — QC e normalização
O que faz:

Para GSE182109 e Synapse: roda QC completo com:

Cálculo de métricas (nFeature, nCount, percent.mt, percent.rb)

Filtros: nFeature_RNA > 300 & < 8000-9000 e percent.mt < 20-25%

Detecção de doublets com scDblFinder (por amostra)

SCTransform v2 com regressão de percent.mt e 3000 HVGs

Para GSE173278: QC leve (apenas filtro de nFeature > 200), mantendo a normalização original do provedor e calculando HVGs no slot data

Gera violinos de QC antes da filtragem para cada dataset

Saída: gse_qc.rds, synapse_qc.rds, gse173278_primary_qc.rds + PDFs de QC

4. 03b_harmony_integration_3datasets.sh — Integração por Harmony
O que faz:

Reconstrói objetos mínimos a partir de gse_qc.rds e synapse_qc.rds, padronizando nomes de genes (uppercase, remove sufixos .1, converte _ para -)

Normaliza cada dataset individualmente e seleciona HVGs

Seleciona features compartilhadas via SelectIntegrationFeatures (fallback: união dos HVGs, depois interseção de genes)

Faz merge dos objetos e define variáveis a regredir (percent.mt, nCount_RNA se presentes)

Cria coluna batch = dataset_Patient e roda Harmony corrigindo por esse batch

Gera múltiplas resoluções de clustering (0.2, 0.4, 0.6, 0.8, 1.0)

Define RNA_snn_res.0.6 como identidade ativa

Roda UMAP com dims = 1:30

Saída: integrated_harmony_v3.rds e pca_elbow_v3.pdf

Observação: o GSE173278 NÃO entra aqui — apenas GSE182109 + Synapse.

5. 04b_downstream_analysis_3datasets.sh — Análise downstream
O que faz:

Carrega o objeto integrado e padroniza metadados

Atribui IDH_status (IDH-mut / IDH-wt) combinando:

Mapa manual para amostras Synapse (SM001...SM019)

Tipo tumoral (LGG -> IDH-mut, GBM -> IDH-wt)

Gera UMAPs coloridos por cluster, dataset, paciente e status IDH

Roda FindAllMarkers por cluster e salva tabela de marcadores

Anotação automática por módulos de genes (Glioma, Myeloid, Tcells, Bcells, Endothelial, Pericytes, Oligodendrocytes) usando AddModuleScore

Define predicted_celltype por célula como o módulo de maior score

Gera tabela de composição IDH por cluster

Cria objetos separados: idh_mut_integrated_v3.rds e idh_wt_integrated_v3.rds

Roda label transfer do GSE173278 (query) para o objeto integrado (reference):

FindTransferAnchors + TransferData

Transfere predicted_celltype e cluster

Salva tabela de predições e UMAP do label transfer

Fallback: UMAP próprio do GSE173278 colorido por label transferido

Saída: annotated_integrated_v3.rds, tabelas de marcadores, UMAPs, objetos IDH-mut/IDH-wt, label_transfer_g173_v3.csv

6. explore_seurat_objects.R — Utilitário de inspeção
O que faz:

Script R independente que recebe dois caminhos de objetos Seurat via argumentos

Para cada objeto, gera um relatório com:

Dimensões (n células, n genes)

Assays e reduções disponíveis

Lista completa de colunas em meta.data

Sumário detalhado de cada coluna (numérica: min/Q1/mediana/média/Q3/max/NA; categórica: tabela de níveis)

Detecção automática de colunas de cluster (res. ou cluster no nome)

Tabela de identidades ativas

Redireciona toda a saída para um arquivo de log e também para o console

Uso:

text
Rscript explore_seurat_objects.R --object1 obj1.rds --object2 obj2.rds --outdir saida/
Saída: exploration_summary.txt com o relatório completo.

7. run_explore_objects.sh — Wrapper SLURM do utilitário
O que faz:

Submete o explore_seurat_objects.R ao SLURM com 64 GB de memória

Ativa o ambiente conda R

Aponta para os objetos annotated_integrated_v3.rds e integrated_harmony_v3.rds

Cria o diretório de saída exploration_output/

Saída: exploration_summary.txt dentro de exploration_output/.

8. scrna_integrated_objects.R — Integração legada (baseada em Patient)
O que faz:

Versão inicial da integração, mais simples, que corrige apenas por Patient (sem dataset)

Usa uma lista de ajuste imunológico (IGKV4-1, IGHV3-30, etc.) e estromal (DCN, FBLN1, LUM, etc.) adicionada aos genes variáveis na PCA

Harmony com group.by.vars = "Patient"

Resolução de clustering fixada em 0.8021

Gera UMAPs por fragmento, paciente, tipo tumoral e subtipo

Saída: vários PDFs de UMAP (UMAP_harmony_*.pdf).

Observação: foi substituído pelo pipeline 01b-04b. Mantido para referência e comparação.

Dependências
Categoria	Pacotes
Single-cell	Seurat, SeuratObject, harmony, scDblFinder, SingleCellExperiment, BiocParallel
Utilitários	dplyr, tidyverse, data.table, Matrix, patchwork, ggplot2
Ambiente	miniconda com ambiente gbm_scrnaseq (jobs) e R (exploração)
Arquivos de entrada esperados
Arquivo	Descrição
gse_raw.rds	Objeto Seurat do GSE182109 já validado (do pipeline anterior)
analysis_scRNAseq_tumor_counts.h5ad	Contagens do Synapse (syn22257780)
41588_2021_926_MOESM2_ESM.xlsx	Metadados clínicos do artigo do Synapse
GSE173278_scRNAseq_filtered_cells_*	Barcodes, genes, metadados e matriz do GSE173278
listagenespaloma.txt	(legado) lista curada de 31 genes
Saídas principais
Arquivo	Conteúdo
gse_raw.rds	Objeto Seurat bruto do GSE182109
synapse_raw.rds	Objeto Seurat bruto do Synapse
gse173278_primary_raw.rds	Objeto Seurat bruto do GSE173278
gse_qc.rds	GSE182109 após QC + SCTransform
synapse_qc.rds	Synapse após QC + SCTransform
gse173278_primary_qc.rds	GSE173278 após QC leve
integrated_harmony_v3.rds	Objeto integrado por Harmony (GSE182109 + Synapse)
annotated_integrated_v3.rds	Objeto com anotação automática e label transfer
idh_mut_integrated_v3.rds	Subset IDH-mutante
idh_wt_integrated_v3.rds	Subset IDH-wild-type
markers_all_clusters_v3.csv	Marcadores diferenciais por cluster
cluster_annotation_auto_v3.csv	Anotação automática por cluster
label_transfer_g173_v3.csv	Predições de tipo celular transferidas do GSE173278
umap_clusters_v3.pdf	UMAP colorido por cluster
umap_dataset_v3.pdf	UMAP colorido por dataset
umap_patient_v3.pdf	UMAP colorido por paciente
umap_idh_v3.pdf	UMAP colorido por status IDH
umap_label_transfer_v3.pdf	UMAP do label transfer
idh_composition_by_cluster_v3.pdf	Composição IDH por cluster
pca_elbow_v3.pdf	Elbow plot da PCA integrada
exploration_summary.txt	Relatório do utilitário de inspeção
Observações importantes
O pipeline está dividido em 4 jobs SLURM encadeados. Sempre submeta via submit_pipeline.sh para garantir a ordem correta.

O GSE173278 NÃO entra na integração Harmony — é usado apenas como referência externa para label transfer no passo 04.

A correção de batch é feita por dataset_Patient (não apenas Patient), o que evita que amostras do mesmo paciente em datasets diferentes sejam colapsadas.

O SCTransform é aplicado por amostra antes do merge, evitando misturar distribuições heterogêneas.

O 03b reconstrói objetos mínimos a partir dos QC, padronizando nomes de genes (uppercase, remove .N, converte _ para -). Isso é essencial para o merge funcionar.

O label transfer usa predicted_celltype do objeto integrado como referência — a qualidade da anotação automática afeta diretamente a transferência.

Caminhos absolutos nos scripts apontam para /home/renanomete/projetos/matdata/met_mat_data/singlecell_novo. Ajuste conforme seu ambiente.

O scrna_integrated_objects.R é legado (integração apenas por Patient) e foi substituído pelo pipeline 01b-04b. Use apenas como referência histórica.

