📊 GSE182109 — Análise scRNA-seq de Gliomas

Pipeline de análise de single-cell RNA-seq do dataset GSE182109 — atlas público de gliomas humanos incluindo GBM primário, GBM recorrente (rGBM) e LGG (astrocytoma e oligodendroglioma).

📖 Visão geral

Este diretório reúne os scripts usados para processar e analisar o dataset GSE182109, cobrindo desde o dado bruto (matrizes 10X do Cell Ranger) até a caracterização de estados celulares e enriquecimento funcional de vias moleculares.

| Categoria | Pacotes |
|---|---|
| Single-cell | Seurat, harmony, dplyr, patchwork, Matrix, tidyverse, RColorBrewer, gridExtra, grid, ggplot2, reshape2 |
| Bulk RNA-seq | DESeq2, circlize, ComplexHeatmap, org.Hs.eg.db |
| Enriquecimento funcional | clusterProfiler, enrichplot, GSEABase, org.Hs.eg.db, pathview, biomaRt, pheatmap, ggupset, tidyverse |

## Arquivos de entrada esperados

| Arquivo | Descrição |
|---|---|
| Meta_Data_GBMatlas.txt | Metadados clínicos das amostras do GSE182109 |
| <GSM_*>/* | Diretórios de saída do Cell Ranger (um por amostra) |
| merged_counts_clean.txt | Matriz de contagens bulk (para DESeq2) |
| samples.txt | Tabela de amostras bulk (para DESeq2) |

## Saídas principais

| Arquivo | Conteúdo |
|---|---|
| merged_objects.rds | Objeto Seurat mesclado e filtrado |
| clustered_harmony_merged_objects.rds | Objeto integrado com clusters e UMAP |
| markers_cluster*_vs_*.rds | Marcadores diferenciais por cluster |
| UMAP_harmony_*.pdf | Visualizações UMAP (cluster, tipo, paciente, fase) |
| featureplots_*.pdf | FeaturePlots de assinaturas e marcadores |
| cell_states_modulescore.pdf | Estados tumorais por célula |
| cluster_pie_charts.pdf | Composição de estados por cluster/tipo |
| enriched_kegg_*.tiff, enriched_GO_*.tiff | Enriquecimento KEGG/GO |
| gseplot_*.tiff, ridgeplot_*.tiff | GSEA |
| filtered_log2FoldChange.txt, foldchanges31.txt | Resultados do DESeq2 |

⚠️ Observações
Caminhos absolutos (C:/Users/... e /data1/projects/...) estão hardcoded nos scripts. Ajuste conforme o ambiente.

📜 Descrição dos scripts
<details open> <summary><b>1. <code>gse182109.R</code> — Carregamento e montagem do objeto Seurat</b></summary>
O que faz:

Lê os metadados das amostras (Meta_Data_GBMatlas.txt) e ajusta os nomes das colunas para casar com os barcodes das células

Percorre todos os diretórios de amostras (um por biblioteca 10X), lê a matriz de contagens com Read10X() e cria um objeto Seurat individual para cada uma, aplicando filtros iniciais de qualidade:

nFeature_RNA > 200 — células com pelo menos 200 genes detectados

nFeature_RNA < 2500 — remove dupletos e células com muita complexidade

percent.mt < 5% — remove células mortas/apoptóticas

Mescla todas as amostras em um único objeto (merged_objects), prefixando os nomes das células com o ID da amostra

Adiciona os metadados clínicos ao objeto (tipo tumoral, paciente, fragmento, etc.)

Salva o resultado como merged_objects.rds

📤 Saída: merged_objects.rds (objeto Seurat não normalizado, mas já filtrado)

</details><details> <summary><b>2. <code>gse182109_harmony.R</code> — Integração, clusterização e anotação</b></summary>
O que faz:

Pré-processamento: normalização (NormalizeData), seleção de genes variáveis (FindVariableFeatures), escalonamento (ScaleData)

PCA com ajuste imunológico: roda RunPCA combinando os 2000 genes variáveis com uma lista personalizada de genes imunes (IGKV4-1, IGHV3-30, IGLC1, etc.). O objetivo é forçar a PCA a capturar também o sinal de células imunes, evitando que elas sejam agrupadas apenas por ruído técnico

Seleção de dimensões: usa JackStraw e ElbowPlot para escolher quantos PCs usar (29 no final)

Integração por Harmony: corrige o efeito de lote entre fragmentos/amostras (group.by.vars = "Fragment")

Clusterização: FindNeighbors + FindClusters (resolução 0.8) e UMAP

Visualizações: gera múltiplos PDFs mostrando os clusters coloridos por fragmento, paciente, tipo tumoral, tipo celular (Assignment), subtipo (SubAssignment) e fase do ciclo celular

Anotação por estados de glioma: define assinaturas gênicas para cinco estados tumorais (AC, OPC, NPC, MES, GSC) e calcula um Module Score para cada célula com AddModuleScore(). Cada célula recebe o estado com maior score (gliomastatescore)

Marcadores diferenciais: roda FindAllMarkers para todos os clusters

Salva o objeto integrado como clustered_harmony_merged_objects.rds

📤 Saída: objeto Seurat integrado com clusters, UMAP, scores de estados e vários PDFs (UMAP_harmony_*.pdf, featureplots_*.pdf, cell_states_modulescore.pdf)

</details><details> <summary><b>3. <code>denovoclustering.R</code> — Re-clusterização focada em glioma + oligodendrócitos</b></summary>
O que faz:

Recarrega o objeto integrado e seleciona apenas as células dos clusters tumorais e oligodendrocíticos (clusters 2, 6, 7, 8, 11, 12, 21 na anotação original)

Reconverte o assay para o formato v3 (compatível com JoinLayers de versões antigas do Seurat)

Refaz todo o pipeline só nessas células:

Normalização, genes variáveis, escalonamento

PCA → Harmony corrigindo por Fragment e Patient simultaneamente

UMAP + clustering (resolução 0.2, mais grosseira)

Anota estados tumorais com as mesmas assinaturas (AC, OPC, NPC, MES, GSC) e calcula o score máximo por célula

Identifica marcadores diferenciais entre clusters tumorais e entre subtipos de tumor (Recurrent GBM vs GBM vs LGG)

Compara marcadores entre pacientes para encontrar genes conservados (FindConservedMarkers)

Gera FeaturePlots dos top 10 e bottom 10 genes para GBM e rGBM, além de heatmaps e plots de módulo

🎯 Objetivo: refinar a heterogeneidade dentro do compartimento tumoral, removendo a influência das células imunes/estromais no agrupamento

📤 Saída: objeto glioma_oligo_cells com anotação refinada e vários PDFs (gliomaoligo_*.pdf)

</details><details> <summary><b>4. <code>cellstates_charts.R</code> — Visualização dos estados celulares</b></summary>
O que faz:

A partir do objeto integrado (harmony_merged_objects), extrai os scores dos cinco estados tumorais (MES, AC, OPC, GSC, NPC) por célula

Calcula, para cada combinação cluster × tipo tumoral, a média dos scores

Gera gráficos de pizza mostrando a composição de estados de cada cluster em cada tipo de tumor (LGG, GBM, rGBM)

Salva tudo em cluster_pie_charts.pdf

🎯 Objetivo: responder visualmente "qual estado de glioma predomina em cada cluster, em cada tipo de tumor?"

</details><details> <summary><b>5. <code>repairgenes_findmarkers_clusterprofiler.Rmd</code> — Genes de reparo de DNA + enriquecimento</b></summary>
O que faz:

Define listas curadas de genes de reparo de DNA por via:

MMR (mismatch repair)

NER (nucleotide excision repair)

BER (base excision repair)

NHEJ (non-homologous end joining)

HR (homologous recombination)

Fanconi anemia pathway

Combina com listas do Hallmarks (GSEA) e de um artigo específico para gerar três conjuntos:

allrepairgenes

allrepairgenes_plus_hallmarksgsea

allrepairgenes_plus_hallmarksgsea_plus_artigo

Calcula Module Score de reparo (Repair_score) e de proliferação (Proliferation_score) por célula

Gera FeaturePlots mostrando a distribuição desses scores no UMAP e por tipo tumoral

DotPlots por via de reparo (NHEJ, MMR, HR, NER, BER, Fanconi), mostrando expressão em clusters tumorais específicos

Reatribui identidades celulares com base nos marcadores (microglia, mac1, mac2, p-mac, DC, MDSC, NK, CD4, prolif T, CD8, B, glioma, oligodendrócitos, endotélio, pericitos)

Roda FindMarkers para todos os pares de clusters relevantes (glioma, T cells, mieloides) e salva cada resultado em .rds

Enriquecimento funcional dos marcadores do cluster 8:

KEGG (enrichKEGG, gseKEGG) — vias enriquecidas e GSEA

GO (enrichGO — BP, CC, MF) — termos enriquecidos

Pathview — mapeia fold changes em vias KEGG visualmente

Gráficos: barplot, dotplot, emapplot, heatplot, ridgeplot, gseaplot2

📤 Saída: centenas de arquivos .rds, .txt, .tiff e .pdf com marcadores, enriquecimentos e visualizações

</details><details> <summary><b>6. <code>DESeq2pipeline_Mateus.R</code> — Análise de bulk RNA-seq com DESeq2</b></summary>
⚠️ Este script NÃO é de single-cell. É uma análise paralela de bulk RNA-seq (contagens agregadas por amostra).

O que faz:

Lê a matriz de contagens (merged_counts_clean.txt) e a tabela de amostras (samples.txt)

Cria um objeto DESeqDataSet com design ~ cell e filtra genes com baixa expressão (rowSums >= 10)

Define ACBRI371 (linha endotelial) como referência e roda DESeq2

Para cada tipo celular, calcula o contraste contra ACBRI371 via results(dds, contrast=...)

Combina todos os contrastes em uma tabela única (combined_results) e exporta apenas as colunas de log2FoldChange

Anota os genes com símbolos via org.Hs.eg.db (chave = ENSEMBL)

Filtra uma lista específica de genes e gera um heatmap com ComplexHeatmap (escala azul-branco-vermelho, −10 a +10)

📤 Saída: filtered_log2FoldChange.txt, foldchanges31.txt e um heatmap de fold changes

</details><details> <summary><b>7. <code>seuratutorial.R</code> — Tutorial Seurat (script de aprendizado)</b></summary>
O que faz:

Reproduz o tutorial oficial do Seurat com o dataset PBMC 3k (células mononucleares de sangue periférico)

Cobre: criação do objeto, QC, normalização, PCA, JackStraw, clustering, UMAP, FindMarkers, FindAllMarkers, ROC test, DoHeatmap e anotação manual

ℹ️ Não faz parte do pipeline de glioma. Está aqui apenas como referência/aprendizado.

<div align="center">
