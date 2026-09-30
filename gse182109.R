library(Seurat)
library(patchwork)
library(dplyr)
library(Matrix)
library(tidyverse)
library(harmony)
library(RColorBrewer)

#setar diretório e ler os metadados
setwd("/data1/projects/LGMB-005/NGS/singlecell/rawdata/GSE182109")
metadata <- read.table("Meta_Data_GBMatlas.txt", header = TRUE, sep = "\t")
metadata$Fragment <- gsub("-", "_", metadata$Fragment)
metadata$combined<-paste(metadata$GSMID, metadata$Fragment, metadata$barcode, sep = "_")
rownames(metadata)<-metadata$combined
metadata$combined<-NULL
View(metadata)

#lendo os arquivos e criando os objetos seurat
parent_dir<-"/data1/projects/LGMB-005/NGS/singlecell/rawdata/GSE182109"
dirs<-list.dirs(path=parent_dir,recursive=F,full.names=F)
for (dataset_dir in dirs) {
  data <- Read10X(data.dir = dataset_dir)
  seurat_obj <- CreateSeuratObject(counts = data, project = basename(dataset_dir), 
                                   min.cells = 5, min.features = 200) %>%
    PercentageFeatureSet(pattern = "^MT-", col.name = "percent.mt") %>% 
    subset(subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
  #Use 'assign' to create separate Seurat objects with names based on directory
  object_name <- basename(dataset_dir)
  assign(object_name, seurat_obj)
}
#combina os objetos em um só
merged_objects<-merge(GSM5518596_rGBM_01_A,y=c(GSM5518597_rGBM_01_B,GSM5518598_rGBM_01_C,GSM5518599_rGBM_01_D,GSM5518600_ndGBM_01_A,GSM5518601_ndGBM_01_C,GSM5518602_ndGBM_01_D,GSM5518603_ndGBM_01_F,GSM5518604_ndGBM_11_A,GSM5518605_ndGBM_11_B,GSM5518606_ndGBM_11_C,GSM5518607_ndGBM_11_D,GSM5518608_ndGBM_02_1,GSM5518609_ndGBM_02_2,GSM5518610_ndGBM_02_4,GSM5518611_ndGBM_02_5,GSM5518612_rGBM_02_2,GSM5518613_rGBM_02_3,GSM5518614_rGBM_02_4,GSM5518615_rGBM_02_5,GSM5518616_rGBM_03_1,GSM5518617_rGBM_03_2,GSM5518618_rGBM_03_3,GSM5518619_rGBM_04_1,GSM5518620_rGBM_04_2,GSM5518621_rGBM_04_3,GSM5518622_rGBM_04_4,GSM5518623_ndGBM_03_1,GSM5518624_ndGBM_03_2,GSM5518625_ndGBM_03_3,GSM5518626_rGBM_05_1,GSM5518627_rGBM_05_2,GSM5518628_rGBM_05_3,GSM5518629_ndGBM_10,GSM5518630_LGG_04_1,GSM5518631_LGG_04_2,GSM5518632_LGG_04_3,GSM5518633_ndGBM_04,GSM5518634_ndGBM_05,GSM5518635_ndGBM_06,GSM5518636_ndGBM_07,GSM5518637_ndGBM_08,GSM5518638_LGG_03,GSM5518639_ndGBM_09), add.cell.ids=ls()[4:47])
rm(list=ls()[4:47])
common_cell_names <- intersect(colnames(merged_objects), rownames(metadata))
merged_objects <- merged_objects[, common_cell_names]

# Adiciona as colunas do metadata ao objeto seurat
merged_objects <- AddMetaData(object = merged_objects, metadata = metadata)

saveRDS(merged_objects, file = "merged_objects.rds")