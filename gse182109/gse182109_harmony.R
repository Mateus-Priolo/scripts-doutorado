library(Seurat)
library(harmony)
library(patchwork)
library(dplyr)
library(Matrix)
library(tidyverse)
library(RColorBrewer)

merged_objects <- subset(merged_objects, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
merged_objects <- NormalizeData(merged_objects)
merged_objects <- FindVariableFeatures(merged_objects, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(merged_objects)
merged_objects <- ScaleData(merged_objects, features = all.genes)
gc()

variable_features <- VariableFeatures(object = merged_objects)

immuneadjustmentlist<-c("IGKV4-1","IGHV3-30","IGLC1","IGLC2","IGLV6-57","IGHG2","IGHA1","IGHV4-61","IGHM","IGLV3-1",
                        "IGHA2","IGHG4","IGHG1","IGKC","IGHG3","IGHGP","IGKV3-20","IGLC3","CD79A","TRDC","TRAC")


merged_objects <- RunPCA(merged_objects,npcs = 50, features=c(variable_features,immuneadjustmentlist),verbose = FALSE)

merged_objects <- JackStraw(merged_objects, num.replicate = 100)
merged_objects <- ScoreJackStraw(merged_objects, dims = 1:40)
pdf(file='elbow.pdf',width =10,height = 8)
ElbowPlot(merged_objects,ndims=40)
dev.off()

#harmony
harmony_merged_objects<-RunHarmony(object=merged_objects,group.by.vars="Fragment",reduction="pca",reduction.save="harmony")

pdf(file='heatmapsPCS_harmony.pdf',width = 8,height = 10)
DimHeatmap(harmony_merged_objects, dims = 1:40, cells = 500, balanced = TRUE)
dev.off()
harmony_merged_objects <- JackStraw(harmony_merged_objects, num.replicate = 100)
merged_objects <- ScoreJackStraw(harmony_merged_objects, dims = 1:40)
JackStrawPlot(harmony_merged_objects, dims = 1:40)
pdf(file='elbow_harmony.pdf',width =10,height = 8)
ElbowPlot(harmony_merged_objects,ndims=40)
dev.off()

harmony_merged_objects <- RunUMAP(harmony_merged_objects,reduction = "harmony",dims = 1:29,reduction.name="UMAP_harmony",verbose = FALSE)
harmony_merged_objects <- FindNeighbors(harmony_merged_objects,reduction = "harmony", dims = 1:29)
harmony_merged_objects <- FindClusters(harmony_merged_objects,reduction = "harmony", resolution = 0.8)

#define axis limits
# Get UMAP coordinates from Seurat object
umap_coords <- as.data.frame(Embeddings(harmony_merged_objects, reduction = "UMAP_harmony"))
x_limits <- c(min(umap_coords$UMAPharmony_1), max(umap_coords$UMAPharmony_1))
y_limits <- c(min(umap_coords$UMAPharmony_2), max(umap_coords$UMAPharmony_2))
#rename clusters
harmony_merged_objects<-RenameIdents(harmony_merged_objects,levels(harmony_merged_objects))
# Create UMAP plot
p <- DimPlot(harmony_merged_objects, reduction = "UMAP_harmony",label = TRUE, pt.size = 0.5)
# Increase axis range
p <- p + scale_x_continuous(limits = x_limits)
p <- p + scale_y_continuous(limits = y_limits)

pdf(file = 'UMAP_harmony_fragment_clusters_immuneadjlist_29pcs_res0.8.pdf', width = 14, height = 12)
# Print the plot
print(p)
dev.off()

# Replace values in the "Type" column
harmony_merged_objects@meta.data$Type <- ifelse(
  harmony_merged_objects@meta.data$Type %in% c("Astrocytoma", "Oligodendroglioma"),
  "LGG", 
  harmony_merged_objects@meta.data$Type
)
# Check the unique values in the "Type" column
unique(harmony_merged_objects@meta.data$Type)
#reorder column level
harmony_merged_objects@meta.data$Type<-factor(x=harmony_merged_objects@meta.data$Type,levels=c("LGG","GBM","Recurrent GBM"))
# Create UMAP plot
p <- DimPlot(harmony_merged_objects, reduction = "UMAP_harmony",label = TRUE, pt.size = 0.5,split.by="Type")
# Increase axis range
p <- p + scale_x_continuous(limits = x_limits)
p <- p + scale_y_continuous(limits = y_limits)

pdf(file = 'UMAP_harmony_fragment_splittedType_immuneadjlist_29pcs_res0.8.pdf', width = 18, height = 10)
# Print the plot
print(p)
dev.off()

p <- DimPlot(harmony_merged_objects, reduction = "UMAP_harmony",group.by = 'Assignment',label = TRUE, pt.size = 0.01)
# Increase axis range
p <- p + scale_x_continuous(limits = x_limits)
p <- p + scale_y_continuous(limits = y_limits)
# Create PDF file with specified dimensions
pdf(file = 'UMAP_harmony_fragment_celltypes_immuneadjlist_29pcs_res0.8.pdf', width = 14, height = 12)
p <- p + ggtitle("Cell Type")
# Print the plot
print(p)
dev.off()

p <- DimPlot(harmony_merged_objects, reduction = "UMAP_harmony",group.by = 'Assignment',label = TRUE, pt.size = 0.01)
# Increase axis range
p <- p + scale_x_continuous(limits = x_limits)
p <- p + scale_y_continuous(limits = y_limits)
# Create PDF file with specified dimensions
pdf(file = 'UMAP_harmony_fragment_celltypes_immuneadjlist_29pcs_res0.8.pdf', width = 14, height = 12)
p <- p + ggtitle("Cell Type")
# Print the plot
print(p)
dev.off()

#phase
p <- DimPlot(harmony_merged_objects, reduction = "UMAP_harmony",group.by = 'Phase')
# Increase axis range
p <- p + scale_x_continuous(limits = x_limits)
p <- p + scale_y_continuous(limits = y_limits)
# Create PDF file with specified dimensions
pdf(file = 'UMAP_harmony_fragment_cellcycle_immuneadjlist_29pcs_res0.8.pdf', width = 14, height = 12)
p <- p + ggtitle("Cell Cycle")
# Print the plot
print(p)
dev.off()

#subcelltype
p <- DimPlot(harmony_merged_objects, reduction = "UMAP_harmony",group.by = 'SubAssignment',label = TRUE, pt.size = 0.01)
# Increase axis range
p <- p + scale_x_continuous(limits = x_limits)
p <- p + scale_y_continuous(limits = y_limits)
# Create PDF file with specified dimensions
pdf(file = 'UMAP_harmony_fragment_subcelltypes_immuneadjlist_29pcs_res0.8.pdf', width = 14, height = 12)
p <- p + ggtitle("Sub-Types")
# Print the plot
print(p)
dev.off()

#rename clusters
#merged_objects<-RenameIdents(merged_objects,levels(merged_objects))
# Create UMAP plot
p <- DimPlot(harmony_merged_objects, reduction = "UMAP_harmony",group.by = 'Fragment',label = TRUE, pt.size = 0.01)
# Increase axis range
p <- p + scale_x_continuous(limits = x_limits)
p <- p + scale_y_continuous(limits = y_limits)
# Create PDF file with specified dimensions
pdf(file = 'UMAP_harmony_fragment_biopsia_immuneadjlist_29pcs_res0.8.pdf', width = 14, height = 12)
p <- p + ggtitle("Biopsy")
# Print the plot
print(p)
dev.off()

#rename clusters
#merged_objects<-RenameIdents(merged_objects,levels(merged_objects))
# Create UMAP plot
p <- DimPlot(harmony_merged_objects, reduction = "UMAP_harmony",group.by = 'Patient',label = TRUE, pt.size = 0.01)
# Increase axis range
p <- p + scale_x_continuous(limits = x_limits)
p <- p + scale_y_continuous(limits = y_limits)
# Create PDF file with specified dimensions
pdf(file = 'UMAP_harmony_fragment_pacientes_immuneadjlist_29pcs_res0.8.pdf', width = 14, height = 12)
p <- p + ggtitle("Patients")
# Print the plot
print(p)
dev.off()

#rename clusters
#merged_objects<-RenameIdents(merged_objects,levels(merged_objects))
# Create UMAP plot
p <- DimPlot(harmony_merged_objects, reduction = "UMAP_harmony",group.by = 'Type',label = TRUE, pt.size = 0.01)
# Increase axis range
p <- p + scale_x_continuous(limits = x_limits)
p <- p + scale_y_continuous(limits = y_limits)
# Create PDF file with specified dimensions
pdf(file = 'UMAP_harmony_fragment_type_immuneadjlist_29pcs_res0.8.pdf', width = 14, height = 12)
p <- p + ggtitle("Tumor Type")
# Print the plot
print(p)
dev.off()

#Findmarkers
harmony_merged_objects<-JoinLayers(harmony_merged_objects)

#all markers of every cluster compared to all remaining cells
all_markers_everycluster<-FindAllMarkers(harmony_merged_objects,only.pos=TRUE) %>%
  dplyr::filter(avg_log2FC>1)

#feature plots for GSC markers
pdf(file = 'featureplots_GSCmarkers.pdf', width = 14, height = 6)
FeaturePlot(harmony_merged_objects, features = c("PROM1", "FUT4", "SOX2", "NANOG", "OLIG2", "NES", "CD9"),min.cutoff='q10' ,ncol = 4)
par(cex.axis = 0.4)
dev.off()

pdf(file = 'featureplots_ACmarkers.pdf', width = 50, height = 30)
FeaturePlot(harmony_merged_objects, features = c("CST3", "S100B", "SLC1A3", "HEPN1", "HOPX", "MT3", "SPARCL1", "MLC1", "GFAP", "FABP7", "BCAN", "PON2", 
                                                 "METTL7B", "SPARC", "GATM", "RAMP1", "PMP2", "AQP4", "DBI", "EDNRB", "PTPRZ1", "CLU", "PMP22", "ATP1A2", 
                                                 "S100A16", "HEY1", "PCDHGC3", "TTYH1", "NDRG2", "PRCP", "ATP1B2", "AGT", "PLTP", "GPM6B", "F3", "RAB31", 
                                                 "PLPP3", "ANXA5", "TSPAN7"),min.cutoff='q10' ,ncol = 10)
par(cex.axis = 0.2)
dev.off()

#Add module Score

AC<-list(c("CST3", "S100B", "SLC1A3", "HEPN1", "HOPX", "MT3", "SPARCL1", "MLC1", "GFAP", "FABP7", "BCAN", "PON2", 
           "METTL7B", "SPARC", "GATM", "RAMP1", "PMP2", "AQP4", "DBI", "EDNRB", "PTPRZ1", "CLU", "PMP22", "ATP1A2", 
           "S100A16", "HEY1", "PCDHGC3", "TTYH1", "NDRG2", "PRCP", "ATP1B2", "AGT", "PLTP", "GPM6B", "F3", "RAB31", 
           "PLPP3", "ANXA5", "TSPAN7"))
OPC<-list(c("BCAN", "PLP1", "GPR17", "FIBIN", "LHFPL3", "OLIG1", "PSAT1", "SCRG1", "OMG", "APOD", "SIRT2", "TNR",
            "THY1", "PHYHIPL", "SOX2-OT", "NKAIN4", "PLPPR1", "PTPRZ1", "VCAN", "DBI", "PMP2", "CNP", "TNS3", "LIMA1",
            "CA10", "PCDHGC3", "CNTN1", "SCD5", "P2RX7", "CADM2", "TTYH1", "FGF12", "TMEM206", "NEU4", "FXYD6", "RNF13",
            "RTKN", "GPM6B", "LMF1", "ALCAM", "PGRMC1", "HRASLS", "BCAS1", "RAB31", "PLLP", "FABP5", "NLGN3", "SERINC5",
            "EPB41L2", "GPR37L1"))
NPC<-list(c("DLL3", "DLL1", "SOX4", "TUBB3", "HES6", "TAGLN3", "NEU4", "MARCKSL1", "CD24", "STMN1", "TCF12", "BEX1",
            "OLIG1", "MAP2", "FXYD6", "PTPRS", "MLLT11", "NPPA", "BCAN", "MEST", "ASCL1", "BTG2", "DCX", "NXPH1",
            "JPT1", "PFN2", "SCG3", "MYT1", "CHD7", "ADGRG1", "TUBA1A", "PCBP4", "ETV1", "SHD", "TNR", "AMOTL2", "DBN1",
            "HIP1", "ABAT", "ELAVL4", "LMF1", "GRIK2", "SERINC5", "TSPAN13", "ELMO1", "GLCCI1", "SEZ6L", "LRRN1", "SEZ6",
            "SOX11"))
MES<-list(c("CHI3L1", "ANXA2", "ANXA1", "CD44", "VIM", "MT2A", "C1S", "NAMPT", "EFEMP1", "C1R", "SOD2", "IFITM3", 
            "TIMP1", "SPP1", "A2M", "S100A11", "MT1X", "S100A10", "FN1", "LGALS1", "S100A16", "CLIC1", "MGST1", 
            "RCAN1", "TAGLN2", "NPC2", "SERPING1", "TCIM", "EMP1", "APOE", "CTSB", "C3", "LGALS3", "MT1E", "EMP3",
            "SERPINA3", "ACTN1", "PRDX6", "IGFBP7", "SERPINE1", "PLP2", "MGP", "CLIC4", "GFPT2", "GSN", "NNMT", "TUBA1C",
            "GJA1", "TNFRSF1A", "WWTR1"))
GSC<-list(c("PROM1", "FUT4", "SOX2", "NANOG", "OLIG2", "NES", "CD9"))

harmony_merged_objects<-JoinLayers(harmony_merged_objects)
harmony_merged_objects <- AddModuleScore(object = harmony_merged_objects, 
                                         features = GSC, 
                                         name = "GSC_score")
harmony_merged_objects <- AddModuleScore(object = harmony_merged_objects, 
                                         features = AC, 
                                         name = "AC_score")
harmony_merged_objects <- AddModuleScore(object = harmony_merged_objects, 
                                         features = OPC, 
                                         name = "OPC_score")
harmony_merged_objects <- AddModuleScore(object = harmony_merged_objects, 
                                         features = NPC, 
                                         name = "NPC_score")
harmony_merged_objects <- AddModuleScore(object = harmony_merged_objects, 
                                         features = MES, 
                                         name = "MES_score")

pdf(file = 'featureplots_AC.pdf')
FeaturePlot(harmony_merged_objects,features = "AC_score1", label = TRUE, repel = TRUE)+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()
pdf(file = 'featureplots_OPC.pdf')
FeaturePlot(harmony_merged_objects,features = "OPC_score1", label = TRUE, repel = TRUE)+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()
pdf(file = 'featureplots_NPC.pdf')
FeaturePlot(harmony_merged_objects,features = "NPC_score1", label = TRUE, repel = TRUE)+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()
pdf(file = 'featureplots_MES.pdf')
FeaturePlot(harmony_merged_objects,features = "MES_score1", label = TRUE, repel = TRUE)+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()
pdf(file = 'featureplots_GSC.pdf')
FeaturePlot(harmony_merged_objects,features = "GSC_score1", label = TRUE, repel = TRUE)+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()

harmony_merged_objects$gliomastatescore <- apply(harmony_merged_objects@meta.data[, c("GSC_score1", "AC_score1", "OPC_score1", "NPC_score1", "MES_score1")], 1, function(x) {
  scores <- c(GSC = x["GSC_score1"], AC = x["AC_score1"], OPC = x["OPC_score1"], NPC = x["NPC_score1"], MES = x["MES_score1"])
  return(names(scores)[which.max(scores)])
})
pdf(file='cell_states_modulescore.pdf')
DimPlot(harmony_merged_objects, reduction = "UMAP_harmony",label=FALSE,group.by ="gliomastatescore")
dev.off()
