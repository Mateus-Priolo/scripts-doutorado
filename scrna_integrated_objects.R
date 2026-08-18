suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(RColorBrewer)
})

merged_objects <- NormalizeData(merged_objects)
merged_objects <- FindVariableFeatures(merged_objects, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(merged_objects)
merged_objects <- ScaleData(merged_objects, features = all.genes)
gc()

variable_features <- VariableFeatures(object = merged_objects)
immuneadjustmentlist<-c("IGKV4-1","IGHV3-30","IGLC1","IGLC2","IGLV6-57","IGHG2","IGHA1","IGHV4-61","IGHM","IGLV3-1",
                        "IGHA2","IGHG4","IGHG1","IGKC","IGHG3","IGHGP","IGKV3-20","IGLC3","CD79A","TRDC","TRAC")
stromaladjustmentlist <- c("DCN", "FBLN1", "LUM","PDGFRB", "RGS5", "ACTA2", "MYH11", "NOTCH3")

merged_objects <- RunPCA(merged_objects,npcs = 50, features=c(variable_features,immuneadjustmentlist,stromaladjustmentlist),verbose = FALSE)
merged_objects <- JackStraw(merged_objects, num.replicate = 100)
merged_objects <- ScoreJackStraw(merged_objects, dims = 1:40)
pdf(file='elbow.pdf',width =10,height = 8)
ElbowPlot(merged_objects,ndims=40)
dev.off()

#harmony
harmony_merged_objects<-RunHarmony(object=merged_objects,group.by.vars="Patient",reduction="pca",reduction.save="harmony")

pdf(file='heatmapsPCS_harmony.pdf',width = 8,height = 10)
DimHeatmap(harmony_merged_objects, dims = 1:40, cells = 500, balanced = TRUE)
dev.off()
harmony_merged_objects <- JackStraw(harmony_merged_objects, num.replicate = 100)
merged_objects <- ScoreJackStraw(harmony_merged_objects, dims = 1:40)
JackStrawPlot(harmony_merged_objects, dims = 1:40)
pdf(file='elbow_harmony.pdf',width =10,height = 8)
ElbowPlot(harmony_merged_objects,ndims=40)
dev.off()

harmony_merged_objects <- RunUMAP(harmony_merged_objects,reduction = "harmony",dims = 1:28,reduction.name="UMAP_harmony",verbose = FALSE)
harmony_merged_objects <- FindNeighbors(harmony_merged_objects,reduction = "harmony", dims = 1:28)
harmony_merged_objects <- FindClusters(harmony_merged_objects,reduction = "harmony", resolution = 0.8021)

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

pdf(file = 'UMAP_harmony_fragment_clusters_immuneadjlist_28pcs_res0.8.pdf', width = 14, height = 12)
# Print the plot
print(p)
dev.off()


harmony_merged_objects@meta.data$Type<-factor(x=harmony_merged_objects@meta.data$Type,levels=c("LGG","GBM","Recurrent GBM"))
# Create UMAP plot
p <- DimPlot(harmony_merged_objects, reduction = "UMAP_harmony",label = TRUE, pt.size = 0.5,split.by="Type")
# Increase axis range
p <- p + scale_x_continuous(limits = x_limits)
p <- p + scale_y_continuous(limits = y_limits)

pdf(file = 'UMAP_harmony_fragment_splittedType_immuneadjlist_28pcs_res0.8.pdf', width = 18, height = 10)
# Print the plot
print(p)
dev.off()

p <- DimPlot(harmony_merged_objects, reduction = "UMAP_harmony",group.by = 'Assignment',label = TRUE, pt.size = 0.01)
# Increase axis range
p <- p + scale_x_continuous(limits = x_limits)
p <- p + scale_y_continuous(limits = y_limits)
# Create PDF file with specified dimensions
pdf(file = 'UMAP_harmony_fragment_celltypes_immuneadjlist_28pcs_res0.8.pdf', width = 14, height = 12)
p <- p + ggtitle("Cell Type")
# Print the plot
print(p)
dev.off()

p <- DimPlot(harmony_merged_objects, reduction = "UMAP_harmony",group.by = 'Assignment',label = TRUE, pt.size = 0.01)
# Increase axis range
p <- p + scale_x_continuous(limits = x_limits)
p <- p + scale_y_continuous(limits = y_limits)
# Create PDF file with specified dimensions
pdf(file = 'UMAP_harmony_fragment_celltypes_immuneadjlist_28pcs_res0.8.pdf', width = 14, height = 12)
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
pdf(file = 'UMAP_harmony_fragment_cellcycle_immuneadjlist_28pcs_res0.8.pdf', width = 14, height = 12)
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
pdf(file = 'UMAP_harmony_fragment_subcelltypes_immuneadjlist_28pcs_res0.8.pdf', width = 14, height = 12)
p <- p + ggtitle("Sub-Types")
# Print the plot
print(p)
dev.off()

#

#


#rename clusters
#merged_objects<-RenameIdents(merged_objects,levels(merged_objects))
# Create UMAP plot
p <- DimPlot(harmony_merged_objects, reduction = "UMAP_harmony",group.by = 'Fragment',label = TRUE, pt.size = 0.01)
# Increase axis range
p <- p + scale_x_continuous(limits = x_limits)
p <- p + scale_y_continuous(limits = y_limits)
# Create PDF file with specified dimensions
pdf(file = 'UMAP_harmony_fragment_biopsia_immuneadjlist_28pcs_res0.8.pdf', width = 14, height = 12)
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
pdf(file = 'UMAP_harmony_fragment_pacientes_immuneadjlist_28pcs_res0.8.pdf', width = 14, height = 12)
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
pdf(file = 'UMAP_harmony_fragment_type_immuneadjlist_28pcs_res0.8.pdf', width = 14, height = 12)
p <- p + ggtitle("Tumor Type")
# Print the plot
print(p)
dev.off()
