merged_objects<-readRDS("merged_objects.rds")
harmony_merged_objects<-readRDS("clustered_harmony_merged_objects.rds")

# convert a v5 assay to a v3 assay
merged_objects[["RNA"]] <- as(object = merged_objects[["RNA"]], Class = "Assay")

specified_clusters <- c(2,6,7,8,11,12,21)
cells_in_specified_clusters <- WhichCells(harmony_merged_objects, idents=specified_clusters)
# Subset merged_objects using the identified cells
glioma_oligo_cells <-subset(merged_objects, cells = cells_in_specified_clusters)

#layersList <-lapply(glioma_oligo_cells@assays$RNA@layers,function(x){dim(x)})
#glioma_oligo_cells@assays$RNA@layers[names(layersList[sapply(layersList, is.null)])] <- NULL


glioma_oligo_cells <- NormalizeData(glioma_oligo_cells)
glioma_oligo_cells <- FindVariableFeatures(glioma_oligo_cells, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(glioma_oligo_cells)
glioma_oligo_cells <- ScaleData(glioma_oligo_cells, features = all.genes)
gc()

variable_features <- VariableFeatures(object = glioma_oligo_cells)

glioma_oligo_cells <- RunPCA(glioma_oligo_cells,npcs = 50)
pdf(file='elbow_gliomaclusters_v3.pdf',width =10,height = 8)
ElbowPlot(glioma_oligo_cells,ndims=40)
dev.off()
#glioma_oligo_cells<-RunHarmony(object=glioma_oligo_cells,group.by.vars="Fragment",reduction="pca",reduction.save="harmony")
glioma_oligo_cells<-RunHarmony(object=glioma_oligo_cells,group.by.vars=c("Fragment","Patient"),reduction="pca",reduction.save="harmony")
glioma_oligo_cells <- RunUMAP(glioma_oligo_cells,reduction = "harmony",dims = 1:10,reduction.name="UMAP_harmony",verbose = FALSE)
glioma_oligo_cells <- FindNeighbors(glioma_oligo_cells,reduction = "harmony", dims = 1:10)
glioma_oligo_cells <- FindClusters(glioma_oligo_cells,reduction = "harmony", resolution = 0.2)

#na_cells<-which(is.na(glioma_oligo_cells$seurat_clusters))
#glioma_oligo_cells<-subset(glioma_oligo_cells,cells=Cells(glioma_oligo_cells)[-na_cells])

#define axis limits
# Get UMAP coordinates from Seurat object
umap_coords <- as.data.frame(Embeddings(glioma_oligo_cells, reduction = "UMAP_harmony"))
x_limits <- c(min(umap_coords$UMAPharmony_1), max(umap_coords$UMAPharmony_1))
y_limits <- c(min(umap_coords$UMAPharmony_2), max(umap_coords$UMAPharmony_2))
#rename clusters
glioma_oligo_cells<-RenameIdents(glioma_oligo_cells,levels(glioma_oligo_cells))
# Create UMAP plot
p <- DimPlot(glioma_oligo_cells, reduction = "UMAP_harmony",label = TRUE, pt.size = 0.5)
# Increase axis range
p <- p + scale_x_continuous(limits = x_limits)
p <- p + scale_y_continuous(limits = y_limits)

pdf(file = 'gliomaoligo_UMAP_harmony_patient_clusters_10pcs_res0.2_v3.pdf', width = 14, height = 12)
# Print the plot
print(p)
dev.off()

# Replace values in the "Type" column
glioma_oligo_cells@meta.data$Type <- ifelse(
  glioma_oligo_cells@meta.data$Type %in% c("Astrocytoma", "Oligodendroglioma"),
  "LGG", 
  glioma_oligo_cells@meta.data$Type
)
# Check the unique values in the "Type" column
unique(glioma_oligo_cells@meta.data$Type)
#reorder column level
glioma_oligo_cells@meta.data$Type<-factor(x=glioma_oligo_cells@meta.data$Type,levels=c("LGG","GBM","Recurrent GBM"))
# Create UMAP plot
p <- DimPlot(glioma_oligo_cells, reduction = "UMAP_harmony",label = TRUE, pt.size = 0.5,split.by="Type")
# Increase axis range
p <- p + scale_x_continuous(limits = x_limits)
p <- p + scale_y_continuous(limits = y_limits)

pdf(file = 'gliomaoligo_UMAP_harmony_patient_splittedType_10pcs_res0.2_v3.pdf', width = 18, height = 10)
# Print the plot
print(p)
dev.off()

p <- DimPlot(glioma_oligo_cells, reduction = "UMAP_harmony",group.by = 'Assignment',label = TRUE, pt.size = 0.01)
# Increase axis range
p <- p + scale_x_continuous(limits = x_limits)
p <- p + scale_y_continuous(limits = y_limits)
# Create PDF file with specified dimensions
pdf(file = 'gliomaoligo_UMAP_harmony_patient_celltypes_10pcs_res0.2_v3.pdf', width = 14, height = 12)
p <- p + ggtitle("Cell Type")
# Print the plot
print(p)
dev.off()

#phase
p <- DimPlot(glioma_oligo_cells, reduction = "UMAP_harmony",group.by = 'Phase')
# Increase axis range
p <- p + scale_x_continuous(limits = x_limits)
p <- p + scale_y_continuous(limits = y_limits)
# Create PDF file with specified dimensions
pdf(file = 'gliomaoligo_UMAP_harmony_patient_cellcycle_10pcs_res0.2_v3.pdf', width = 14, height = 12)
p <- p + ggtitle("Cell Cycle")
# Print the plot
print(p)
dev.off()

#subcelltype
p <- DimPlot(glioma_oligo_cells, reduction = "UMAP_harmony",group.by = 'SubAssignment',label = TRUE, pt.size = 0.01)
# Increase axis range
p <- p + scale_x_continuous(limits = x_limits)
p <- p + scale_y_continuous(limits = y_limits)
# Create PDF file with specified dimensions
pdf(file = 'gliomaoligo_UMAP_harmony_patient_subcelltypes_10pcs_res0.2_v3.pdf', width = 14, height = 12)
p <- p + ggtitle("Sub-Types")
# Print the plot
print(p)
dev.off()

#

#
#glioma_oligo_cells<-JoinLayers(glioma_oligo_cells) #não roda no objeto assay v3
all_markers_cluster0<-FindMarkers(glioma_oligo_cells,ident.1=0)
saveRDS(all_markers_cluster0,file="all_markers_cluster0.rds")
all_markers_cluster1<-FindMarkers(glioma_oligo_cells,ident.1=1)
saveRDS(all_markers_cluster1,file="all_markers_cluster1.rds")
all_markers_cluster2<-FindMarkers(glioma_oligo_cells,ident.1=2)
saveRDS(all_markers_cluster2,file="all_markers_cluster2.rds")
all_markers_cluster3<-FindMarkers(glioma_oligo_cells,ident.1=3)
saveRDS(all_markers_cluster3,file="all_markers_cluster3.rds")
all_markers_cluster4<-FindMarkers(glioma_oligo_cells,ident.1=4)
saveRDS(all_markers_cluster4,file="all_markers_cluster4.rds")
all_markers_cluster5<-FindMarkers(glioma_oligo_cells,ident.1=5)
saveRDS(all_markers_cluster5,file="all_markers_cluster5.rds")



#rename clusters
#glioma_oligo_cells<-RenameIdents(glioma_oligo_cells,levels(glioma_oligo_cells))
# Create UMAP plot
p <- DimPlot(glioma_oligo_cells, reduction = "UMAP_harmony",group.by = 'Fragment',label = TRUE, pt.size = 0.01)
# Increase axis range
p <- p + scale_x_continuous(limits = x_limits)
p <- p + scale_y_continuous(limits = y_limits)
# Create PDF file with specified dimensions
pdf(file = 'gliomaoligo_UMAP_harmony_patient_biopsia_10pcs_res0.2_v3.pdf', width = 14, height = 12)
p <- p + ggtitle("Biopsy")
# Print the plot
print(p)
dev.off()

#rename clusters
#glioma_oligo_cells<-RenameIdents(glioma_oligo_cells,levels(glioma_oligo_cells))
# Create UMAP plot
p <- DimPlot(glioma_oligo_cells, reduction = "UMAP_harmony",group.by = 'Patient',label = TRUE, pt.size = 0.01)
# Increase axis range
p <- p + scale_x_continuous(limits = x_limits)
p <- p + scale_y_continuous(limits = y_limits)
# Create PDF file with specified dimensions
pdf(file = 'gliomaoligo_UMAP_harmony_patient_pacientes_10pcs_res0.2_v3.pdf', width = 14, height = 12)
p <- p + ggtitle("Patients")
# Print the plot
print(p)
dev.off()

#rename clusters
#glioma_oligo_cells<-RenameIdents(glioma_oligo_cells,levels(glioma_oligo_cells))
# Create UMAP plot
p <- DimPlot(glioma_oligo_cells, reduction = "UMAP_harmony",group.by = 'Type',label = TRUE, pt.size = 0.01)
# Increase axis range
p <- p + scale_x_continuous(limits = x_limits)
p <- p + scale_y_continuous(limits = y_limits)
# Create PDF file with specified dimensions
pdf(file = 'gliomaoligo_UMAP_harmony_patient_type_10pcs_res0.2_v3.pdf', width = 14, height = 12)
p <- p + ggtitle("Tumor Type")
# Print the plot
print(p)
dev.off()

glioma_oligo_cells<-JoinLayers(glioma_oligo_cells)
all_markers_cluster3<-FindMarkers(glioma_oligo_cells,ident.1=3)
saveRDS(all_markers_cluster3,file="all_markers_cluster3.rds")
all_markers_cluster19<-FindMarkers(glioma_oligo_cells,ident.1=19)
saveRDS(all_markers_cluster19,file="all_markers_cluster19.rds")



markers_cluster4_between_tumortypes<-FindConservedMarkers(glioma_oligo_cells,ident.1=4,grouping.var="Type")#sera que é util?
markers_cluster4_between_patients<-FindConservedMarkers(glioma_oligo_cells,ident.1=4,grouping.var="Patient")

#markers cluster4 vs 0,1,8,10,11 e 12;
markers_cluster4_vs_othersfrommyeloids<-FindMarkers(glioma_oligo_cells,ident.1=4,ident.2=c(0,1,8,10,11,12))
saveRDS(markers_cluster4_vs_othersfrommyeloids,file="markers_cluster4_vs_othersfrommyeloids.rds")

#all markers of every cluster compared to all remaining cells
all_markers_everycluster<-FindAllMarkers(glioma_oligo_cells,only.pos=TRUE) %>%
  dplyr::filter(avg_log2FC>1)

#all markers cluster 4 vs all other
all_markers_cluster3<-FindMarkers(glioma_oligo_cells,ident.1=3)
saveRDS(all_markers_cluster3,file="all_markers_cluster3.rds")
all_markers_cluster19<-FindMarkers(glioma_oligo_cells,ident.1=19)
saveRDS(all_markers_cluster19,file="all_markers_cluster19.rds")

top10C4<-top_n(all_markers_cluster4,10,avg_log2FC)
bottom10C4<-top_n(all_markers_cluster4,-10,avg_log2FC)
pdf(file='top10_dfgenes_cluster4_vs_all_spltype_v3.pdf',width = 16, height = 52)
FeaturePlot(glioma_oligo_cells, features =rownames(top10C4),split.by="Type",min.cutoff='q10',ncol = 5)& 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
par(cex.axis = 0.2)
dev.off()
pdf(file='bottom10_dfgenes_cluster4_vs_all_spltype_v3.pdf',width = 16, height = 52)
FeaturePlot(glioma_oligo_cells, features =rownames(bottom10C4),split.by="Type",min.cutoff='q10',ncol = 5)& 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
par(cex.axis = 0.2)
dev.off()

#top10 cluster 4
glioma_oligo_cells@meta.data$Type <- factor(glioma_oligo_cells@meta.data$Type,levels=c("LGG", "GBM", "Recurrent GBM"))
pdf(file='top10_cluster4features_type_v3.pdf',width = 26, height = 28)
FeaturePlot(glioma_oligo_cells, features =c(rownames(markers_cluster4)[1:10],"CXCL5",
                                                "MT1H","SLAMF9","SPINK1","ATP6V0D2","TM4SF19","BICDL2"),
            min.cutoff='q10')& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
par(cex.axis = 0.2)
dev.off()
#bottom10 cluster 4
pdf(file='bottom10_cluster4features_type_v3.pdf',width = 22, height = 30)
FeaturePlot(glioma_oligo_cells, features =c("ADGRL3","PHYHIPL","COL20A1","DCX","LRRC4C","CNTN1","LRRTM2","FGF14","PLEKHH2","KLRC2","GRIA2","SCG3","NRXN1","DLL3","NKAIN4","SMOC1","STMN4","FOXG1","MEG3","TNR")
            ,min.cutoff='q10')& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
par(cex.axis = 0.2)
dev.off()

#azulvermelho
library(RColorBrewer)
brewer.pal(11, "RdBu")
pdf(file='top10_cluster4features_patiens_v3.pdf',width = 18, height = 30)
FeaturePlot(glioma_oligo_cells, features =rownames(markers_cluster4_between_patients)[1:10],
            split.by="Type",min.cutoff='q10',ncol = 5)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
par(cex.axis = 0.2)
dev.off()

#top 10 e bottom 10 GBM, rGBM
top10c4GBM<-top_n(markers_cluster4,10,GBM_avg_log2FC)
top10c4rGBM<-top_n(markers_cluster4,10,`Recurrent GBM_avg_log2FC`)
bottom10c4GBM<-top_n(markers_cluster4,-10,GBM_avg_log2FC)
bottom10c4rGBM<-top_n(markers_cluster4,-10,`Recurrent GBM_avg_log2FC`)

#top GBM diff exp genes of cluster 4 markers
pdf(file='top10_dfgenes_cluster4_GBM_v3.pdf',width = 16, height = 6)
FeaturePlot(glioma_oligo_cells, features =rownames(top10c4GBM),min.cutoff='q10',ncol = 5)
par(cex.axis = 0.2)
dev.off()
#top rGBM diff exp genes of cluster 4
pdf(file='top10_dfgenes_cluster4_rGBM_v3.pdf',width = 16, height = 6)
FeaturePlot(glioma_oligo_cells, features =rownames(top10c4rGBM),min.cutoff='q10',ncol = 5)
par(cex.axis = 0.2)
dev.off()
#bottom GBM diff exp genes of cluster 4 markers
pdf(file='bottom10_dfgenes_cluster4_GBM_v3.pdf',width = 16, height = 6)
FeaturePlot(glioma_oligo_cells, features =rownames(bottom10c4GBM),min.cutoff='q10',ncol = 5)
par(cex.axis = 0.2)
dev.off()
#top rGBM diff exp genes of cluster 4
pdf(file='bottom10_dfgenes_cluster4_rGBM_v3.pdf',width = 16, height = 6)
FeaturePlot(glioma_oligo_cells, features =rownames(bottom10c4rGBM),min.cutoff='q10',ncol = 5)
par(cex.axis = 0.2)
dev.off()


#feature plots for GSC markers
pdf(file = 'gliomaoligo_featureplots_GSCmarkers_v3.pdf', width = 14, height = 6)
FeaturePlot(glioma_oligo_cells, features = c("PROM1", "FUT4", "SOX2", "NANOG", "OLIG2", "NES", "CD9"),min.cutoff='q10' ,ncol = 4)
par(cex.axis = 0.4)
dev.off()

pdf(file = 'gliomaoligo_featureplots_ACmarkers_v3.pdf', width = 50, height = 30)
FeaturePlot(glioma_oligo_cells, features = c("CST3", "S100B", "SLC1A3", "HEPN1", "HOPX", "MT3", "SPARCL1", "MLC1", "GFAP", "FABP7", "BCAN", "PON2", 
                                                 "METTL7B", "SPARC", "GATM", "RAMP1", "PMP2", "AQP4", "DBI", "EDNRB", "PTPRZ1", "CLU", "PMP22", "ATP1A2", 
                                                 "S100A16", "HEY1", "PCDHGC3", "TTYH1", "NDRG2", "PRCP", "ATP1B2", "AGT", "PLTP", "GPM6B", "F3", "RAB31", 
                                                 "PLPP3", "ANXA5", "TSPAN7"),min.cutoff='q10' ,ncol = 10)
par(cex.axis = 0.2)
dev.off()
#Add module Score
topC4<-list(c(rownames(top10c4GBM),rownames(top10c4rGBM)))
bottomC4<-list(c(rownames(bottom10c4GBM),rownames(bottom10c4rGBM)))
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
glioma_oligo_cells<-JoinLayers(glioma_oligo_cells)
glioma_oligo_cells <- AddModuleScore(object = glioma_oligo_cells, 
                                         features = GSC, 
                                         name = "GSC_score")
glioma_oligo_cells <- AddModuleScore(object = glioma_oligo_cells, 
                                         features = AC, 
                                         name = "AC_score")
glioma_oligo_cells <- AddModuleScore(object = glioma_oligo_cells, 
                                         features = OPC, 
                                         name = "OPC_score")
glioma_oligo_cells <- AddModuleScore(object = glioma_oligo_cells, 
                                         features = NPC, 
                                         name = "NPC_score")
glioma_oligo_cells <- AddModuleScore(object = glioma_oligo_cells, 
                                         features = MES, 
                                         name = "MES_score")


pdf(file = 'featureplots_topC4_v3.pdf')
FeaturePlot(glioma_oligo_cells,features = "C4_score1", label = TRUE, repel = TRUE)+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()
pdf(file = 'featureplots_bottomC4_v3.pdf')
FeaturePlot(glioma_oligo_cells,features = "bottomC4_score1", label = TRUE, repel = TRUE)+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()

pdf(file = 'gliomaoligo_featureplots_AC_v3.pdf')
FeaturePlot(glioma_oligo_cells,features = "AC_score1", label = TRUE, repel = TRUE)+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()
pdf(file = 'gliomaoligo_featureplots_OPC_v3.pdf')
FeaturePlot(glioma_oligo_cells,features = "OPC_score1", label = TRUE, repel = TRUE)+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()
pdf(file = 'gliomaoligo_featureplots_NPC_v3.pdf')
FeaturePlot(glioma_oligo_cells,features = "NPC_score1", label = TRUE, repel = TRUE)+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()
pdf(file = 'gliomaoligo_featureplots_MES_v3.pdf')
FeaturePlot(glioma_oligo_cells,features = "MES_score1", label = TRUE, repel = TRUE)+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()
pdf(file = 'gliomaoligo_featureplots_GSC_v3.pdf')
FeaturePlot(glioma_oligo_cells,features = "GSC_score1", label = TRUE, repel = TRUE)+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()


glioma_oligo_cells$gliomastatescore <- apply(glioma_oligo_cells@meta.data[, c("GSC_score1", "AC_score1", "OPC_score1", "NPC_score1", "MES_score1")], 1, function(x) {
  scores <- c(GSC = x["GSC_score1"], AC = x["AC_score1"], OPC = x["OPC_score1"], NPC = x["NPC_score1"], MES = x["MES_score1"])
  return(names(scores)[which.max(scores)])
})
pdf(file='gliomaoligo_cell_states_modulescore_v3.pdf')
DimPlot(glioma_oligo_cells, reduction = "UMAP_harmony",label=FALSE,group.by ="gliomastatescore")
dev.off()
