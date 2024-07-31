library(dplyr)
library(DESeq2)
library(org.Hs.eg.db)
library(circlize)
library(ComplexHeatmap)

setwd("C:/Users/...")

features<-read.table("merged_counts_clean.txt",header = TRUE)
samples<-read.table("samples.txt", header = TRUE)

#colapsar technical replicates
#dds<-collapseReplicates(dds, dds$cell)
#dds<-collapseReplicates(dds, groupby = dds$cell, run = dds$run)

#chatgpt:
# Assuming 'samples' contains the sample information with 'cell' as a factor

# Create a DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = features, colData = samples, design = ~ cell)
keep<-rowSums(counts(dds))>=10
dds<-dds[keep,]
# Set the reference level
dds$cell <- relevel(dds$cell, ref = "ACBRI371")
# Run DESeq
dds <- DESeq(dds)
res_adjp0.01<-results(dds,alpha = 0.01)

# List of cells (excluding ACBRI371)
other_cells <- levels(dds$cell)[levels(dds$cell) != "ACBRI371"]

# Initialize an empty list to store results
results_list <- list()

# Loop through each cell and perform a contrast
for (cell in other_cells) {
  contrast_name <- paste(cell, "vsACBRI371", sep = "_")
  res <- results(dds, contrast = c("cell", cell, "ACBRI371"))
  results_list[[contrast_name]] <- res
}

# Combine all results into a single DataFrame
combined_results <- do.call(cbind, results_list)
View(combined_results)

# Select columns containing '.log2FoldChange'
log2FoldChange_columns <- grep("\\.log2FoldChange", colnames(combined_results), value = TRUE)

# Subset combined_results to keep only columns with log2FoldChange
log2FoldChange_data <- combined_results[, log2FoldChange_columns]

# Write the filtered data to a new file
write.table(log2FoldChange_data, file = "filtered_log2FoldChange.txt", sep = "\t", quote = FALSE, row.names = TRUE)


res_fc<-as.data.frame(log2FoldChange_data)
res_fc$symbol<-mapIds(org.Hs.eg.db,keys = rownames(res_fc),keytype = "ENSEMBL",column = "SYMBOL")
View(res_fc)

library(org.Hs.eg.db)
res_df<-as.data.frame(res_adjp0.01)
View(res_df)
res_df$symbol<-mapIds(org.Hs.eg.db,keys = rownames(res_df),keytype = "ENSEMBL",column = "SYMBOL")
View(res_df)
normal<-counts(dds, normalized=TRUE)[rownames(res_df),] #normalized counts from dds object, rownames: pega so os genes do res_df

#norm_zscore<-t(apply(normal,1,scale)) #transpose because its aplied to the columns
View(normal)
#View(norm_zscore)

normal$symbol<-mapIds(org.Hs.eg.db,keys = rownames(normal),keytype = "ENSEMBL",column = "SYMBOL")

# colnames nome das amostras )
colnames(normal)<-samples$cell
?vst
vsdata<-vst(dds,blind = FALSE)
plotPCA(vsdata,intgroup="cell")
plotDispEsts(dds) #variabilidade das replicatas em fç de normalized read counts

genespaloma<-read.table("listagenespaloma.txt", header = FALSE)
genespaloma <- scan("listagenespaloma.txt", what = "character", sep = "\n")
res_fc_31genes <- res_fc[res_fc$symbol %in% genespaloma, ]
write.table(res_fc_31genes, file = "foldchanges31.txt", sep = "\t", quote = FALSE, row.names = TRUE)
rownames(res_fc_31genes) <- res_fc_31genes[, ncol(res_fc_31genes)]

# Remove the last column (as it's now the row names)
res_fc_31genes <- res_fc_31genes[, -ncol(res_fc_31genes)]
View(res_fc_31genes)
library(circlize)
fold_changes <- as.matrix(res_fc_31genes) # Extract fold change values for each gene
res_fc_31genes_matrix <-as.matrix(res_fc_31genes)

fold_change_range <- range(res_fc_31genes_matrix, na.rm = TRUE)
fold_change_range
my_color_mapping <- colorRamp2(c(-10, 0, 10), c("blue", "white", "red"))
Heatmap(res_fc_31genes_matrix, cluster_rows = FALSE, column_labels = colnames(res_fc_31genes_matrix),name="FC", row_labels =rownames(res_fc_31genes_matrix), row_names_gp = gpar(fontsize = 5),column_names_gp = gpar(fontsize = 6),heatmap_height = unit(40, "npc"))
