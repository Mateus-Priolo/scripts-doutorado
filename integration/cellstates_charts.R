library(gridExtra)
library(grid)
library(ggplot2)
library(Seurat)
types <- c("LGG", types[types != "LGG"])
# Calculate the percentage of each score within each cluster
scores <- harmony_merged_objects@meta.data[, c("seurat_clusters", "Type", "MES_score1", "AC_score1", "OPC_score1", "GSC_score1", "NPC_score1")]
scores <- scores[rowSums(is.na(scores)) == 0, ]  # Remove rows with NA values

# Reshape data to long format
scores_long <- reshape2::melt(scores, id.vars = c("seurat_clusters", "Type"), variable.name = "Score")

# Group by Type, seurat_clusters, and Score, calculate percentages
score_percentages <- scores_long %>%
  group_by(Type, seurat_clusters, Score) %>%
  summarise(Percentage = mean(value, na.rm = TRUE))

# Get unique clusters and types
clusters <- sort(unique(score_percentages$seurat_clusters))
types <- sort(unique(score_percentages$Type))

# Create a grid of plots
n_clusters <- length(clusters)
n_types <- length(types)
plot_list <- list()
# Create plots
legend_drawn <- FALSE
for (i in 1:n_clusters) {
  cluster <- clusters[i]
  for (j in 1:n_types) {
    type <- types[j]
    data <- score_percentages[score_percentages$Type == type & score_percentages$seurat_clusters == cluster, ]
    if (nrow(data) > 0) {
      plot_obj <- ggplot(data, aes(x = "", y = Percentage, fill = Score)) +
        geom_bar(stat = "identity", width = 1) +
        coord_polar(theta = "y") +
        labs(title = paste("Cluster:", cluster, "| Type:", type)) +
        theme_void() + scale_fill_brewer(palette = "Set1")+
        theme(legend.position = ifelse(!legend_drawn, "right", "none"),  # Place legend on the first plot only
              plot.title = element_text(size = 8),                    # Adjust title size
              legend.key.size = unit(0.6, "cm"),                        # Adjust legend size
              legend.margin = margin(0, 0, 0, 0, "cm"),               # Add margin to legend
              plot.margin = margin(0, 0, 0, 0, "cm"))                 # Adjust the plot margins
      plot_list[[length(plot_list) + 1]] <- ggplotGrob(plot_obj)
      legend_drawn <- TRUE  # Set legend_drawn to TRUE after drawing the first legend
    } else {
      # If the data is empty for the current cluster and type, add a blank plot
      plot_list[[length(plot_list) + 1]] <- ggplotGrob(ggplot() + theme_void())
    }
  }
}
# Remove NULL plots
plot_list <- plot_list[!sapply(plot_list, is.null)]
# Arrange plots in a grid
grid_arranged_plots <- grid.arrange(grobs = plot_list, ncol = 3)
# Save as PDF
pdf("cluster_pie_charts.pdf", width = 12, height = 30)  # Adjust the height for two pages
grid.draw(grid_arranged_plots)
dev.off()

#library(RColorBrewer)
#display.brewer.all()