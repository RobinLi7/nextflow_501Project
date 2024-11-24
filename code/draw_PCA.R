# Load necessary libraries
library(Seurat)
library(ggplot2)
library(dplyr)

# Load the preprocessed Seurat object
args <- commandArgs(trailingOnly = TRUE)
seurat_file <- args[1]
output_dir <- args[2]

# seurat_file <- "/projects/rli_prj/CPSC501/project/work/e3/d7bb67e92222ee62b9806e671affbd/seurat_obj.rds"
# output_dir <- "."

seurat_obj <- readRDS(seurat_file)

# Find neighbors and clusters
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.8)

# Find markers for each cluster
markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.15, logfc.threshold = 0.15)

# Check if markers were found
if (nrow(markers) > 0) {
  # Get top 5 markers for each cluster
  top_markers <- markers %>%
    group_by(cluster) %>%
    top_n(n = 5, wt = avg_log2FC)

  # Save top markers to a CSV file
  write.csv(top_markers, file = file.path(output_dir, "top_markers_pca.csv"))
} else {
  message("No markers found. Please check the preprocessing steps and cluster sizes.")
}

# Define new cluster identities
new_cluster_ids <- c("Neutrophils", "NK Cells", "B cells")
# new_cluster_ids <- c("T Cells", "Monocytes", "B cells", "Neutrophils", "NK Cells", "NK Cells", "Macrophages")


# Assign new identities to clusters
names(new_cluster_ids) <- levels(seurat_obj)
seurat_obj <- RenameIdents(seurat_obj, new_cluster_ids)

# Run PCA if not already done
if (!"pca" %in% names(seurat_obj@reductions)) {
  seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
}

# Plot PCA
pca_plot <- DimPlot(seurat_obj, reduction = "pca", label = TRUE) + NoLegend()

# Save the PCA plot
output_pca_file <- file.path(output_dir, "pca_plot.png")
ggsave(output_pca_file, plot = pca_plot, width = 8, height = 6)
