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
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.8)

# Find markers for each cluster
markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.15, logfc.threshold = 0.15)

# Check if markers were found
if (nrow(markers) > 0) {
  # Get top 5 markers for each group
  top_markers <- markers %>%
    group_by(cluster) %>%
    top_n(n = 5, wt = avg_log2FC)

  # Print top markers to help with cell type identification
  write.csv(top_markers, file = file.path(output_dir, "top_markers_umap.csv"))
} else {
  message("No markers found. Please check the preprocessing steps and cluster sizes.")
}

# Define new cluster identities
new_cluster_ids <- c("Neutrophils", "NK Cells", "B cells")
# new_cluster_ids <- c("T Cells", "Monocytes", "B cells", "Neutrophils", "NK Cells", "NK Cells", "Macrophages")

# Assign new identities to clusters
names(new_cluster_ids) <- levels(seurat_obj)
seurat_obj <- RenameIdents(seurat_obj, new_cluster_ids)

seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)

# Plot UMAP
umap_plot <- DimPlot(seurat_obj, reduction = "umap", label = TRUE) + NoLegend()

# Save the UMAP plot
output_umap_file <- file.path(output_dir, "umap_plot.png")
ggsave(output_umap_file, plot = umap_plot, width = 8, height = 6)
