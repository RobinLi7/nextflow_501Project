# Load necessary libraries
library(Seurat)
library(ggplot2)

# Load the data
args <- commandArgs(trailingOnly = TRUE)

raw_file <- args[1]
output_dir <- args[2]
# raw_file <- "/projects/rli_prj/CPSC501/project/work/aa/ee6e864d335bcf73f34d315c400183/result/Solo.out/Gene/raw"
# print(raw_file)
data <- Read10X(data.dir = raw_file)
# data <- Read10X(data.dir = file.path(output_dir, "Solo.out", "Gene", "raw"))

# Create Seurat object
seurat_obj <- CreateSeuratObject(counts = data, project = "scRNAseq", min.cells = 3, min.features = 200)

# Preprocessing
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)

# Determine the number of features and cells
num_features <- nrow(seurat_obj)
num_cells <- ncol(seurat_obj)

# Set npcs to a value less than the smaller dimension
npcs_to_use <- min(num_features, num_cells) - 1

# Run PCA with the adjusted npcs value
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj), npcs = npcs_to_use)
# seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))

# Save the Seurat object
output_file <- file.path(output_dir, "seurat_obj.rds")
saveRDS(seurat_obj, file = output_file)