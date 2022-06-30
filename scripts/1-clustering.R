# Load library
library(Seurat) # v4.0.1

# Make results directories if they do not exist
output_dirs <- c("results",
                 "results/objects")

for (i in output_dirs) {
  if (!dir.exists(i)) {
    dir.create(i)
    print(paste0("made ", i, " directory"))
  } else {
    print(paste0(i, " directory already exists."))
  }
}

# Load in data and add meta data
a1 <- Load10X_Spatial(
  data.dir = "data/spaceranger/643-1_A1_manualAlignment2/outs",
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "KO",
  filter.matrix = T,
  to.upper = F
)
b1 <- Load10X_Spatial(
  data.dir = "data/spaceranger/643-2_B1/outs",
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "Control",
  filter.matrix = T,
  to.upper = F
)

a1@meta.data$mouse_line <- "1955"
a1@meta.data$genotype_long <- "iAdipoQCre-het_ATGL-flox-homo"
a1@meta.data$genotype <- "KO"
a1@meta.data$animal_number <- "283075"
a1@meta.data$surgery <- "IR-24h"

b1@meta.data$mouse_line <- "1955"
b1@meta.data$genotype_long <- "iAdipoQCre-wt_ATGL-flox-homo"
b1@meta.data$genotype <- "Control"
b1@meta.data$animal_number <- "283076"
b1@meta.data$surgery <- "IR-24h"

# SCTransform
a1 <- SCTransform(a1, assay = "Spatial")
b1 <- SCTransform(b1, assay = "Spatial")

# Using direct list instead of merge,
# otherwise not all genes are present in the SCT model:
# https://github.com/satijalab/seurat/issues/3198
hearts_list <- list(b1, a1)

# Integrate
features <- SelectIntegrationFeatures(
  object.list = hearts_list,
  nfeatures = 3000
)
hearts_list <- PrepSCTIntegration(
  object.list = hearts_list,
  anchor.features = features
)
hearts_anchors <- FindIntegrationAnchors(
  object.list = hearts_list,
  normalization.method = "SCT",
  anchor.features = features
)
hearts <- IntegrateData(
  anchorset = hearts_anchors,
  normalization.method = "SCT"
)

# Cluster using integrated assay
DefaultAssay(hearts) <- "integrated"
hearts <- RunPCA(hearts, verbose = T)
ElbowPlot(hearts, ndims = 50)
hearts <- FindNeighbors(hearts, dims = 1:20)
hearts <- FindClusters(hearts, verbose = T, res = 0.4)
hearts <- RunUMAP(hearts, reduction = "pca", dims = 1:20)

# Normalizing and scaling the Spatial assay for DEG testing
DefaultAssay(hearts) <- "Spatial"
hearts <- NormalizeData(hearts)
all_genes <- rownames(hearts)
hearts <- ScaleData(hearts, features = all_genes)

# Set variable features in SCT assay, because merge wasn't used
DefaultAssay(hearts) <- "SCT"
VariableFeatures(hearts) <- c(VariableFeatures(a1), VariableFeatures(b1))

# Factor surgery variable 
hearts@meta.data$genotype <- factor(hearts@meta.data$genotype,
  levels = c("Control", "KO"))

# Save clustered object
save(hearts, file = "results/objects/hearts.Rdata")
