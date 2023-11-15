# In this script we will run the pre-processing for each sample
#if (!"BiocManager" %in% installed.packages()) install.packages("BiocManager")

x <- c("Seurat", "ggplot2", "patchwork", "clustree", "glmGamPoi")
#BiocManager::install(x)

# Load libraries
invisible(lapply(x, library, character.only = TRUE))
#-----------------------------------

# In this script we will perform integrated clustering of both genotypes

# Set up output dirs
output_dirs <- c("results",
                 "results/02_integrated-clustering_mine")
for (i in output_dirs) {
  if (!dir.exists(i)) {
    dir.create(i)
    print(paste0("made ", i, " directory"))
  } else {
    print(paste0(i, " directory already exists."))
  }
}

# Load obj list
load(file = "results/01_objects_mine/01_obj_list.Rdata")

# SCTransform
obj_list <- lapply(X = obj_list, FUN = function(x){
  x <- SCTransform(x, assay = "Spatial", vst.flavor = "v2")
  })

# Integration
features <- SelectIntegrationFeatures(obj_list, nfeatures = 3000)
obj_list <- PrepSCTIntegration(obj_list, anchor.features = features)
anchors <- FindIntegrationAnchors(obj_list, normalization.method = "SCT",
                                  anchor.features = features)
obj <- IntegrateData(anchors, normalization.method = "SCT")

# Clustering
obj <- RunPCA(obj, verbose = T)
ElbowPlot(obj, ndims = 50)
obj <- FindNeighbors(obj, dims = 1:20)
res_range <- c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2)
obj <- FindClusters(obj, verbose = T, res = res_range)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:20)
DefaultAssay(obj) <- "Spatial"

# Find the optimum clustering resolution with clustree
# https://cran.r-project.org/web/packages/clustree/vignettes/clustree.html
clust <- clustree(obj, prefix = "integrated_snn_res.", layout = "sugiyama")
pdf(file ="results/02_integrated-clustering_mine/clustree.pdf",
    useDingbats = F)
print(clust)
dev.off()

# Normalizing and scaling Spatial assay for DEG testing
obj <- NormalizeData(obj)
all_genes <- rownames(obj)
obj <- ScaleData(obj, features = all_genes)

# Factor surgery variable 
obj@meta.data$genotype <- factor(obj@meta.data$genotype,
                                    levels = c("Control", "KO"))

# Save clustered object
save(obj, file = "results/01_objects_mine/02_obj_integrated.Rdata")
