# Differential gene expression, each cluster, Control v KO

# Make results directories if they do not exist
output_dirs <- c("results",
                 "results/differential-gene-expression")

for (i in output_dirs) {
  if (!dir.exists(i)) {
    dir.create(i)
    print(paste0("made ", i, " directory"))
  } else {
    print(paste0(i, " directory already exists."))
  }
}

# Load libraries, functions and objects
library(Seurat) # v4.0.1
library(ggplot2)
library(patchwork)
library(scales)
library(dplyr)
library(yulab.utils)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(ggrepel)
source("scripts/SpatialFeaturePlotScaled.R")
source("scripts/VolcanoPlot.R")
load("results/objects/hearts.Rdata")
default_colors <- (hue_pal()(7))

# No significance threshold
de_genes <- list()
for (i in levels(Idents(hearts))) {
  results <- FindMarkers(hearts,
    group.by = "genotype",
    subset.ident = i,
    ident.1 = "KO",
    assay = "Spatial",
    logfc.threshold = 0
  )
  results$cluster <- i
  results$gene <- row.names(results)
  row.names(results) <- NULL
  results$regulation <- ifelse(results$avg_log2FC > 0, "Up", "Down")
  de_genes[[i]] <- results
}
deg_clusters <- do.call(rbind, de_genes)
rownames(deg_clusters) <- NULL
write.csv(deg_clusters,
  file = "results/differential-gene-expression/deg_clusters.csv",
  row.names = F
)
remove(de_genes)

# With significance threshold
de_genes <- list()
for (i in levels(Idents(hearts))) {
  results <- FindMarkers(hearts,
                         group.by = "genotype",
                         subset.ident = i,
                         ident.1 = "KO",
                         assay = "Spatial",
                         logfc.threshold = 0.25,
                         
  )
  results$cluster <- i
  results$gene <- row.names(results)
  row.names(results) <- NULL
  results$regulation <- ifelse(results$avg_log2FC > 0, "Up", "Down")
  de_genes[[i]] <- results
}
deg_clusters_sig <- do.call(rbind, de_genes)
deg_clusters_sig <- deg_clusters_sig %>% filter(p_val_adj <= 0.01)
rownames(deg_clusters_sig) <- NULL
write.csv(deg_clusters_sig,
          file = "results/differential-gene-expression/deg_clusters_sig.csv",
          row.names = F
)
remove(de_genes)

# Volcano plot
for (i in unique(deg_clusters$cluster)) {
  pdf(
    file = paste0(
      "results/differential-gene-expression/VolcanoPlot_cluster_",
      i,
      "_Control_vs_KO.pdf"
    ),
    height = 6,
    width = 8,
    useDingbats = F
  )
  VolcanoPlot(
    df = deg_clusters,
    identity = i,
    title = paste0("DEG cluster ", i, " - KO/Control")
  )
  dev.off()
}
