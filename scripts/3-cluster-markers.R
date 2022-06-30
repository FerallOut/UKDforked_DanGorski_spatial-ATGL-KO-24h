# Identify cluster markers

# Make results directories if they do not exist
output_dirs <- c("results",
                 "results/cluster-markers")

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
source("scripts/SpatialFeaturePlotScaled.R")
source("scripts/GOBP_fold_enrich.R")
load("results/objects/hearts.Rdata")
default_colors <- (hue_pal()(7))

# Find cluster markers
cluster_markers <- FindAllMarkers(hearts, assay = "SCT", only.pos = T)
top_10_cluster_markers <- cluster_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)
top_100_cluster_markers <- cluster_markers %>%
  group_by(cluster) %>%
  top_n(n = 100, wt = avg_log2FC)
write.csv(cluster_markers,
  file = "results/cluster-markers/cluster_markers.csv",
  row.names = F
)

# SpatialFeaturePlots of top 10 cluster markers
genes <- top_10_cluster_markers$gene
for (i in genes) {
  pdf(
    file = paste0("results/cluster-markers/SpatialFeaturePlot_", i, ".pdf"),
    height = 6,
    width = 12.75,
    useDingbats = F
  )
  SpatialFeaturePlotScaled(
    object = hearts,
    group = "genotype",
    group.1 = "Control",
    group.2 = "KO",
    feature_of_interest = i,
    from.meta.data = FALSE,
    group.1.title = "Control",
    group.2.title = "KO"
  )
  dev.off()
}

# DotPlot
top_5_cluster_markers <- cluster_markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)
pdf(
  file = "results/cluster-markers/DotPlot_top_5_cluster_markers.pdf",
  height = 6,
  width = 12,
  useDingbats = F
)
DotPlot(hearts, features = unique(top_5_cluster_markers$gene)) +
  RotatedAxis() +
  ylab("Cluster") +
  xlab("Marker genes") +
  theme(axis.text.x = element_text(face = "italic"))
dev.off()

# Heatmap
pdf(
  file = "results/cluster-markers/Heatmap_top_30_cluster_markers.pdf",
  height = 8,
  width = 6,
  useDingbats = F
)
top_30_cluster_markers <- cluster_markers %>%
  group_by(cluster) %>%
  top_n(n = 30, wt = avg_log2FC)
heatmap_genes <- unique(top_30_cluster_markers$gene)
DoHeatmap(hearts,
  features = heatmap_genes,
  assay = "Spatial",
  angle = 0,
  size = 3.5
) +
  ylab("Marker genes") +
  scale_color_manual(values = default_colors) +
  theme(axis.text.y = element_blank()) +
  labs(color = "Cluster")
dev.off()

# Gene ontology biological processes over representation test
for (i in unique(cluster_markers$cluster)) {
  genes <- cluster_markers[cluster_markers$cluster == i, ]$gene
  p <- GOBP_fold_enrich(
    rownames(hearts),
    genes,
    paste0(
      "GO Biological Processes Overrepresentation: Cluster ",
      i,
      " Markers"
    )
  )
  pdf(
    file = paste0("results/cluster-markers/GOBP_Overrep_Cluster_", i, ".pdf"),
    height = 6,
    width = 10,
    useDingbats = F
  )
  print(p)
  dev.off()
}
