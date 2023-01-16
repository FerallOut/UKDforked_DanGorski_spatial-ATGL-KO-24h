# Identify niche markers

# Make results directories if they do not exist
output_dirs <- c("results",
                 "results/niche-markers")

for (i in output_dirs) {
  if (!dir.exists(i)) {
    dir.create(i)
    print(paste0("made ", i, " directory"))
  } else {
    print(paste0(i, " directory already exists."))
  }
}

# Load libraries, functions and objects
library(Seurat)
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
load("results/objects/obj_annotated.Rdata")
default_colors <- (hue_pal()(7))

# Find niche markers
obj <- PrepSCTFindMarkers(obj, assay = "SCT", verbose = TRUE)
niche_markers <- FindAllMarkers(obj, assay = "SCT", only.pos = T)
colnames(niche_markers)[6] <- "niche" # rename cluster
top_10_niche_markers <- niche_markers %>%
  group_by(niche) %>%
  top_n(n = 10, wt = avg_log2FC)
top_100_niche_markers <- niche_markers %>%
  group_by(niche) %>%
  top_n(n = 100, wt = avg_log2FC)
write.csv(niche_markers,
  file = "results/niche-markers/niche_markers.csv",
  row.names = F
)

# SpatialFeaturePlots of top 10 niche markers
genes <- top_10_niche_markers$gene
for (i in genes) {
  pdf(
    file = paste0("results/niche-markers/SpatialFeaturePlot_", i, ".pdf"),
    height = 6,
    width = 12.75,
    useDingbats = F
  )
  SpatialFeaturePlotScaled(
    object = obj,
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
top_5_niche_markers <- niche_markers %>%
  group_by(niche) %>%
  top_n(n = 5, wt = avg_log2FC)
pdf(
  file = "results/niche-markers/DotPlot_top_5_niche_markers.pdf",
  height = 6,
  width = 12,
  useDingbats = F
)
DotPlot(obj, features = unique(top_5_niche_markers$gene)) +
  RotatedAxis() +
  ylab("Niche") +
  xlab("Marker genes") +
  theme(axis.text.x = element_text(face = "italic"))
dev.off()

# Heatmap
pdf(
  file = "results/niche-markers/Heatmap_top_30_niche_markers.pdf",
  height = 8,
  width = 6,
  useDingbats = F
)
top_30_niche_markers <- niche_markers %>%
  group_by(niche) %>%
  top_n(n = 30, wt = avg_log2FC)
heatmap_genes <- unique(top_30_niche_markers$gene)
DoHeatmap(obj,
  features = heatmap_genes,
  assay = "Spatial",
  angle = 0,
  size = 3.5
) +
  ylab("Marker genes") +
  scale_color_manual(values = default_colors) +
  theme(axis.text.y = element_blank()) +
  labs(color = "niche")
dev.off()

# Gene ontology biological processes over representation test
for (i in unique(niche_markers$niche)) {
  genes <- niche_markers[niche_markers$niche == i, ]$gene
  p <- GOBP_fold_enrich(
    rownames(obj),
    genes,
    paste0(
      "GO Biological Processes Overrepresentation: niche ",
      i,
      " Markers"
    )
  )
  pdf(
    file = paste0("results/niche-markers/GOBP_Overrep_Niche_", i, ".pdf"),
    height = 6,
    width = 10,
    useDingbats = F
  )
  print(p)
  dev.off()
}
