# Export basic figures

# Make results directories if they do not exist
output_dirs <- c("results",
                 "results/basic-figures")

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
source("scripts/SpatialFeaturePlotScaled.R")
load("results/objects/obj_annotated.Rdata")
default_colors <- (hue_pal()(7))

# VlnPlot of UMI counts
pdf(
  file = "results/basic-figures/VlnPlot_UMI_count.pdf",
  height = 6,
  width = 6,
  useDingbats = F
)
VlnPlot(obj,
  features = "nCount_Spatial",
  pt.size = 0.25,
  group.by = "genotype"
) +
  NoLegend() +
  ggtitle("") +
  ylab("UMI count") +
  xlab("") +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.title = element_text(size = 16)
  )
dev.off()

# SpatialFeaturePlot of UMI counts
pdf(
  file = "results/basic-figures/SpatialFeaturePlot_UMI_count.pdf",
  height = 6,
  width = 12.75,
  useDingbats = F
)
SpatialFeaturePlotScaled(
  object = obj,
  group = "genotype",
  group.1 = "Control",
  group.2 = "KO",
  feature_of_interest = "nCount_Spatial",
  from.meta.data = TRUE,
  group.1.title = "Control",
  group.2.title = "KO",
  legend.title = "UMI count"
)
dev.off()

# SpatialDimPlot
a <- SpatialDimPlot(obj, images = "Control", pt.size.factor = 1.6) + 
  labs(fill = "Niche") +
  ggtitle("Control") +
  theme(
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.justification = "top",
    legend.key.size = unit(3, "point"),
    title = element_text(size = 16)
  ) +
  guides(fill = guide_legend(override.aes = list(size = 7))) + NoLegend()
b <- SpatialDimPlot(obj, images = "KO", pt.size.factor = 1.6) +
  labs(fill = "Niche") +
  ggtitle("KO") +
  theme( 
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.justification = "top",
    legend.key.size = unit(3, "point"),
    title = element_text(size = 16)
  ) +
  guides(fill = guide_legend(override.aes = list(size = 5)))
pdf(
  file = "results/basic-figures/SpatialDimPlot.pdf",
  height = 6,
  width = 12.75,
  useDingbats = F
)
print(a + b + plot_layout(guides = "collect"))
dev.off()

# Normal DimPlot
p <- DimPlot(obj, pt.size = 1.5) +
  xlab("UMAP-1") +
  ylab("UMAP-2") +
  labs(color = "Niche") +
  theme(
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.justification = "top",
    legend.key.size = unit(3, "point"),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_text(size = 16)
  ) +
  guides(color = guide_legend(override.aes = list(size = 5)))
p2 <- LabelClusters(
  plot = p,
  id = "ident",
  repel = F,
  box = T,
  fill = alpha("white", 0.25),
  size = 8,
  label.r = unit(0.75, "lines"),
  label.size = 0
)
pdf(
  file = "results/basic-figures/DimPlot.pdf",
  height = 6,
  width = 8,
  useDingbats = F
)
print(p2)
dev.off()

# DimPlot grouped by genotype
pdf(
  file = "results/basic-figures/DimPlot_genotype.pdf",
  height = 6,
  width = 8,
  useDingbats = F
)
DimPlot(obj, pt.size = 1.5, group.by = "genotype") +
  xlab("UMAP-1") +
  ylab("UMAP-2") +
  labs(color = "Genotype") +
  theme(
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.justification = "top",
    legend.key.size = unit(3, "point"),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_text(size = 16),
    plot.title = element_blank()
  ) +
  guides(color = guide_legend(override.aes = list(size = 5)))
dev.off()

# Single cluster highlight figures
for (i in levels(Idents(obj))) {
  coi <- default_colors[(as.numeric(i) + 1)]
  a <- SpatialDimPlot(obj,
    images = "Control",
    cells.highlight = WhichCells(obj, idents = i),
    stroke = NA) +
    NoLegend() +
    labs(title = paste0("Niche ", i), subtitle = "Control") +
    theme(
      plot.title = element_text(color = coi, size = 16, face = "bold"),
      plot.subtitle = element_text(size = 14)
    )
  b <- SpatialDimPlot(obj,
    images = "KO",
    cells.highlight = WhichCells(obj, idents = i),
    stroke = NA) +
    NoLegend() +
    labs(title = "", subtitle = "KO") +
    theme(
      plot.title = element_text(color = coi, size = 16, face = "bold"),
      plot.subtitle = element_text(size = 14)
    )
  pdf(
    file = paste0("results/basic-figures/SpatialDimPlot_Cluster_", i, ".pdf"),
    height = 6,
    width = 11.25,
    useDingbats = F
  )
  print(a + b)
  dev.off()
}
