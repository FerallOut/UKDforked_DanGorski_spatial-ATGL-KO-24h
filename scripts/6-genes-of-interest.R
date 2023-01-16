# Visualize genes of interest

# Make results directories if they do not exist
output_dirs <- c("results",
                 "results/genes-of-interest")

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
source("scripts/SpatialFeaturePlotScaled.R")
load("results/objects/obj_annotated.Rdata")
default_colors <- (hue_pal()(7))

# Genes of interest
genes <- c(
  "Plin5", "Plin2", "Pdk4", "Pdp2", "Cpt1a", "Acacb",
  "Angptl4", "Slc2a4", "Cd36", "Lpl", "Acadm",
  "Ppargc1a", "Has1", "Has2", "Has3", "Vcan", "Acan",
  "Ncan", "Bcan", "Dcn", "Bgn", "Col1a1", "Fn1", "Il2",
  "Cemip", "Hyal1", "Hyal2", "Hmmr", "Adamts5", "Dmkn",
  "Tbx18", "Wt1", "Aqp1", "Lcn2", "Serpinb2", "Krt19",
  "Krt8", "Upk3b", "Clu", "Gpm6a", "Bnc1", "Saa3",
  "Tcf21", "Pdgfra", "Postn", "Cthrc1", "Pecam1",
  "Adgre1", "Cx3cr1", "S100a8", "S100a9", "Ccr1", "Csf3r",
  "Cd3e", "Cd79a", "Ppard", "Pparg", "Lipa", "Plin3",
  "Slc2a1", "Myh7", "Map1lc3b", "Ccl6", "Ccl9", "Thbs1",
  "Ccl2", "Il1b", "Ccl7", "Pf4", "Cxcl2", "Ppbp", "Ednrb",
  "Ccl7", "Ccl8", "Adipor1", "Adipor2", "Lepr", "Leprot",
  "Leprotl1", "Lep", "Snx1", "Snx2", "Snx4", "Il4ra", "Tmsb4",
  "Fabp3", "Pnpla2", "Tgfb2"
  )

# Prune genes of interest, making sure they are expressed in the object
expressed_genes <- list()
for (i in genes) {
  
  # Error handling
  possibleError <- tryCatch(
  FetchData(obj, vars = i),
  error=function(e) e)
  
  # If no error (exists in data set), add it to expressed genes list
  if(!inherits(possibleError, "error")) {
  expressed_genes[[i]] <- i
  } else (
    print(paste0(i, " - not expressed, removing from list."))
  )
}
expressed_genes <- unlist(expressed_genes)

# VlnPlots
for (i in expressed_genes) {
  pdf(
    file = paste0("results/genes-of-interest/VlnPlot_", i, ".pdf"),
    height = 6,
    width = 8,
    useDingbats = F
  )
  p <- VlnPlot(obj,
               features = i,
               split.by = "genotype",
               assay = "Spatial",
               split.plot = TRUE
  ) +
    theme(
      plot.title = element_text(face = "italic", size = 20),
      axis.text.x = element_text(size = 14, angle = 0, hjust = 0.5),
      axis.title.x = element_blank(),
      legend.position = "right"
    )
  print(p)
  dev.off()
  print(paste0(i, " - done."))
  }

# SpatialFeaturePlots
for (i in expressed_genes) {
  pdf(
    file = paste0("results/genes-of-interest/SpatialFeaturePlot_", i, ".pdf"),
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
