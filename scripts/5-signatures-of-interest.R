# Visualize gene signatures of interest

# Make results directories if they do not exist
output_dirs <- c("results",
                 "results/signatures-of-interest")

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
library(readr)
source("scripts/SpatialFeaturePlotScaledSig.R")
load("results/objects/hearts.Rdata")
default_colors <- (hue_pal()(7))

# Reading in gene sets, score cells
gene_ontology <- read.csv("data/gene_signatures.csv")
for (i in colnames(gene_ontology)) {
  genes <- gene_ontology[i]
  colnames(genes) <- "gene"
  genes <- genes$gene
  genes <- unique(genes)
  genes <- genes[!is.na(genes)]
  hearts <- AddModuleScore(hearts,
    features = list(genes),
    name = i,
    assay = "Spatial"
  )
}
R
# removing extra "1" added by Seurat
meta_data_columns <- colnames(hearts@meta.data)
all <- as.numeric(length(meta_data_columns))
new <- as.numeric(length(colnames(gene_ontology)))
original <- all - new
left <- meta_data_columns[1:original]
right <- gsub(".{1}$", "", meta_data_columns)[(original + 1):all]
new_col_names <- c(left, right)
colnames(hearts@meta.data) <- new_col_names
signatures <- tail(colnames(hearts@meta.data), length(colnames(gene_ontology)))

# SpatialFeaturePlots
for (i in signatures) {
  pdf(
    file = paste0(
      "results/signatures-of-interest/SpatialFeaturePlot_",
      i,
      ".pdf"
    ),
    height = 6,
    width = 11,
    useDingbats = F
  )
  SpatialFeaturePlotScaledSig(
    object = hearts,
    group = "genotype",
    group.1 = "Control",
    group.2 = "KO",
    feature_of_interest = i,
    group.1.title = "Control",
    group.2.title = "KO",
    legend.title = i
  )
  dev.off()
}
