# Over-representation analysis of differentially expressed genes

# Make results directories if they do not exist
output_dirs <- c("results",
                 "results/ora")

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
library(knitr)
library(kableExtra)
source("scripts/ORA_up.R")
source("scripts/ORA_down.R")
load("results/objects/obj_annotated.Rdata")
default_colors <- (hue_pal()(7))

# Load dge data
deg_niches_sig <- read.csv(file = "results/differential-gene-expression/deg_niches_sig.csv")

# Background genes
detected <- rownames(obj)

# Loop through niches and export upregulated ORA analyses
for (i in unique(deg_niches_sig$cluster)) {
  pdf(file = paste0("results/ora/ORA_up_", i, ".pdf"),
      #height = 4,
      width = 7,
      useDingbats = F)
  ORA_up(deg_niches_sig, i, detected, "Upregulated in KO")
  dev.off()
}

# Loop through clusters and export downregulated ORA analyses
for (i in unique(deg_niches_sig$cluster)) {
  pdf(file = paste0("results/ora/ORA_down_", i, ".pdf"),
      #height = 4,
      width = 7,
      useDingbats = F)
  ORA_down(deg_niches_sig, i, detected, "Downregulated in KO")
  dev.off()
}
