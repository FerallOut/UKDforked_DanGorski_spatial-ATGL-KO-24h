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
library(Seurat)
library(ggplot2)
library(patchwork)
library(scales)
library(dplyr)
library(readr)
library(rstatix)
library(ggpubr)
library(orthogene)
source("scripts/SpatialFeaturePlotScaledSig.R")
load("results/objects/obj_annotated.Rdata")
default_colors <- (hue_pal()(7))

# Read in vCM marker genes from Kruppe et al. 2022, convert to mouse gene
# orthologs and score cells
Kruppe_et_al_vCM_markers <- read_csv("data/Kruppe_et_al_vCM_markers.csv")
for (i in colnames(Kruppe_et_al_vCM_markers)) {
  df <- convert_orthologs(gene_df = Kruppe_et_al_vCM_markers,
                                 gene_input = i,
                                 input_species = "human",
                                 output_species = "mouse",
                                 method = "gprofiler",
                                 drop_nonorths = T,
                                 non121_strategy = "drop_both_species")
  obj <- AddModuleScore(obj,
                        features = list(row.names(df)),
                        name = i,
                        assay = "Spatial")
}

# Read in stromal cell population marker genes from Forte et al. 2020, 
# scoring cells
Forte_et_al_stromal_markers <- read_csv("data/Forte_et_al_stromal_markers.csv")
for (i in colnames(Forte_et_al_stromal_markers)) {
  genes <- Forte_et_al_stromal_markers[i]
  colnames(genes) <- "gene"
  genes <- genes$gene
  genes <- unique(genes)
  genes <- genes[!is.na(genes)]
  obj <- AddModuleScore(obj,
                        features = list(genes),
                        name = i,
                        assay = "Spatial")
}

# Reading in gene sets, score cells
gene_ontology <- read.csv("data/gene_signatures.csv")
for (i in colnames(gene_ontology)) {
  genes <- gene_ontology[i]
  colnames(genes) <- "gene"
  genes <- genes$gene
  genes <- unique(genes)
  genes <- genes[!is.na(genes)]
  obj <- AddModuleScore(obj,
    features = list(genes),
    name = i,
    assay = "Spatial"
  )
}

# removing extra "1" added by Seurat
meta_data_columns <- colnames(obj@meta.data)
all <- as.numeric(length(meta_data_columns))
new <- as.numeric(length(colnames(gene_ontology))) +
  as.numeric(length(colnames(Kruppe_et_al_vCM_markers))) +
  as.numeric(length(colnames(Forte_et_al_stromal_markers))) 
original <- all - new
left <- meta_data_columns[1:original]
right <- gsub(".{1}$", "", meta_data_columns)[(original + 1):all]
new_col_names <- c(left, right)
colnames(obj@meta.data) <- new_col_names
signatures <- tail(colnames(obj@meta.data), new)

# VlnPlots
for (i in signatures) {
  pdf(
    file = paste0("results/signatures-of-interest/VlnPlot_", i, ".pdf"),
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
      plot.title = element_text(size = 20),
      axis.text.x = element_text(size = 14, angle = 0, hjust = 0.5),
      axis.title.x = element_blank(),
      legend.position = "right"
    )
  print(p)
  dev.off()
  print(paste0(i, " - done."))
}

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
    object = obj,
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

# Wilcoxon sum rank test with bonferroni correction of gene signatures

# Initiate list
signature_stats <- list()

# Loop through each signature
for (i in signatures) {
  
  # as.formula for looping through left hand side of formula
  formula <- as.formula(paste0(i, "~ genotype"))
  print(formula)
  
  # wilcox test
  stat_test <- obj@meta.data %>%
    group_by(basic_annotation) %>% 
    wilcox_test(formula) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")
  
  # append to list
  signature_stats[[i]] <- stat_test
  
  # add xy coordinates to stats 
  stat_test <- stat_test %>% add_xy_position(x = "basic_annotation")
  
  # boxplot with statistical sig level
  p <- ggplot(data = obj@meta.data,
              aes_string(x = "basic_annotation", y = i)) +
    geom_boxplot(aes(fill=genotype),
                 outlier.size = .1) + 
    ggtitle(i) +
    theme(
      plot.title = element_text(size = 16, hjust = 0.5),
      axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
      axis.title = element_blank(),
      axis.line = element_line(linewidth = .25),
      panel.background = element_blank(),
      legend.key = element_rect(fill = "transparent"),
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(size = 10)) +
    stat_pvalue_manual(stat_test, label = "p.adj.signif", tip.length = 0)
  
  # output boxplot
  pdf(
    file = paste0("results/signatures-of-interest/BoxPlot_", i, ".pdf"),
    height = 6,
    width = 8,
    useDingbats = F
  )
  print(p)
  dev.off()
}
signature_stats <- do.call(rbind, signature_stats)
row.names(signature_stats) <- NULL
write_csv(signature_stats, file = "results/signatures-of-interest/signature_stats.csv")
