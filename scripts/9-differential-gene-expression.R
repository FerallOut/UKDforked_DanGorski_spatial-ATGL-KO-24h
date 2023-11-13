#if (!"BiocManager" %in% installed.packages()) install.packages("BiocManager")
# Load libraries, functions and objects

x <- c("Seurat", "ggplot2", "patchwork", "scales", "dplyr", "yulab.utils", "clusterProfiler", "org.Mm.eg.db", "enrichplot", "ggrepel")
#BiocManager::install(x)

# Load libraries
invisible(lapply(x, library, character.only = TRUE))
#---------------------------------------------------------

# Differential gene expression, each niche, Control v KO

# Make results directories if they do not exist
output_dirs <- c("results",
                 "results/09_differential-gene-expression")

for (i in output_dirs) {
  if (!dir.exists(i)) {
    dir.create(i)
    print(paste0("made ", i, " directory"))
  } else {
    print(paste0(i, " directory already exists."))
  }
}

# Load libraries, functions and objects
source("scripts/SpatialFeaturePlotScaled.R")
source("scripts/VolcanoPlot.R")
load("results/01_objects/03_obj_annotated.Rdata")
default_colors <- (hue_pal()(7))

# No significance threshold
de_genes <- list()
for (i in levels(Idents(obj))) {
  results <- FindMarkers(obj,
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
deg_niches <- do.call(rbind, de_genes)
rownames(deg_niches) <- NULL
write.csv(deg_niches,
  file = "results/09_differential-gene-expression/deg_niches.csv",
  row.names = F
)

# With significance threshold
de_genes <- list()
for (i in levels(Idents(obj))) {
  results <- FindMarkers(obj,
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
deg_niches_sig <- do.call(rbind, de_genes)
deg_niches_sig <- deg_niches_sig %>% filter(p_val_adj <= 0.01)
rownames(deg_niches_sig) <- NULL
write.csv(deg_niches_sig,
          file = "results/09_differential-gene-expression/deg_niches_sig.csv",
          row.names = F
)

# Volcano plot
source("scripts/VolcanoPlot.R")
for (i in unique(deg_niches$cluster)) {
  pdf(
    file = paste0(
      "results/09_differential-gene-expression/VolcanoPlot_Niche_",
      i,
      "_Control_vs_KO.pdf"
    ),
    height = 6,
    width = 8,
    useDingbats = F
  )
  VolcanoPlot(
    df = deg_niches,
    identity = i,
    title = paste0("DEG Niche: ", i, " - KO/Control")
  )
  dev.off()
}
