# In this script we will annotate the molecular niches

# Load libraries
library(Seurat)
library(ggplot2)
library(ape)

# Set up output dirs
output_dirs <- c("results",
                 "results/annotation")
for (i in output_dirs) {
  if (!dir.exists(i)) {
    dir.create(i)
    print(paste0("made ", i, " directory"))
  } else {
    print(paste0(i, " directory already exists."))
  }
}

# Load object
load(file = "results/objects/obj_integrated.Rdata")

# Set Idents to the optimum resolution
Idents(obj) <- "integrated_snn_res.0.4"

# Calculate niche markers
niche_markers <- FindAllMarkers(obj,
                                assay = "Spatial",
                                only.pos = T)
write.csv(niche_markers,
          file = "results/annotation/niche_markers.csv",
          row.names = F)
top5 <- niche_markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

# Dendrogram
obj <- BuildClusterTree(obj, dims = 1:25)
PlotClusterTreeDJG <- function(object, ...) {
  if (is.null(x = Tool(object = object, slot = "BuildClusterTree"))) {
    stop("Phylogenetic tree does not exist, build using BuildClusterTree")
  }
  data.tree <- Tool(object = object, slot = "BuildClusterTree")
  plot.phylo(x = data.tree, font = 2, direction = "rightwards", label.offset = 2, edge.width = 1)
}
pdf(file = "results/annotation/Dendrogram.pdf",
    useDingbats = F)
PlotClusterTreeDJG(obj)
dev.off()

###############################################################################
# Annotation notes
###############################################################################

#-0-Border zone (BZ-1): Nppa, Nppb, Ankrd1, Similar to Calcagno et al. 2022 BZ-1

#-1-Border zone (BZ-2): Flnc, Xirp2, Similar to Calcagno et al. 2022 BZ-2,
#   very thing layer, similar to Alcagno et al. 2022 BZ-2

#-6-Border zone (BZ-3): Nppa, Nppb, Ankrd1, Similar to Calcagno et al. 2022 BZ-1
#   Unique markers  Mgp, Tpm1, Myh7

#-3-Remote zone (RZ): Fabp3, Myh6, Lpl, clearly surviving myocardium

#-4-Ischemic zone (IZ-1): Ccl6, Ccl9, Arg1, Lyz2, classic inflammatory response
#-2-Ischemic zone (IZ-3): Col1a1, Col1a2, fibrotic ischemic zone, dendrogram
#   shows close relationship to ischemic clusters 4,5. Likely fibrotic response
#   from suture.
#-5-Ischemic zone RBC (IZ-2): Hbb-bs, Hba-a2, Hba-a1, hemoglobins from erythrocytes
#   secondary markers are S100a8, S100a9, Mmp9, Ccl6

###############################################################################

# Rename Idents to annotations
obj <- RenameIdents(obj,
                    "0" = "BZ-1",
                    "1" = "BZ-2",
                    "2" = "IZ-3",
                    "3" = "RZ",
                    "4" = "IZ-1",
                    "5" = "IZ-2",
                    "6" = "BZ-3")

# Store renamed idents as a new meta data column, set as Idents
obj@meta.data$basic_annotation <- Idents(obj)

# Refactor annotation levels
source("scripts/dimplotlevels.R")
obj@meta.data$basic_annotation <- factor(obj@meta.data$basic_annotation,
                                         levels = dimplotlevels)
DimPlot(obj,
        group.by = "basic_annotation",
        label = T,
        repel = T)

# Set Idents as re-factored basic_annotation identities
Idents(obj) <- "basic_annotation"

# Save object with basic annotations
save(obj, file = "results/objects/obj_annotated.Rdata")
