# In this script we will run the pre-processing for each sample
#if (!"BiocManager" %in% installed.packages()) install.packages("BiocManager")

x <- c("Seurat", "cowplot", "tidyverse", "hdf5r"#,"ggpointdensity", "viridis"
       )
#BiocManager::install(x)

# Load libraries
invisible(lapply(x, library, character.only = TRUE))
#---------------------------------------------------------------

# Set up output dirs
output_dirs <- c("results",
                 "results/01_preprocessing",
                 "results/01_objects")

for (i in output_dirs) {
  if (!dir.exists(i)) {
    dir.create(i)
    print(paste0("made ", i, " directory"))
  } else {
    print(paste0(i, " directory already exists."))
  }
}
#---------------------------------------------------------------

# Load metadata
metadata <- read.csv("data/sample-metadata.csv")
#metadata <- metadata %>% filter(mouse_line == "1955")

# Load and process each spatial object
obj_list <- list()

for (i in metadata$sample) {
  
  # Isolate metadata for each sample
  df <- metadata[metadata$sample==i,]
  
  # Create the Seurat object
  obj <- Load10X_Spatial(
    data.dir = paste0("data/spaceranger/", df$bmfz_sample, "/outs"),
    filename = "filtered_feature_bc_matrix.h5",
    assay = "Spatial",
    slice = df$sample,
    filter.matrix = T)
  
  # Add in metadata
  obj[["bmfz_sample"]] <- df$bmfz_sample
  obj[["genotype"]] <- df$sample
  obj[["genotype_long"]] <- df$genotype_long
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
  
  # Plot pre-filter QC
  qc_p1 <- obj@meta.data %>%
    ggplot(aes(x = nCount_Spatial, y = nFeature_Spatial)) +
    geom_point() +
    #geom_pointdensity() +
    #scale_color_viridis() +
    theme_classic() +
    ggtitle(paste0("nspots ", ncol(obj)))
  
  qc_p2 <- obj@meta.data %>%
    ggplot(aes(x = nCount_Spatial, y = percent.mt)) +
    geom_point() +
    theme(aspect.ratio=1) + 
    #geom_pointdensity() +
    #scale_color_viridis() +
    theme_classic() +
    ggtitle(paste0("nspots ", ncol(obj)))
  
  qc_panel <- cowplot::plot_grid(qc_p1, qc_p2, ncol = 2, align = "hv")
  
  slide_qc_p <- SpatialFeaturePlot(obj,
                                   features = c("nCount_Spatial", 
                                                "nFeature_Spatial", 
                                                "percent.mt"),
                                   ncol = 3) &
    theme(legend.text=element_text(size=7),
          legend.title=element_text(size=7) )
  
  
  qc_panel_pre <- cowplot::plot_grid(qc_panel, slide_qc_p, 
                                   nrow = 2, ncol = 1, 
                                   rel_heights = c(0.5, 0.5))
  #------------------------------
  
  # Filter genes expressed in less than 10 spots
  obj <- obj[rowSums(GetAssayData(obj, assay = "Spatial") > 0) > 10, ]
  
  # Re-calculate genes and reads
  obj$nFeature_Spatial_filt <- colSums(GetAssayData(obj, assay = "Spatial") > 0)
  obj$nCount_Spatial_filt <- colSums(GetAssayData(obj, assay = "Spatial"))
  
  # Plot filtered QC
  qc_p1_filt <- obj@meta.data %>%
    ggplot(aes(x = nCount_Spatial_filt, y = nFeature_Spatial_filt)) +
    geom_point() +
    #geom_pointdensity() +
    theme(aspect.ratio=1) +
    #scale_color_viridis() +
    theme_classic() +
    geom_vline(xintercept = 300) +
    geom_hline(yintercept = 500) +
    ggtitle(paste0("nspots ", ncol(obj)))
  
  qc_p2_filt <- obj@meta.data %>%
    ggplot(aes(x = nFeature_Spatial_filt, y = percent.mt)) +
    geom_point() +
    #geom_pointdensity() +
    theme(aspect.ratio=1) +
    #scale_color_viridis() +
    theme_classic() +
    geom_vline(xintercept = 300) +
    ggtitle(paste0("nspots ", ncol(obj)))

  qc_panel_filt <- cowplot::plot_grid(qc_p1_filt,
                                      qc_p2_filt,
                                      ncol = 2,
                                      align = "hv")
  
  # Filter spots that should have at least 1 cell (Same as HCA and Kuppe et al.)
  obj <- subset(obj, subset = nFeature_Spatial_filt > 300 &
                  nCount_Spatial_filt > 500)
  #--------------------------------------------------
  
  # Export plots
  pdf(file = paste0("results/01_preprocessing/qc_panel_pre_",
                    i,
                    ".pdf"),
      useDingbats = F)
  print(qc_panel_pre)
  dev.off()
  
  pdf(file = paste0("results/01_preprocessing/qc_panel_filt_",
                    i,
                    ".pdf"),
      useDingbats = F)
  print(qc_panel_filt)
  dev.off()
  
  # Append to list
  obj_list[i] <- obj
  
}

# Save list
save(obj_list, file = "results/01_objects/01_obj_list.Rdata")
