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
                 "results/01_preprocessing_mine",
                 "results/01_objects_mine")

for (i in output_dirs) {
  if (!dir.exists(i)) {
    dir.create(i)
    print(paste0("made ", i, " directory"))
  } else {
    print(paste0(i, " directory already exists."))
  }
}
#---------------------------------------------------------------

### FUNCTIONS
## 2023-11-10 exit script without error
f.exit_no_error <- function() { invokeRestart("abort") }

## 2023-11-10 check if the metadata has the necessary column names
# tested with:
##  [v] all necessary columns present among others
##  [v] only some of the columns present
##  [v] only some columns present, and some columns have no name or named NA
##  [v] all columns present, and some columns have no name or named NA

# to still do:
##  [] exit the run gracefully with message
##  [v] print the vector within the function
f.check_metadata_cols <- function(df.metadata) {
  v.mandatory_metadata_cols <- c("sample", "mouse_line", "genotype_long", "bmfz_sample")
  
  message_true <- "All necessary columns exist in metadata"
  message_false <- paste(c("Check metadata. The necessary column names are: ", v.mandatory_metadata_cols), collapse = "  ")
  
  ifelse(all(v.mandatory_metadata_cols %in% names(df.metadata)), 
         message_true,
         #message_false
         #f.exit_no_error()}
         stop(message_false)
  )
}

## 2023-11-15 process string; 
## it doesn't work, same issue as when using rlang::erquo directly into the plot
# f.str_short <- function(func_output, delimiter) {
#   defuse_func <- rlang::enquo(func_output)
#   make_string <- toString(defuse_func)
#   final_x <- gsub(".*\\$", "", make_string) 
# }

## 2023-11-10 dotplots for qc
f.dotplot_qc <- function(seurat_obj, x_data, y_data, x_intercept=0, y_intercept=0) {
  xto_str <- toString(rlang::enquo(x_data))
  yto_str <- toString(rlang::enquo(y_data))
  x_ax <- gsub(".*\\$","",xto_str)
  y_ax <- gsub(".*\\$","",yto_str)
  
  seurat_obj@meta.data %>%
    ggplot(aes(x=x_data, y=y_data)) +
    geom_point() +
    #geom_pointdensity() +
    #scale_color_viridis() +
    theme(aspect.ratio=1) +
    theme_classic() +
    geom_vline(xintercept = x_intercept) +
    geom_hline(yintercept = y_intercept) +
    labs(title=seurat_obj@meta.data$bmfz_sample, x=x_ax, y=y_ax, subtitle = paste0("nspots ", ncol(seurat_obj))) 
}
#############################################################################

## 2023-11-10 preprocess a spatial transcriptomics dataset
## from loading the data to simple filtering: 
## nCount and nFeature filtering
## remove genes expressed in less than N spots (nr_spots)
## with qc plots before and after filtering

# to add:
## [] obj
## [] i

f.preprocess_spatial_ds <- function(str.spaceranger_dir, df.metadata, i, nr_spots=10, 
                                    nFeat_threshold=300, nCount_threshold=500) {
  # Isolate metadata for each sample
  df <- df.metadata[df.metadata$bmfz_sample==i,]
  
  # Create the Seurat object
  obj <- Load10X_Spatial(
    data.dir = paste0(str.spaceranger_dir, df$bmfz_sample, "/outs"),
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
  qc_p1 <- f.dotplot_qc(obj, obj@meta.data$nCount_Spatial, obj@meta.data$nFeature_Spatial, 0, 0)
  qc_p2 <- f.dotplot_qc(obj, obj@meta.data$nCount_Spatial, obj@meta.data$percent.mt, 0, 0)
  qc_panel <- cowplot::plot_grid(qc_p1, qc_p2, 
                                 ncol = 2, align = "hv")
  
  slide_qc_p <- SpatialFeaturePlot(obj, 
                                   ncol = 3, 
                                   features = c("nCount_Spatial", 
                                                "nFeature_Spatial", 
                                                "percent.mt")
  ) &
    theme(legend.text=element_text(size=7),
          legend.title=element_text(size=7) 
    )
  
  qc_panel_pre <- cowplot::plot_grid(qc_panel, slide_qc_p, 
                                     nrow = 2, ncol = 1, 
                                     rel_heights = c(0.5, 0.5))
  
  # Filter genes expressed in less than 10 spots
  obj <- obj[rowSums(GetAssayData(obj, assay = "Spatial", layer="counts") > 0) > as.integer(nr_spots), ]
  
  # Re-calculate genes and reads
  obj$nFeature_Spatial_filt <- colSums(GetAssayData(obj, assay = "Spatial", layer="counts") > 0)
  obj$nCount_Spatial_filt <- colSums(GetAssayData(obj, assay = "Spatial", layer="counts"))
  
  # Plot filtered QC
  qc_p1_filt <- f.dotplot_qc(obj, obj@meta.data$nCount_Spatial_filt, obj@meta.data$nFeature_Spatial_filt, 
                             nFeat_threshold, nCount_threshold)
  qc_p2_filt <- f.dotplot_qc(obj, obj@meta.data$nFeature_Spatial_filt, obj@meta.data$percent.mt, 
                             nFeat_threshold, 0)
  qc_panel_filt <- cowplot::plot_grid(qc_p1_filt,
                                      qc_p2_filt,
                                      ncol = 2,
                                      align = "hv")
  
  # Filter spots that should have at least 1 cell (Same as HCA and Kuppe et al.)
  obj <- subset(obj, 
                subset = nFeature_Spatial_filt > nFeat_threshold &
                         nCount_Spatial_filt > nCount_threshold)

  # Export plots
  pdf(file = paste0("results/01_preprocessing_mine/qc_panel_pre_",
                    i,
                    ".pdf"),
      useDingbats = F)
  print(qc_panel_pre)
  dev.off()
  
  pdf(file = paste0("results/01_preprocessing_mine/qc_panel_filt_",
                    i,
                    ".pdf"),
      useDingbats = F)
  print(qc_panel_filt)
  dev.off()
  
  return(obj) 
}
#--------------------------------------------------
  
#------------------------------------
#------------------------------------
#------------------------------------
# Load metadata
metadata <- read.csv("metadata/2023.11.13_metadata_all_samples.csv")
#---------------------------------------------------------------

## CALLS
# check metadata
f.check_metadata_cols(metadata)

# Load and process each spatial object
obj_list <- list()
str.spaceranger_dir <- "data/spaceranger/"

for (i in metadata$bmfz_sample) {
  obj <- f.preprocess_spatial_ds(str.spaceranger_dir, metadata, i, nr_spots=10,
                                 nFeat_threshold=300, nCount_threshold=500)

  # Append to list
  obj_list[i] <- obj
 }

# obj_list <- lapply(list(metadata$bmfz_sample), FUN = function(x){
#   x <- f.preprocess_spatial_ds(str.spaceranger_dir, metadata, i, nr_spots=10, 
#                                nFeat_threshold=300, nCount_threshold=500)
# })


# Save list
save(obj_list, file = "results/01_objects_mine/01_obj_list.Rdata")
