# spatial-Atgl-KO-24h

This is a spatial transcriptomics project led by [Dr. Katharina Botterman](mailto:katharina.bottermann@hhu.de). It includes 10x Genomics Visium data of cardiac tissue 24 h after ischemia/reperfusion injury from adipose specific (iAdipoQCre) *Atgl*-KO mice.

## Sequencing data
Sequencing data, including fastq files and count matrices will be available upon publication or request.

## Analysis
To recreate the full analysis you can follow the steps below.

### Libraries
The following transcriptomics-specific libraries were used:

* [Seurat](https://satijalab.org/seurat/index.html) v4.0.1
* [clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html) v4.2.1
* [enrichplot](https://bioconductor.org/packages/release/bioc/html/enrichplot.html) v1.14.1
* [org.Mm.eg.db](https://bioconductor.org/packages/release/data/annotation/html/org.Mm.eg.db.html) v3.14.0
* [yulab.utils](https://cran.r-project.org/package=yulab.utils) v0.0.4

They can be installed with the following commands:
```R
# Seurat
remotes::install_version("Seurat", version = "4.0.1")

# clusterProfiler etc.
BiocManager::install("clusterProfiler")
BiocManager::install("enrichplot")
BiocManager::install("org.Mm.eg.db")
install.packages("yulab.utils")
```

Additionally, the following standard R libraries were used:
* dplyr
* ggplot2
* patchwork
* scales
* readr
* ggrepel
* rstudioapi
* shiny

### Instructions
To reproduce the analysis, clone this repository and place the spaceranger outputs inside the project folder under a sub-folder titled data. i.e. `spatial-Atgl-KO-24h/data/spaceranger`

By starting your R session with the R project file, `spatial-Atgl-KO-24h.Rproj`, your working directory will be set to project folder, no matter the location on your machine. This will allow easy reading/writing of data/results using relative paths.

`0-full-analysis.R` will create all necessary directories and run the full analysis in the appropriate order. Each analysis step can also be run individually for better interactivity, starting from `1-clustering.R`.

## Examples
<p align="center">
  <img src="/examples/SpatialFeaturePlot_UMI_count.png" width="1000">
</p>
<p align="center">
  <img src="/examples/SpatialDimPlot.png" width="1000">
</p>

## To-do
* Add abstract and author list at submission
* Make input data available at submission, i.e. gene signatures, cluster markers, reference maps, slide images and spaceranger outputs.
