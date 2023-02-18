# spatial-ATGL-KO-24h

This is a spatial transcriptomics project led by [Dr. Katharina Botterman](mailto:katharina.bottermann@hhu.de). It includes 10x Genomics Visium data of cardiac tissue 24 h after ischemia/reperfusion injury from adipocyte-specific (iAdipoQCre) ATGL-KO mice.

## Sequencing data
Sequencing data, including fastq files and spaceranger outputs can be found at ArrayExpress accession E-MTAB-12700, a public link will be available upon publication or request.

## Instructions

To recreate the full analysis you can follow the steps below.

### Data

Clone this repository and place the extracted `spaceranger` outputs from ArrayExpress under a sub-folder titled 'spaceranger' inside 'data'. i.e. `spatial-ATGL-KO-24h/data/spaceranger`

### R & Libraries

R version 4.2.2 was used, packages versions and external sources were recorded with `renv` v0.16.0. To create a library with matching package versions from the ```renv.lock``` file, start an R session in the project directory and run:

```r
renv::restore()
```

More information on ```renv``` can be found [here](https://rstudio.github.io/renv/articles/renv.html).

### Analysis

By starting your R session with the R project file, `spatial-ATGL-KO-24h.Rproj`, your working directory will be set to project folder, no matter the location on your machine. This will allow easy reading/writing of data/results using relative paths.

`0-full-analysis.R` will run the full analysis in the appropriate order. Each analysis step can also be run individually for better interactivity, starting from `1-preprocessing.R`.

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
