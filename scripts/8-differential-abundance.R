# Differential abundance

# Make results directories if they do not exist
output_dirs <- c("results",
                 "results/differential-abundance")

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
load("results/objects/obj_annotated.Rdata")

# Proportion
differential_abundance <- (prop.table(table(
  Idents(obj),
  obj$genotype
),
margin = 2
))
differential_abundance <- as.data.frame(differential_abundance)
colnames(differential_abundance) <- c("Cluster", "Sample", "Fraction")
write.csv(differential_abundance,
  file = "results/differential-abundance/differential_abundance.csv",
  row.names = F
)
