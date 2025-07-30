args <- commandArgs(trailingOnly = TRUE)

output_file <- args[1]  # First argument: output file path
dataset_name <- args[2]  # Second argument: dataset name

cat("Output File:", output_file, "\n")
cat("Dataset Name:", dataset_name, "\n")

suppressPackageStartupMessages({
    library(scRNAseq)
    library(matrixStats)
    library(AnnotationHub)
    library(BiocFileCache)
    library(dbplyr)
})

library(ExperimentHub)
library(AnnotationHub)

ah <- AnnotationHub(cache = tempfile())

eh <- ExperimentHub(cache = tempfile())

sce <- do.call(dataset_name, list(ensembl = TRUE))

save(sce, file = output_file)
cat("Dataset downloaded and saved to:", output_file, "\n")