# Export raw ATAC peak counts in original SE cell order for the historical LSI.
source("scripts/_paths.R")
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(Matrix)
  library(data.table)
})
se <- readRDS(file.path(sctarf_blood_root(), "scATAC-Healthy-Hematopoiesis-191120.rds"))
out <- file.path(sctarf_blood_root(), "multiomics_module1_pairing/method_B_scSHEFT/scsheft_input")
if (file.exists(file.path(out, "blood_lsi.npy"))) {
  stop("Existing blood_lsi.npy: reuse the matching cell CSV and skip both LSI preparation steps.")
}
counts <- SummarizedExperiment::assay(se)
cells <- colnames(se)
stopifnot(length(cells) == ncol(counts), !anyNA(cells), !anyDuplicated(cells))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
Matrix::writeMM(counts, file.path(out, "atac_counts.mtx"))
data.table::fwrite(data.frame(cell = cells), file.path(out, "atac_cells.csv"))
cat("Exported", nrow(counts), "peaks and", length(cells), "cells in SE order\n")
