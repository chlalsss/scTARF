# module4: rerun neighborhood-level scDRS FDR with selected RNA grouping
# -----------------------------------------------------------------------------
# Why this wrapper exists:
#   The first neighborhood-FDR run used a clustering resolution that was too high,
#   so most RNA "groups" became single cells. The granularity sweep showed that
#   k = 30 and Leiden resolution = 0.005 gives a more reasonable grouping:
#     - about 94 unsupervised RNA neighborhoods
#     - median group size about 1035 cells
#     - only about 6.4% of groups have fewer than 30 cells
#
# Method boundary:
#   - Cell type is NOT used to build RNA neighborhoods.
#   - Cell type is NOT used for scDRS FDR, mapping, or candidate seed selection.
#   - Cell type is only used inside the sourced diagnostic script for post hoc
#     checking after the unsupervised/statistical steps are done.
#
# Run in RStudio:
#   setwd("D:/RStudio/work/module4_traitaware_pseudopair")
#   source("04_evaluation/run_neighborhood_fdr_k30_res0p005.R")
# -----------------------------------------------------------------------------

ROOT <- "D:/RStudio/work/module4_traitaware_pseudopair"

# These two values come from:
#   04_evaluation/tune_rna_neighborhood_granularity.R
RNA_KNN_K <- 30L
LEIDEN_RESOLUTION <- 0.005

# Keep the same scDRS group-FDR cutoffs as the diagnostic script.
GROUP_FDR_THRESHOLDS <- c(0.05, 0.10, 0.20)

# Keep the downstream diagnostic comparable to the old top-1% seed size:
# 350 = 1% of 35,000 ATAC cells.
TOP_N_ATAC <- 350L

cat(">>> Run neighborhood-level scDRS FDR with selected RNA grouping\n\n")
cat("RNA_KNN_K = ", RNA_KNN_K, "\n", sep = "")
cat("LEIDEN_RESOLUTION = ", LEIDEN_RESOLUTION, "\n", sep = "")
cat("GROUP_FDR_THRESHOLDS = ", paste(GROUP_FDR_THRESHOLDS, collapse = ", "), "\n", sep = "")
cat("TOP_N_ATAC = ", TOP_N_ATAC, "\n\n", sep = "")

source(file.path(ROOT, "04_evaluation/diagnose_scdrs_neighborhood_fdr.R"))
