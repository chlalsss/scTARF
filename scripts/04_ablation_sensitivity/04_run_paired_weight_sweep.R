# Run after the paired stage. Reuse the original graph, RWR and evaluation code.
# Existing graph files are reused; choose a fresh data/output root for a clean run.
ATAC_W_LIST <- c(0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
source("scripts/01_paired_analysis/15_build_weighted_knn.R")
source("scripts/01_paired_analysis/17_rwr_trs.R")
source("scripts/01_paired_analysis/18_evaluate_auroc.R")
