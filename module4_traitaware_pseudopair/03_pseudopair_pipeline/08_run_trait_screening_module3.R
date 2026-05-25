# module4 trait screening: run module3 v2 KNN + RWR/TRS
# -----------------------------------------------------------------------------
# Thin wrapper around:
#   03_pseudopair_pipeline/03_run_module3_v2_knn_rwr.R
#
# Run in RStudio:
#   setwd("D:/RStudio/work/module4_traitaware_pseudopair")
#   source("03_pseudopair_pipeline/08_run_trait_screening_module3.R")
# -----------------------------------------------------------------------------

ROOT <- "D:/RStudio/work/module4_traitaware_pseudopair"

if (!exists("MODULE3_COMPAT_OUT_DIR")) {
  MODULE3_COMPAT_OUT_DIR <- file.path(ROOT, "03_pseudopair_pipeline/output/module3_v2_trait_screening_ge0p80")
}
if (!exists("TRAITS_TO_RUN")) TRAITS_TO_RUN <- c("MCH", "RDW", "RBC")
if (!exists("ATAC_W_LIST")) ATAC_W_LIST <- c(0.7, 1.0)

cat(">>> Run module3 v2 for trait screening\n\n")
cat("Traits: ", paste(TRAITS_TO_RUN, collapse = ", "), "\n", sep = "")
cat("Output dir: ", MODULE3_COMPAT_OUT_DIR, "\n\n", sep = "")

source(file.path(ROOT, "03_pseudopair_pipeline/03_run_module3_v2_knn_rwr.R"), local = TRUE)
