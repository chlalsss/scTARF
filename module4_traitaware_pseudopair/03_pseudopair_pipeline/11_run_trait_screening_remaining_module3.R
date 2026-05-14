# module4 trait screening: run module3 v2 KNN + RWR/TRS for EOS/WBC
# -----------------------------------------------------------------------------
# Thin wrapper around:
#   03_pseudopair_pipeline/08_run_trait_screening_module3.R
#
# Both weights are retained because we need to record atac_w=0.7 and atac_w=1.0.
#
# Run in RStudio:
#   setwd("D:/RStudio/work/module4_traitaware_pseudopair")
#   source("03_pseudopair_pipeline/11_run_trait_screening_remaining_module3.R")
# -----------------------------------------------------------------------------

ROOT <- "D:/RStudio/work/module4_traitaware_pseudopair"

MODULE3_COMPAT_OUT_DIR <- file.path(ROOT, "03_pseudopair_pipeline/output/module3_v2_trait_screening_remaining_ge0p80")
TRAITS_TO_RUN <- c("EOS", "WBC")
ATAC_W_LIST <- c(0.7, 1.0)

cat(">>> Run module3 v2 for remaining trait screening\n\n")
cat("Traits: ", paste(TRAITS_TO_RUN, collapse = ", "), "\n", sep = "")
cat("Output dir: ", MODULE3_COMPAT_OUT_DIR, "\n\n", sep = "")

source(file.path(ROOT, "03_pseudopair_pipeline/08_run_trait_screening_module3.R"))
