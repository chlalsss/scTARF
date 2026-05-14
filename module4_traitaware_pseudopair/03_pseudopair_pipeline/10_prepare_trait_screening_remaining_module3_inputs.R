# module4 trait screening: prepare module3 inputs for EOS/WBC
# -----------------------------------------------------------------------------
# Thin wrapper around:
#   03_pseudopair_pipeline/07_prepare_trait_screening_module3_inputs.R
#
# It reuses transferred pseudo RNA/zRNA outputs and replaces the default top1%
# seeds with neighborhood-FDR mapping seeds:
#   candidate_group_weight_ge0p80
#
# Run in RStudio:
#   setwd("D:/RStudio/work/module4_traitaware_pseudopair")
#   source("03_pseudopair_pipeline/10_prepare_trait_screening_remaining_module3_inputs.R")
# -----------------------------------------------------------------------------

ROOT <- "D:/RStudio/work/module4_traitaware_pseudopair"

MODULE3_TRANSFER_DIR <- file.path(ROOT, "03_pseudopair_pipeline/output/transfer_trait_screening_remaining")
MODULE3_COMPAT_OUT_DIR <- file.path(ROOT, "03_pseudopair_pipeline/output/module3_v2_trait_screening_remaining_ge0p80")
TRAITS_TO_RUN <- c("EOS", "WBC")

SEED_CELL_FILE <- file.path(
  ROOT,
  "02_trait_mapping/output/trait_screening_remaining/diagnostics/trait_screening_atac_candidate_cell_level.csv.gz"
)
SEED_SUMMARY_FILE <- file.path(MODULE3_COMPAT_OUT_DIR, "seed_summary_neighborhood_fdr_ge0p80.csv")

cat(">>> Prepare module3 inputs for remaining trait screening\n\n")
cat("Traits: ", paste(TRAITS_TO_RUN, collapse = ", "), "\n", sep = "")
cat("Transfer dir: ", MODULE3_TRANSFER_DIR, "\n", sep = "")
cat("Output dir  : ", MODULE3_COMPAT_OUT_DIR, "\n\n", sep = "")

source(file.path(ROOT, "03_pseudopair_pipeline/07_prepare_trait_screening_module3_inputs.R"))
