# module4: run the complete EOS/WBC trait-aware pseudopair flow
# -----------------------------------------------------------------------------
# This is a convenience driver. It only sources thin wrappers; all heavy logic is
# reused from the existing mapping, transfer, module3, and evaluation scripts.
#
# Run in RStudio:
#   setwd("D:/RStudio/work/module4_traitaware_pseudopair")
#   source("RUN_EOS_WBC_FULL_FLOW.R")
# -----------------------------------------------------------------------------

ROOT <- "D:/RStudio/work/module4_traitaware_pseudopair"

source(file.path(ROOT, "01_data_prep_external_rna/11_build_trait_screening_remaining_score_rds.R"))
source(file.path(ROOT, "02_trait_mapping/05_run_trait_screening_remaining_mapping.R"))
source(file.path(ROOT, "02_trait_mapping/06_screen_trait_mapping_remaining_atac_candidates.R"))
source(file.path(ROOT, "03_pseudopair_pipeline/09_transfer_trait_screening_remaining.R"))
source(file.path(ROOT, "03_pseudopair_pipeline/10_prepare_trait_screening_remaining_module3_inputs.R"))
source(file.path(ROOT, "03_pseudopair_pipeline/11_run_trait_screening_remaining_module3.R"))
source(file.path(ROOT, "04_evaluation/run_trait_screening_remaining_auroc.R"))
