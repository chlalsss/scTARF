# module4: run top3 raw ATAC read drop-out robustness
# -----------------------------------------------------------------------------
# This driver runs the full robustness workflow for the selected top3 traits:
#   MONO, RDW, RBC
#
# It reuses module3 robustness code through:
#   05_robustness/02_run_top3_raw_atac_dropout.R
#
# Run in RStudio:
#   setwd("D:/RStudio/work/module4_traitaware_pseudopair")
#   source("RUN_TOP3_RAW_ATAC_DROPOUT_ROBUSTNESS.R")
#
# For a quick sanity run, run the step scripts manually and set:
#   P_DROP_LIST <- c(0.0)
#   SEED_LIST <- c(1)
#
# This driver always runs the planned full grid:
#   p_drop = 0.0, 0.1, ..., 0.9
#   seeds  = 1, 2, 3
# -----------------------------------------------------------------------------

ROOT <- "D:/RStudio/work/module4_traitaware_pseudopair"

P_DROP_LIST <- seq(0.0, 0.9, by = 0.1)
SEED_LIST <- c(1L, 2L, 3L)

source(file.path(ROOT, "05_robustness/00_select_top3_traits.R"))
source(file.path(ROOT, "05_robustness/01_prepare_top3_module3_inputs.R"))
source(file.path(ROOT, "05_robustness/02_run_top3_raw_atac_dropout.R"))
source(file.path(ROOT, "05_robustness/03_summarize_top3_raw_atac_dropout.R"))
