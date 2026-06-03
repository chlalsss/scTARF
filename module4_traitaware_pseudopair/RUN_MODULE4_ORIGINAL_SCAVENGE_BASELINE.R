# module4: original pure-ATAC SCAVENGE robustness baseline
# -----------------------------------------------------------------------------
# This is only a thin module4 driver. The actual SCAVENGE/drop-out implementation
# is reused from:
#   D:/RStudio/work/scavenge_atac_baseline_robustness/01_run_scavenge_atac_dropout.R
#
# Run in RStudio:
#   setwd("D:/RStudio/work/module4_traitaware_pseudopair")
#   source("RUN_MODULE4_ORIGINAL_SCAVENGE_BASELINE.R")
# -----------------------------------------------------------------------------

ROOT <- "D:/RStudio/work/module4_traitaware_pseudopair"
WORK_ROOT <- "D:/RStudio/work"

ROBUST_SCRIPT_DIR <- file.path(WORK_ROOT, "scavenge_atac_baseline_robustness")
COMMON_BASELINE_RUNNER <- file.path(ROBUST_SCRIPT_DIR, "01_run_scavenge_atac_dropout.R")

# Module4 traits that have been run/evaluated in this experiment.
# MONO has already been completed in:
#   04_evaluation/baseline_original_scavenge_robustness/
#   module4_original_scavenge_robust_MONO_originalflow_seed9527
# This driver only runs the remaining traits to avoid recomputing MONO.
# To rerun only a subset, edit TRAITS and keep the matching TARGETS entries.
# The previous tag below is kept on disk but should not be reused:
#   module4_original_scavenge_robust_remaining_traits_seed9527
# It was generated before the common runner fixed a data.table name-collision
# that caused every trait to import all selected GWAS rows.
RUN_TAG <- "module4_original_scavenge_robust_remaining_traits_seed9527_traitfilter_fixed"
OUTPUT_DIR <- file.path(ROOT, "04_evaluation/baseline_original_scavenge_robustness", RUN_TAG)

RAW_ATAC_SE_FILE <- "D:/SCAVENGE-reproducibility-main/blood_data/scATAC-Healthy-Hematopoiesis-191120.rds"
GWAS_FILE <- "D:/SCAVENGE-reproducibility-main/blood_data/CS.trans_all.PP01.strong_atac_overlap.tsv"
META_FILE <- file.path(ROOT, "04_evaluation/inputs/granja_cell_meta.rds")

TRAITS <- c("PLT", "LYM", "MCH", "RDW", "RBC", "EOS", "WBC")
TARGETS <- list(
  PLT  = c("02_Early.Eryth", "03_Late.Eryth", "04_Early.Baso"),
  LYM  = c("16_Pre.B", "17_B", "18_Plasma",
           "19_CD8.N", "23_CD8.EM", "24_CD8.CM",
           "20_CD4.N1", "21_CD4.N2", "22_CD4.M",
           "25_NK"),
  MCH  = c("02_Early.Eryth", "03_Late.Eryth", "04_Early.Baso"),
  RDW  = c("02_Early.Eryth", "03_Late.Eryth", "04_Early.Baso"),
  RBC  = c("02_Early.Eryth", "03_Late.Eryth", "04_Early.Baso"),
  EOS  = c("04_Early.Baso", "07_GMP", "08_GMP.Neut"),
  WBC  = c("04_Early.Baso",
           "06_CLP.1", "07_GMP", "08_GMP.Neut", "09_pDC", "10_cDC",
           "11_CD14.Mono.1", "12_CD14.Mono.2",
           "15_CLP.2", "16_Pre.B", "17_B", "18_Plasma",
           "19_CD8.N", "20_CD4.N1", "21_CD4.N2", "22_CD4.M",
           "23_CD8.EM", "24_CD8.CM", "25_NK")
)

P_DROP_LIST <- round(seq(0.0, 0.9, by = 0.1), 2)
SEED_LIST <- c(1L, 2L, 3L)

# Granja blood data used by blood_small.R is hg19.
GENOME_PACKAGE <- "BSgenome.Hsapiens.UCSC.hg19"
GENOME_OBJECT <- "BSgenome.Hsapiens.UCSC.hg19"
CLEAN_RNG_SEED <- 9527L

SKIP_IF_EXISTS <- TRUE
APPEND_EXISTING_EVAL <- FALSE

stopifnot(file.exists(COMMON_BASELINE_RUNNER))
source(COMMON_BASELINE_RUNNER, local = environment())
