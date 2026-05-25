# module4: MONO original pure-ATAC SCAVENGE robustness baseline
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

# Only MONO for now.
RUN_TAG <- "module4_original_scavenge_robust_MONO_originalflow"
OUTPUT_DIR <- file.path(ROOT, "04_evaluation/baseline_original_scavenge_robustness", RUN_TAG)

RAW_ATAC_SE_FILE <- "D:/SCAVENGE-reproducibility-main/blood_data/scATAC-Healthy-Hematopoiesis-191120.rds"
GWAS_FILE <- "D:/SCAVENGE-reproducibility-main/blood_data/CS.trans_all.PP01.strong_atac_overlap.tsv"
META_FILE <- file.path(ROOT, "04_evaluation/inputs/granja_cell_meta.rds")

TRAITS <- c("MONO")
TARGETS <- list(
  MONO = c("11_CD14.Mono.1", "12_CD14.Mono.2")
)

P_DROP_LIST <- round(seq(0.0, 0.9, by = 0.1), 2)
SEED_LIST <- c(1L, 2L, 3L)

# Granja blood data used by blood_small.R is hg19.
GENOME_PACKAGE <- "BSgenome.Hsapiens.UCSC.hg19"
GENOME_OBJECT <- "BSgenome.Hsapiens.UCSC.hg19"

SKIP_IF_EXISTS <- TRUE
APPEND_EXISTING_EVAL <- FALSE

stopifnot(file.exists(COMMON_BASELINE_RUNNER))
source(COMMON_BASELINE_RUNNER, local = environment())
