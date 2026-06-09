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

# Module4 traits needed for the original-vs-improved robustness plot.
# Keep this narrow because the original SCAVENGE baseline is expensive; the
# improved method currently has robustness curves for MCH/RDW/RBC/LYM only.
# The previous tag below is kept on disk but should not be reused:
#   module4_original_scavenge_robust_remaining_traits_seed9527
# It was generated before the common runner fixed a data.table name-collision
# that caused every trait to import all selected GWAS rows.
RUN_TAG <- "module4_original_scavenge_robust_remaining_traits_seed9527_traitfilter_fixed"
OUTPUT_DIR <- file.path(ROOT, "04_evaluation/baseline_original_scavenge_robustness", RUN_TAG)

RAW_ATAC_SE_FILE <- "D:/SCAVENGE-reproducibility-main/blood_data/scATAC-Healthy-Hematopoiesis-191120.rds"
GWAS_FILE <- "D:/SCAVENGE-reproducibility-main/blood_data/CS.trans_all.PP01.strong_atac_overlap.tsv"
META_FILE <- file.path(ROOT, "04_evaluation/inputs/granja_cell_meta.rds")

# Only the three traits we need to display: MONO, RDW, LYM.
# RDW and LYM are already done in this output tag and will be skipped via
# SKIP_IF_EXISTS; MONO is new and reuses the per-p_drop LSI/KNN/background-peak
# caches (trait-independent), so it only computes the MONO z-score and TRS.
TRAITS <- c("MONO", "RDW", "LYM")
TARGETS <- list(
  MONO = c("11_CD14.Mono.1", "12_CD14.Mono.2"),
  RDW  = c("02_Early.Eryth", "03_Late.Eryth", "04_Early.Baso"),
  LYM  = c("16_Pre.B", "17_B", "18_Plasma",
           "19_CD8.N", "23_CD8.EM", "24_CD8.CM",
           "20_CD4.N1", "21_CD4.N2", "22_CD4.M",
           "25_NK")
)

# Allow a quick-test override: set P_DROP_LIST / SEED_LIST before sourcing this
# driver and they will be respected; otherwise the full grid runs.
if (!exists("P_DROP_LIST")) P_DROP_LIST <- round(seq(0.0, 0.9, by = 0.1), 2)
if (!exists("SEED_LIST")) SEED_LIST <- c(1L, 2L, 3L)

# Pin input schema and SCAVENGE parameters so stale objects in the RStudio
# global environment cannot change this expensive rerun.
ATAC_ASSAY <- "counts"
GWAS_TRAIT_COL <- "trait"
GWAS_SCORE_COL <- "PP"
GWAS_VARIANT_COL <- "var"
GWAS_CHR_COL <- "chr"
GWAS_START_COL <- "start"
GWAS_END_COL <- "end"
META_CELL_COL <- "barcode"
META_CELLTYPE_COL <- "cell_type"

# Granja blood data used by blood_small.R is hg19.
GENOME_PACKAGE <- "BSgenome.Hsapiens.UCSC.hg19"
GENOME_OBJECT <- "BSgenome.Hsapiens.UCSC.hg19"
CLEAN_RNG_SEED <- 9527L

BG_NITER <- 200L
LSI_DIMS <- 30L
MUTUAL_K <- 30L
GAMMA <- 0.05
SEED_PERCENT_CUT <- 0.05
SCALE_PERCENT_CUT <- 0.01
CAP_Q <- 0.95
DROP_ALL_ZERO_PEAKS <- TRUE

SKIP_IF_EXISTS <- TRUE
APPEND_EXISTING_EVAL <- FALSE

stopifnot(file.exists(COMMON_BASELINE_RUNNER))
source(COMMON_BASELINE_RUNNER, local = environment())
