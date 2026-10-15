source("scripts/_paths.R")
# Configuration for pure-ATAC original SCAVENGE drop-out robustness.
#
# For a new dataset or new traits, edit this file first. The main runner and
# plotting script should not need changes.

if (!exists("ROBUST_SCRIPT_DIR")) {
  ROBUST_SCRIPT_DIR <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

# ---------------------------------------------------------------------------
# User-editable block
# ---------------------------------------------------------------------------

if (!exists("RUN_TAG")) {
  RUN_TAG <- "blood_small_original_scavenge"
}

if (!exists("RAW_ATAC_SE_FILE")) {
  RAW_ATAC_SE_FILE <- file.path(sctarf_blood_root(), "scATAC-Healthy-Hematopoiesis-191120.rds")
}

if (!exists("GWAS_FILE")) {
  GWAS_FILE <- file.path(sctarf_blood_root(), "CS.trans_all.PP01.strong_atac_overlap.tsv")
}

default_meta_file <- file.path(sctarf_nonpaired_root(), "04_evaluation/inputs/granja_cell_meta.rds")
if (!exists("META_FILE")) {
  META_FILE <- if (file.exists(default_meta_file)) default_meta_file else ""
}

if (!exists("OUTPUT_DIR")) {
  OUTPUT_DIR <- file.path(sctarf_output_root(), "scavenge_robustness", RUN_TAG)
}

# Traits must match the GWAS trait column. TARGETS maps each trait to the cell
# types used only for AUROC evaluation.
if (!exists("TRAITS")) {
  TRAITS <- c("MONO", "RDW", "RBC")
}

if (!exists("TARGETS")) {
  TARGETS <- list(
    MONO = c("11_CD14.Mono.1", "12_CD14.Mono.2"),
    RDW  = c("02_Early.Eryth", "03_Late.Eryth", "04_Early.Baso"),
    RBC  = c("02_Early.Eryth", "03_Late.Eryth", "04_Early.Baso")
  )
}

# Drop-out grid. p_drop is the probability that each ATAC count is removed.
if (!exists("P_DROP_LIST")) {
  P_DROP_LIST <- round(seq(0.0, 0.9, by = 0.1), 2)
}
if (!exists("SEED_LIST")) {
  SEED_LIST <- c(1L, 2L, 3L)
}

# Input column names.
if (!exists("ATAC_ASSAY")) ATAC_ASSAY <- "counts"
if (!exists("GWAS_TRAIT_COL")) GWAS_TRAIT_COL <- "trait"
if (!exists("GWAS_SCORE_COL")) GWAS_SCORE_COL <- "PP"
if (!exists("GWAS_VARIANT_COL")) GWAS_VARIANT_COL <- "var"
if (!exists("GWAS_CHR_COL")) GWAS_CHR_COL <- "chr"
if (!exists("GWAS_START_COL")) GWAS_START_COL <- "start"
if (!exists("GWAS_END_COL")) GWAS_END_COL <- "end"
if (!exists("META_CELL_COL")) META_CELL_COL <- "barcode"
if (!exists("META_CELLTYPE_COL")) META_CELLTYPE_COL <- "cell_type"

# Genome used by chromVAR::addGCBias if the SE has no GC/bias column yet.
if (!exists("GENOME_PACKAGE")) GENOME_PACKAGE <- "BSgenome.Hsapiens.UCSC.hg19"
if (!exists("GENOME_OBJECT")) GENOME_OBJECT <- "BSgenome.Hsapiens.UCSC.hg19"

# Original SCAVENGE-style parameters.
if (!exists("BG_NITER")) BG_NITER <- 200L
if (!exists("LSI_DIMS")) LSI_DIMS <- 30L
if (!exists("MUTUAL_K")) MUTUAL_K <- 30L
if (!exists("GAMMA")) GAMMA <- 0.05
if (!exists("SEED_PERCENT_CUT")) SEED_PERCENT_CUT <- 0.05
if (!exists("SCALE_PERCENT_CUT")) SCALE_PERCENT_CUT <- 0.01
if (!exists("CAP_Q")) CAP_Q <- 0.95
if (!exists("CLEAN_RNG_SEED")) CLEAN_RNG_SEED <- 9527L

# Keep this TRUE for high drop rates. All-zero peaks after downsampling are
# removed and background peaks are recomputed for the remaining peak universe.
if (!exists("DROP_ALL_ZERO_PEAKS")) DROP_ALL_ZERO_PEAKS <- TRUE

# Resume behavior.
if (!exists("SKIP_IF_EXISTS")) SKIP_IF_EXISTS <- TRUE
if (!exists("APPEND_EXISTING_EVAL")) APPEND_EXISTING_EVAL <- TRUE

# ---------------------------------------------------------------------------
# Derived paths
# ---------------------------------------------------------------------------

CACHE_DIR <- file.path(OUTPUT_DIR, "cache")
TRS_DIR <- file.path(OUTPUT_DIR, "trs")
BED_DIR <- file.path(OUTPUT_DIR, "trait_beds")
FIG_DIR <- file.path(OUTPUT_DIR, "figures")

EVAL_FILE <- file.path(OUTPUT_DIR, "scavenge_atac_dropout_eval.csv")
SUMMARY_FILE <- file.path(OUTPUT_DIR, "scavenge_atac_dropout_summary.csv")
RUN_MANIFEST_FILE <- file.path(OUTPUT_DIR, "scavenge_atac_dropout_run_manifest.csv")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TRS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(BED_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
