# module4 v3.2: top5 + random-background significant RNA completion
# -----------------------------------------------------------------------------
# Purpose:
#   Discard the plain top1 strategy and test the teacher-suggested correction:
#   keep top5 RNA references, but retain only genes whose top5 average is
#   significantly higher than random top5 RNA averages.
#
# Run in RStudio:
#   setwd("D:/RStudio/work/module4_traitaware_pseudopair")
#   source("RUN_V3_TOP5_SIGNIFICANT.R")
#
# Optional partial run:
#   RUN_STEPS <- c("top5sig_mono_transfer", "module3_blocks_mono")
#   source("RUN_V3_TOP5_SIGNIFICANT.R")
#
# Optional non-overwriting diagnostic run:
#   V3_OUTPUT_TAG <- "v3_top5_significant_p005"
#   V3_EMPIRICAL_P_THRESHOLD <- 0.05
#   source("RUN_V3_TOP5_SIGNIFICANT.R")
# -----------------------------------------------------------------------------

ROOT <- "D:/RStudio/work/module4_traitaware_pseudopair"

if (!exists("RUN_STEPS")) {
  RUN_STEPS <- c(
    "top5sig_mono_transfer",
    "module3_blocks_mono",
    "knn_mono",
    "neighborhood_seed_mono",
    "rwr_mono",
    "auroc_mono",
    "top5sig_rdw_rbc_transfer",
    "module3_rdw_rbc",
    "run_rdw_rbc",
    "auroc_rdw_rbc"
  )
}

if (!exists("V3_ATAC_W_LIST")) V3_ATAC_W_LIST <- c(0.7, 1.0)
if (!exists("V3_B_NULL_TOP5")) V3_B_NULL_TOP5 <- 1000L
if (!exists("V3_SIGNIFICANCE_METHOD")) V3_SIGNIFICANCE_METHOD <- "empirical_p"
if (!exists("V3_EMPIRICAL_P_THRESHOLD")) V3_EMPIRICAL_P_THRESHOLD <- 0.01
if (!exists("V3_FILTER_Q_THRESHOLD")) V3_FILTER_Q_THRESHOLD <- 0.10
if (!exists("V3_MIN_DELTA_OVER_NULL_MEDIAN")) V3_MIN_DELTA_OVER_NULL_MEDIAN <- 0
if (!exists("V3_MAX_RETAINED_GENES_PER_CELL")) V3_MAX_RETAINED_GENES_PER_CELL <- Inf
if (!exists("V3_FILTER_BATCH_SIZE")) V3_FILTER_BATCH_SIZE <- 512L
if (!exists("V3_RNG_SEED")) V3_RNG_SEED <- 20260518L
if (!exists("V3_OUTPUT_TAG")) V3_OUTPUT_TAG <- "v3_top5_significant"
if (!exists("V3_METHOD_LABEL")) V3_METHOD_LABEL <- V3_OUTPUT_TAG

if (!grepl("^[A-Za-z0-9_.-]+$", V3_OUTPUT_TAG)) {
  stop("V3_OUTPUT_TAG may only contain letters, numbers, underscore, dot, or hyphen.")
}

run_script <- function(script, vars = list()) {
  e <- new.env(parent = globalenv())
  for (nm in names(vars)) assign(nm, vars[[nm]], envir = e)
  cat("\n>>> source(", script, ")\n", sep = "")
  source(file.path(ROOT, script), local = e)
  invisible(e)
}

step_enabled <- function(x) x %in% RUN_STEPS

V3_ROOT <- file.path(ROOT, V3_OUTPUT_TAG)

MONO_MAP_DIR <- file.path(ROOT, "02_trait_mapping/output")
MONO_TRANSFER_DIR <- file.path(ROOT, "03_pseudopair_pipeline/output", V3_OUTPUT_TAG, "transfer_mono")
MONO_DIAG_DIR <- file.path(ROOT, "02_trait_mapping/output", V3_OUTPUT_TAG, "mono_completion_diagnostics")
MONO_COMPAT_DIR <- file.path(ROOT, "03_pseudopair_pipeline/output", V3_OUTPUT_TAG, "module3_mono_top5significant")
MONO_NBD_COMPAT_DIR <- file.path(ROOT, "03_pseudopair_pipeline/output", V3_OUTPUT_TAG, "module3_mono_neighborhood_fdr_ge0p80")
MONO_NBD_DIR <- file.path(ROOT, "04_evaluation/module4_v2_diagnostics/neighborhood_fdr")
MONO_AUROC_DIR <- file.path(ROOT, "04_evaluation", V3_OUTPUT_TAG, "mono_neighborhood_fdr_ge0p80")
MONO_AUROC_CSV <- paste0(V3_OUTPUT_TAG, "_mono_neighborhood_fdr_ge0p80_auroc.csv")

TR_MAP_DIR <- file.path(ROOT, "02_trait_mapping/output/trait_screening")
TR_TRANSFER_DIR <- file.path(ROOT, "03_pseudopair_pipeline/output", V3_OUTPUT_TAG, "transfer_trait_screening")
TR_DIAG_DIR <- file.path(ROOT, "02_trait_mapping/output", V3_OUTPUT_TAG, "trait_screening_completion_diagnostics")
TR_COMPAT_DIR <- file.path(ROOT, "03_pseudopair_pipeline/output", V3_OUTPUT_TAG, "module3_trait_screening_ge0p80")
TR_SEED_CELL_FILE <- file.path(
  ROOT,
  "02_trait_mapping/output/trait_screening/diagnostics/trait_screening_atac_candidate_cell_level.csv.gz"
)
TR_AUROC_DIR <- file.path(ROOT, "04_evaluation", V3_OUTPUT_TAG, "trait_screening_ge0p80")
TR_AUROC_CSV <- paste0(V3_OUTPUT_TAG, "_trait_screening_ge0p80_auroc.csv")

cat(">>> module4 v3.2 top5-significant RNA completion\n\n")
cat("V3_OUTPUT_TAG = ", V3_OUTPUT_TAG, "; METHOD_LABEL = ", V3_METHOD_LABEL, "\n", sep = "")
cat("RUN_STEPS:\n")
cat(paste0("  - ", RUN_STEPS, collapse = "\n"), "\n\n", sep = "")
cat("B_NULL_TOP5 = ", V3_B_NULL_TOP5,
    "; SIGNIFICANCE = ", V3_SIGNIFICANCE_METHOD,
    "; EMPIRICAL_P = ", V3_EMPIRICAL_P_THRESHOLD,
    "; FILTER_Q = ", V3_FILTER_Q_THRESHOLD,
    "; MIN_DELTA = ", V3_MIN_DELTA_OVER_NULL_MEDIAN, "\n\n", sep = "")

if (step_enabled("top5sig_mono_transfer")) {
  run_script("02_trait_mapping/07_build_top5_significant_completion.R", list(
    MAP_DIR = MONO_MAP_DIR,
    TRANSFER_OUT_DIR = MONO_TRANSFER_DIR,
    DIAG_DIR = MONO_DIAG_DIR,
    CITE_MTX = file.path(MONO_MAP_DIR, "cite_extract/cite_trait_counts.mtx"),
    CITE_GENES = file.path(MONO_MAP_DIR, "cite_extract/cite_trait_counts_genes.txt"),
    CITE_CELLS = file.path(MONO_MAP_DIR, "cite_extract/cite_trait_counts_cells.txt"),
    CITE_SCORE_RDS = file.path(ROOT, "01_data_prep_external_rna/output/cite_rna_scdrs_scores.rds"),
    TRAITS_TO_RUN = c("MONO"),
    B_NULL_TOP5 = V3_B_NULL_TOP5,
    SIGNIFICANCE_METHOD = V3_SIGNIFICANCE_METHOD,
    EMPIRICAL_P_THRESHOLD = V3_EMPIRICAL_P_THRESHOLD,
    FILTER_Q_THRESHOLD = V3_FILTER_Q_THRESHOLD,
    MIN_DELTA_OVER_NULL_MEDIAN = V3_MIN_DELTA_OVER_NULL_MEDIAN,
    MAX_RETAINED_GENES_PER_CELL = V3_MAX_RETAINED_GENES_PER_CELL,
    FILTER_BATCH_SIZE = V3_FILTER_BATCH_SIZE,
    RNG_SEED = V3_RNG_SEED,
    NULL_WEIGHT_MODE = "sample_observed"
  ))
}

if (step_enabled("module3_blocks_mono")) {
  run_script("03_pseudopair_pipeline/02_prepare_module3_v2_inputs.R", list(
    MODULE3_TRANSFER_DIR = MONO_TRANSFER_DIR,
    MODULE3_COMPAT_OUT_DIR = MONO_COMPAT_DIR,
    TRAITS_TO_RUN = c("MONO")
  ))
}

if (step_enabled("knn_mono")) {
  run_script("03_pseudopair_pipeline/03_run_module3_v2_knn_rwr.R", list(
    MODULE3_COMPAT_OUT_DIR = MONO_COMPAT_DIR,
    TRAITS_TO_RUN = c("MONO"),
    ATAC_W_LIST = V3_ATAC_W_LIST
  ))
}

if (step_enabled("neighborhood_seed_mono")) {
  run_script("03_pseudopair_pipeline/04_prepare_neighborhood_fdr_seed_inputs.R", list(
    BASE_COMPAT_DIR = MONO_COMPAT_DIR,
    NBD_DIR = MONO_NBD_DIR,
    OUT_DIR = MONO_NBD_COMPAT_DIR,
    NBD_TRAITS_TO_PREPARE = c("MONO"),
    NBD_ATAC_W_LIST = V3_ATAC_W_LIST,
    NBD_GROUP_FDR_THRESHOLD = 0.10,
    NBD_SEED_CANDIDATE_COLUMN = "candidate_group_weight_ge0p80"
  ))
}

if (step_enabled("rwr_mono")) {
  run_script("03_pseudopair_pipeline/05_run_module3_v2_neighborhood_fdr_rwr.R", list(
    OUT_DIR = MONO_NBD_COMPAT_DIR,
    NBD_TRAITS_TO_RUN = c("MONO"),
    NBD_ATAC_W_LIST = V3_ATAC_W_LIST
  ))
}

if (step_enabled("auroc_mono")) {
  run_script("04_evaluation/run_neighborhood_fdr_seed_auroc.R", list(
    V2_FILES_DIR_NEW = MONO_NBD_COMPAT_DIR,
    OUT_DIR_NEW = MONO_AUROC_DIR,
    OUT_CSV_NAME_NEW = MONO_AUROC_CSV,
    PATCHED_SCRIPT = file.path(MONO_AUROC_DIR, paste0("eval_auroc.", V3_OUTPUT_TAG, "_mono_patched.R")),
    METHOD_LABEL_NEW = V3_METHOD_LABEL
  ))
}

if (step_enabled("top5sig_rdw_rbc_transfer")) {
  run_script("02_trait_mapping/07_build_top5_significant_completion.R", list(
    MAP_DIR = TR_MAP_DIR,
    TRANSFER_OUT_DIR = TR_TRANSFER_DIR,
    DIAG_DIR = TR_DIAG_DIR,
    CITE_MTX = file.path(TR_MAP_DIR, "cite_extract/cite_trait_counts.mtx"),
    CITE_GENES = file.path(TR_MAP_DIR, "cite_extract/cite_trait_counts_genes.txt"),
    CITE_CELLS = file.path(TR_MAP_DIR, "cite_extract/cite_trait_counts_cells.txt"),
    CITE_SCORE_RDS = file.path(
      ROOT,
      "01_data_prep_external_rna/output/trait_screening/cite_rna_scdrs_scores_trait_screening.rds"
    ),
    TRAITS_TO_RUN = c("RDW", "RBC"),
    B_NULL_TOP5 = V3_B_NULL_TOP5,
    SIGNIFICANCE_METHOD = V3_SIGNIFICANCE_METHOD,
    EMPIRICAL_P_THRESHOLD = V3_EMPIRICAL_P_THRESHOLD,
    FILTER_Q_THRESHOLD = V3_FILTER_Q_THRESHOLD,
    MIN_DELTA_OVER_NULL_MEDIAN = V3_MIN_DELTA_OVER_NULL_MEDIAN,
    MAX_RETAINED_GENES_PER_CELL = V3_MAX_RETAINED_GENES_PER_CELL,
    FILTER_BATCH_SIZE = V3_FILTER_BATCH_SIZE,
    RNG_SEED = V3_RNG_SEED + 100L,
    NULL_WEIGHT_MODE = "sample_observed"
  ))
}

if (step_enabled("module3_rdw_rbc")) {
  run_script("03_pseudopair_pipeline/07_prepare_trait_screening_module3_inputs.R", list(
    MODULE3_TRANSFER_DIR = TR_TRANSFER_DIR,
    MODULE3_COMPAT_OUT_DIR = TR_COMPAT_DIR,
    TRAITS_TO_RUN = c("RDW", "RBC"),
    SEED_CELL_FILE = TR_SEED_CELL_FILE,
    SEED_SUMMARY_FILE = file.path(TR_COMPAT_DIR, "seed_summary_neighborhood_fdr_ge0p80.csv")
  ))
}

if (step_enabled("run_rdw_rbc")) {
  run_script("03_pseudopair_pipeline/08_run_trait_screening_module3.R", list(
    MODULE3_COMPAT_OUT_DIR = TR_COMPAT_DIR,
    TRAITS_TO_RUN = c("RDW", "RBC"),
    ATAC_W_LIST = V3_ATAC_W_LIST
  ))
}

if (step_enabled("auroc_rdw_rbc")) {
  run_script("04_evaluation/run_trait_screening_auroc.R", list(
    V2_FILES_DIR_NEW = TR_COMPAT_DIR,
    OUT_DIR_NEW = TR_AUROC_DIR,
    OUT_CSV_NAME_NEW = TR_AUROC_CSV,
    PATCHED_SCRIPT = file.path(TR_AUROC_DIR, paste0("eval_auroc.", V3_OUTPUT_TAG, "_trait_screening_patched.R")),
    METHOD_LABEL_NEW = V3_METHOD_LABEL,
    TARGET_REPLACEMENT = c(
      'TARGETS <- list(',
      '  RDW = c("02_Early.Eryth", "03_Late.Eryth", "04_Early.Baso"),',
      '  RBC = c("02_Early.Eryth", "03_Late.Eryth", "04_Early.Baso")',
      ')'
    ),
    PRED_REPLACEMENT = c(
      'PRED_TRAIT <- list(RDW = "RDW", RBC = "RBC")',
      ''
    )
  ))
}

cat("\n=== v3.2 top5-significant driver finished selected steps ===\n")
cat("Key result locations after a full run:\n")
cat("  MONO completion diagnostics: ", MONO_DIAG_DIR, "\n", sep = "")
cat("  RDW/RBC completion diagnostics: ", TR_DIAG_DIR, "\n", sep = "")
cat("  MONO AUROC: ", file.path(MONO_AUROC_DIR, MONO_AUROC_CSV), "\n", sep = "")
cat("  RDW/RBC AUROC: ", file.path(TR_AUROC_DIR, TR_AUROC_CSV), "\n", sep = "")
