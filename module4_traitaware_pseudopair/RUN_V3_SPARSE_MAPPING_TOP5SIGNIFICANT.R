# module4 v3.3: sparse mapping + top5-significant RNA completion
# -----------------------------------------------------------------------------
# Purpose:
#   First move sparse trait features into the ATAC->RNA mapping stage, then
#   reuse the existing top5-significant completion filter on the new top5 RNA
#   references.
#
# Default run is MONO mapping + completion diagnostics only. Downstream KNN/RWR
# can be run later after mapping diagnostics look reasonable.
#
# Run in RStudio:
#   setwd("D:/RStudio/work/module4_traitaware_pseudopair")
#   source("RUN_V3_SPARSE_MAPPING_TOP5SIGNIFICANT.R")
# -----------------------------------------------------------------------------

ROOT <- "D:/RStudio/work/module4_traitaware_pseudopair"

if (!exists("RUN_STEPS")) {
  RUN_STEPS <- c(
    "sparse_mapping_mono",
    "top5sig_mono_transfer"
  )
}

if (!exists("V3_SPARSE_MAP_TAG")) V3_SPARSE_MAP_TAG <- "v3_sparse_mapping"
if (!exists("V3_SPARSE_COMPLETION_TAG")) {
  V3_SPARSE_COMPLETION_TAG <- "v3_sparse_mapping_top5significant"
}
if (!exists("V3_METHOD_LABEL")) V3_METHOD_LABEL <- V3_SPARSE_COMPLETION_TAG

if (!exists("V3_B_NULL_CELL")) V3_B_NULL_CELL <- 1000L
if (!exists("V3_RNA_EMPIRICAL_P_THRESHOLD")) V3_RNA_EMPIRICAL_P_THRESHOLD <- 0.05
if (!exists("V3_ATAC_EMPIRICAL_P_THRESHOLD")) V3_ATAC_EMPIRICAL_P_THRESHOLD <- 0.05
if (!exists("V3_MIN_DELTA_OVER_NULL_MEDIAN")) V3_MIN_DELTA_OVER_NULL_MEDIAN <- 0
if (!exists("V3_SPARSE_RNG_SEED")) V3_SPARSE_RNG_SEED <- 20260519L

if (!exists("V3_SPARSE_COMPLETION_B_NULL_TOP5")) V3_SPARSE_COMPLETION_B_NULL_TOP5 <- 200L
if (!exists("V3_SPARSE_COMPLETION_SIGNIFICANCE_METHOD")) {
  V3_SPARSE_COMPLETION_SIGNIFICANCE_METHOD <- "empirical_p"
}
if (!exists("V3_SPARSE_COMPLETION_EMPIRICAL_P_THRESHOLD")) {
  V3_SPARSE_COMPLETION_EMPIRICAL_P_THRESHOLD <- 0.05
}
if (!exists("V3_SPARSE_COMPLETION_FILTER_Q_THRESHOLD")) {
  V3_SPARSE_COMPLETION_FILTER_Q_THRESHOLD <- 0.10
}
if (!exists("V3_SPARSE_COMPLETION_MAX_RETAINED_GENES_PER_CELL")) {
  V3_SPARSE_COMPLETION_MAX_RETAINED_GENES_PER_CELL <- Inf
}
if (!exists("V3_SPARSE_COMPLETION_FILTER_BATCH_SIZE")) {
  V3_SPARSE_COMPLETION_FILTER_BATCH_SIZE <- 512L
}
if (!exists("V3_TOP5_RNG_SEED")) V3_TOP5_RNG_SEED <- 20260520L

if (!exists("V3_ATAC_W_LIST")) V3_ATAC_W_LIST <- c(0.7, 1.0)

validate_tag <- function(x, name) {
  if (!grepl("^[A-Za-z0-9_.-]+$", x)) {
    stop(name, " may only contain letters, numbers, underscore, dot, or hyphen.",
         call. = FALSE)
  }
}
validate_tag(V3_SPARSE_MAP_TAG, "V3_SPARSE_MAP_TAG")
validate_tag(V3_SPARSE_COMPLETION_TAG, "V3_SPARSE_COMPLETION_TAG")

run_script <- function(script, vars = list()) {
  e <- new.env(parent = globalenv())
  for (nm in names(vars)) assign(nm, vars[[nm]], envir = e)
  cat("\n>>> source(", script, ")\n", sep = "")
  source(file.path(ROOT, script), local = e)
  invisible(e)
}

step_enabled <- function(x) x %in% RUN_STEPS

SPARSE_MAP_DIR <- file.path(ROOT, "02_trait_mapping/output", V3_SPARSE_MAP_TAG)
SPARSE_MAP_DIAG_DIR <- file.path(SPARSE_MAP_DIR, "diagnostics")

MONO_TRANSFER_DIR <- file.path(ROOT, "03_pseudopair_pipeline/output", V3_SPARSE_COMPLETION_TAG, "transfer_mono")
MONO_COMPLETION_DIAG_DIR <- file.path(ROOT, "02_trait_mapping/output", V3_SPARSE_COMPLETION_TAG, "mono_completion_diagnostics")
MONO_COMPAT_DIR <- file.path(ROOT, "03_pseudopair_pipeline/output", V3_SPARSE_COMPLETION_TAG, "module3_mono_top5significant")
MONO_NBD_COMPAT_DIR <- file.path(ROOT, "03_pseudopair_pipeline/output", V3_SPARSE_COMPLETION_TAG, "module3_mono_neighborhood_fdr_ge0p80")
MONO_NBD_DIR <- file.path(ROOT, "04_evaluation/module4_v2_diagnostics/neighborhood_fdr")
MONO_AUROC_DIR <- file.path(ROOT, "04_evaluation", V3_SPARSE_COMPLETION_TAG, "mono_neighborhood_fdr_ge0p80")
MONO_AUROC_CSV <- paste0(V3_SPARSE_COMPLETION_TAG, "_mono_neighborhood_fdr_ge0p80_auroc.csv")

BASE_CITE_EXTRACT_DIR <- file.path(ROOT, "02_trait_mapping/output/cite_extract")

cat(">>> module4 v3.3 sparse mapping + top5-significant completion\n\n")
cat("RUN_STEPS:\n")
cat(paste0("  - ", RUN_STEPS, collapse = "\n"), "\n\n", sep = "")
cat("V3_SPARSE_MAP_TAG = ", V3_SPARSE_MAP_TAG, "\n", sep = "")
cat("V3_SPARSE_COMPLETION_TAG = ", V3_SPARSE_COMPLETION_TAG, "\n", sep = "")
cat("Sparse mapping: B_NULL_CELL=", V3_B_NULL_CELL,
    "; RNA_P=", V3_RNA_EMPIRICAL_P_THRESHOLD,
    "; ATAC_P=", V3_ATAC_EMPIRICAL_P_THRESHOLD, "\n", sep = "")
cat("Completion: B_NULL_TOP5=", V3_SPARSE_COMPLETION_B_NULL_TOP5,
    "; EMPIRICAL_P=", V3_SPARSE_COMPLETION_EMPIRICAL_P_THRESHOLD, "\n\n", sep = "")

if (step_enabled("sparse_mapping_mono")) {
  run_script("02_trait_mapping/08_build_sparse_trait_mapping.R", list(
    SPARSE_MAP_OUT_DIR = SPARSE_MAP_DIR,
    DIAG_DIR = SPARSE_MAP_DIAG_DIR,
    CITE_EXTRACT_DIR = BASE_CITE_EXTRACT_DIR,
    TRAITS_TO_RUN = c("MONO"),
    B_NULL_CELL = V3_B_NULL_CELL,
    RNA_EMPIRICAL_P_THRESHOLD = V3_RNA_EMPIRICAL_P_THRESHOLD,
    ATAC_EMPIRICAL_P_THRESHOLD = V3_ATAC_EMPIRICAL_P_THRESHOLD,
    MIN_DELTA_OVER_NULL_MEDIAN = V3_MIN_DELTA_OVER_NULL_MEDIAN,
    RNG_SEED = V3_SPARSE_RNG_SEED
  ))
}

if (step_enabled("top5sig_mono_transfer")) {
  run_script("02_trait_mapping/07_build_top5_significant_completion.R", list(
    MAP_DIR = SPARSE_MAP_DIR,
    TRANSFER_OUT_DIR = MONO_TRANSFER_DIR,
    DIAG_DIR = MONO_COMPLETION_DIAG_DIR,
    CITE_MTX = file.path(BASE_CITE_EXTRACT_DIR, "cite_trait_counts.mtx"),
    CITE_GENES = file.path(BASE_CITE_EXTRACT_DIR, "cite_trait_counts_genes.txt"),
    CITE_CELLS = file.path(BASE_CITE_EXTRACT_DIR, "cite_trait_counts_cells.txt"),
    CITE_SCORE_RDS = file.path(ROOT, "01_data_prep_external_rna/output/cite_rna_scdrs_scores.rds"),
    TRAITS_TO_RUN = c("MONO"),
    B_NULL_TOP5 = V3_SPARSE_COMPLETION_B_NULL_TOP5,
    SIGNIFICANCE_METHOD = V3_SPARSE_COMPLETION_SIGNIFICANCE_METHOD,
    EMPIRICAL_P_THRESHOLD = V3_SPARSE_COMPLETION_EMPIRICAL_P_THRESHOLD,
    FILTER_Q_THRESHOLD = V3_SPARSE_COMPLETION_FILTER_Q_THRESHOLD,
    MIN_DELTA_OVER_NULL_MEDIAN = V3_MIN_DELTA_OVER_NULL_MEDIAN,
    MAX_RETAINED_GENES_PER_CELL = V3_SPARSE_COMPLETION_MAX_RETAINED_GENES_PER_CELL,
    FILTER_BATCH_SIZE = V3_SPARSE_COMPLETION_FILTER_BATCH_SIZE,
    RNG_SEED = V3_TOP5_RNG_SEED,
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
    PATCHED_SCRIPT = file.path(MONO_AUROC_DIR, paste0("eval_auroc.", V3_SPARSE_COMPLETION_TAG, "_mono_patched.R")),
    METHOD_LABEL_NEW = V3_METHOD_LABEL
  ))
}

cat("\n=== v3.3 sparse mapping driver finished selected steps ===\n")
cat("Key result locations:\n")
cat("  Sparse mapping: ", file.path(SPARSE_MAP_DIR, "granja_cite_mapping_MONO.rds"), "\n", sep = "")
cat("  Sparse mapping diagnostics: ", SPARSE_MAP_DIAG_DIR, "\n", sep = "")
cat("  Completion diagnostics: ", MONO_COMPLETION_DIAG_DIR, "\n", sep = "")
cat("  Completion transfer: ", MONO_TRANSFER_DIR, "\n", sep = "")
cat("  Optional MONO AUROC: ", file.path(MONO_AUROC_DIR, MONO_AUROC_CSV), "\n", sep = "")
