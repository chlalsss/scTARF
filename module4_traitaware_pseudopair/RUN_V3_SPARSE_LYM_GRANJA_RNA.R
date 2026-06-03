# module4 v3 sparse LYM run using Granja RNA reference
# -----------------------------------------------------------------------------
# This is a thin driver that reuses the fixed v3 sparse workflow:
#   08 sparse mapping
#   07 top5-significant pseudo RNA completion
#   diagnose_scdrs_neighborhood_fdr for RNA groups
#   09 neighborhood-FDR seed candidates
#   07 module3 input preparation with ge0p80 seed replacement
#   03 KNN/RWR
#   run_trait_screening_auroc
#
# It only changes the RNA reference for LYM from external CITE RNA to Granja RNA.
# Run in RStudio:
#   setwd("D:/RStudio/work/module4_traitaware_pseudopair")
#   source("RUN_V3_SPARSE_LYM_GRANJA_RNA.R")
# -----------------------------------------------------------------------------

ROOT <- "D:/RStudio/work/module4_traitaware_pseudopair"

if (!exists("RUN_STEPS")) {
  RUN_STEPS <- c(
    "prepare_granja_scores",
    "sparse_mapping",
    "top5sig_transfer",
    "rna_neighborhood_fdr",
    "seed_candidates",
    "module3_blocks_and_seed",
    "knn_rwr",
    "auroc"
  )
}

if (!exists("V3_LYM_GRANJA_TAG")) V3_LYM_GRANJA_TAG <- "v3_sparse_lym_granja_rna"
if (!exists("V3_METHOD_LABEL")) V3_METHOD_LABEL <- V3_LYM_GRANJA_TAG

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
if (!exists("V3_GROUP_FDR_THRESHOLD")) V3_GROUP_FDR_THRESHOLD <- 0.10
if (!exists("V3_RNA_KNN_K")) V3_RNA_KNN_K <- 30L
if (!exists("V3_LEIDEN_RESOLUTION")) V3_LEIDEN_RESOLUTION <- 0.005

TRAITS_TO_RUN <- c("LYM")
TRAIT_IDS <- c(LYM = "UKB_460K.blood_LYMPHOCYTE_COUNT")
LYM_TARGETS <- c(
  "16_Pre.B", "17_B", "18_Plasma",
  "19_CD8.N", "20_CD4.N1", "21_CD4.N2", "22_CD4.M",
  "23_CD8.EM", "24_CD8.CM", "25_NK"
)
GRANJA_RNA_MYELOID <- c(
  "07_GMP", "08_GMP.Neut", "09_pDC", "10_cDC",
  "11_CD14.Mono.1", "12_CD14.Mono.2", "13_CD16.Mono"
)

BASE_OUT_DIR <- file.path(ROOT, "02_trait_mapping/output", V3_LYM_GRANJA_TAG)
INPUT_DIR <- file.path(BASE_OUT_DIR, "_inputs")
SPARSE_MAP_DIR <- BASE_OUT_DIR
SPARSE_MAP_DIAG_DIR <- file.path(SPARSE_MAP_DIR, "diagnostics")
RNA_EXTRACT_DIR <- file.path(SPARSE_MAP_DIR, "granja_rna_extract")
TRANSFER_DIR <- file.path(ROOT, "03_pseudopair_pipeline/output", V3_LYM_GRANJA_TAG, "transfer")
COMPLETION_DIAG_DIR <- file.path(BASE_OUT_DIR, "top5significant_diagnostics")
NBD_DIR <- file.path(BASE_OUT_DIR, "rna_neighborhood_fdr")
SEED_DIAG_DIR <- file.path(BASE_OUT_DIR, "neighborhood_seed_diagnostics")
COMPAT_DIR <- file.path(ROOT, "03_pseudopair_pipeline/output", V3_LYM_GRANJA_TAG, "module3_neighborhood_fdr_ge0p80")
AUROC_DIR <- file.path(ROOT, "04_evaluation", V3_LYM_GRANJA_TAG, "neighborhood_fdr_ge0p80")

GS_FILE <- file.path(ROOT, "01_data_prep/03_rna_scdrs_python/granja_3traits.gs")
GRANJA_RNA_H5AD <- file.path(ROOT, "01_data_prep/03_rna_scdrs_python/input/granja_rna_raw_counts.h5ad")
GRANJA_RNA_META_RDS <- file.path(ROOT, "01_data_prep/output/granja_rna_cell_meta.rds")
GRANJA_RNA_SCORE_RDS_RAW <- file.path(ROOT, "01_data_prep/output/granja_rna_scdrs_scores.rds")
GRANJA_FULL_SCORE_DIR <- file.path(ROOT, "01_data_prep/03_rna_scdrs_python/output")
GRANJA_SCORE_RDS <- file.path(INPUT_DIR, "granja_rna_scores_LYM_for_v3_sparse.rds")
GRANJA_SCORE_CSV <- file.path(INPUT_DIR, "granja_rna_scores_LYM_for_v3_sparse.csv.gz")
SEED_CELL_FILE <- file.path(SEED_DIAG_DIR, "sparse_neighborhood_seed_cell_level.csv.gz")

dir.create(INPUT_DIR, recursive = TRUE, showWarnings = FALSE)

step_enabled <- function(x) x %in% RUN_STEPS

stop_missing <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0L) {
    stop("Missing required file(s):\n  ", paste(missing, collapse = "\n  "),
         call. = FALSE)
  }
}

run_script <- function(script, vars = list()) {
  e <- new.env(parent = globalenv())
  for (nm in names(vars)) assign(nm, vars[[nm]], envir = e)
  cat("\n>>> source(", script, ")\n", sep = "")
  source(file.path(ROOT, script), local = e)
  invisible(e)
}

r_literal_chr <- function(x) {
  paste(sprintf('"%s"', x), collapse = ", ")
}

make_lym_target_replacement <- function() {
  c(
    "TARGETS <- list(",
    paste0("  LYM = c(", r_literal_chr(LYM_TARGETS), ")"),
    ")"
  )
}

make_lym_pred_replacement <- function() {
  c('PRED_TRAIT <- list(LYM = "LYM")', "")
}

write_granja_score_rds <- function(out_rds, out_csv) {
  suppressPackageStartupMessages({ library(data.table) })
  stop_missing(c(GRANJA_RNA_META_RDS, GRANJA_RNA_SCORE_RDS_RAW))

  meta <- as.data.table(readRDS(GRANJA_RNA_META_RDS))
  score <- as.data.table(readRDS(GRANJA_RNA_SCORE_RDS_RAW))

  if (!"barcode" %in% names(meta) || !"barcode" %in% names(score)) {
    stop("Granja RNA metadata and scores must contain barcode.", call. = FALSE)
  }
  if (!"LYM_norm_score" %in% names(score)) {
    stop("Granja RNA score RDS is missing LYM_norm_score.", call. = FALSE)
  }
  if (!"cell_type" %in% names(meta)) {
    stop("Granja RNA metadata is missing cell_type.", call. = FALSE)
  }

  lib_col <- intersect(c("nUMI", "nUMI_pre", "scdrs_input_total_counts"), names(meta))[1]
  if (is.na(lib_col)) {
    stop("Granja RNA metadata has no nUMI/nUMI_pre/scdrs_input_total_counts column.",
         call. = FALSE)
  }

  out <- merge(
    meta[, .(
      cell = as.character(barcode),
      cell_type = as.character(cell_type),
      scdrs_input_total_counts = as.numeric(get(lib_col))
    )],
    score[, .(
      cell = as.character(barcode),
      LYM_norm_score = as.numeric(LYM_norm_score),
      LYM_pval = if ("LYM_pval" %in% names(score)) as.numeric(LYM_pval) else NA_real_
    )],
    by = "cell",
    all = FALSE,
    sort = FALSE
  )
  if (nrow(out) == 0L) stop("No overlap between Granja RNA metadata and scores.", call. = FALSE)
  bad_lib <- !is.finite(out$scdrs_input_total_counts) | out$scdrs_input_total_counts < 1
  if (any(bad_lib)) {
    warning("Replacing ", sum(bad_lib), " invalid Granja RNA library sizes with 1.")
    out$scdrs_input_total_counts[bad_lib] <- 1
  }

  saveRDS(as.data.frame(out), out_rds)
  fwrite(out, out_csv)
  out_rds
}

ensure_granja_score_rds <- function() {
  if (!file.exists(GRANJA_SCORE_RDS)) {
    cat("[prepare] Granja RNA score RDS is missing; writing it now ...\n")
    write_granja_score_rds(GRANJA_SCORE_RDS, GRANJA_SCORE_CSV)
    cat("  wrote: ", GRANJA_SCORE_RDS, "\n", sep = "")
  }
  invisible(GRANJA_SCORE_RDS)
}

ensure_rna_neighborhood_fdr <- function() {
  group_file <- file.path(NBD_DIR, "rna_unsupervised_group_membership.csv.gz")
  test_file <- file.path(NBD_DIR, "neighborhood_fdr_group_tests.csv")
  if (!file.exists(group_file) || !file.exists(test_file)) {
    cat("[rna_neighborhood_fdr] Required Granja RNA neighborhood-FDR files are missing; generating them now ...\n")
    ensure_granja_score_rds()
    run_script("04_evaluation/diagnose_scdrs_neighborhood_fdr.R", list(
      MAP_DIR = SPARSE_MAP_DIR,
      CITE_SCORE_RDS = GRANJA_SCORE_RDS,
      CITE_MTX = file.path(RNA_EXTRACT_DIR, "cite_trait_counts.mtx"),
      CITE_GENES = file.path(RNA_EXTRACT_DIR, "cite_trait_counts_genes.txt"),
      CITE_CELLS = file.path(RNA_EXTRACT_DIR, "cite_trait_counts_cells.txt"),
      FULL_SCORE_DIR = GRANJA_FULL_SCORE_DIR,
      OUT_DIR = NBD_DIR,
      TRAITS = TRAITS_TO_RUN,
      TRAIT_IDS = TRAIT_IDS,
      RNA_REFERENCE_LABEL = "Granja RNA",
      RNA_KNN_K = V3_RNA_KNN_K,
      LEIDEN_RESOLUTION = V3_LEIDEN_RESOLUTION,
      CITE_TARGETS = list(LYM = LYM_TARGETS),
      CITE_MYELOID = GRANJA_RNA_MYELOID,
      GROUP_FDR_THRESHOLDS = c(0.05, 0.10, 0.20)
    ))
  }
  invisible(list(group_file = group_file, test_file = test_file))
}

cat(">>> v3 sparse LYM with Granja RNA reference\n\n")
cat("RUN_STEPS:\n")
cat(paste0("  - ", RUN_STEPS, collapse = "\n"), "\n\n", sep = "")
cat("Tag: ", V3_LYM_GRANJA_TAG, "\n", sep = "")
cat("RNA reference: Granja scRNA Healthy Hematopoiesis\n")
cat("Trait: LYM / ", TRAIT_IDS[["LYM"]], "\n\n", sep = "")
cat("RNA neighborhood: k=", V3_RNA_KNN_K,
    ", Leiden resolution=", V3_LEIDEN_RESOLUTION, "\n\n", sep = "")

if (step_enabled("prepare_granja_scores")) {
  cat("[prepare] Writing Granja RNA score RDS for LYM ...\n")
  write_granja_score_rds(GRANJA_SCORE_RDS, GRANJA_SCORE_CSV)
  cat("  wrote: ", GRANJA_SCORE_RDS, "\n", sep = "")
}

if (step_enabled("sparse_mapping")) {
  ensure_granja_score_rds()
  run_script("02_trait_mapping/08_build_sparse_trait_mapping.R", list(
    SPARSE_MAP_OUT_DIR = SPARSE_MAP_DIR,
    DIAG_DIR = SPARSE_MAP_DIAG_DIR,
    CITE_EXTRACT_DIR = RNA_EXTRACT_DIR,
    GS_FILE = GS_FILE,
    CITE_H5AD = GRANJA_RNA_H5AD,
    CITE_SCORE_RDS = GRANJA_SCORE_RDS,
    TRAITS_ALL = TRAIT_IDS,
    TRAITS_TO_RUN = TRAITS_TO_RUN,
    B_NULL_CELL = V3_B_NULL_CELL,
    RNA_EMPIRICAL_P_THRESHOLD = V3_RNA_EMPIRICAL_P_THRESHOLD,
    ATAC_EMPIRICAL_P_THRESHOLD = V3_ATAC_EMPIRICAL_P_THRESHOLD,
    MIN_DELTA_OVER_NULL_MEDIAN = V3_MIN_DELTA_OVER_NULL_MEDIAN,
    RNG_SEED = V3_SPARSE_RNG_SEED
  ))
}

if (step_enabled("top5sig_transfer")) {
  ensure_granja_score_rds()
  run_script("02_trait_mapping/07_build_top5_significant_completion.R", list(
    MAP_DIR = SPARSE_MAP_DIR,
    TRANSFER_OUT_DIR = TRANSFER_DIR,
    DIAG_DIR = COMPLETION_DIAG_DIR,
    CITE_MTX = file.path(RNA_EXTRACT_DIR, "cite_trait_counts.mtx"),
    CITE_GENES = file.path(RNA_EXTRACT_DIR, "cite_trait_counts_genes.txt"),
    CITE_CELLS = file.path(RNA_EXTRACT_DIR, "cite_trait_counts_cells.txt"),
    CITE_SCORE_RDS = GRANJA_SCORE_RDS,
    TRAITS_TO_RUN = TRAITS_TO_RUN,
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

if (step_enabled("rna_neighborhood_fdr")) {
  ensure_granja_score_rds()
  run_script("04_evaluation/diagnose_scdrs_neighborhood_fdr.R", list(
    MAP_DIR = SPARSE_MAP_DIR,
    CITE_SCORE_RDS = GRANJA_SCORE_RDS,
    CITE_MTX = file.path(RNA_EXTRACT_DIR, "cite_trait_counts.mtx"),
    CITE_GENES = file.path(RNA_EXTRACT_DIR, "cite_trait_counts_genes.txt"),
    CITE_CELLS = file.path(RNA_EXTRACT_DIR, "cite_trait_counts_cells.txt"),
    FULL_SCORE_DIR = GRANJA_FULL_SCORE_DIR,
    OUT_DIR = NBD_DIR,
    TRAITS = TRAITS_TO_RUN,
    TRAIT_IDS = TRAIT_IDS,
    RNA_REFERENCE_LABEL = "Granja RNA",
    RNA_KNN_K = V3_RNA_KNN_K,
    LEIDEN_RESOLUTION = V3_LEIDEN_RESOLUTION,
    CITE_TARGETS = list(LYM = LYM_TARGETS),
    CITE_MYELOID = GRANJA_RNA_MYELOID,
    GROUP_FDR_THRESHOLDS = c(0.05, 0.10, 0.20)
  ))
}

if (step_enabled("seed_candidates")) {
  ensure_rna_neighborhood_fdr()
  run_script("02_trait_mapping/09_screen_sparse_mapping_neighborhood_seeds.R", list(
    MAP_DIR = SPARSE_MAP_DIR,
    OUT_DIR = SEED_DIAG_DIR,
    TRAITS_TO_SCREEN = TRAITS_TO_RUN,
    RNA_GROUP_FILE = file.path(NBD_DIR, "rna_unsupervised_group_membership.csv.gz"),
    RNA_GROUP_TEST_FILES = file.path(NBD_DIR, "neighborhood_fdr_group_tests.csv"),
    GROUP_TEST_TRAIT_MAP = c(),
    GROUP_FDR_THRESHOLD = V3_GROUP_FDR_THRESHOLD,
    POSTHOC_TARGETS_LIST = list(LYM = LYM_TARGETS),
    POSTHOC_TARGET_NAMES = c(LYM = "lymph_proxy")
  ))
}

if (step_enabled("module3_blocks_and_seed")) {
  run_script("03_pseudopair_pipeline/07_prepare_trait_screening_module3_inputs.R", list(
    MODULE3_TRANSFER_DIR = TRANSFER_DIR,
    MODULE3_COMPAT_OUT_DIR = COMPAT_DIR,
    TRAITS_TO_RUN = TRAITS_TO_RUN,
    SEED_CELL_FILE = SEED_CELL_FILE,
    SEED_SUMMARY_FILE = file.path(COMPAT_DIR, "seed_summary_neighborhood_fdr_ge0p80.csv")
  ))
}

if (step_enabled("knn_rwr")) {
  run_script("03_pseudopair_pipeline/03_run_module3_v2_knn_rwr.R", list(
    MODULE3_COMPAT_OUT_DIR = COMPAT_DIR,
    TRAITS_TO_RUN = TRAITS_TO_RUN,
    ATAC_W_LIST = V3_ATAC_W_LIST
  ))
}

if (step_enabled("auroc")) {
  run_script("04_evaluation/run_trait_screening_auroc.R", list(
    V2_FILES_DIR_NEW = COMPAT_DIR,
    OUT_DIR_NEW = AUROC_DIR,
    OUT_CSV_NAME_NEW = paste0(V3_LYM_GRANJA_TAG, "_neighborhood_fdr_ge0p80_auroc.csv"),
    PATCHED_SCRIPT = file.path(AUROC_DIR, paste0("eval_auroc.", V3_LYM_GRANJA_TAG, "_patched.R")),
    METHOD_LABEL_NEW = V3_METHOD_LABEL,
    TARGET_REPLACEMENT = make_lym_target_replacement(),
    PRED_REPLACEMENT = make_lym_pred_replacement()
  ))
}

cat("\n=== v3 sparse LYM Granja RNA driver finished selected steps ===\n")
cat("Key result locations:\n")
cat("  Granja score RDS: ", GRANJA_SCORE_RDS, "\n", sep = "")
cat("  Sparse mapping: ", SPARSE_MAP_DIR, "\n", sep = "")
cat("  Transfer: ", TRANSFER_DIR, "\n", sep = "")
cat("  RNA neighborhood-FDR: ", NBD_DIR, "\n", sep = "")
cat("  Seed diagnostics: ", SEED_DIAG_DIR, "\n", sep = "")
cat("  Module3 KNN/RWR/TRS: ", COMPAT_DIR, "\n", sep = "")
cat("  AUROC: ", file.path(AUROC_DIR, paste0(V3_LYM_GRANJA_TAG, "_neighborhood_fdr_ge0p80_auroc.csv")), "\n", sep = "")
