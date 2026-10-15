# module4 v3 sparse pipeline for non-MONO traits
# -----------------------------------------------------------------------------
# Purpose:
#   Run the trait-aware v3 sparse workflow for the final non-MONO non-paired
#   trait reported in the main table.
#
# This script intentionally does not write into the MONO v3 result folders.
#
# Run from the repository root:
#   source("scripts/03_nonpaired_analysis/05_run_v3_sparse_other_traits.R")
#
# Optional partial run:
#   OTHER_TRAITS_TO_RUN <- c("RDW", "RBC")
#   RUN_STEPS <- c("prepare_unified_inputs", "sparse_mapping", "top5sig_transfer")
#   source("scripts/03_nonpaired_analysis/05_run_v3_sparse_other_traits.R")
# -----------------------------------------------------------------------------

current_script_dir <- function(default) {
  path <- tryCatch(normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = FALSE),
                   error = function(e) "")
  if (nzchar(path)) dirname(path) else normalizePath(default, winslash = "/", mustWork = FALSE)
}

CODE_ROOT <- current_script_dir("scripts/03_nonpaired_analysis")
PROJECT_ROOT <- normalizePath(file.path(CODE_ROOT, "..", ".."), winslash = "/", mustWork = FALSE)
ROOT <- Sys.getenv("SCTARF_NONPAIRED_ROOT", file.path(PROJECT_ROOT, "outputs", "nonpaired"))

if (!exists("OTHER_TRAITS_TO_RUN")) {
  OTHER_TRAITS_TO_RUN <- c("RDW")
}

if (!exists("RUN_STEPS")) {
  RUN_STEPS <- c(
    "prepare_unified_inputs",
    "sparse_mapping",
    "top5sig_transfer",
    "seed_candidates",
    "module3_blocks_and_seed",
    "knn_rwr",
    "auroc"
  )
}

if (!exists("V3_SPARSE_OTHER_TAG")) V3_SPARSE_OTHER_TAG <- "v3_sparse_other_traits"
if (!exists("V3_METHOD_LABEL")) V3_METHOD_LABEL <- V3_SPARSE_OTHER_TAG

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

TRAIT_IDS_ALL <- c(
  PLT = "UKB_460K.blood_PLATELET_COUNT",
  LYM = "UKB_460K.blood_LYMPHOCYTE_COUNT",
  MCH = "UKB_460K.blood_MEAN_CORPUSCULAR_HEMOGLOBIN",
  RDW = "UKB_460K.blood_RBC_DISTRIB_WIDTH",
  RBC = "UKB_460K.blood_RED_COUNT",
  EOS = "UKB_460K.blood_EOSINOPHIL_COUNT",
  WBC = "UKB_460K.blood_WHITE_COUNT"
)

red_proxy <- c("02_Early.Eryth", "03_Late.Eryth", "04_Early.Baso")
lymph_proxy <- c(
  "16_Pre.B", "17_B", "18_Plasma",
  "19_CD8.N", "20_CD4.N1", "21_CD4.N2", "22_CD4.M",
  "23_CD8.EM", "24_CD8.CM", "25_NK"
)
granulocyte_proxy <- c("04_Early.Baso", "07_GMP", "08_GMP.Neut")
leukocyte_proxy <- c(
  "04_Early.Baso",
  "06_CLP.1", "07_GMP", "08_GMP.Neut", "09_pDC", "10_cDC",
  "11_CD14.Mono.1", "12_CD14.Mono.2",
  "15_CLP.2", "16_Pre.B", "17_B", "18_Plasma",
  "19_CD8.N", "20_CD4.N1", "21_CD4.N2", "22_CD4.M",
  "23_CD8.EM", "24_CD8.CM", "25_NK"
)
POSTHOC_TARGETS_LIST <- list(
  PLT = red_proxy,
  LYM = lymph_proxy,
  MCH = red_proxy,
  RDW = red_proxy,
  RBC = red_proxy,
  EOS = granulocyte_proxy,
  WBC = leukocyte_proxy
)
POSTHOC_TARGET_NAMES <- c(
  PLT = "red_lineage_proxy",
  LYM = "lymph_proxy",
  MCH = "red_lineage_proxy",
  RDW = "red_lineage_proxy",
  RBC = "red_lineage_proxy",
  EOS = "granulocyte_like_proxy",
  WBC = "leukocyte_non_erythroid_proxy"
)

validate_tag <- function(x, name) {
  if (!grepl("^[A-Za-z0-9_.-]+$", x)) {
    stop(name, " may only contain letters, numbers, underscore, dot, or hyphen.",
         call. = FALSE)
  }
}

stop_missing <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0L) {
    stop("Missing required file(s):\n  ", paste(missing, collapse = "\n  "),
         call. = FALSE)
  }
}

step_enabled <- function(x) x %in% RUN_STEPS

run_script <- function(script, vars = list()) {
  script_alias <- c(
    "02_trait_mapping/08_build_sparse_trait_mapping.R" = "39_build_sparse_trait_mapping.R",
    "02_trait_mapping/07_build_top5_significant_completion.R" = "38_build_top5_significant_completion.R",
    "02_trait_mapping/09_screen_sparse_mapping_neighborhood_seeds.R" = "40_screen_sparse_mapping_neighborhood_seeds.R",
    "03_pseudopair_pipeline/07_prepare_trait_screening_module3_inputs.R" = "49_prepare_trait_screening_module3_inputs.R",
    "03_pseudopair_pipeline/03_run_module3_v2_knn_rwr.R" = "45_run_module3_v2_knn_rwr.R",
    "04_evaluation/run_trait_screening_auroc.R" = "58_run_trait_screening_auroc.R"
  )
  script_path <- file.path(CODE_ROOT, if (script %in% names(script_alias)) script_alias[[script]] else script)
  e <- new.env(parent = globalenv())
  assign("ROOT", ROOT, envir = e)
  assign("SCTARF_SCRIPT_DIR", CODE_ROOT, envir = e)
  for (nm in names(vars)) assign(nm, vars[[nm]], envir = e)
  cat("\n>>> source(", script_path, ")\n", sep = "")
  source(script_path, local = e)
  invisible(e)
}

r_literal_chr <- function(x) {
  paste(sprintf('"%s"', x), collapse = ", ")
}

make_target_replacement <- function(traits, target_list) {
  lines <- character(length(traits) + 2L)
  lines[1] <- "TARGETS <- list("
  for (i in seq_along(traits)) {
    trait <- traits[i]
    suffix <- if (i < length(traits)) "," else ""
    lines[i + 1L] <- paste0("  ", trait, " = c(",
                            r_literal_chr(target_list[[trait]]), ")", suffix)
  }
  lines[length(lines)] <- ")"
  lines
}

make_pred_replacement <- function(traits) {
  pairs <- paste0(traits, ' = "', traits, '"')
  c(paste0("PRED_TRAIT <- list(", paste(pairs, collapse = ", "), ")"), "")
}

read_gs_table <- function(path) {
  stop_missing(path)
  dt <- data.table::fread(path, sep = "\t", header = TRUE)
  if (!all(c("TRAIT", "GENESET") %in% names(dt))) {
    stop("Gene-set file must contain TRAIT and GENESET columns: ", path,
         call. = FALSE)
  }
  dt[, .(TRAIT, GENESET)]
}

write_unified_gene_set <- function(out_file, trait_ids) {
  gs_files <- c(
    file.path(ROOT, "01_data_prep/03_rna_scdrs_python/granja_3traits.gs"),
    file.path(ROOT, "01_data_prep_external_rna/03_cite_scdrs_python/cite_trait_screening.gs")
  )
  gs <- data.table::rbindlist(lapply(gs_files, read_gs_table), use.names = TRUE, fill = TRUE)
  keep <- gs[TRAIT %in% unname(trait_ids)]
  missing_ids <- setdiff(unname(trait_ids), keep$TRAIT)
  if (length(missing_ids) > 0L) {
    stop("Missing trait IDs in gene-set files: ",
         paste(missing_ids, collapse = ", "), call. = FALSE)
  }
  keep <- keep[match(unname(trait_ids), TRAIT)]
  data.table::fwrite(keep, out_file, sep = "\t")
  out_file
}

merge_score_rds <- function(out_rds, out_csv, traits) {
  base_rds <- file.path(ROOT, "01_data_prep_external_rna/output/cite_rna_scdrs_scores.rds")
  screening_rds <- file.path(
    ROOT,
    "01_data_prep_external_rna/output/trait_screening/cite_rna_scdrs_scores_trait_screening.rds"
  )
  remaining_rds <- file.path(
    ROOT,
    "01_data_prep_external_rna/output/trait_screening/cite_rna_scdrs_scores_trait_screening_remaining.rds"
  )
  stop_missing(base_rds)
  out <- data.table::as.data.table(readRDS(base_rds))
  if (!all(c("cell", "scdrs_input_total_counts") %in% names(out))) {
    stop("Base CITE score RDS must contain cell and scdrs_input_total_counts.",
         call. = FALSE)
  }

  add_scores <- function(path) {
    if (!file.exists(path)) return(invisible(NULL))
    dt <- data.table::as.data.table(readRDS(path))
    if (!"cell" %in% names(dt)) stop("Score RDS missing cell column: ", path, call. = FALSE)
    score_cols <- grep("_norm_score$", names(dt), value = TRUE)
    score_cols <- setdiff(score_cols, names(out))
    if (length(score_cols) == 0L) return(invisible(NULL))
    out <<- merge(
      out,
      dt[, c("cell", score_cols), with = FALSE],
      by = "cell",
      all.x = TRUE,
      sort = FALSE
    )
    invisible(NULL)
  }

  add_scores(screening_rds)
  add_scores(remaining_rds)

  needed <- paste0(traits, "_norm_score")
  missing_cols <- setdiff(needed, names(out))
  if (length(missing_cols) > 0L) {
    stop("Unified CITE score RDS is missing required score columns: ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  for (col in needed) {
    bad <- !is.finite(out[[col]]) | is.na(out[[col]])
    if (any(bad)) {
      warning("Replacing ", sum(bad), " missing values in ", col, " with 0.")
      data.table::set(out, i = which(bad), j = col, value = 0)
    }
  }

  saveRDS(as.data.frame(out), out_rds)
  data.table::fwrite(out, out_csv)
  out_rds
}

validate_tag(V3_SPARSE_OTHER_TAG, "V3_SPARSE_OTHER_TAG")

if ("MONO" %in% OTHER_TRAITS_TO_RUN) {
  stop("This driver is for non-MONO traits only. Remove MONO from OTHER_TRAITS_TO_RUN.",
       call. = FALSE)
}
unknown_traits <- setdiff(OTHER_TRAITS_TO_RUN, names(TRAIT_IDS_ALL))
if (length(unknown_traits) > 0L) {
  stop("Unknown trait(s): ", paste(unknown_traits, collapse = ", "), call. = FALSE)
}

TRAITS_TO_RUN <- OTHER_TRAITS_TO_RUN
TRAIT_IDS <- TRAIT_IDS_ALL[TRAITS_TO_RUN]

BASE_OUT_DIR <- file.path(ROOT, "02_trait_mapping/output", V3_SPARSE_OTHER_TAG)
UNIFIED_INPUT_DIR <- file.path(BASE_OUT_DIR, "_inputs")
SPARSE_MAP_DIR <- BASE_OUT_DIR
SPARSE_MAP_DIAG_DIR <- file.path(SPARSE_MAP_DIR, "diagnostics")
CITE_EXTRACT_DIR <- file.path(SPARSE_MAP_DIR, "cite_extract")
TRANSFER_DIR <- file.path(ROOT, "03_pseudopair_pipeline/output", V3_SPARSE_OTHER_TAG, "transfer")
COMPLETION_DIAG_DIR <- file.path(BASE_OUT_DIR, "top5significant_diagnostics")
SEED_DIAG_DIR <- file.path(BASE_OUT_DIR, "neighborhood_seed_diagnostics")
COMPAT_DIR <- file.path(ROOT, "03_pseudopair_pipeline/output", V3_SPARSE_OTHER_TAG, "module3_neighborhood_fdr_ge0p80")
AUROC_DIR <- file.path(ROOT, "04_evaluation", V3_SPARSE_OTHER_TAG, "neighborhood_fdr_ge0p80")

UNIFIED_GS_FILE <- file.path(UNIFIED_INPUT_DIR, "v3_sparse_other_traits.gs")
UNIFIED_SCORE_RDS <- file.path(UNIFIED_INPUT_DIR, "cite_rna_scdrs_scores_v3_sparse_other_traits.rds")
UNIFIED_SCORE_CSV <- file.path(UNIFIED_INPUT_DIR, "cite_rna_scdrs_scores_v3_sparse_other_traits.csv.gz")
SEED_CELL_FILE <- file.path(SEED_DIAG_DIR, "sparse_neighborhood_seed_cell_level.csv.gz")

dir.create(UNIFIED_INPUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat(">>> module4 v3 sparse pipeline for non-MONO traits\n\n")
cat("RUN_STEPS:\n")
cat(paste0("  - ", RUN_STEPS, collapse = "\n"), "\n\n", sep = "")
cat("Traits: ", paste(TRAITS_TO_RUN, collapse = ", "), "\n", sep = "")
cat("Tag   : ", V3_SPARSE_OTHER_TAG, "\n\n", sep = "")
cat("Sparse mapping: B_NULL_CELL=", V3_B_NULL_CELL,
    "; RNA_P=", V3_RNA_EMPIRICAL_P_THRESHOLD,
    "; ATAC_P=", V3_ATAC_EMPIRICAL_P_THRESHOLD, "\n", sep = "")
cat("Completion: B_NULL_TOP5=", V3_SPARSE_COMPLETION_B_NULL_TOP5,
    "; EMPIRICAL_P=", V3_SPARSE_COMPLETION_EMPIRICAL_P_THRESHOLD, "\n", sep = "")
cat("Seed: RNA neighborhood FDR=", V3_GROUP_FDR_THRESHOLD,
    "; candidate_group_weight_ge0p80\n\n", sep = "")

if (step_enabled("prepare_unified_inputs")) {
  suppressPackageStartupMessages({ library(data.table) })
  cat("[prepare] Writing unified non-MONO gene-set file ...\n")
  write_unified_gene_set(UNIFIED_GS_FILE, TRAIT_IDS)
  cat("  wrote: ", UNIFIED_GS_FILE, "\n", sep = "")

  cat("[prepare] Writing unified CITE score RDS ...\n")
  merge_score_rds(UNIFIED_SCORE_RDS, UNIFIED_SCORE_CSV, TRAITS_TO_RUN)
  cat("  wrote: ", UNIFIED_SCORE_RDS, "\n", sep = "")
}

if (step_enabled("sparse_mapping")) {
  run_script("02_trait_mapping/08_build_sparse_trait_mapping.R", list(
    SPARSE_MAP_OUT_DIR = SPARSE_MAP_DIR,
    DIAG_DIR = SPARSE_MAP_DIAG_DIR,
    CITE_EXTRACT_DIR = CITE_EXTRACT_DIR,
    GS_FILE = UNIFIED_GS_FILE,
    CITE_SCORE_RDS = UNIFIED_SCORE_RDS,
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
  run_script("02_trait_mapping/07_build_top5_significant_completion.R", list(
    MAP_DIR = SPARSE_MAP_DIR,
    TRANSFER_OUT_DIR = TRANSFER_DIR,
    DIAG_DIR = COMPLETION_DIAG_DIR,
    CITE_MTX = file.path(CITE_EXTRACT_DIR, "cite_trait_counts.mtx"),
    CITE_GENES = file.path(CITE_EXTRACT_DIR, "cite_trait_counts_genes.txt"),
    CITE_CELLS = file.path(CITE_EXTRACT_DIR, "cite_trait_counts_cells.txt"),
    CITE_SCORE_RDS = UNIFIED_SCORE_RDS,
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

if (step_enabled("seed_candidates")) {
  run_script("02_trait_mapping/09_screen_sparse_mapping_neighborhood_seeds.R", list(
    MAP_DIR = SPARSE_MAP_DIR,
    OUT_DIR = SEED_DIAG_DIR,
    TRAITS_TO_SCREEN = TRAITS_TO_RUN,
    RNA_GROUP_FILE = file.path(
      ROOT,
      "04_evaluation/module4_v2_diagnostics/neighborhood_fdr/rna_unsupervised_group_membership.csv.gz"
    ),
    RNA_GROUP_TEST_FILES = c(
      file.path(ROOT, "04_evaluation/module4_v2_diagnostics/neighborhood_fdr/neighborhood_fdr_group_tests.csv"),
      file.path(ROOT, "01_data_prep_external_rna/output/trait_screening/cite_trait_screening_neighborhood_fdr_group_tests.csv")
    ),
    GROUP_TEST_TRAIT_MAP = c(EOS = "EOSINOPHIL", WBC = "WHITE"),
    GROUP_FDR_THRESHOLD = V3_GROUP_FDR_THRESHOLD,
    POSTHOC_TARGETS_LIST = POSTHOC_TARGETS_LIST,
    POSTHOC_TARGET_NAMES = POSTHOC_TARGET_NAMES
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
    OUT_CSV_NAME_NEW = paste0(V3_SPARSE_OTHER_TAG, "_neighborhood_fdr_ge0p80_auroc.csv"),
    BASE_EVAL_SCRIPT = file.path(PROJECT_ROOT, "scripts", "03_nonpaired_analysis", "55_eval_auroc.R"),
    PATCHED_SCRIPT = file.path(AUROC_DIR, paste0("eval_auroc.", V3_SPARSE_OTHER_TAG, "_patched.R")),
    METHOD_LABEL_NEW = V3_METHOD_LABEL,
    TARGET_REPLACEMENT = make_target_replacement(TRAITS_TO_RUN, POSTHOC_TARGETS_LIST),
    PRED_REPLACEMENT = make_pred_replacement(TRAITS_TO_RUN)
  ))
}

cat("\n=== v3 sparse non-MONO driver finished selected steps ===\n")
cat("Key result locations:\n")
cat("  Unified gene set: ", UNIFIED_GS_FILE, "\n", sep = "")
cat("  Unified CITE score RDS: ", UNIFIED_SCORE_RDS, "\n", sep = "")
cat("  Sparse mappings: ", SPARSE_MAP_DIR, "\n", sep = "")
cat("  Sparse mapping diagnostics: ", SPARSE_MAP_DIAG_DIR, "\n", sep = "")
cat("  Completion transfer: ", TRANSFER_DIR, "\n", sep = "")
cat("  Completion diagnostics: ", COMPLETION_DIAG_DIR, "\n", sep = "")
cat("  Seed diagnostics: ", SEED_DIAG_DIR, "\n", sep = "")
cat("  Module3 KNN/RWR/TRS: ", COMPAT_DIR, "\n", sep = "")
cat("  AUROC: ", file.path(AUROC_DIR, paste0(V3_SPARSE_OTHER_TAG, "_neighborhood_fdr_ge0p80_auroc.csv")), "\n", sep = "")
