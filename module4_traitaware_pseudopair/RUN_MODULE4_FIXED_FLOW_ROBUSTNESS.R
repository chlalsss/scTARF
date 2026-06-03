# module4 fixed-flow raw ATAC drop-out robustness
# -----------------------------------------------------------------------------
# Current module4 status after RNA-reference screening:
#   - LYM uses Granja RNA reference and is now a validated positive example.
#   - MCH/RDW/RBC use the external CITE RNA reference; RDW improves, MCH/RBC are
#     near-ceiling stress cases.
#   - PLT/EOS/WBC are not run by default because the current reference-screening
#     result says they need a better/split RNA reference first.
#
# This driver dispatches the full raw ATAC drop-out pipeline by RNA-reference
# family. It reuses:
#   05_robustness/10_run_v3_sparse_mono_raw_atac_full_pipeline_dropout.R
#
# Run in RStudio:
#   setwd("D:/RStudio/work/module4_traitaware_pseudopair")
#   source("RUN_MODULE4_FIXED_FLOW_ROBUSTNESS.R")
#
# Optional one-shot examples:
#   MODULE4_FIXED_TRAITS_TO_RUN <- c("LYM")
#   MODULE4_ROBUST_P_DROP_LIST <- c(0)
#   MODULE4_ROBUST_SEED_LIST <- c(1)
#   MODULE4_ROBUST_SPARSE_MAP_BATCH_SIZE <- 128L
#   MODULE4_DRY_RUN <- TRUE
#   source("RUN_MODULE4_FIXED_FLOW_ROBUSTNESS.R")
#
#   MODULE4_FIXED_TRAITS_TO_RUN <- c("RDW", "RBC")
#   source("RUN_MODULE4_FIXED_FLOW_ROBUSTNESS.R")
# -----------------------------------------------------------------------------

ROOT <- "D:/RStudio/work/module4_traitaware_pseudopair"
ROBUST_RUNNER <- file.path(
  ROOT,
  "05_robustness/10_run_v3_sparse_mono_raw_atac_full_pipeline_dropout.R"
)

SCRIPT_ENV <- environment()

consume_one_shot <- function(name, default) {
  if (!exists(name, envir = SCRIPT_ENV, inherits = FALSE)) return(default)
  value <- get(name, envir = SCRIPT_ENV, inherits = FALSE)
  rm(list = name, envir = SCRIPT_ENV)
  cat("Using one-shot setting ", name, ": ",
      paste(value, collapse = ", "), "\n", sep = "")
  value
}

drop_legacy_override <- function(name) {
  if (!exists(name, envir = SCRIPT_ENV, inherits = FALSE)) return(invisible(FALSE))
  rm(list = name, envir = SCRIPT_ENV)
  warning(
    "Ignored and removed legacy setting ", name,
    ". Use MODULE4_* one-shot settings shown at the top of this file.",
    call. = FALSE
  )
  invisible(TRUE)
}

invisible(lapply(c(
  "FIXED_TRAITS_TO_RUN_OVERRIDE",
  "ROBUST_P_DROP_LIST_OVERRIDE",
  "ROBUST_SEED_LIST_OVERRIDE",
  "ROBUST_ATAC_W_LIST_OVERRIDE"
), drop_legacy_override))

FIXED_TRAITS_TO_RUN <- consume_one_shot(
  "MODULE4_FIXED_TRAITS_TO_RUN",
  c("LYM", "MCH", "RDW", "RBC")
)
ROBUST_P_DROP_LIST <- consume_one_shot(
  "MODULE4_ROBUST_P_DROP_LIST",
  round(seq(0.0, 0.9, by = 0.1), 2)
)
ROBUST_SEED_LIST <- consume_one_shot(
  "MODULE4_ROBUST_SEED_LIST",
  c(1L, 2L, 3L)
)
ROBUST_ATAC_W_LIST <- consume_one_shot(
  "MODULE4_ROBUST_ATAC_W_LIST",
  c(0.7, 1.0)
)
ROBUST_SPARSE_MAP_BATCH_SIZE <- consume_one_shot(
  "MODULE4_ROBUST_SPARSE_MAP_BATCH_SIZE",
  128L
)
MODULE4_DRY_RUN <- consume_one_shot("MODULE4_DRY_RUN", FALSE)

unsupported <- setdiff(FIXED_TRAITS_TO_RUN, c("LYM", "MCH", "RDW", "RBC"))
if (length(unsupported) > 0L) {
  stop(
    "Unsupported in current fixed-flow robustness default: ",
    paste(unsupported, collapse = ", "),
    "\nCurrent reference-screening status: PLT/EOS need new RNA references; ",
    "WBC needs component-wise targets; MONO has its own MONO robustness route.",
    call. = FALSE
  )
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

stop_missing <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0L) {
    stop("Missing required file(s):\n  ", paste(missing, collapse = "\n  "),
         call. = FALSE)
  }
}

run_reference_group <- function(group_name, traits, cfg) {
  traits <- intersect(traits, FIXED_TRAITS_TO_RUN)
  if (length(traits) == 0L) return(invisible(NULL))

  e <- new.env(parent = globalenv())
  assign("ROOT", ROOT, envir = e)
  assign("ROBUST_OUTPUT_TAG", cfg$output_tag, envir = e)
  assign("ROBUST_TRAITS_TO_RUN", traits, envir = e)
  assign("ROBUST_RUN_LABEL", cfg$run_label, envir = e)
  assign("ROBUST_METHOD_LABEL", cfg$method_label, envir = e)
  assign("ROBUST_SEED_MODE", "sparse_neighborhood_screen", envir = e)
  assign("ROBUST_NBD_MODE", "recompute", envir = e)
  assign("ROBUST_P0_USE_CLEAN_RAW_PRODUCTS", TRUE, envir = e)
  assign("ROBUST_CACHE_DEGRADED_COUNTS", FALSE, envir = e)
  assign("ROBUST_SPARSE_MAP_BATCH_SIZE", ROBUST_SPARSE_MAP_BATCH_SIZE, envir = e)
  assign("ROBUST_P_DROP_LIST", ROBUST_P_DROP_LIST, envir = e)
  assign("ROBUST_SEED_LIST", ROBUST_SEED_LIST, envir = e)
  assign("ROBUST_ATAC_W_LIST", ROBUST_ATAC_W_LIST, envir = e)
  assign("ROBUST_RUN_STEPS", c(
    "raw_products",
    "sparse_mapping",
    "top5sig_transfer",
    "neighborhood_fdr",
    "module3_blocks",
    "knn_rwr",
    "auroc",
    "summarize"
  ), envir = e)

  assign("ROBUST_AUROC_RUNNER", "04_evaluation/run_trait_screening_auroc.R", envir = e)
  assign("ROBUST_AUROC_DIR_LABEL", paste0(cfg$run_label, "_neighborhood_fdr_ge0p80"), envir = e)
  assign("ROBUST_TRAITS_ALL", cfg$trait_ids[traits], envir = e)
  assign("ROBUST_GS_FILE", cfg$gs_file, envir = e)
  assign("ROBUST_CITE_H5AD", cfg$h5ad, envir = e)
  assign("BASE_CITE_EXTRACT_DIR", cfg$extract_dir, envir = e)
  assign("CITE_SCORE_RDS", cfg$score_rds, envir = e)
  assign("CLEAN_COMPAT_DIR", cfg$clean_compat_dir, envir = e)
  assign("CLEAN_ATAC_Z_PATH", file.path(cfg$clean_compat_dir, "atac_z_block_v2.rds"), envir = e)
  assign("CLEAN_AUROC_CSV", cfg$clean_auroc_csv, envir = e)
  assign("ROBUST_RNA_GROUP_FILE", cfg$rna_group_file, envir = e)
  assign("ROBUST_RNA_GROUP_TEST_FILES", cfg$rna_group_test_files, envir = e)
  assign("ROBUST_GROUP_TEST_TRAIT_MAP", cfg$group_test_trait_map, envir = e)
  assign("ROBUST_GROUP_FDR_THRESHOLD", 0.10, envir = e)
  assign("ROBUST_POSTHOC_TARGETS_LIST", cfg$targets[traits], envir = e)
  assign("ROBUST_POSTHOC_TARGET_NAMES", cfg$target_names[traits], envir = e)
  assign("ROBUST_TARGET_REPLACEMENT", make_target_replacement(traits, cfg$targets), envir = e)
  assign("ROBUST_PRED_REPLACEMENT", make_pred_replacement(traits), envir = e)

  required <- c(
    ROBUST_RUNNER,
    cfg$gs_file,
    cfg$score_rds,
    cfg$extract_dir,
    file.path(cfg$extract_dir, "cite_trait_counts.mtx"),
    file.path(cfg$extract_dir, "cite_trait_counts_genes.txt"),
    file.path(cfg$extract_dir, "cite_trait_counts_cells.txt"),
    cfg$clean_compat_dir,
    file.path(cfg$clean_compat_dir, "atac_z_block_v2.rds"),
    cfg$clean_auroc_csv,
    cfg$rna_group_file,
    cfg$rna_group_test_files
  )
  if (!is.null(cfg$h5ad) && nzchar(cfg$h5ad)) required <- c(required, cfg$h5ad)
  stop_missing(required)

  cat("\n", strrep("=", 78), "\n", sep = "")
  cat("Reference group: ", group_name, "\n", sep = "")
  cat("Traits         : ", paste(traits, collapse = ", "), "\n", sep = "")
  cat("Output tag     : ", cfg$output_tag, "\n", sep = "")
  cat(strrep("=", 78), "\n", sep = "")

  source(ROBUST_RUNNER, local = e)
  invisible(TRUE)
}

red_proxy <- c("02_Early.Eryth", "03_Late.Eryth", "04_Early.Baso")
lymph_proxy <- c(
  "16_Pre.B", "17_B", "18_Plasma",
  "19_CD8.N", "20_CD4.N1", "21_CD4.N2", "22_CD4.M",
  "23_CD8.EM", "24_CD8.CM", "25_NK"
)

group_configs <- list(
  lym_granja = list(
    traits = c("LYM"),
    output_tag = "v3_sparse_lym_granja_rna_raw_atac_full_pipeline_dropout",
    run_label = "lym_granja_rna",
    method_label = "v3_sparse_lym_granja_rna_raw_full_pipeline",
    trait_ids = c(LYM = "UKB_460K.blood_LYMPHOCYTE_COUNT"),
    targets = list(LYM = lymph_proxy),
    target_names = c(LYM = "lymph_proxy"),
    gs_file = file.path(ROOT, "01_data_prep/03_rna_scdrs_python/granja_3traits.gs"),
    h5ad = file.path(ROOT, "01_data_prep/03_rna_scdrs_python/input/granja_rna_raw_counts.h5ad"),
    extract_dir = file.path(ROOT, "02_trait_mapping/output/v3_sparse_lym_granja_rna/granja_rna_extract"),
    score_rds = file.path(
      ROOT,
      "02_trait_mapping/output/v3_sparse_lym_granja_rna/_inputs/granja_rna_scores_LYM_for_v3_sparse.rds"
    ),
    clean_compat_dir = file.path(
      ROOT,
      "03_pseudopair_pipeline/output/v3_sparse_lym_granja_rna/module3_neighborhood_fdr_ge0p80"
    ),
    clean_auroc_csv = file.path(
      ROOT,
      "04_evaluation/v3_sparse_lym_granja_rna/neighborhood_fdr_ge0p80",
      "v3_sparse_lym_granja_rna_neighborhood_fdr_ge0p80_auroc.csv"
    ),
    rna_group_file = file.path(
      ROOT,
      "02_trait_mapping/output/v3_sparse_lym_granja_rna/rna_neighborhood_fdr",
      "rna_unsupervised_group_membership.csv.gz"
    ),
    rna_group_test_files = file.path(
      ROOT,
      "02_trait_mapping/output/v3_sparse_lym_granja_rna/rna_neighborhood_fdr",
      "neighborhood_fdr_group_tests.csv"
    ),
    group_test_trait_map = c()
  ),
  cite_erythroid = list(
    traits = c("MCH", "RDW", "RBC"),
    output_tag = "v3_sparse_cite_erythroid_raw_atac_full_pipeline_dropout",
    run_label = "cite_erythroid",
    method_label = "v3_sparse_cite_erythroid_raw_full_pipeline",
    trait_ids = c(
      MCH = "UKB_460K.blood_MEAN_CORPUSCULAR_HEMOGLOBIN",
      RDW = "UKB_460K.blood_RBC_DISTRIB_WIDTH",
      RBC = "UKB_460K.blood_RED_COUNT"
    ),
    targets = list(MCH = red_proxy, RDW = red_proxy, RBC = red_proxy),
    target_names = c(
      MCH = "red_lineage_proxy",
      RDW = "red_lineage_proxy",
      RBC = "red_lineage_proxy"
    ),
    gs_file = file.path(
      ROOT,
      "02_trait_mapping/output/v3_sparse_other_traits/_inputs/v3_sparse_other_traits.gs"
    ),
    h5ad = file.path(
      ROOT,
      "01_data_prep_external_rna/03_cite_scdrs_python/input/cite_rna_raw_counts_gex_only.h5ad"
    ),
    extract_dir = file.path(ROOT, "02_trait_mapping/output/v3_sparse_other_traits/cite_extract"),
    score_rds = file.path(
      ROOT,
      "02_trait_mapping/output/v3_sparse_other_traits/_inputs/cite_rna_scdrs_scores_v3_sparse_other_traits.rds"
    ),
    clean_compat_dir = file.path(
      ROOT,
      "03_pseudopair_pipeline/output/v3_sparse_other_traits/module3_neighborhood_fdr_ge0p80"
    ),
    clean_auroc_csv = file.path(
      ROOT,
      "04_evaluation/v3_sparse_other_traits/neighborhood_fdr_ge0p80",
      "v3_sparse_other_traits_neighborhood_fdr_ge0p80_auroc.csv"
    ),
    rna_group_file = file.path(
      ROOT,
      "04_evaluation/module4_v2_diagnostics/neighborhood_fdr/rna_unsupervised_group_membership.csv.gz"
    ),
    rna_group_test_files = c(
      file.path(ROOT, "04_evaluation/module4_v2_diagnostics/neighborhood_fdr/neighborhood_fdr_group_tests.csv"),
      file.path(ROOT, "01_data_prep_external_rna/output/trait_screening/cite_trait_screening_neighborhood_fdr_group_tests.csv")
    ),
    group_test_trait_map = c()
  )
)

cat(">>> module4 fixed-flow robustness driver\n\n")
cat("Current supported traits: LYM, MCH, RDW, RBC\n")
cat("Requested traits        : ", paste(FIXED_TRAITS_TO_RUN, collapse = ", "), "\n", sep = "")
cat("p_drop                  : ", paste(ROBUST_P_DROP_LIST, collapse = ", "), "\n", sep = "")
cat("seeds                   : ", paste(ROBUST_SEED_LIST, collapse = ", "), "\n", sep = "")
cat("ATAC weights            : ", paste(ROBUST_ATAC_W_LIST, collapse = ", "), "\n", sep = "")
cat("sparse map batch size   : ", ROBUST_SPARSE_MAP_BATCH_SIZE, "\n", sep = "")
cat("\nDispatch plan:\n")
for (nm in names(group_configs)) {
  cfg <- group_configs[[nm]]
  group_traits <- intersect(cfg$traits, FIXED_TRAITS_TO_RUN)
  if (length(group_traits) == 0L) next
  cat("  - ", cfg$run_label, ": ", paste(group_traits, collapse = ", "),
      " -> ", cfg$output_tag, "\n", sep = "")
}

if (isTRUE(MODULE4_DRY_RUN)) {
  cat("\nDry run only; no robustness jobs launched.\n")
} else {
  run_reference_group("LYM / Granja RNA", group_configs$lym_granja$traits, group_configs$lym_granja)
  run_reference_group("MCH-RDW-RBC / CITE RNA", group_configs$cite_erythroid$traits, group_configs$cite_erythroid)

  cat("\n=== module4 fixed-flow robustness dispatch finished ===\n")
}
