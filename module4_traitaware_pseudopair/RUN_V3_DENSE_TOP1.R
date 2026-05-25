# module4 v3.1: dense G_trait + top1 RNA mapping
# -----------------------------------------------------------------------------
# Purpose:
#   Test whether the current top5 RNA-neighbor weighted average is over-smoothing
#   cell-specific trait signal.
#
# What changes:
#   - Keep the original dense G_trait feature space.
#   - Set K_NEIGHBOR = 1 in trait-aware mapping.
#   - Recompute transfer, RNA-neighborhood seed evidence, module3 KNN/RWR/TRS,
#     and AUROC in v3-specific output directories.
#
# What does not change:
#   - No cell type is used in mapping, seed selection, KNN, RWR, or TRS.
#   - AUROC still uses the shared 04_evaluation/eval_auroc.R evaluator.
#
# Run in RStudio:
#   setwd("D:/RStudio/work/module4_traitaware_pseudopair")
#   source("RUN_V3_DENSE_TOP1.R")
#
# Optional partial run:
#   RUN_STEPS <- c("map_mono", "transfer_mono")
#   source("RUN_V3_DENSE_TOP1.R")
# -----------------------------------------------------------------------------

ROOT <- "D:/RStudio/work/module4_traitaware_pseudopair"

if (!exists("RUN_STEPS")) {
  RUN_STEPS <- c(
    "map_mono",
    "transfer_mono",
    "module3_blocks_mono",
    "knn_mono",
    "neighborhood_fdr_mono",
    "neighborhood_seed_mono",
    "rwr_mono",
    "auroc_mono",
    "map_rdw_rbc",
    "transfer_rdw_rbc",
    "screen_rdw_rbc",
    "module3_rdw_rbc",
    "run_rdw_rbc",
    "auroc_rdw_rbc"
  )
}

if (!exists("V3_ATAC_W_LIST")) V3_ATAC_W_LIST <- c(0.7, 1.0)

run_script <- function(script, vars = list()) {
  e <- new.env(parent = globalenv())
  for (nm in names(vars)) assign(nm, vars[[nm]], envir = e)
  cat("\n>>> source(", script, ")\n", sep = "")
  source(file.path(ROOT, script), local = e)
  invisible(e)
}

step_enabled <- function(x) x %in% RUN_STEPS

V3_ROOT <- file.path(ROOT, "v3_dense_top1")

MONO_MAP_DIR <- file.path(ROOT, "02_trait_mapping/output/v3_dense_top1/mono")
MONO_TRANSFER_DIR <- file.path(ROOT, "03_pseudopair_pipeline/output/v3_dense_top1/transfer_mono")
MONO_COMPAT_DIR <- file.path(ROOT, "03_pseudopair_pipeline/output/v3_dense_top1/module3_mono_dense_seed")
MONO_NBD_DIR <- file.path(ROOT, "04_evaluation/v3_dense_top1/neighborhood_fdr_mono")
MONO_NBD_COMPAT_DIR <- file.path(ROOT, "03_pseudopair_pipeline/output/v3_dense_top1/module3_mono_neighborhood_fdr_ge0p80")
MONO_AUROC_DIR <- file.path(ROOT, "04_evaluation/v3_dense_top1/mono_neighborhood_fdr_ge0p80")

TR_MAP_DIR <- file.path(ROOT, "02_trait_mapping/output/v3_dense_top1/trait_screening")
TR_TRANSFER_DIR <- file.path(ROOT, "03_pseudopair_pipeline/output/v3_dense_top1/transfer_trait_screening")
TR_DIAG_DIR <- file.path(ROOT, "02_trait_mapping/output/v3_dense_top1/trait_screening/diagnostics")
TR_COMPAT_DIR <- file.path(ROOT, "03_pseudopair_pipeline/output/v3_dense_top1/module3_trait_screening_ge0p80")
TR_AUROC_DIR <- file.path(ROOT, "04_evaluation/v3_dense_top1/trait_screening_ge0p80")

cat(">>> module4 v3.1 dense G_trait + top1 mapping\n\n")
cat("RUN_STEPS:\n")
cat(paste0("  - ", RUN_STEPS, collapse = "\n"), "\n\n", sep = "")
cat("Outputs are isolated under v3_dense_top1-specific directories.\n\n")

if (step_enabled("map_mono")) {
  run_script("02_trait_mapping/01_align_granja_cite_by_trait.R", list(
    ALIGN_OUT_DIR = MONO_MAP_DIR,
    TRAITS_ALL = c(MONO = "UKB_460K.blood_MONOCYTE_COUNT"),
    TRAITS_TO_RUN = c("MONO"),
    K_NEIGHBOR = 1L,
    BATCH_SIZE = 512L,
    USE_GWAS_WEIGHT = TRUE,
    FORCE_REEXTRACT_CITE = FALSE
  ))
}

if (step_enabled("transfer_mono")) {
  run_script("03_pseudopair_pipeline/01_transfer_dual.R", list(
    TRANSFER_MAP_DIR = MONO_MAP_DIR,
    TRANSFER_OUT_DIR = MONO_TRANSFER_DIR,
    CITE_MTX = file.path(MONO_MAP_DIR, "cite_extract/cite_trait_counts.mtx"),
    CITE_GENES = file.path(MONO_MAP_DIR, "cite_extract/cite_trait_counts_genes.txt"),
    CITE_CELLS = file.path(MONO_MAP_DIR, "cite_extract/cite_trait_counts_cells.txt"),
    TRAITS_ALL = c("MONO"),
    TRAITS_TO_RUN = c("MONO")
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

if (step_enabled("neighborhood_fdr_mono")) {
  run_script("04_evaluation/diagnose_scdrs_neighborhood_fdr.R", list(
    MAP_DIR = MONO_MAP_DIR,
    # Use the original v2 CITE RNA feature matrix for RNA neighborhoods, so this
    # step isolates the mapping change (top5 -> top1) rather than changing RNA
    # neighborhood construction at the same time.
    CITE_MTX = file.path(ROOT, "02_trait_mapping/output/cite_extract/cite_trait_counts.mtx"),
    CITE_GENES = file.path(ROOT, "02_trait_mapping/output/cite_extract/cite_trait_counts_genes.txt"),
    CITE_CELLS = file.path(ROOT, "02_trait_mapping/output/cite_extract/cite_trait_counts_cells.txt"),
    OUT_DIR = MONO_NBD_DIR,
    TRAITS = c("MONO"),
    RNA_KNN_K = 30L,
    LEIDEN_RESOLUTION = 0.005,
    GROUP_FDR_THRESHOLDS = c(0.10),
    TOP_N_ATAC = 350L
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
    OUT_CSV_NAME_NEW = "v3_dense_top1_mono_neighborhood_fdr_ge0p80_auroc.csv",
    PATCHED_SCRIPT = file.path(MONO_AUROC_DIR, "eval_auroc.v3_dense_top1_mono_patched.R"),
    METHOD_LABEL_NEW = "v3_dense_top1"
  ))
}

if (step_enabled("map_rdw_rbc")) {
  run_script("02_trait_mapping/01_align_granja_cite_by_trait.R", list(
    GS_FILE = file.path(ROOT, "01_data_prep_external_rna/03_cite_scdrs_python/cite_trait_screening.gs"),
    ALIGN_OUT_DIR = TR_MAP_DIR,
    TRAITS_ALL = c(
      RDW = "UKB_460K.blood_RBC_DISTRIB_WIDTH",
      RBC = "UKB_460K.blood_RED_COUNT"
    ),
    TRAITS_TO_RUN = c("RDW", "RBC"),
    K_NEIGHBOR = 1L,
    BATCH_SIZE = 512L,
    USE_GWAS_WEIGHT = TRUE,
    FORCE_REEXTRACT_CITE = FALSE
  ))
}

if (step_enabled("transfer_rdw_rbc")) {
  run_script("03_pseudopair_pipeline/01_transfer_dual.R", list(
    TRANSFER_MAP_DIR = TR_MAP_DIR,
    TRANSFER_OUT_DIR = TR_TRANSFER_DIR,
    CITE_MTX = file.path(TR_MAP_DIR, "cite_extract/cite_trait_counts.mtx"),
    CITE_GENES = file.path(TR_MAP_DIR, "cite_extract/cite_trait_counts_genes.txt"),
    CITE_CELLS = file.path(TR_MAP_DIR, "cite_extract/cite_trait_counts_cells.txt"),
    CITE_SCORE_RDS = file.path(
      ROOT,
      "01_data_prep_external_rna/output/trait_screening/cite_rna_scdrs_scores_trait_screening.rds"
    ),
    TRAITS_ALL = c("RDW", "RBC"),
    TRAITS_TO_RUN = c("RDW", "RBC")
  ))
}

if (step_enabled("screen_rdw_rbc")) {
  run_script("02_trait_mapping/04_screen_trait_mapping_atac_candidates.R", list(
    MAP_DIR = TR_MAP_DIR,
    OUT_DIR = TR_DIAG_DIR,
    TRAITS_TO_SCREEN = c("RDW", "RBC"),
    GROUP_FDR_THRESHOLD = 0.10,
    TOP_N_ATAC = 350L
  ))
}

if (step_enabled("module3_rdw_rbc")) {
  run_script("03_pseudopair_pipeline/07_prepare_trait_screening_module3_inputs.R", list(
    MODULE3_TRANSFER_DIR = TR_TRANSFER_DIR,
    MODULE3_COMPAT_OUT_DIR = TR_COMPAT_DIR,
    TRAITS_TO_RUN = c("RDW", "RBC"),
    SEED_CELL_FILE = file.path(TR_DIAG_DIR, "trait_screening_atac_candidate_cell_level.csv.gz"),
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
    OUT_CSV_NAME_NEW = "v3_dense_top1_trait_screening_ge0p80_auroc.csv",
    PATCHED_SCRIPT = file.path(TR_AUROC_DIR, "eval_auroc.v3_dense_top1_trait_screening_patched.R"),
    METHOD_LABEL_NEW = "v3_dense_top1",
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

cat("\n=== v3.1 dense top1 driver finished selected steps ===\n")
cat("Key result locations after a full run:\n")
cat("  MONO AUROC: ", file.path(MONO_AUROC_DIR, "v3_dense_top1_mono_neighborhood_fdr_ge0p80_auroc.csv"), "\n", sep = "")
cat("  RDW/RBC AUROC: ", file.path(TR_AUROC_DIR, "v3_dense_top1_trait_screening_ge0p80_auroc.csv"), "\n", sep = "")
cat("  MONO seed diagnostics: ", file.path(MONO_NBD_DIR, "neighborhood_fdr_seed_summary.csv"), "\n", sep = "")
cat("  RDW/RBC seed diagnostics: ", file.path(TR_DIAG_DIR, "trait_screening_atac_candidate_summary.csv"), "\n", sep = "")
