# RStudio entry point for reverse ATAC -> RNA experiment after Step 0.
#
# Usage from the repository root:
#   source("scripts/06_reverse_analysis/06_run_reverse_step1_to_step4_from_rstudio.R")
#
# To resume, edit RUN_STEPS before source(), for example:
#   RUN_STEPS <- c("prepare_inputs", "knn_rwr", "auroc")
#   source("scripts/06_reverse_analysis/06_run_reverse_step1_to_step4_from_rstudio.R")

current_script_dir <- function(default) {
  path <- tryCatch(normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = FALSE),
                   error = function(e) "")
  if (nzchar(path)) dirname(path) else normalizePath(default, winslash = "/", mustWork = FALSE)
}

CODE_ROOT <- current_script_dir("scripts/06_reverse_analysis")
PROJECT_ROOT <- normalizePath(file.path(CODE_ROOT, "..", ".."), winslash = "/", mustWork = FALSE)
ROOT <- Sys.getenv("SCTARF_WORK_ROOT", PROJECT_ROOT)
REV_ROOT <- Sys.getenv("SCTARF_REVERSE_ROOT", file.path(PROJECT_ROOT, "outputs", "reverse"))

DRIVER_ENV <- environment()
has_driver <- function(x) exists(x, envir = DRIVER_ENV, inherits = FALSE)
get_driver <- function(x) get(x, envir = DRIVER_ENV, inherits = FALSE)

if (!has_driver("REVERSE_TAG")) REVERSE_TAG <- "reverse_granja_rna_main_plt_mch_rbc"
if (!has_driver("TRAITS_TO_RUN")) TRAITS_TO_RUN <- c("PLT", "MCH", "RBC")
if (!has_driver("ATAC_W_LIST")) ATAC_W_LIST <- c(0.7, 1.0)
if (!has_driver("ALLOW_OUTSIDE_REV_ROOT")) ALLOW_OUTSIDE_REV_ROOT <- FALSE

if (!has_driver("RUN_STEPS")) {
  RUN_STEPS <- c(
    "align",
    "transfer",
    "rna_main_z",
    "prepare_inputs",
    "knn_rwr",
    "auroc"
  )
}

step_enabled <- function(x) x %in% RUN_STEPS

OPTIONAL_PASSTHROUGH <- c(
  "K_NEIGHBOR",
  "BATCH_SIZE",
  "USE_GWAS_WEIGHT",
  "FORCE_REEXTRACT_RNA",
  "FORCE_RETRANSFER",
  "FORCE_REBUILD_RNA_MAIN_Z",
  "FORCE_REEXTRACT_RNA_MAIN",
  "FORCE_REMAP",
  "RNA_N_HVG",
  "RNA_N_PCS",
  "SEED_TOP_PCT"
)

run_local_script <- function(script_name, vars = list()) {
  script_alias <- c(
    "02_reverse_align_granja_rna_to_atac.R" = "07_reverse_align_granja_rna_to_atac.R",
    "03_reverse_transfer_atac_to_granja_rna.R" = "08_reverse_transfer_atac_to_granja_rna.R",
    "04_reverse_build_granja_rna_main_z_block.R" = "09_reverse_build_granja_rna_main_z_block.R",
    "05_reverse_prepare_module3_inputs.R" = "10_reverse_prepare_module3_inputs.R",
    "06_reverse_run_module3_knn_rwr.R" = "11_reverse_run_module3_knn_rwr.R",
    "07_reverse_eval_auroc.R" = "12_reverse_eval_auroc.R"
  )
  e <- new.env(parent = globalenv())
  assign("ROOT", ROOT, envir = e)
  assign("REV_ROOT", REV_ROOT, envir = e)
  assign("SCTARF_REVERSE_SCRIPT_DIR", CODE_ROOT, envir = e)
  assign("SCTARF_NONPAIRED_SCRIPT_DIR", file.path(PROJECT_ROOT, "scripts", "03_nonpaired_analysis"), envir = e)
  assign("BASE_EVAL_SCRIPT", file.path(PROJECT_ROOT, "scripts", "03_nonpaired_analysis", "55_eval_auroc.R"), envir = e)
  assign("REVERSE_TAG", REVERSE_TAG, envir = e)
  assign("TRAITS_TO_RUN", TRAITS_TO_RUN, envir = e)
  assign("ATAC_W_LIST", ATAC_W_LIST, envir = e)
  assign("ALLOW_OUTSIDE_REV_ROOT", ALLOW_OUTSIDE_REV_ROOT, envir = e)
  for (nm in OPTIONAL_PASSTHROUGH) {
    if (has_driver(nm)) assign(nm, get_driver(nm), envir = e)
  }
  for (nm in names(vars)) assign(nm, vars[[nm]], envir = e)
  script_path <- file.path(CODE_ROOT, if (script_name %in% names(script_alias)) script_alias[[script_name]] else script_name)
  cat("\n", strrep("=", 78), "\n", sep = "")
  cat("source: ", script_path, "\n", sep = "")
  cat(strrep("=", 78), "\n", sep = "")
  source(script_path, local = e)
  invisible(e)
}

cat(">>> Reverse ATAC -> RNA pipeline driver\n\n")
cat("Tag    : ", REVERSE_TAG, "\n", sep = "")
cat("Traits : ", paste(TRAITS_TO_RUN, collapse = ", "), "\n", sep = "")
cat("Weights: ", paste(ATAC_W_LIST, collapse = ", "), "  (weight on main RNA block)\n", sep = "")
cat("Steps  : ", paste(RUN_STEPS, collapse = " -> "), "\n\n", sep = "")

if (step_enabled("align")) {
  run_local_script("02_reverse_align_granja_rna_to_atac.R")
}

if (step_enabled("transfer")) {
  run_local_script("03_reverse_transfer_atac_to_granja_rna.R")
}

if (step_enabled("rna_main_z")) {
  run_local_script("04_reverse_build_granja_rna_main_z_block.R")
}

if (step_enabled("prepare_inputs")) {
  run_local_script("05_reverse_prepare_module3_inputs.R")
}

if (step_enabled("knn_rwr")) {
  run_local_script("06_reverse_run_module3_knn_rwr.R")
}

if (step_enabled("auroc")) {
  run_local_script("07_reverse_eval_auroc.R")
}

cat("\n=== Reverse driver finished selected steps ===\n")
cat("Key output locations:\n")
cat("  Mapping : ", file.path(REV_ROOT, "outputs", REVERSE_TAG, "mapping"), "\n", sep = "")
cat("  Transfer: ", file.path(REV_ROOT, "outputs", REVERSE_TAG, "transfer"), "\n", sep = "")
cat("  Compat  : ", file.path(REV_ROOT, "outputs", REVERSE_TAG, "module3_v2_compat"), "\n", sep = "")
cat("  AUROC   : ", file.path(REV_ROOT, "outputs", REVERSE_TAG, "auroc"), "\n", sep = "")
