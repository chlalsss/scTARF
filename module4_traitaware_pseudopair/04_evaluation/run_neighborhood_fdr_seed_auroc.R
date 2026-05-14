# module4: AUROC evaluation for neighborhood-FDR seed TRS
# -----------------------------------------------------------------------------
# Purpose:
#   Reuse the existing AUROC evaluator, but point it to the TRS files generated
#   from neighborhood-FDR seeds.
#
# Run in RStudio after the RWR script finishes:
#   setwd("D:/RStudio/work/module4_traitaware_pseudopair")
#   source("04_evaluation/run_neighborhood_fdr_seed_auroc.R")
# -----------------------------------------------------------------------------

ROOT <- "D:/RStudio/work/module4_traitaware_pseudopair"
BASE_EVAL_SCRIPT <- file.path(ROOT, "04_evaluation/eval_auroc.R")
V2_FILES_DIR_NEW <- file.path(ROOT, "03_pseudopair_pipeline/output/module3_v2_neighborhood_fdr_ge0p80")
OUT_DIR_NEW <- file.path(ROOT, "04_evaluation/module4_v2_neighborhood_fdr_ge0p80")
OUT_CSV_NAME_NEW <- "module4_v2_neighborhood_fdr_ge0p80_auroc.csv"
PATCHED_SCRIPT <- file.path(OUT_DIR_NEW, "eval_auroc.neighborhood_fdr_ge0p80_patched.R")

dir.create(OUT_DIR_NEW, recursive = TRUE, showWarnings = FALSE)

stop_missing <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0L) {
    stop("Missing required file(s):\n  ", paste(missing, collapse = "\n  "),
         call. = FALSE)
  }
}

patch_eval_script <- function(in_file, out_file) {
  txt <- readLines(in_file, warn = FALSE)
  txt <- sub(
    "^V2_FILES_DIR\\s*<-.*$",
    paste0('V2_FILES_DIR <- "', normalizePath(V2_FILES_DIR_NEW, winslash = "/", mustWork = TRUE), '"'),
    txt
  )
  txt <- sub(
    "^OUT_DIR\\s*<-.*$",
    paste0('OUT_DIR <- "', normalizePath(OUT_DIR_NEW, winslash = "/", mustWork = TRUE), '"'),
    txt
  )
  txt <- sub(
    "^OUT_CSV_NAME\\s*<-.*$",
    paste0('OUT_CSV_NAME <- "', OUT_CSV_NAME_NEW, '"'),
    txt
  )
  writeLines(txt, out_file, useBytes = TRUE)
  out_file
}

cat(">>> Run AUROC evaluation for neighborhood-FDR seed TRS\n\n")
cat("TRS files dir: ", V2_FILES_DIR_NEW, "\n", sep = "")
cat("Output dir   : ", OUT_DIR_NEW, "\n\n", sep = "")

trs_files <- list.files(
  V2_FILES_DIR_NEW,
  pattern = "^TRS_v2_[A-Z_]+_w[0-9.]+\\.tsv\\.gz$",
  full.names = TRUE
)
if (length(trs_files) == 0L) {
  stop("No TRS_v2_*.tsv.gz files found in: ", V2_FILES_DIR_NEW,
       call. = FALSE)
}
cat("Found TRS files:\n")
cat(paste0("  ", basename(trs_files), collapse = "\n"), "\n\n", sep = "")

stop_missing(BASE_EVAL_SCRIPT)
patch_eval_script(BASE_EVAL_SCRIPT, PATCHED_SCRIPT)
cat("Patched evaluator: ", PATCHED_SCRIPT, "\n\n", sep = "")

source(PATCHED_SCRIPT, local = new.env(parent = globalenv()))

cat("\n=== neighborhood-FDR seed AUROC evaluation finished ===\n")
cat("Result CSV: ", file.path(OUT_DIR_NEW, OUT_CSV_NAME_NEW), "\n", sep = "")
