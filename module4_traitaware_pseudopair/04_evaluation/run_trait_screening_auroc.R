# module4 trait screening: AUROC wrapper for configurable trait proxy targets
# -----------------------------------------------------------------------------
# Reuses 04_evaluation/eval_auroc.R by writing a patched evaluator.
#
# Target definition:
#   By default this script evaluates the original red-lineage trait-screening
#   set. External drivers can replace TARGETS and PRED_TRAIT at runtime, e.g.
#   RUN_V3_SPARSE_OTHER_TRAITS.R supplies proxy targets for PLT, LYM, MCH, RDW,
#   RBC, EOS, and WBC.
#
# Run in RStudio after module3 finishes:
#   setwd("D:/RStudio/work/module4_traitaware_pseudopair")
#   source("04_evaluation/run_trait_screening_auroc.R")
# -----------------------------------------------------------------------------

ROOT <- "D:/RStudio/work/module4_traitaware_pseudopair"
if (!exists("BASE_EVAL_SCRIPT")) BASE_EVAL_SCRIPT <- file.path(ROOT, "04_evaluation/eval_auroc.R")
if (!exists("V2_FILES_DIR_NEW")) {
  V2_FILES_DIR_NEW <- file.path(ROOT, "03_pseudopair_pipeline/output/module3_v2_trait_screening_ge0p80")
}
if (!exists("OUT_DIR_NEW")) OUT_DIR_NEW <- file.path(ROOT, "04_evaluation/trait_screening_ge0p80")
if (!exists("OUT_CSV_NAME_NEW")) OUT_CSV_NAME_NEW <- "trait_screening_ge0p80_auroc.csv"
if (!exists("PATCHED_SCRIPT")) {
  PATCHED_SCRIPT <- file.path(OUT_DIR_NEW, "eval_auroc.trait_screening_ge0p80_patched.R")
}
if (!exists("METHOD_LABEL_NEW")) METHOD_LABEL_NEW <- "v2"
dir.create(OUT_DIR_NEW, recursive = TRUE, showWarnings = FALSE)

if (!exists("TARGET_REPLACEMENT")) {
  TARGET_REPLACEMENT <- c(
    'TARGETS <- list(',
    '  MCH = c("02_Early.Eryth", "03_Late.Eryth", "04_Early.Baso"),',
    '  RDW = c("02_Early.Eryth", "03_Late.Eryth", "04_Early.Baso"),',
    '  RBC = c("02_Early.Eryth", "03_Late.Eryth", "04_Early.Baso")',
    ')'
  )
}
if (!exists("PRED_REPLACEMENT")) {
  PRED_REPLACEMENT <- c(
    'PRED_TRAIT <- list(MCH = "MCH", RDW = "RDW", RBC = "RBC")',
    ''
  )
}

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

  target_start <- grep("^TARGETS\\s*<-\\s*list\\(", txt)[1]
  pred_start <- grep("^PRED_TRAIT\\s*<-\\s*list\\(", txt)[1]
  if (is.na(target_start) || is.na(pred_start) || pred_start <= target_start) {
    stop("Could not locate TARGETS/PRED_TRAIT blocks in eval script.", call. = FALSE)
  }
  txt <- c(
    txt[seq_len(target_start - 1L)],
    TARGET_REPLACEMENT,
    txt[pred_start:length(txt)]
  )

  pred_start <- grep("^PRED_TRAIT\\s*<-\\s*list\\(", txt)[1]
  boot_start <- grep("^BOOTSTRAP_CI\\s*<-", txt)[1]
  if (is.na(pred_start) || is.na(boot_start) || boot_start <= pred_start) {
    stop("Could not locate PRED_TRAIT block in eval script.", call. = FALSE)
  }
  txt <- c(
    txt[seq_len(pred_start - 1L)],
    PRED_REPLACEMENT,
    txt[boot_start:length(txt)]
  )

  writeLines(txt, out_file, useBytes = TRUE)
  out_file
}

relabel_auroc_method <- function(path, label) {
  if (identical(label, "v2")) return(invisible(FALSE))
  if (!file.exists(path)) return(invisible(FALSE))
  out <- read.csv(path, check.names = FALSE)
  if (!"method" %in% colnames(out)) return(invisible(FALSE))
  out$method[out$method == "v2"] <- label
  write.csv(out, path, row.names = FALSE, quote = TRUE)
  invisible(TRUE)
}

cat(">>> Run AUROC evaluation for trait screening\n\n")
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

result_csv <- file.path(OUT_DIR_NEW, OUT_CSV_NAME_NEW)
relabel_auroc_method(result_csv, METHOD_LABEL_NEW)

cat("\n=== trait-screening AUROC evaluation finished ===\n")
cat("Result CSV: ", result_csv, "\n", sep = "")
