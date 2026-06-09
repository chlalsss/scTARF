# Reverse Step 4: AUROC on Granja RNA cells.
#
# Reuses module4's eval_auroc.R by writing a patched evaluator with:
#   META_FILE    = Granja RNA cell meta
#   BASELINE_TRS = Granja RNA scDRS baseline matrix
#   V2_FILES_DIR = reverse module3 outputs
#   TARGETS      = RNA-side Granja cell_type targets

suppressPackageStartupMessages({
  library(data.table)
})

ROOT <- "D:/RStudio/work"
MODULE4_ROOT <- file.path(ROOT, "module4_traitaware_pseudopair")
REV_ROOT <- file.path(ROOT, "reverse_atac_to_rna_experiment_20260608")

SCRIPT_ENV <- environment()
has_local <- function(x) exists(x, envir = SCRIPT_ENV, inherits = FALSE)

if (!has_local("TRAITS_TO_RUN")) TRAITS_TO_RUN <- c("PLT", "MCH", "RBC")
if (!has_local("REVERSE_TAG")) REVERSE_TAG <- "reverse_granja_rna_main_plt_mch_rbc"

if (!has_local("GRANJA_RNA_META_RDS")) {
  GRANJA_RNA_META_RDS <- file.path(MODULE4_ROOT, "01_data_prep/output/granja_rna_cell_meta.rds")
}
if (!has_local("GRANJA_RNA_SCORE_RDS")) {
  GRANJA_RNA_SCORE_RDS <- file.path(MODULE4_ROOT, "01_data_prep/output/granja_rna_scdrs_scores_with_missing_traits.rds")
}
if (!has_local("V2_FILES_DIR_NEW")) {
  V2_FILES_DIR_NEW <- file.path(REV_ROOT, "outputs", REVERSE_TAG, "module3_v2_compat")
}
if (!has_local("OUT_DIR_NEW")) {
  OUT_DIR_NEW <- file.path(REV_ROOT, "outputs", REVERSE_TAG, "auroc")
}
if (!has_local("OUT_CSV_NAME_NEW")) OUT_CSV_NAME_NEW <- paste0(REVERSE_TAG, "_auroc.csv")
if (!has_local("BASELINE_RDS")) {
  BASELINE_RDS <- file.path(REV_ROOT, "outputs", REVERSE_TAG, "_inputs", "granja_rna_scdrs_baseline_matrix.rds")
}
if (!has_local("PATCHED_SCRIPT")) {
  PATCHED_SCRIPT <- file.path(OUT_DIR_NEW, paste0("eval_auroc.", REVERSE_TAG, "_patched.R"))
}
if (!has_local("BASE_EVAL_SCRIPT")) BASE_EVAL_SCRIPT <- file.path(MODULE4_ROOT, "04_evaluation/eval_auroc.R")
if (!has_local("ALLOW_OUTSIDE_REV_ROOT")) ALLOW_OUTSIDE_REV_ROOT <- FALSE

RNA_TARGETS <- list(
  PLT = c("02_Early.Eryth", "03_Late.Eryth", "04_Early.Baso"),
  MCH = c("02_Early.Eryth", "03_Late.Eryth", "04_Early.Baso"),
  RBC = c("02_Early.Eryth", "03_Late.Eryth", "04_Early.Baso")
)

stop_missing <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0L) {
    stop("Missing required file(s):\n  ", paste(missing, collapse = "\n  "), call. = FALSE)
  }
}

assert_write_under_rev <- function(paths) {
  if (isTRUE(ALLOW_OUTSIDE_REV_ROOT)) return(invisible(TRUE))
  root <- normalizePath(REV_ROOT, winslash = "/", mustWork = TRUE)
  for (path in paths) {
    resolved <- normalizePath(path, winslash = "/", mustWork = FALSE)
    if (!startsWith(resolved, paste0(root, "/")) && !identical(resolved, root)) {
      stop("Refusing to write outside reverse workspace: ", path,
           "\nSet ALLOW_OUTSIDE_REV_ROOT <- TRUE only if this is intentional.",
           call. = FALSE)
    }
  }
  invisible(TRUE)
}

r_literal_chr <- function(x) paste(sprintf('"%s"', x), collapse = ", ")

make_target_replacement <- function(traits) {
  lines <- "TARGETS <- list("
  for (i in seq_along(traits)) {
    tr <- traits[[i]]
    suffix <- if (i < length(traits)) "," else ""
    lines <- c(lines, paste0("  ", tr, " = c(", r_literal_chr(RNA_TARGETS[[tr]]), ")", suffix))
  }
  c(lines, ")")
}

make_pred_replacement <- function(traits) {
  parts <- paste0(traits, ' = "', traits, '"')
  c(paste0("PRED_TRAIT <- list(", paste(parts, collapse = ", "), ")"), "")
}

write_baseline_rds <- function() {
  score <- as.data.frame(readRDS(GRANJA_RNA_SCORE_RDS))
  if (!"barcode" %in% colnames(score)) {
    stop("Granja RNA scDRS score RDS must contain barcode.", call. = FALSE)
  }
  out <- data.frame(row.names = as.character(score$barcode))
  for (trait in TRAITS_TO_RUN) {
    col <- paste0(trait, "_norm_score")
    if (!col %in% colnames(score)) {
      stop("Missing Granja RNA scDRS column: ", col, call. = FALSE)
    }
    out[[trait]] <- as.numeric(score[[col]])
  }
  saveRDS(out, BASELINE_RDS)
  fwrite(
    data.table(barcode = rownames(out), out),
    sub("\\.rds$", ".csv.gz", BASELINE_RDS)
  )
  BASELINE_RDS
}

patch_eval_script <- function(in_file, out_file) {
  txt <- readLines(in_file, warn = FALSE)
  txt <- sub(
    "^META_FILE\\s*<-.*$",
    paste0('META_FILE <- "', normalizePath(GRANJA_RNA_META_RDS, winslash = "/", mustWork = TRUE), '"'),
    txt
  )
  txt <- sub(
    "^BASELINE_TRS\\s*<-.*$",
    paste0('BASELINE_TRS <- "', normalizePath(BASELINE_RDS, winslash = "/", mustWork = TRUE), '"'),
    txt
  )
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
    make_target_replacement(TRAITS_TO_RUN),
    txt[pred_start:length(txt)]
  )

  pred_start <- grep("^PRED_TRAIT\\s*<-\\s*list\\(", txt)[1]
  boot_start <- grep("^BOOTSTRAP_CI\\s*<-", txt)[1]
  if (is.na(pred_start) || is.na(boot_start) || boot_start <= pred_start) {
    stop("Could not locate PRED_TRAIT block in eval script.", call. = FALSE)
  }
  txt <- c(
    txt[seq_len(pred_start - 1L)],
    make_pred_replacement(TRAITS_TO_RUN),
    txt[boot_start:length(txt)]
  )

  writeLines(txt, out_file, useBytes = TRUE)
  out_file
}

relabel_methods <- function(path) {
  out <- read.csv(path, check.names = FALSE)
  if ("method" %in% colnames(out)) {
    out$method[out$method == "baseline_atac"] <- "baseline_rna_scdrs"
    out$method[out$method == "v2"] <- "reverse_atac_to_rna"
  }
  write.csv(out, path, row.names = FALSE, quote = TRUE)
  out
}

cat(">>> Reverse Step 4: AUROC on Granja RNA cells\n\n")
cat("Traits      : ", paste(TRAITS_TO_RUN, collapse = ", "), "\n", sep = "")
cat("TRS dir     : ", V2_FILES_DIR_NEW, "\n", sep = "")
cat("Output dir  : ", OUT_DIR_NEW, "\n\n", sep = "")

stop_missing(c(GRANJA_RNA_META_RDS, GRANJA_RNA_SCORE_RDS, BASE_EVAL_SCRIPT))
assert_write_under_rev(c(OUT_DIR_NEW, BASELINE_RDS, PATCHED_SCRIPT))
dir.create(dirname(BASELINE_RDS), recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DIR_NEW, recursive = TRUE, showWarnings = FALSE)

trs_files <- list.files(
  V2_FILES_DIR_NEW,
  pattern = "^TRS_v2_[A-Z_]+_w[0-9.]+\\.tsv\\.gz$",
  full.names = TRUE
)
if (length(trs_files) == 0L) {
  stop("No TRS_v2_*.tsv.gz files found in: ", V2_FILES_DIR_NEW, call. = FALSE)
}

cat("[1/3] Writing RNA scDRS baseline matrix ...\n")
write_baseline_rds()
cat("  wrote: ", BASELINE_RDS, "\n", sep = "")

cat("[2/3] Writing patched evaluator ...\n")
patch_eval_script(BASE_EVAL_SCRIPT, PATCHED_SCRIPT)
cat("  wrote: ", PATCHED_SCRIPT, "\n", sep = "")

cat("[3/3] Running AUROC evaluator ...\n")
source(PATCHED_SCRIPT, local = new.env(parent = globalenv()))

result_csv <- file.path(OUT_DIR_NEW, OUT_CSV_NAME_NEW)
result <- relabel_methods(result_csv)
cat("\n=== Reverse AUROC finished ===\n")
cat("Result CSV: ", result_csv, "\n\n", sep = "")
print(result)
