# module4 robustness: v3 sparse MONO raw ATAC drop-out
# -----------------------------------------------------------------------------
# Core computation is reused from module3:
#   D:/RStudio/work/module3_multiome/05_robustness/01_atac_sparsify.R
#
# This is the single entry point for the v3 sparse MONO robustness figure.
# It runs the degraded-data sweep for w=0.7 and w=1.0, recomputes the
# dynamic ATAC-only baseline on the same degraded ATAC data, then draws the
# three-line figure.
#
# Run in RStudio:
#   setwd("D:/RStudio/work/module4_traitaware_pseudopair")
#   source("05_robustness/07_run_v3_sparse_mono_raw_atac_dropout.R")
#
# Optional partial resume:
#   ROBUST_RESUME_STEPS <- c("dynamic_baseline", "plot")
#   source("05_robustness/07_run_v3_sparse_mono_raw_atac_dropout.R")
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

ROOT <- "D:/RStudio/work/module4_traitaware_pseudopair"
MODULE3_ATAC_SPARSIFY <- "D:/RStudio/work/module3_multiome/05_robustness/01_atac_sparsify.R"

V3_TAG <- "v3_sparse_mapping_top5significant"
TRAITS <- c("MONO")
P_DROP_LIST <- seq(0.0, 0.9, by = 0.1)
SEED_LIST <- c(1L, 2L, 3L)
ATAC_W_LIST <- c(0.7, 1.0)

VALID_ROBUST_RUN_STEPS <- c(
  "prepare_inputs",
  "patch_module3",
  "run_module3",
  "dynamic_baseline",
  "plot"
)
if (exists("ROBUST_RESUME_STEPS")) {
  ROBUST_RUN_STEPS <- ROBUST_RESUME_STEPS
} else {
  ROBUST_RUN_STEPS <- VALID_ROBUST_RUN_STEPS
}
bad_steps <- setdiff(ROBUST_RUN_STEPS, VALID_ROBUST_RUN_STEPS)
if (length(bad_steps) > 0L) {
  stop(
    "Unknown ROBUST_RUN_STEPS: ", paste(bad_steps, collapse = ", "),
    "\nValid steps: ", paste(VALID_ROBUST_RUN_STEPS, collapse = ", "),
    call. = FALSE
  )
}

RAW_ATAC_SE_FILE <- "D:/SCAVENGE-reproducibility-main/blood_data/scATAC-Healthy-Hematopoiesis-191120.rds"
BG_PEAKS_FILE <- "D:/SCAVENGE-reproducibility-main/blood_data/scavenge_out_small_all_traits/bg_peaks_n200.rds"
GWAS_FILE <- "D:/SCAVENGE-reproducibility-main/blood_data/CS.trans_all.PP01.strong_atac_overlap.tsv"
META_FILE <- file.path(ROOT, "04_evaluation/inputs/granja_cell_meta.rds")
MODULE3_INPUT_DIR <- file.path(
  ROOT,
  "03_pseudopair_pipeline/output",
  V3_TAG,
  "module3_mono_neighborhood_fdr_ge0p80"
)

OUT_DIR <- file.path(ROOT, "05_robustness/output/v3_sparse_mono_raw_atac_dropout")
DUAL_DIR <- file.path(OUT_DIR, "module3_inputs")
PATCHED_DIR <- file.path(OUT_DIR, "_patched_module3_scripts")
dir.create(DUAL_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PATCHED_DIR, recursive = TRUE, showWarnings = FALSE)

TARGETS <- list(
  MONO = c("11_CD14.Mono.1", "12_CD14.Mono.2")
)
PRED_TRAIT <- list(MONO = "MONO")

stop_missing <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0L) {
    stop("Missing required file(s):\n  ", paste(missing, collapse = "\n  "),
         call. = FALSE)
  }
}

copy_required <- function(from, to) {
  stop_missing(from)
  ok <- file.copy(from, to, overwrite = TRUE)
  if (!ok) stop("Failed to copy:\n  ", from, "\n  -> ", to, call. = FALSE)
}

r_chr <- function(x) paste(sprintf('"%s"', x), collapse = ", ")
r_num <- function(x) paste(format(x, trim = TRUE, scientific = FALSE), collapse = ", ")

replace_line <- function(txt, pattern, replacement) {
  hit <- grep(pattern, txt)
  if (length(hit) == 0L) stop("Pattern not found: ", pattern, call. = FALSE)
  txt[hit[1]] <- replacement
  txt
}

replace_block <- function(txt, start_pattern, end_pattern, replacement) {
  start <- grep(start_pattern, txt)[1]
  if (is.na(start)) {
    stop("Could not replace block: ", start_pattern, " -> ", end_pattern,
         call. = FALSE)
  }
  end_candidates <- grep(end_pattern, txt)
  end_candidates <- end_candidates[end_candidates > start]
  end <- end_candidates[1]
  if (is.na(end)) {
    stop("Could not replace block: ", start_pattern, " -> ", end_pattern,
         call. = FALSE)
  }
  c(txt[seq_len(start - 1L)], replacement, txt[end:length(txt)])
}

patch_module3_script <- function(in_file, out_file) {
  txt <- readLines(in_file, warn = FALSE)

  txt <- replace_line(txt, "^P_DROP_LIST\\s*<-", paste0("P_DROP_LIST <- c(", r_num(P_DROP_LIST), ")"))
  txt <- replace_line(txt, "^SEED_LIST\\s*<-", paste0("SEED_LIST <- c(", r_num(SEED_LIST), ")"))
  txt <- replace_line(txt, "^ATAC_W_LIST\\s*<-", paste0("ATAC_W_LIST <- c(", r_num(ATAC_W_LIST), ")"))
  txt <- replace_line(txt, "^TRAITS\\s*<-", paste0("TRAITS <- c(", r_chr(TRAITS), ")"))
  txt <- replace_line(txt, "^SE_FILE\\s*<-", paste0('SE_FILE <- "', normalizePath(RAW_ATAC_SE_FILE, winslash = "/", mustWork = TRUE), '"'))
  txt <- replace_line(txt, "^DUAL_DIR\\s*<-", paste0('DUAL_DIR <- "', normalizePath(DUAL_DIR, winslash = "/", mustWork = TRUE), '"'))
  txt <- replace_line(txt, "^META_FILE\\s*<-", paste0('META_FILE <- "', normalizePath(META_FILE, winslash = "/", mustWork = TRUE), '"'))
  txt <- replace_line(txt, "^OUT_DIR\\s*<-", paste0('OUT_DIR <- "', normalizePath(OUT_DIR, winslash = "/", mustWork = FALSE), '"'))

  target_block <- c(
    "TARGETS <- list(",
    paste0("  MONO = c(", r_chr(TARGETS$MONO), ")"),
    ")"
  )
  txt <- replace_block(txt, "^TARGETS\\s*<-\\s*list\\(", "^PRED_TRAIT\\s*<-\\s*list\\(", target_block)
  txt <- replace_block(
    txt,
    "^PRED_TRAIT\\s*<-\\s*list\\(",
    "^# ============================================================================$",
    'PRED_TRAIT <- list(MONO = "MONO")'
  )

  auc_block <- c(
    "compute_auroc <- function(predictor, labels) {",
    "  keep <- is.finite(predictor) & !is.na(labels)",
    "  predictor <- predictor[keep]",
    "  labels <- labels[keep]",
    "  if (length(predictor) < 20L) return(NA_real_)",
    "  if (sum(labels == 1L) < 10L || sum(labels == 0L) < 10L) return(NA_real_)",
    "  if (!is.finite(stats::var(predictor)) || stats::var(predictor) < 1e-12) return(NA_real_)",
    "  roc_obj <- suppressMessages(pROC::roc(labels, predictor, quiet = TRUE,",
    "                                         direction = \"<\", levels = c(0L, 1L)))",
    "  as.numeric(pROC::auc(roc_obj))",
    "}",
    ""
  )
  txt <- replace_block(txt, "^compute_auroc\\s*<-\\s*function", "^# 0\\. meta", auc_block)

  headline_start <- grep("^# 4\\. sanity headline", txt)[1]
  if (!is.na(headline_start)) {
    txt <- c(
      txt[seq_len(headline_start - 1L)],
      "cat(\"\\nmodule3 01_atac_sparsify.R reused for module4 v3 sparse MONO OK.\\n\")"
    )
  }

  writeLines(txt, out_file, useBytes = TRUE)
  out_file
}

prepare_inputs <- function() {
  copy_required(
    file.path(MODULE3_INPUT_DIR, "atac_z_block_v2.rds"),
    file.path(DUAL_DIR, "atac_z_block_v2.rds")
  )
  for (trait in TRAITS) {
    copy_required(
      file.path(MODULE3_INPUT_DIR, paste0("rna_z_block_v2_", trait, ".rds")),
      file.path(DUAL_DIR, paste0("rna_z_block_v2_", trait, ".rds"))
    )
    copy_required(
      file.path(MODULE3_INPUT_DIR, paste0("seed_v2_", trait, ".rds")),
      file.path(DUAL_DIR, paste0("seed_v2_", trait, ".rds"))
    )
  }

  atac <- readRDS(file.path(DUAL_DIR, "atac_z_block_v2.rds"))
  rows <- lapply(TRAITS, function(trait) {
    rna <- readRDS(file.path(DUAL_DIR, paste0("rna_z_block_v2_", trait, ".rds")))
    seed <- readRDS(file.path(DUAL_DIR, paste0("seed_v2_", trait, ".rds")))
    if (!identical(rownames(atac), rownames(rna))) {
      stop("ATAC/RNA rownames differ for ", trait, call. = FALSE)
    }
    if (!identical(rownames(atac), names(seed))) {
      stop("ATAC/seed names differ for ", trait, call. = FALSE)
    }
    data.table(
      trait = trait,
      n_cells = nrow(atac),
      atac_dim = ncol(atac),
      rna_dim = ncol(rna),
      n_seed_positive = sum(seed > 0),
      seed_min_positive = min(seed[seed > 0]),
      seed_median_positive = median(seed[seed > 0]),
      seed_max = max(seed)
    )
  })
  manifest <- rbindlist(rows)
  fwrite(manifest, file.path(DUAL_DIR, "input_manifest.csv"))
  manifest
}

step_enabled <- function(x) x %in% ROBUST_RUN_STEPS

cat(">>> v3 sparse MONO raw ATAC drop-out robustness\n\n")
if (exists("RUN_STEPS")) {
  cat("Ignoring RUN_STEPS from current R session.\n")
}
if (exists("ROBUST_RESUME_STEPS")) {
  cat("Using ROBUST_RESUME_STEPS for partial resume.\n")
} else if (exists("ROBUST_RUN_STEPS")) {
  cat("Using default full workflow; set ROBUST_RESUME_STEPS only for partial resume.\n")
}
cat("ROBUST_RUN_STEPS: ", paste(ROBUST_RUN_STEPS, collapse = ", "), "\n", sep = "")
cat("P_DROP_LIST: ", paste(P_DROP_LIST, collapse = ", "), "\n", sep = "")
cat("SEED_LIST: ", paste(SEED_LIST, collapse = ", "), "\n", sep = "")
cat("ATAC_W_LIST: ", paste(ATAC_W_LIST, collapse = ", "), "\n", sep = "")
cat("Baseline: dynamic ATAC-only baseline recomputed at each p_drop/seed\n")
cat("Output: ", OUT_DIR, "\n\n", sep = "")

stop_missing(c(MODULE3_ATAC_SPARSIFY, RAW_ATAC_SE_FILE, META_FILE, MODULE3_INPUT_DIR))

if (step_enabled("prepare_inputs")) {
  manifest <- prepare_inputs()
  print(manifest)
}

patched <- file.path(PATCHED_DIR, "01_atac_sparsify.module4_v3_sparse_mono.R")
if (step_enabled("patch_module3")) {
  patch_module3_script(MODULE3_ATAC_SPARSIFY, patched)
  cat("Patched module3 script:\n  ", patched, "\n", sep = "")
}

if (step_enabled("run_module3")) {
  stop_missing(patched)
  source(patched, local = new.env(parent = globalenv()))
}

if (step_enabled("dynamic_baseline")) {
  ROBUST_DIR <- OUT_DIR
  BASELINE_DIR <- file.path(OUT_DIR, "dynamic_baseline")
  source(file.path(ROOT, "05_robustness/05_run_top3_raw_atac_dropout_dynamic_baseline.R"))
}

if (step_enabled("plot")) {
  source(file.path(ROOT, "05_robustness/08_plot_v3_sparse_mono_raw_atac_dropout.R"))
}

cat("\n=== v3 sparse MONO raw ATAC drop-out driver finished ===\n")
cat("Output: ", OUT_DIR, "\n", sep = "")
