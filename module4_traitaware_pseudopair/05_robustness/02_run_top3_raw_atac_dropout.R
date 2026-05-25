# module4 robustness: raw ATAC read drop-out for MONO/RDW/RBC
# -----------------------------------------------------------------------------
# Purpose:
#   Reuse module3's robustness implementation as much as possible. The script
#   below patches and runs:
#
#     D:/RStudio/work/module3_multiome/05_robustness/01_atac_sparsify.R
#
#   The patched module3 script does the same core operations:
#     raw peak-by-cell ATAC counts
#       -> binomial read drop-out
#       -> TF-IDF
#       -> LSI
#       -> weighted mutual KNN
#       -> RWR/TRS
#       -> AUROC
#
#   Only paths, traits, targets, and output directory are changed for module4.
#
# Inputs:
#   Raw Granja ATAC:
#     D:/SCAVENGE-reproducibility-main/blood_data/scATAC-Healthy-Hematopoiesis-191120.rds
#
#   Module4 top3 module3-compatible inputs prepared by:
#     05_robustness/01_prepare_top3_module3_inputs.R
#
# Run in RStudio:
#   setwd("D:/RStudio/work/module4_traitaware_pseudopair")
#   source("05_robustness/02_run_top3_raw_atac_dropout.R")
#
# Recommended first sanity run:
#   P_DROP_LIST <- c(0.0)
#   SEED_LIST <- c(1)
#   source("05_robustness/02_run_top3_raw_atac_dropout.R")
#
# Full run:
#   P_DROP_LIST <- seq(0.0, 0.9, by = 0.1)
#   SEED_LIST <- c(1, 2, 3)
#   source("05_robustness/02_run_top3_raw_atac_dropout.R")
# -----------------------------------------------------------------------------

if (!exists("ROOT")) ROOT <- "D:/RStudio/work/module4_traitaware_pseudopair"

if (!exists("MODULE3_ROBUST_SCRIPT")) {
  MODULE3_ROBUST_SCRIPT <- "D:/RStudio/work/module3_multiome/05_robustness/01_atac_sparsify.R"
}
if (!exists("RAW_ATAC_SE_FILE")) {
  RAW_ATAC_SE_FILE <- "D:/SCAVENGE-reproducibility-main/blood_data/scATAC-Healthy-Hematopoiesis-191120.rds"
}
if (!exists("META_FILE_NEW")) {
  META_FILE_NEW <- file.path(ROOT, "04_evaluation/inputs/granja_cell_meta.rds")
}

if (!exists("ROBUST_DIR")) {
  ROBUST_DIR <- file.path(ROOT, "05_robustness/output/top3_raw_atac_dropout")
}
if (!exists("BASE_INPUT_DIR")) {
  BASE_INPUT_DIR <- file.path(ROOT, "05_robustness/output/top3_embedding_stress/base_module3_inputs")
}
if (!exists("PATCHED_DIR")) {
  PATCHED_DIR <- file.path(ROBUST_DIR, "_patched_module3_scripts")
}
dir.create(PATCHED_DIR, recursive = TRUE, showWarnings = FALSE)

if (!exists("P_DROP_LIST")) P_DROP_LIST <- seq(0.0, 0.9, by = 0.1)
if (!exists("SEED_LIST")) SEED_LIST <- c(1L, 2L, 3L)
if (!exists("ATAC_W_LIST")) ATAC_W_LIST <- c(0.7, 1.0)
if (!exists("ROBUST_TRAITS")) ROBUST_TRAITS <- c("MONO", "RDW", "RBC")
if (!exists("ROBUST_PATCH_LABEL")) ROBUST_PATCH_LABEL <- "module4_top3"

TARGETS_ALL <- list(
  MONO = c("11_CD14.Mono.1", "12_CD14.Mono.2"),
  RDW = c("02_Early.Eryth", "03_Late.Eryth", "04_Early.Baso"),
  RBC = c("02_Early.Eryth", "03_Late.Eryth", "04_Early.Baso")
)

stop_missing <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0L) {
    stop("Missing required file(s):\n  ", paste(missing, collapse = "\n  "),
         call. = FALSE)
  }
}

r_literal_chr <- function(x) {
  paste(sprintf('"%s"', x), collapse = ", ")
}

r_literal_num <- function(x) {
  paste(format(x, trim = TRUE, scientific = FALSE), collapse = ", ")
}

r_named_target_block <- function(traits) {
  bad <- setdiff(traits, names(TARGETS_ALL))
  if (length(bad) > 0L) {
    stop("No target definition for trait(s): ", paste(bad, collapse = ", "),
         call. = FALSE)
  }
  lines <- character(length(traits))
  for (i in seq_along(traits)) {
    trait <- traits[[i]]
    comma <- if (i < length(traits)) "," else ""
    lines[[i]] <- paste0("  ", trait, " = c(", r_literal_chr(TARGETS_ALL[[trait]]), ")", comma)
  }
  c("TARGETS <- list(", lines, ")")
}

r_named_pred_trait <- function(traits) {
  entries <- paste0(traits, ' = "', traits, '"')
  paste0("PRED_TRAIT <- list(", paste(entries, collapse = ", "), ")")
}

patch_module3_robustness <- function(in_file, out_file) {
  txt <- readLines(in_file, warn = FALSE)

  txt <- sub(
    "^P_DROP_LIST\\s*<-.*$",
    paste0("P_DROP_LIST <- c(", r_literal_num(P_DROP_LIST), ")"),
    txt
  )
  txt <- sub(
    "^SEED_LIST\\s*<-.*$",
    paste0("SEED_LIST <- c(", r_literal_num(SEED_LIST), ")"),
    txt
  )
  txt <- sub(
    "^ATAC_W_LIST\\s*<-.*$",
    paste0("ATAC_W_LIST <- c(", r_literal_num(ATAC_W_LIST), ")"),
    txt
  )
  txt <- sub(
    "^TRAITS\\s*<-.*$",
    paste0("TRAITS <- c(", r_literal_chr(ROBUST_TRAITS), ")"),
    txt
  )
  txt <- sub(
    "^SE_FILE\\s*<-.*$",
    paste0('SE_FILE <- "', normalizePath(RAW_ATAC_SE_FILE, winslash = "/", mustWork = TRUE), '"'),
    txt
  )
  txt <- sub(
    "^DUAL_DIR\\s*<-.*$",
    paste0('DUAL_DIR <- "', normalizePath(BASE_INPUT_DIR, winslash = "/", mustWork = TRUE), '"'),
    txt
  )
  txt <- sub(
    "^META_FILE\\s*<-.*$",
    paste0('META_FILE <- "', normalizePath(META_FILE_NEW, winslash = "/", mustWork = TRUE), '"'),
    txt
  )
  txt <- sub(
    "^OUT_DIR\\s*<-.*$",
    paste0('OUT_DIR <- "', normalizePath(ROBUST_DIR, winslash = "/", mustWork = FALSE), '"'),
    txt
  )

  target_start <- grep("^TARGETS\\s*<-\\s*list\\(", txt)[1]
  pred_start <- grep("^PRED_TRAIT\\s*<-\\s*list\\(", txt)[1]
  if (is.na(target_start) || is.na(pred_start) || pred_start <= target_start) {
    stop("Could not locate TARGETS/PRED_TRAIT blocks in module3 robustness script.",
         call. = FALSE)
  }
  target_block <- r_named_target_block(ROBUST_TRAITS)
  txt <- c(
    txt[seq_len(target_start - 1L)],
    target_block,
    txt[pred_start:length(txt)]
  )

  pred_line <- grep("^PRED_TRAIT\\s*<-", txt)[1]
  if (is.na(pred_line)) {
    stop("Could not locate PRED_TRAIT in patched module3 robustness script.",
         call. = FALSE)
  }
  pred_after <- seq.int(pred_line + 1L, length(txt))
  pred_stop <- pred_line
  if (length(pred_after) > 0L) {
    pred_end_markers <- pred_after[grepl("^\\s*$|^#\\s*=+", txt[pred_after])]
    if (length(pred_end_markers) > 0L) pred_stop <- pred_end_markers[1L] - 1L
  }
  pred_block <- r_named_pred_trait(ROBUST_TRAITS)
  txt <- c(
    if (pred_line > 1L) txt[seq_len(pred_line - 1L)] else character(),
    pred_block,
    if (pred_stop < length(txt)) txt[(pred_stop + 1L):length(txt)] else character()
  )

  # Keep module3's robustness pipeline, but align AUROC direction with module4
  # evaluation so low-scoring positives are not silently flipped to high AUROC.
  auc_start <- grep("^compute_auroc\\s*<-\\s*function", txt)[1]
  zero_start <- grep("^# 0\\. meta", txt)[1]
  if (is.na(auc_start) || is.na(zero_start) || zero_start <= auc_start) {
    stop("Could not locate compute_auroc block in module3 robustness script.",
         call. = FALSE)
  }
  auc_block <- c(
    "compute_auroc <- function(predictor, labels) {",
    "  keep <- is.finite(predictor) & !is.na(labels)",
    "  predictor <- predictor[keep]",
    "  labels <- labels[keep]",
    "  if (length(predictor) < 20L) return(NA_real_)",
    "  if (sum(labels == 1L) < 10L || sum(labels == 0L) < 10L) return(NA_real_)",
    "  v <- stats::var(predictor)",
    "  if (!is.finite(v) || v < 1e-12) return(NA_real_)",
    "  roc_obj <- suppressMessages(pROC::roc(labels, predictor, quiet = TRUE,",
    "                                         direction = \"<\", levels = c(0L, 1L)))",
    "  as.numeric(pROC::auc(roc_obj))",
    "}",
    ""
  )
  txt <- c(
    txt[seq_len(auc_start - 1L)],
    auc_block,
    txt[zero_start:length(txt)]
  )

  sanity_start <- grep("^# 4\\. sanity headline", txt)[1]
  if (!is.na(sanity_start)) {
    sanity_block <- c(
      "# ============================================================================",
      "# 4. module4 headline at p_drop = 0",
      "# ============================================================================",
      "cat(\"\\n\", strrep(\"=\", 70), \"\\n\", sep = \"\")",
      "cat(\"MODULE4 RAW ATAC DROP-OUT HEADLINE (p_drop = 0)\\n\")",
      "cat(strrep(\"=\", 70), \"\\n\", sep = \"\")",
      "sanity <- eval_results[eval_results$p_drop == 0, ]",
      "if (nrow(sanity) > 0L) {",
      "  print(sanity[order(sanity$trait, sanity$target, sanity$seed, sanity$atac_w), ])",
      "} else {",
      "  cat(\"No p_drop=0 rows were generated in this run.\\n\")",
      "}",
      "",
      "cat(\"\\nmodule3 01_atac_sparsify.R patched for module4 raw ATAC drop-out OK.\\n\")"
    )
    txt <- c(txt[seq_len(sanity_start - 1L)], sanity_block)
  }

  writeLines(txt, out_file, useBytes = TRUE)
  out_file
}

cat(">>> Run module4 top3 raw ATAC read drop-out robustness\n\n")
cat("Traits     : ", paste(ROBUST_TRAITS, collapse = ", "), "\n", sep = "")
cat("P_DROP_LIST: ", paste(P_DROP_LIST, collapse = ", "), "\n", sep = "")
cat("SEED_LIST  : ", paste(SEED_LIST, collapse = ", "), "\n", sep = "")
cat("ATAC_W_LIST: ", paste(ATAC_W_LIST, collapse = ", "), "\n", sep = "")
cat("Raw ATAC SE: ", RAW_ATAC_SE_FILE, "\n", sep = "")
cat("Module4 inputs: ", BASE_INPUT_DIR, "\n", sep = "")
cat("Output dir: ", ROBUST_DIR, "\n\n", sep = "")

required_inputs <- c(
  MODULE3_ROBUST_SCRIPT,
  RAW_ATAC_SE_FILE,
  META_FILE_NEW,
  file.path(BASE_INPUT_DIR, "atac_z_block_v2.rds"),
  file.path(BASE_INPUT_DIR, paste0("rna_z_block_v2_", ROBUST_TRAITS, ".rds")),
  file.path(BASE_INPUT_DIR, paste0("seed_v2_", ROBUST_TRAITS, ".rds"))
)
stop_missing(required_inputs)

dir.create(ROBUST_DIR, recursive = TRUE, showWarnings = FALSE)
patched <- file.path(PATCHED_DIR, paste0("module3_01_atac_sparsify.", ROBUST_PATCH_LABEL, "_patched.R"))
patch_module3_robustness(MODULE3_ROBUST_SCRIPT, patched)
cat("Patched module3 robustness script:\n  ", patched, "\n\n", sep = "")

source(patched, local = new.env(parent = globalenv()))

cat("\n=== module4 top3 raw ATAC read drop-out robustness finished ===\n")
cat("Outputs are in: ", ROBUST_DIR, "\n", sep = "")
