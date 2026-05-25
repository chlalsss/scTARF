# module4 robustness: current-code ATAC feature drop-out curves
# -----------------------------------------------------------------------------
# Purpose:
#   Generate plot-ready AUROC data for the current v3 pipeline only.
#
# Important boundary:
#   This script starts from the exact current clean module3-compatible directory:
#     atac_z_block_v2.rds
#     rna_z_block_v2_<trait>.rds
#     seed_v2_<trait>.rds
#     clean knn/TRS files
#
#   Therefore p_drop = 0 reuses the clean TRS outputs directly and must match
#   the clean AUROC CSV. For p_drop > 0, only the current ATAC feature block is
#   degraded. No raw-count LSI is recomputed here.
#
# Outputs:
#   05_robustness/output/v3_current_code_atac_dropout/
#     v3_current_code_atac_dropout_auroc_all.csv
#     v3_current_code_atac_dropout_plot_data.csv
#     v3_current_code_atac_dropout_run_manifest.csv
#
# Run in RStudio:
#   setwd("D:/RStudio/work/module4_traitaware_pseudopair")
#   source("05_robustness/09_run_v3_current_code_atac_dropout.R")
#
# Future trait example:
#   ROBUST_TRAITS_TO_RUN <- c("RDW", "RBC")
#   ROBUST_BASE_COMPAT_DIR <- "path/to/current/module3-compatible/dir"
#   ROBUST_CLEAN_AUROC_CSV <- "path/to/current/clean_auroc.csv"
#   ROBUST_TARGETS <- list(
#     RDW = c("02_Early.Eryth", "03_Late.Eryth", "04_Early.Baso"),
#     RBC = c("02_Early.Eryth", "03_Late.Eryth", "04_Early.Baso")
#   )
#   ROBUST_PRED_TRAIT <- list(RDW = "RDW", RBC = "RBC")
#   source("05_robustness/09_run_v3_current_code_atac_dropout.R")
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

ROOT <- "D:/RStudio/work/module4_traitaware_pseudopair"

if (!exists("ROBUST_OUTPUT_TAG")) ROBUST_OUTPUT_TAG <- "v3_current_code_atac_dropout"

if (exists("ROBUST_BASE_COMPAT_DIR")) {
  BASE_COMPAT_DIR <- ROBUST_BASE_COMPAT_DIR
} else {
  BASE_COMPAT_DIR <- file.path(
    ROOT,
    "03_pseudopair_pipeline/output/v3_sparse_mapping_top5significant",
    "module3_mono_neighborhood_fdr_ge0p80"
  )
}
if (exists("ROBUST_CLEAN_AUROC_CSV")) {
  CLEAN_AUROC_CSV <- ROBUST_CLEAN_AUROC_CSV
} else {
  CLEAN_AUROC_CSV <- file.path(
    ROOT,
    "04_evaluation/v3_sparse_mapping_top5significant/mono_neighborhood_fdr_ge0p80",
    "v3_sparse_mapping_top5significant_mono_neighborhood_fdr_ge0p80_auroc.csv"
  )
}

if (exists("ROBUST_TRAITS_TO_RUN")) {
  TRAITS_TO_RUN <- ROBUST_TRAITS_TO_RUN
} else {
  TRAITS_TO_RUN <- c("MONO")
}
if (exists("ROBUST_TARGETS")) {
  TARGETS <- ROBUST_TARGETS
} else {
  TARGETS <- list(MONO = c("11_CD14.Mono.1", "12_CD14.Mono.2"))
}
if (exists("ROBUST_PRED_TRAIT")) {
  PRED_TRAIT <- ROBUST_PRED_TRAIT
} else {
  PRED_TRAIT <- list(MONO = "MONO")
}
if (exists("ROBUST_P_DROP_LIST")) {
  P_DROP_LIST <- ROBUST_P_DROP_LIST
} else {
  P_DROP_LIST <- round(seq(0.0, 0.9, by = 0.1), 2)
}
if (exists("ROBUST_SEED_LIST")) {
  SEED_LIST <- ROBUST_SEED_LIST
} else {
  SEED_LIST <- c(1L, 2L, 3L)
}
if (exists("ROBUST_ATAC_W_LIST")) {
  ATAC_W_LIST <- ROBUST_ATAC_W_LIST
} else {
  ATAC_W_LIST <- c(0.7, 1.0)
}
if (exists("ROBUST_BOOTSTRAP_CI")) {
  BOOTSTRAP_CI_ROBUST <- ROBUST_BOOTSTRAP_CI
} else {
  BOOTSTRAP_CI_ROBUST <- FALSE
}
if (exists("ROBUST_ENFORCE_P0_CLEAN_MATCH")) {
  ENFORCE_P0_CLEAN_MATCH <- ROBUST_ENFORCE_P0_CLEAN_MATCH
} else {
  ENFORCE_P0_CLEAN_MATCH <- TRUE
}

MODULE3_RUN_SCRIPT <- file.path(ROOT, "03_pseudopair_pipeline/03_run_module3_v2_knn_rwr.R")
BASE_EVAL_SCRIPT <- file.path(ROOT, "04_evaluation/eval_auroc.R")

ROBUST_DIR <- file.path(ROOT, "05_robustness/output", ROBUST_OUTPUT_TAG)
RUNS_DIR <- file.path(ROBUST_DIR, "runs")
EVAL_DIR <- file.path(ROBUST_DIR, "eval")
dir.create(RUNS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(EVAL_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_ALL <- file.path(ROBUST_DIR, paste0(ROBUST_OUTPUT_TAG, "_auroc_all.csv"))
PLOT_DATA_FILE <- file.path(ROBUST_DIR, paste0(ROBUST_OUTPUT_TAG, "_plot_data.csv"))
RUN_MANIFEST_FILE <- file.path(ROBUST_DIR, paste0(ROBUST_OUTPUT_TAG, "_run_manifest.csv"))

stop_missing <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0L) {
    stop("Missing required file(s):\n  ", paste(missing, collapse = "\n  "),
         call. = FALSE)
  }
}

copy_required <- function(from, to, overwrite = TRUE) {
  stop_missing(from)
  dir.create(dirname(to), recursive = TRUE, showWarnings = FALSE)
  ok <- file.copy(from, to, overwrite = overwrite)
  if (!ok) stop("Failed to copy:\n  ", from, "\n  -> ", to, call. = FALSE)
}

r_literal_chr <- function(x) {
  paste(sprintf('"%s"', x), collapse = ", ")
}

r_literal_num <- function(x) {
  paste(format(x, trim = TRUE, scientific = FALSE), collapse = ", ")
}

target_block_lines <- function(targets) {
  nms <- names(targets)
  lines <- "TARGETS <- list("
  for (i in seq_along(nms)) {
    comma <- if (i < length(nms)) "," else ""
    lines <- c(
      lines,
      paste0("  ", nms[[i]], " = c(", r_literal_chr(targets[[i]]), ")", comma)
    )
  }
  c(lines, ")")
}

pred_trait_line <- function(pred_trait) {
  entries <- paste(
    sprintf('%s = "%s"', names(pred_trait), unlist(pred_trait, use.names = FALSE)),
    collapse = ", "
  )
  paste0("PRED_TRAIT <- list(", entries, ")")
}

replace_block <- function(txt, start_pattern, end_pattern, replacement) {
  start <- grep(start_pattern, txt)[1]
  if (is.na(start)) stop("Could not find block start: ", start_pattern, call. = FALSE)
  end_candidates <- grep(end_pattern, txt)
  end_candidates <- end_candidates[end_candidates > start]
  end <- end_candidates[1]
  if (is.na(end)) stop("Could not find block end: ", end_pattern, call. = FALSE)
  c(txt[seq_len(start - 1L)], replacement, txt[end:length(txt)])
}

patch_eval_script <- function(in_file, out_file, v2_dir, out_dir, out_csv_name) {
  txt <- readLines(in_file, warn = FALSE)

  no_baseline_file <- file.path(out_dir, "no_baseline_for_current_code_dropout.rds")
  txt <- sub(
    "^BASELINE_TRS\\s*<-.*$",
    paste0('BASELINE_TRS <- "', normalizePath(no_baseline_file, winslash = "/", mustWork = FALSE), '"'),
    txt
  )
  txt <- sub(
    "^V2_FILES_DIR\\s*<-.*$",
    paste0('V2_FILES_DIR <- "', normalizePath(v2_dir, winslash = "/", mustWork = TRUE), '"'),
    txt
  )
  txt <- sub(
    "^OUT_DIR\\s*<-.*$",
    paste0('OUT_DIR <- "', normalizePath(out_dir, winslash = "/", mustWork = TRUE), '"'),
    txt
  )
  txt <- sub("^OUT_CSV_NAME\\s*<-.*$", paste0('OUT_CSV_NAME <- "', out_csv_name, '"'), txt)
  txt <- sub("^BOOTSTRAP_CI\\s*<-.*$", paste0("BOOTSTRAP_CI <- ", ifelse(BOOTSTRAP_CI_ROBUST, "TRUE", "FALSE")), txt)

  txt <- replace_block(
    txt,
    "^TARGETS\\s*<-\\s*list\\(",
    "^PRED_TRAIT\\s*<-\\s*list\\(",
    target_block_lines(TARGETS)
  )
  txt <- replace_block(
    txt,
    "^PRED_TRAIT\\s*<-\\s*list\\(",
    "^BOOTSTRAP_CI\\s*<-",
    c(pred_trait_line(PRED_TRAIT), "")
  )

  writeLines(txt, out_file, useBytes = TRUE)
  out_file
}

degrade_atac_block <- function(atac_z, p_drop, seed) {
  if (p_drop <= 0) return(atac_z)
  set.seed(seed)
  out <- as.matrix(atac_z)
  mask <- matrix(
    stats::runif(length(out)) < p_drop,
    nrow = nrow(out),
    ncol = ncol(out)
  )
  out[mask] <- 0
  rownames(out) <- rownames(atac_z)
  colnames(out) <- colnames(atac_z)
  out
}

format_p <- function(p_drop) sprintf("%.2f", p_drop)

run_id_for <- function(p_drop, seed) {
  paste0("p", format_p(p_drop), "_seed", seed)
}

prepare_run_dir <- function(p_drop, seed) {
  run_id <- run_id_for(p_drop, seed)
  run_dir <- file.path(RUNS_DIR, run_id)
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

  copy_required(
    file.path(BASE_COMPAT_DIR, "atac_z_block_v2.rds"),
    file.path(run_dir, "atac_z_block_v2.rds")
  )
  for (trait in TRAITS_TO_RUN) {
    copy_required(
      file.path(BASE_COMPAT_DIR, paste0("rna_z_block_v2_", trait, ".rds")),
      file.path(run_dir, paste0("rna_z_block_v2_", trait, ".rds"))
    )
    copy_required(
      file.path(BASE_COMPAT_DIR, paste0("seed_v2_", trait, ".rds")),
      file.path(run_dir, paste0("seed_v2_", trait, ".rds"))
    )
  }

  if (p_drop == 0) {
    for (trait in TRAITS_TO_RUN) {
      for (w in ATAC_W_LIST) {
        copy_required(
          file.path(BASE_COMPAT_DIR, sprintf("knn_v2_%s_w%.1f.rds", trait, w)),
          file.path(run_dir, sprintf("knn_v2_%s_w%.1f.rds", trait, w))
        )
        copy_required(
          file.path(BASE_COMPAT_DIR, sprintf("TRS_v2_%s_w%.1f.tsv.gz", trait, w)),
          file.path(run_dir, sprintf("TRS_v2_%s_w%.1f.tsv.gz", trait, w))
        )
      }
    }
  } else {
    base_atac <- readRDS(file.path(BASE_COMPAT_DIR, "atac_z_block_v2.rds"))
    stressed_atac <- degrade_atac_block(base_atac, p_drop = p_drop, seed = seed)
    saveRDS(stressed_atac, file.path(run_dir, "atac_z_block_v2.rds"))
  }

  data.table(
    run_id = run_id,
    p_drop = p_drop,
    stress_seed = seed,
    run_dir = run_dir,
    base_compat_dir = BASE_COMPAT_DIR,
    n_traits = length(TRAITS_TO_RUN),
    atac_dropout_unit = "ATAC z-block entries",
    p0_reuses_clean_trs = p_drop == 0
  )
}

run_module3_if_needed <- function(run_dir, p_drop) {
  if (p_drop == 0) {
    cat("  p_drop=0: reusing clean KNN/TRS from BASE_COMPAT_DIR\n")
    return(invisible(NULL))
  }

  e <- new.env(parent = globalenv())
  assign("MODULE3_COMPAT_OUT_DIR", run_dir, envir = e)
  assign("TRAITS_TO_RUN", TRAITS_TO_RUN, envir = e)
  assign("ATAC_W_LIST", ATAC_W_LIST, envir = e)
  source(MODULE3_RUN_SCRIPT, local = e)
  invisible(NULL)
}

evaluate_one_run <- function(run_id, run_dir, p_drop, seed) {
  run_eval_dir <- file.path(EVAL_DIR, run_id)
  dir.create(run_eval_dir, recursive = TRUE, showWarnings = FALSE)
  out_csv_name <- paste0(run_id, "_auroc.csv")
  patched_eval <- file.path(run_eval_dir, paste0("eval_auroc.", run_id, ".patched.R"))

  patch_eval_script(
    in_file = BASE_EVAL_SCRIPT,
    out_file = patched_eval,
    v2_dir = run_dir,
    out_dir = run_eval_dir,
    out_csv_name = out_csv_name
  )

  source(patched_eval, local = new.env(parent = globalenv()))

  one <- fread(file.path(run_eval_dir, out_csv_name))
  one <- one[method != "baseline_atac"]
  one[, `:=`(
    run_id = run_id,
    p_drop = p_drop,
    stress_seed = seed
  )]
  one
}

cat(">>> v3 current-code ATAC drop-out robustness data\n\n")
cat("BASE_COMPAT_DIR : ", BASE_COMPAT_DIR, "\n", sep = "")
cat("TRAITS_TO_RUN   : ", paste(TRAITS_TO_RUN, collapse = ", "), "\n", sep = "")
cat("ATAC_W_LIST     : ", paste(ATAC_W_LIST, collapse = ", "), "\n", sep = "")
cat("P_DROP_LIST     : ", paste(P_DROP_LIST, collapse = ", "), "\n", sep = "")
cat("SEED_LIST       : ", paste(SEED_LIST, collapse = ", "), "\n", sep = "")
cat("Output          : ", ROBUST_DIR, "\n\n", sep = "")

stop_missing(c(
  MODULE3_RUN_SCRIPT,
  BASE_EVAL_SCRIPT,
  file.path(BASE_COMPAT_DIR, "atac_z_block_v2.rds"),
  file.path(BASE_COMPAT_DIR, paste0("rna_z_block_v2_", TRAITS_TO_RUN, ".rds")),
  file.path(BASE_COMPAT_DIR, paste0("seed_v2_", TRAITS_TO_RUN, ".rds"))
))

run_rows <- list()
eval_rows <- list()
idx <- 0L

for (p_drop in P_DROP_LIST) {
  for (seed in SEED_LIST) {
    run_id <- run_id_for(p_drop, seed)
    cat("\n", strrep("=", 72), "\n", sep = "")
    cat("Run: ", run_id, "\n", sep = "")
    cat(strrep("=", 72), "\n", sep = "")

    manifest_row <- prepare_run_dir(p_drop, seed)
    run_rows[[run_id]] <- manifest_row

    run_module3_if_needed(manifest_row$run_dir, p_drop)

    idx <- idx + 1L
    eval_rows[[idx]] <- evaluate_one_run(
      run_id = run_id,
      run_dir = manifest_row$run_dir,
      p_drop = p_drop,
      seed = seed
    )
  }
}

run_manifest <- rbindlist(run_rows, use.names = TRUE, fill = TRUE)
fwrite(run_manifest, RUN_MANIFEST_FILE)

all_dt <- rbindlist(eval_rows, use.names = TRUE, fill = TRUE)
setcolorder(all_dt, c("run_id", "p_drop", "stress_seed"))
setorder(all_dt, target, atac_w, p_drop, stress_seed)
fwrite(all_dt, OUT_ALL)

summary_dt <- all_dt[, .(
  auroc = mean(auroc, na.rm = TRUE),
  sd_auroc = stats::sd(auroc, na.rm = TRUE),
  n_runs = .N,
  n_positive = n_positive[1],
  n_total = n_total[1]
), by = .(target, atac_w, p_drop)]
summary_dt[, se_auroc := sd_auroc / sqrt(n_runs)]
summary_dt[, lo95 := auroc - 1.96 * se_auroc]
summary_dt[, hi95 := auroc + 1.96 * se_auroc]
summary_dt[, method := sprintf("w=%.1f", atac_w)]

pred_map <- data.table(
  target = names(PRED_TRAIT),
  trait = unlist(PRED_TRAIT, use.names = FALSE)
)
summary_dt <- merge(summary_dt, pred_map, by = "target", all.x = TRUE)
plot_dt <- summary_dt[, .(
  trait,
  target,
  p_drop,
  method,
  atac_w,
  auroc,
  lo95,
  hi95,
  sd_auroc,
  n_runs,
  n_positive,
  n_total
)]
setorder(plot_dt, trait, target, atac_w, p_drop)

if (ENFORCE_P0_CLEAN_MATCH && file.exists(CLEAN_AUROC_CSV)) {
  clean <- fread(CLEAN_AUROC_CSV)
  clean <- clean[!is.na(atac_w)]
  clean <- clean[target %in% plot_dt$target & atac_w %in% ATAC_W_LIST]
  p0 <- plot_dt[p_drop == 0, .(target, atac_w, p0_auroc = auroc)]
  cmp <- merge(
    p0,
    clean[, .(target, atac_w, clean_auroc = auroc)],
    by = c("target", "atac_w"),
    all.x = TRUE
  )
  cmp[, abs_delta := abs(p0_auroc - clean_auroc)]
  bad <- cmp[is.na(clean_auroc) | abs_delta > 1e-8]
  if (nrow(bad) > 0L) {
    stop(
      "p_drop=0 does not match clean AUROC:\n",
      paste(capture.output(print(bad)), collapse = "\n"),
      call. = FALSE
    )
  }
  cat("\np_drop=0 clean AUROC check passed.\n")
} else if (ENFORCE_P0_CLEAN_MATCH) {
  cat("\nClean AUROC CSV not found; skipped p_drop=0 match check:\n  ",
      CLEAN_AUROC_CSV, "\n", sep = "")
}

fwrite(plot_dt, PLOT_DATA_FILE)

cat("\nPlot-ready data:\n")
print(plot_dt[, .(
  trait,
  target,
  p_drop,
  method,
  auroc = round(auroc, 4),
  lo95 = round(lo95, 4),
  hi95 = round(hi95, 4),
  n_runs
)])

cat("\nWrote:\n")
cat("  ", OUT_ALL, "\n", sep = "")
cat("  ", PLOT_DATA_FILE, "\n", sep = "")
cat("  ", RUN_MANIFEST_FILE, "\n", sep = "")
