# Pure scDRS robustness for the reverse Granja RNA traits.
#
# Usage in RStudio:
#   setwd("D:/RStudio/work")
#   source("reverse_atac_to_rna_experiment_20260608/reverse_robustness/scripts/01_run_pure_scdrs_robustness_from_rstudio.R")
#
# Quick one-condition sanity:
#   ROBUST_P_DROP_LIST <- c(0.1)
#   ROBUST_SEED_LIST <- c(1L)
#   source("reverse_atac_to_rna_experiment_20260608/reverse_robustness/scripts/01_run_pure_scdrs_robustness_from_rstudio.R")

suppressPackageStartupMessages({
  library(data.table)
  library(pROC)
})

ROOT <- "D:/RStudio/work"
REV_ROOT <- file.path(ROOT, "reverse_atac_to_rna_experiment_20260608")
ROBUST_ROOT <- file.path(REV_ROOT, "reverse_robustness")
MODULE4_ROOT <- file.path(ROOT, "module4_traitaware_pseudopair")

RUN_ENV <- environment()
has_local <- function(x) exists(x, envir = RUN_ENV, inherits = FALSE)

if (!has_local("ROBUST_P_DROP_LIST")) ROBUST_P_DROP_LIST <- round(seq(0.1, 0.9, by = 0.1), 2)
if (!has_local("ROBUST_SEED_LIST")) ROBUST_SEED_LIST <- c(1L, 2L, 3L)
if (!has_local("ROBUST_TRAITS_TO_RUN")) ROBUST_TRAITS_TO_RUN <- c("PLT", "MCH", "RBC")
if (!has_local("N_CTRL")) N_CTRL <- 1000L
if (!has_local("RUN_DOWNSTREAM")) RUN_DOWNSTREAM <- TRUE
if (!has_local("RESUME_SKIP_COMPLETED")) RESUME_SKIP_COMPLETED <- TRUE
if (!has_local("FORCE_DROPOUT_H5AD")) FORCE_DROPOUT_H5AD <- FALSE
if (!has_local("FORCE_SCDRS")) FORCE_SCDRS <- FALSE
if (!has_local("FORCE_DOWNSTREAM")) FORCE_DOWNSTREAM <- FALSE

PYTHON_EXE <- "D:/python/python.exe"
SCDRS_CLI <- "D:/python/Scripts/scdrs"
BASE_H5AD <- file.path(
  MODULE4_ROOT,
  "01_data_prep/03_rna_scdrs_python/input/granja_rna_raw_counts.h5ad"
)
META_RDS <- file.path(MODULE4_ROOT, "01_data_prep/output/granja_rna_cell_meta.rds")
BASELINE_LONG <- file.path(REV_ROOT, "outputs/step0b_extended_rna_only_scdrs_baseline_auroc_long.csv")
DROP_HELPER <- file.path(ROBUST_ROOT, "scripts/make_granja_rna_dropout_h5ad.py")

GS_CANDIDATES <- c(
  file.path(MODULE4_ROOT, "01_data_prep/03_rna_scdrs_python/granja_3traits.gs"),
  file.path(MODULE4_ROOT, "01_data_prep_external_rna/03_cite_scdrs_python/cite_trait_screening.gs")
)

TRAIT_IDS <- c(
  PLT = "UKB_460K.blood_PLATELET_COUNT",
  MCH = "UKB_460K.blood_MEAN_CORPUSCULAR_HEMOGLOBIN",
  RBC = "UKB_460K.blood_RED_COUNT"
)

GRANJA_TARGETS <- list(
  PLT = c("02_Early.Eryth", "03_Late.Eryth", "04_Early.Baso"),
  MCH = c("02_Early.Eryth", "03_Late.Eryth", "04_Early.Baso"),
  RBC = c("02_Early.Eryth", "03_Late.Eryth", "04_Early.Baso")
)
GRANJA_TARGET_NAMES <- c(
  PLT = "erythroid_baso",
  MCH = "erythroid_baso",
  RBC = "erythroid_baso"
)

PURE_OUT_DIR <- file.path(ROBUST_ROOT, "outputs/pure_scdrs")
RUNS_DIR <- file.path(PURE_OUT_DIR, "runs")
LOG_DIR <- file.path(ROBUST_ROOT, "logs/pure_scdrs")
GS_SPLIT_DIR <- file.path(ROBUST_ROOT, "cache/gs_split")
NUMBA_CACHE_DIR <- file.path(ROBUST_ROOT, "cache/numba_cache")
dir.create(RUNS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(GS_SPLIT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(NUMBA_CACHE_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_EVAL <- file.path(PURE_OUT_DIR, "pure_scdrs_robustness_eval.csv")
OUT_SUMMARY <- file.path(PURE_OUT_DIR, "pure_scdrs_robustness_summary.csv")
OUT_PLOT_DATA <- file.path(PURE_OUT_DIR, "pure_scdrs_robustness_plot_data.csv")
OUT_MANIFEST <- file.path(PURE_OUT_DIR, "pure_scdrs_robustness_run_manifest.csv")

check_file <- function(path, min_size = 1L) {
  file.exists(path) && !is.na(file.info(path)$size) && file.info(path)$size >= min_size
}

win_path <- function(path, mustWork = TRUE) {
  if (mustWork || file.exists(path)) {
    return(normalizePath(path, winslash = "\\", mustWork = mustWork))
  }
  parent <- dirname(path)
  if (dir.exists(parent)) {
    return(paste0(normalizePath(parent, winslash = "\\", mustWork = TRUE),
                  "\\", basename(path)))
  }
  gsub("/", "\\\\", path, fixed = TRUE)
}

run_or_stop <- function(cmd, args, log_prefix) {
  stdout_log <- file.path(LOG_DIR, paste0(log_prefix, ".stdout.log"))
  stderr_log <- file.path(LOG_DIR, paste0(log_prefix, ".stderr.log"))
  cat("Running:\n  ", cmd, " ", paste(args, collapse = " "), "\n", sep = "")
  status <- system2(cmd, args = args, stdout = stdout_log, stderr = stderr_log)
  if (!identical(status, 0L)) {
    cat("\nSTDOUT log: ", stdout_log, "\n", sep = "")
    cat("STDERR log: ", stderr_log, "\n", sep = "")
    stop("Command failed with status ", status, call. = FALSE)
  }
  cat("  OK. Logs:\n")
  cat("    ", stdout_log, "\n", sep = "")
  cat("    ", stderr_log, "\n", sep = "")
}

read_trait_row <- function(trait_id) {
  for (gs_file in GS_CANDIDATES[file.exists(GS_CANDIDATES)]) {
    lines <- readLines(gs_file, warn = FALSE)
    header <- lines[1]
    hit <- grep(paste0("^", trait_id, "\t"), lines[-1], value = TRUE)
    if (length(hit) == 1L) {
      return(list(header = header, row = hit, source = gs_file))
    }
  }
  stop("Could not find gene-set row for trait: ", trait_id, call. = FALSE)
}

write_one_trait_gs <- function(short_trait) {
  trait_id <- TRAIT_IDS[[short_trait]]
  x <- read_trait_row(trait_id)
  out_gs <- file.path(GS_SPLIT_DIR, paste0(trait_id, ".gs"))
  writeLines(c(x$header, x$row), out_gs, useBytes = TRUE)
  data.table(
    trait = short_trait,
    trait_id = trait_id,
    gs_file = out_gs,
    source_gs = x$source
  )
}

score_col <- function(dt) {
  hit <- intersect(c("norm_score", "zscore", "score"), names(dt))[1]
  if (is.na(hit)) stop("No scDRS score column found in: ", paste(names(dt), collapse = ", "))
  hit
}

compute_auc <- function(score, label) {
  keep <- is.finite(score) & !is.na(label)
  score <- score[keep]
  label <- label[keep]
  if (length(score) < 20 || sum(label == 1L) < 10 || sum(label == 0L) < 10) {
    return(NA_real_)
  }
  as.numeric(pROC::auc(pROC::roc(
    label, score,
    quiet = TRUE,
    direction = "<",
    levels = c(0L, 1L)
  )))
}

load_meta <- function() {
  meta <- readRDS(META_RDS)
  if (is.null(rownames(meta))) {
    if (!("barcode" %in% colnames(meta))) {
      stop("Granja metadata lacks rownames and barcode column.", call. = FALSE)
    }
    rownames(meta) <- as.character(meta$barcode)
  }
  data.table(
    barcode = rownames(meta),
    cell_type = as.character(meta$cell_type)
  )
}

condition_eval <- function(scdrs_out_dir, p_drop, seed, meta_dt) {
  rows <- lapply(ROBUST_TRAITS_TO_RUN, function(short) {
    trait_id <- TRAIT_IDS[[short]]
    score_file <- file.path(scdrs_out_dir, paste0(trait_id, ".score.gz"))
    if (!check_file(score_file, 1000L)) {
      stop("Missing scDRS score file: ", score_file, call. = FALSE)
    }
    dt <- fread(score_file)
    bc_col <- intersect(c("cell", "cell_id", "barcode", "V1"), names(dt))[1]
    if (is.na(bc_col)) bc_col <- names(dt)[1]
    sc_col <- score_col(dt)
    score_dt <- dt[, .(
      barcode = as.character(get(bc_col)),
      score = as.numeric(get(sc_col))
    )]
    x <- merge(meta_dt, score_dt, by = "barcode", all.x = FALSE, sort = FALSE)
    label <- as.integer(x$cell_type %in% GRANJA_TARGETS[[short]])
    keep <- is.finite(x$score) & !is.na(label)
    data.table(
      method = "pure_scdrs",
      dataset = "GRANJA_RNA",
      trait = short,
      trait_id = trait_id,
      target_definition = GRANJA_TARGET_NAMES[[short]],
      p_drop = as.numeric(p_drop),
      seed = as.integer(seed),
      target = GRANJA_TARGET_NAMES[[short]],
      n_positive = sum(label[keep] == 1L),
      n_total = sum(keep),
      auroc = compute_auc(x$score, label),
      source_file = score_file,
      note = "official scDRS on degraded Granja RNA"
    )
  })
  rbindlist(rows, use.names = TRUE, fill = TRUE)
}

clean_anchor <- function() {
  if (!file.exists(BASELINE_LONG)) {
    warning("Missing clean anchor file: ", BASELINE_LONG)
    return(data.table())
  }
  dt <- fread(BASELINE_LONG)
  dt <- dt[dataset == "GRANJA_RNA" & trait %in% ROBUST_TRAITS_TO_RUN]
  if (nrow(dt) == 0L) return(data.table())
  data.table(
    method = "pure_scdrs",
    dataset = "GRANJA_RNA",
    trait = dt$trait,
    trait_id = dt$trait_id,
    target_definition = dt$target_definition,
    p_drop = 0,
    seed = 0L,
    target = dt$target_definition,
    n_positive = dt$n_positive_target,
    n_total = dt$n_scored_cells,
    auroc = dt$baseline_scdrs_auroc,
    source_file = BASELINE_LONG,
    note = "clean p_drop=0 anchor from existing Step0b table"
  )
}

summarize_eval <- function(all_eval) {
  summary_dt <- all_eval[is.finite(auroc), .(
    mean_auroc = mean(auroc, na.rm = TRUE),
    sd_auroc = if (.N > 1L) stats::sd(auroc, na.rm = TRUE) else NA_real_,
    n_runs = .N,
    n_positive = n_positive[1],
    n_total = n_total[1]
  ), by = .(method, dataset, trait, trait_id, target_definition, p_drop)]
  summary_dt[, se_auroc := fifelse(n_runs > 1L, sd_auroc / sqrt(n_runs), NA_real_)]
  summary_dt[, lo95 := fifelse(n_runs > 1L, mean_auroc - 1.96 * se_auroc, mean_auroc)]
  summary_dt[, hi95 := fifelse(n_runs > 1L, mean_auroc + 1.96 * se_auroc, mean_auroc)]
  setorder(summary_dt, trait, p_drop)
  summary_dt
}

make_dropout_h5ad <- function(run_dir, p_drop, seed) {
  out_h5ad <- file.path(run_dir, "granja_rna_dropout.h5ad")
  summary_csv <- file.path(run_dir, "granja_rna_dropout_summary.csv")
  if (!FORCE_DROPOUT_H5AD && check_file(out_h5ad, 1000L) && check_file(summary_csv, 1L)) {
    cat("  [dropout] skip existing h5ad: ", out_h5ad, "\n", sep = "")
    return(out_h5ad)
  }
  run_or_stop(
    win_path(PYTHON_EXE),
    args = c(
      win_path(DROP_HELPER),
      "--input-h5ad", win_path(BASE_H5AD),
      "--output-h5ad", win_path(out_h5ad, mustWork = FALSE),
      "--p-drop", as.character(p_drop),
      "--seed", as.character(seed),
      "--summary-csv", win_path(summary_csv, mustWork = FALSE)
    ),
    log_prefix = sprintf("p%.2f_seed%d_make_dropout_h5ad", p_drop, seed)
  )
  out_h5ad
}

condition_paths <- function(p_drop, seed) {
  run_id <- sprintf("p%.2f_seed%d", p_drop, seed)
  run_dir <- file.path(RUNS_DIR, run_id)
  list(
    run_id = run_id,
    run_dir = run_dir,
    h5ad_file = file.path(run_dir, "granja_rna_dropout.h5ad"),
    scdrs_out_dir = file.path(run_dir, "scdrs"),
    condition_eval_file = file.path(run_dir, "pure_scdrs_condition_eval.csv")
  )
}

expected_score_files <- function(scdrs_out_dir, suffix) {
  file.path(scdrs_out_dir, paste0(unname(TRAIT_IDS[ROBUST_TRAITS_TO_RUN]), suffix))
}

condition_missing_items <- function(paths) {
  missing <- character()
  if (!check_file(paths$h5ad_file, 1000L)) {
    missing <- c(missing, "dropout_h5ad")
  }
  score_files <- expected_score_files(paths$scdrs_out_dir, ".score.gz")
  full_files <- expected_score_files(paths$scdrs_out_dir, ".full_score.gz")
  if (!all(vapply(score_files, check_file, logical(1), min_size = 1000L))) {
    missing <- c(missing, "score_files")
  }
  if (!all(vapply(full_files, check_file, logical(1), min_size = 1000L))) {
    missing <- c(missing, "full_score_files")
  }
  if (isTRUE(RUN_DOWNSTREAM)) {
    group_files <- expected_score_files(paths$scdrs_out_dir, ".scdrs_group.cell_type")
    if (!all(vapply(group_files, check_file, logical(1), min_size = 1L))) {
      missing <- c(missing, "group_analysis_files")
    }
  }
  if (!check_file(paths$condition_eval_file, 1L)) {
    missing <- c(missing, "condition_eval")
  }
  missing
}

condition_complete <- function(paths) {
  length(condition_missing_items(paths)) == 0L
}

manifest_row <- function(paths, p_drop, seed, status, missing_items = character()) {
  data.table(
    run_id = paths$run_id,
    p_drop = p_drop,
    seed = seed,
    h5ad_file = paths$h5ad_file,
    scdrs_out_dir = paths$scdrs_out_dir,
    condition_eval_file = paths$condition_eval_file,
    status = status,
    missing_items = if (length(missing_items) == 0L) "" else paste(missing_items, collapse = ";")
  )
}

run_scdrs_condition <- function(run_id, h5ad_file, scdrs_out_dir) {
  dir.create(scdrs_out_dir, recursive = TRUE, showWarnings = FALSE)
  gs_manifest <- rbindlist(lapply(ROBUST_TRAITS_TO_RUN, write_one_trait_gs))
  fwrite(gs_manifest, file.path(scdrs_out_dir, "trait_gs_manifest.csv"))

  for (short in ROBUST_TRAITS_TO_RUN) {
    trait_id <- TRAIT_IDS[[short]]
    score_file <- file.path(scdrs_out_dir, paste0(trait_id, ".score.gz"))
    full_file <- file.path(scdrs_out_dir, paste0(trait_id, ".full_score.gz"))
    if (!FORCE_SCDRS && check_file(score_file, 1000L) && check_file(full_file, 1000L)) {
      cat("  [scDRS] skip existing score files for ", short, "\n", sep = "")
      next
    }
    one_gs <- file.path(GS_SPLIT_DIR, paste0(trait_id, ".gs"))
    run_or_stop(
      win_path(PYTHON_EXE),
      args = c(
        win_path(SCDRS_CLI),
        "compute-score",
        "--h5ad-file", win_path(h5ad_file),
        "--h5ad-species", "human",
        "--gs-file", win_path(one_gs),
        "--gs-species", "human",
        "--out-folder", win_path(scdrs_out_dir),
        "--n-ctrl", as.character(N_CTRL)
      ),
      log_prefix = paste0(run_id, "_", short, "_compute_score")
    )
  }

  if (isTRUE(RUN_DOWNSTREAM)) {
    trait_ids <- unname(TRAIT_IDS[ROBUST_TRAITS_TO_RUN])
    expected_group <- file.path(scdrs_out_dir, paste0(trait_ids, ".scdrs_group.cell_type"))
    if (!FORCE_DOWNSTREAM && all(vapply(expected_group, check_file, logical(1), min_size = 1L))) {
      cat("  [scDRS] skip perform-downstream; group files already exist.\n")
    } else {
      full_files <- file.path(scdrs_out_dir, paste0(trait_ids, ".full_score.gz"))
      missing_full <- full_files[!vapply(full_files, check_file, logical(1), min_size = 1000L)]
      if (length(missing_full) > 0L) {
        stop("Cannot run perform-downstream because full_score files are missing:\n  ",
             paste(missing_full, collapse = "\n  "), call. = FALSE)
      }
      run_or_stop(
        win_path(PYTHON_EXE),
        args = c(
          win_path(SCDRS_CLI),
          "perform-downstream",
          "--h5ad-file", win_path(h5ad_file),
          "--score-file", win_path(file.path(scdrs_out_dir, "@.full_score.gz"), mustWork = FALSE),
          "--out-folder", win_path(scdrs_out_dir),
          "--group-analysis", "cell_type"
        ),
        log_prefix = paste0(run_id, "_perform_downstream")
      )
    }
  }

  invisible(TRUE)
}

cat(">>> Pure scDRS robustness for reverse Granja RNA traits\n\n")
cat("Traits      : ", paste(ROBUST_TRAITS_TO_RUN, collapse = ", "), "\n", sep = "")
cat("P_DROP_LIST : ", paste(ROBUST_P_DROP_LIST, collapse = ", "), "\n", sep = "")
cat("SEED_LIST   : ", paste(ROBUST_SEED_LIST, collapse = ", "), "\n", sep = "")
cat("Output      : ", PURE_OUT_DIR, "\n\n", sep = "")

required <- c(PYTHON_EXE, SCDRS_CLI, BASE_H5AD, META_RDS, DROP_HELPER)
missing <- required[!file.exists(required)]
if (length(missing) > 0L) {
  stop("Missing required file(s):\n  ", paste(missing, collapse = "\n  "), call. = FALSE)
}
if (any(!ROBUST_TRAITS_TO_RUN %in% names(TRAIT_IDS))) {
  stop("Unknown trait(s): ", paste(setdiff(ROBUST_TRAITS_TO_RUN, names(TRAIT_IDS)), collapse = ", "),
       call. = FALSE)
}

Sys.setenv(NUMBA_CACHE_DIR = NUMBA_CACHE_DIR)
meta_dt <- load_meta()
manifest_rows <- list()

for (p_drop in ROBUST_P_DROP_LIST) {
  for (seed in ROBUST_SEED_LIST) {
    paths <- condition_paths(p_drop, seed)
    run_id <- paths$run_id
    run_dir <- paths$run_dir
    scdrs_out_dir <- paths$scdrs_out_dir
    dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

    cat("\n", strrep("=", 78), "\n", sep = "")
    cat("[condition] ", run_id, "\n", sep = "")
    cat(strrep("=", 78), "\n", sep = "")

    missing_before <- condition_missing_items(paths)
    if (isTRUE(RESUME_SKIP_COMPLETED) && !FORCE_DROPOUT_H5AD && !FORCE_SCDRS &&
        !FORCE_DOWNSTREAM && length(missing_before) == 0L) {
      cat("  [skip] condition already complete.\n")
      manifest_rows[[run_id]] <- manifest_row(
        paths, p_drop, seed,
        status = "skipped_complete",
        missing_items = character()
      )
      fwrite(rbindlist(manifest_rows, use.names = TRUE, fill = TRUE), OUT_MANIFEST)
      next
    }
    if (length(missing_before) > 0L) {
      cat("  [resume] missing: ", paste(missing_before, collapse = ", "), "\n", sep = "")
    }

    h5ad_file <- make_dropout_h5ad(run_dir, p_drop = p_drop, seed = seed)
    run_scdrs_condition(run_id, h5ad_file, scdrs_out_dir)

    eval_dt <- condition_eval(scdrs_out_dir, p_drop = p_drop, seed = seed, meta_dt = meta_dt)
    fwrite(eval_dt, paths$condition_eval_file)

    missing_after <- condition_missing_items(paths)
    manifest_rows[[run_id]] <- manifest_row(
      paths, p_drop, seed,
      status = if (length(missing_after) == 0L) "completed" else "incomplete",
      missing_items = missing_after
    )
    fwrite(rbindlist(manifest_rows, use.names = TRUE, fill = TRUE), OUT_MANIFEST)
  }
}

condition_eval_files <- list.files(
  RUNS_DIR,
  pattern = "^pure_scdrs_condition_eval[.]csv$",
  recursive = TRUE,
  full.names = TRUE
)
stress_eval <- rbindlist(lapply(condition_eval_files, fread), use.names = TRUE, fill = TRUE)
all_eval <- rbindlist(list(clean_anchor(), stress_eval), use.names = TRUE, fill = TRUE)
setorder(all_eval, trait, p_drop, seed)

summary_dt <- summarize_eval(all_eval)
plot_dt <- copy(summary_dt)

fwrite(all_eval, OUT_EVAL)
fwrite(summary_dt, OUT_SUMMARY)
fwrite(plot_dt, OUT_PLOT_DATA)

cat("\n=== Pure scDRS robustness finished ===\n")
cat("Wrote:\n")
cat("  ", OUT_EVAL, "\n", sep = "")
cat("  ", OUT_SUMMARY, "\n", sep = "")
cat("  ", OUT_PLOT_DATA, "\n", sep = "")
cat("  ", OUT_MANIFEST, "\n", sep = "")
cat("\nSummary:\n")
print(summary_dt[, .(
  trait,
  p_drop,
  mean_auroc = round(mean_auroc, 6),
  lo95 = round(lo95, 6),
  hi95 = round(hi95, 6),
  n_runs
)])
