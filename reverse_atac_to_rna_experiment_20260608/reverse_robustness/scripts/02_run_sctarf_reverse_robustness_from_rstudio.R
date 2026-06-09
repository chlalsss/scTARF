# scTARF reverse robustness for Granja RNA dropout.
#
# Usage in RStudio:
#   setwd("D:/RStudio/work")
#   source("reverse_atac_to_rna_experiment_20260608/reverse_robustness/scripts/02_run_sctarf_reverse_robustness_from_rstudio.R")
#
# Quick one-condition sanity:
#   ROBUST_P_DROP_LIST <- c(0.1)
#   ROBUST_SEED_LIST <- c(1L)
#   source("reverse_atac_to_rna_experiment_20260608/reverse_robustness/scripts/02_run_sctarf_reverse_robustness_from_rstudio.R")

suppressPackageStartupMessages({
  library(data.table)
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
if (!has_local("ROBUST_ATAC_W_LIST")) ROBUST_ATAC_W_LIST <- c(0.7, 1.0)
if (!has_local("RESUME_SKIP_COMPLETED")) RESUME_SKIP_COMPLETED <- TRUE
if (!has_local("FORCE_DROPOUT_H5AD")) FORCE_DROPOUT_H5AD <- FALSE
if (!has_local("FORCE_CONDITION_META")) FORCE_CONDITION_META <- FALSE
if (!has_local("FORCE_REMAP")) FORCE_REMAP <- FALSE
if (!has_local("FORCE_RETRANSFER")) FORCE_RETRANSFER <- FALSE
if (!has_local("FORCE_REBUILD_RNA_MAIN_Z")) FORCE_REBUILD_RNA_MAIN_Z <- FALSE
if (!has_local("FORCE_REEXTRACT_RNA")) FORCE_REEXTRACT_RNA <- FALSE
if (!has_local("FORCE_REEXTRACT_RNA_MAIN")) FORCE_REEXTRACT_RNA_MAIN <- FALSE

PYTHON_EXE <- "D:/python/python.exe"
BASE_H5AD <- file.path(
  MODULE4_ROOT,
  "01_data_prep/03_rna_scdrs_python/input/granja_rna_raw_counts.h5ad"
)
BASE_META_RDS <- file.path(MODULE4_ROOT, "01_data_prep/output/granja_rna_cell_meta.rds")
DROP_HELPER <- file.path(ROBUST_ROOT, "scripts/make_granja_rna_dropout_h5ad.py")
COUNT_HELPER <- file.path(ROBUST_ROOT, "scripts/summarize_h5ad_cell_counts.py")

CLEAN_REVERSE_TAG <- "reverse_granja_rna_main_plt_mch_rbc"
CLEAN_AUROC_CSV <- file.path(
  REV_ROOT,
  "outputs",
  CLEAN_REVERSE_TAG,
  "auroc",
  paste0(CLEAN_REVERSE_TAG, "_auroc.csv")
)
PURE_SCDRS_SUMMARY <- file.path(
  ROBUST_ROOT,
  "outputs/pure_scdrs/pure_scdrs_robustness_summary.csv"
)

SCTARF_OUT_DIR <- file.path(ROBUST_ROOT, "outputs/sctarf_reverse")
RUNS_DIR <- file.path(SCTARF_OUT_DIR, "runs")
LOG_DIR <- file.path(ROBUST_ROOT, "logs/sctarf_reverse")
COMBINED_DIR <- file.path(ROBUST_ROOT, "outputs/combined")
dir.create(RUNS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(COMBINED_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_EVAL <- file.path(SCTARF_OUT_DIR, "sctarf_reverse_robustness_eval.csv")
OUT_SUMMARY <- file.path(SCTARF_OUT_DIR, "sctarf_reverse_robustness_summary.csv")
OUT_PLOT_DATA <- file.path(SCTARF_OUT_DIR, "sctarf_reverse_robustness_plot_data.csv")
OUT_MANIFEST <- file.path(SCTARF_OUT_DIR, "sctarf_reverse_robustness_run_manifest.csv")

COMBINED_SUMMARY <- file.path(COMBINED_DIR, "reverse_robustness_combined_summary.csv")
COMBINED_PLOT_DATA <- file.path(COMBINED_DIR, "reverse_robustness_combined_plot_data.csv")
COMBINED_PNG <- file.path(COMBINED_DIR, "reverse_robustness_auroc_curves.png")
COMBINED_PDF <- file.path(COMBINED_DIR, "reverse_robustness_auroc_curves.pdf")

TRAIT_IDS <- c(
  PLT = "UKB_460K.blood_PLATELET_COUNT",
  MCH = "UKB_460K.blood_MEAN_CORPUSCULAR_HEMOGLOBIN",
  RBC = "UKB_460K.blood_RED_COUNT"
)

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

weight_label <- function(w) format(as.numeric(w), trim = TRUE, scientific = FALSE, nsmall = 1L)

condition_paths <- function(p_drop, seed) {
  run_id <- sprintf("p%.2f_seed%d", p_drop, seed)
  run_dir <- file.path(RUNS_DIR, run_id)
  inputs_dir <- file.path(run_dir, "_inputs")
  pure_run_dir <- file.path(ROBUST_ROOT, "outputs/pure_scdrs/runs", run_id)
  list(
    run_id = run_id,
    run_dir = run_dir,
    inputs_dir = inputs_dir,
    pure_h5ad_file = file.path(pure_run_dir, "granja_rna_dropout.h5ad"),
    local_h5ad_file = file.path(inputs_dir, "granja_rna_dropout.h5ad"),
    local_h5ad_summary = file.path(inputs_dir, "granja_rna_dropout_summary.csv"),
    cell_counts_csv = file.path(inputs_dir, "granja_rna_dropout_cell_counts.csv"),
    meta_rds = file.path(inputs_dir, "granja_rna_cell_meta_dropout.rds"),
    align_dir = file.path(run_dir, "mapping"),
    rna_extract_dir = file.path(run_dir, "mapping/granja_rna_extract"),
    transfer_dir = file.path(run_dir, "transfer"),
    rna_main_z_rds = file.path(inputs_dir, "granja_rna_main_z_block_v2.rds"),
    rna_main_extract_dir = file.path(inputs_dir, "granja_rna_main_extract"),
    module3_dir = file.path(run_dir, "module3_v2_compat"),
    auroc_dir = file.path(run_dir, "auroc"),
    auroc_csv = file.path(run_dir, "auroc", paste0(run_id, "_sctarf_reverse_auroc.csv")),
    baseline_rds = file.path(inputs_dir, "granja_rna_scdrs_baseline_matrix.rds"),
    patched_eval_script = file.path(run_dir, "auroc", paste0("eval_auroc.", run_id, "_patched.R")),
    condition_eval_file = file.path(run_dir, "sctarf_reverse_condition_eval.csv")
  )
}

available_h5ad <- function(paths) {
  if (check_file(paths$pure_h5ad_file, 1000L) && !isTRUE(FORCE_DROPOUT_H5AD)) {
    return(paths$pure_h5ad_file)
  }
  if (check_file(paths$local_h5ad_file, 1000L) && !isTRUE(FORCE_DROPOUT_H5AD)) {
    return(paths$local_h5ad_file)
  }
  ""
}

expected_mapping_files <- function(paths) {
  file.path(paths$align_dir, paste0("granja_rna_to_atac_mapping_", ROBUST_TRAITS_TO_RUN, ".rds"))
}

expected_transfer_files <- function(paths) {
  c(
    file.path(paths$transfer_dir, paste0("granja_pseudo_rna_profile_", ROBUST_TRAITS_TO_RUN, ".rds")),
    file.path(paths$transfer_dir, paste0("granja_pseudo_zrna_", ROBUST_TRAITS_TO_RUN, ".rds")),
    file.path(paths$transfer_dir, paste0("granja_pseudo_zrna_", ROBUST_TRAITS_TO_RUN, ".csv"))
  )
}

expected_trs_files <- function(paths) {
  unlist(lapply(ROBUST_TRAITS_TO_RUN, function(trait) {
    file.path(paths$module3_dir, paste0("TRS_v2_", trait, "_w", weight_label(ROBUST_ATAC_W_LIST), ".tsv.gz"))
  }), use.names = FALSE)
}

condition_missing_items <- function(paths) {
  missing <- character()
  if (!check_file(available_h5ad(paths), 1000L)) missing <- c(missing, "dropout_h5ad")
  if (!check_file(paths$meta_rds, 1000L)) missing <- c(missing, "dropout_meta_rds")
  if (!all(vapply(expected_mapping_files(paths), check_file, logical(1), min_size = 1000L))) {
    missing <- c(missing, "mapping")
  }
  if (!all(vapply(expected_transfer_files(paths), check_file, logical(1), min_size = 1000L))) {
    missing <- c(missing, "transfer")
  }
  if (!check_file(paths$rna_main_z_rds, 1000L)) missing <- c(missing, "rna_main_z")
  if (!all(vapply(expected_trs_files(paths), check_file, logical(1), min_size = 1000L))) {
    missing <- c(missing, "module3_trs")
  }
  if (!check_file(paths$auroc_csv, 1L)) missing <- c(missing, "auroc_csv")
  if (!check_file(paths$condition_eval_file, 1L)) missing <- c(missing, "condition_eval")
  missing
}

manifest_row <- function(paths, p_drop, seed, h5ad_file, status, missing_items = character()) {
  data.table(
    run_id = paths$run_id,
    p_drop = p_drop,
    seed = seed,
    h5ad_file = h5ad_file,
    run_dir = paths$run_dir,
    auroc_csv = paths$auroc_csv,
    condition_eval_file = paths$condition_eval_file,
    status = status,
    missing_items = if (length(missing_items) == 0L) "" else paste(missing_items, collapse = ";")
  )
}

ensure_dropout_h5ad <- function(paths, p_drop, seed) {
  h5ad <- available_h5ad(paths)
  if (nzchar(h5ad)) {
    cat("  [dropout] reuse h5ad: ", h5ad, "\n", sep = "")
    return(h5ad)
  }
  dir.create(paths$inputs_dir, recursive = TRUE, showWarnings = FALSE)
  run_or_stop(
    win_path(PYTHON_EXE),
    args = c(
      win_path(DROP_HELPER),
      "--input-h5ad", win_path(BASE_H5AD),
      "--output-h5ad", win_path(paths$local_h5ad_file, mustWork = FALSE),
      "--p-drop", as.character(p_drop),
      "--seed", as.character(seed),
      "--summary-csv", win_path(paths$local_h5ad_summary, mustWork = FALSE)
    ),
    log_prefix = sprintf("%s_make_dropout_h5ad", paths$run_id)
  )
  paths$local_h5ad_file
}

ensure_condition_meta <- function(paths, h5ad_file) {
  if (!FORCE_CONDITION_META && check_file(paths$meta_rds, 1000L)) {
    cat("  [meta] skip existing dropout metadata: ", paths$meta_rds, "\n", sep = "")
    return(paths$meta_rds)
  }
  dir.create(paths$inputs_dir, recursive = TRUE, showWarnings = FALSE)
  run_or_stop(
    win_path(PYTHON_EXE),
    args = c(
      win_path(COUNT_HELPER),
      "--h5ad", win_path(h5ad_file),
      "--out-csv", win_path(paths$cell_counts_csv, mustWork = FALSE)
    ),
    log_prefix = paste0(paths$run_id, "_summarize_h5ad_counts")
  )

  meta <- readRDS(BASE_META_RDS)
  counts <- fread(paths$cell_counts_csv)
  if (is.null(rownames(meta)) || identical(rownames(meta), as.character(seq_len(nrow(meta))))) {
    if (!"barcode" %in% colnames(meta)) {
      stop("Base Granja RNA metadata has no rownames or barcode column.", call. = FALSE)
    }
    rownames(meta) <- as.character(meta$barcode)
  }
  counts <- counts[barcode %in% rownames(meta)]
  if (nrow(counts) == 0L) stop("No h5ad barcodes overlap metadata.", call. = FALSE)
  for (col in intersect(c("scdrs_input_total_counts", "nUMI", "nUMI_pre"), colnames(meta))) {
    meta[counts$barcode, col] <- counts$dropout_total_counts
  }
  meta[counts$barcode, "dropout_detected_genes"] <- counts$dropout_detected_genes
  saveRDS(meta, paths$meta_rds)
  paths$meta_rds
}

run_local_script <- function(script_name, vars) {
  script_path <- file.path(REV_ROOT, "scripts", script_name)
  cat("\n", strrep("-", 78), "\n", sep = "")
  cat("source: ", script_path, "\n", sep = "")
  cat(strrep("-", 78), "\n", sep = "")

  # The reused reverse scripts check configuration with exists(..., inherits=FALSE).
  # Sourcing them in a child environment can therefore miss values supplied by
  # this driver. Temporarily inject only the requested configuration into the
  # global environment, source the script, then restore the user's session.
  ge <- .GlobalEnv
  before_names <- ls(envir = ge, all.names = TRUE)
  var_names <- names(vars)
  old_exists <- var_names %in% before_names
  old_values <- vector("list", length(var_names))
  names(old_values) <- var_names
  for (nm in var_names[old_exists]) old_values[[nm]] <- get(nm, envir = ge, inherits = FALSE)

  cleanup <- function() {
    after_names <- ls(envir = ge, all.names = TRUE)
    new_names <- setdiff(after_names, before_names)
    if (length(new_names) > 0L) rm(list = new_names, envir = ge)
    for (nm in var_names) {
      if (old_exists[match(nm, var_names)]) {
        assign(nm, old_values[[nm]], envir = ge)
      } else if (exists(nm, envir = ge, inherits = FALSE)) {
        rm(list = nm, envir = ge)
      }
    }
  }

  for (nm in var_names) assign(nm, vars[[nm]], envir = ge)
  on.exit(cleanup(), add = TRUE)
  source(script_path, local = ge)
  invisible(TRUE)
}

run_condition_pipeline <- function(paths, h5ad_file, meta_rds) {
  common_vars <- list(
    TRAITS_TO_RUN = ROBUST_TRAITS_TO_RUN,
    REVERSE_TAG = paths$run_id,
    GRANJA_RNA_H5AD = h5ad_file,
    GRANJA_RNA_META_RDS = meta_rds,
    ALLOW_OUTSIDE_REV_ROOT = FALSE,
    FORCE_REMAP = FORCE_REMAP,
    FORCE_RETRANSFER = FORCE_RETRANSFER,
    FORCE_REBUILD_RNA_MAIN_Z = FORCE_REBUILD_RNA_MAIN_Z,
    FORCE_REEXTRACT_RNA = FORCE_REEXTRACT_RNA,
    FORCE_REEXTRACT_RNA_MAIN = FORCE_REEXTRACT_RNA_MAIN
  )

  run_local_script("02_reverse_align_granja_rna_to_atac.R", c(common_vars, list(
    ALIGN_OUT_DIR = paths$align_dir,
    RNA_EXTRACT_DIR = paths$rna_extract_dir
  )))

  run_local_script("03_reverse_transfer_atac_to_granja_rna.R", c(common_vars, list(
    MAP_DIR = paths$align_dir,
    TRANSFER_OUT_DIR = paths$transfer_dir
  )))

  run_local_script("04_reverse_build_granja_rna_main_z_block.R", c(common_vars, list(
    OUT_DIR = paths$inputs_dir,
    RNA_MAIN_Z_RDS = paths$rna_main_z_rds,
    RNA_MAIN_EXTRACT_DIR = paths$rna_main_extract_dir
  )))

  run_local_script("05_reverse_prepare_module3_inputs.R", c(common_vars, list(
    MODULE3_TRANSFER_DIR = paths$transfer_dir,
    MODULE3_COMPAT_OUT_DIR = paths$module3_dir,
    RNA_MAIN_Z_RDS = paths$rna_main_z_rds,
    ATAC_Z_BLOCK_OVERRIDE = paths$rna_main_z_rds
  )))

  run_local_script("06_reverse_run_module3_knn_rwr.R", c(common_vars, list(
    MODULE3_COMPAT_OUT_DIR = paths$module3_dir,
    ATAC_W_LIST = ROBUST_ATAC_W_LIST,
    PATCHED_DIR_OVERRIDE = file.path(paths$module3_dir, "_patched_module3_scripts")
  )))

  run_local_script("07_reverse_eval_auroc.R", c(common_vars, list(
    V2_FILES_DIR_NEW = paths$module3_dir,
    OUT_DIR_NEW = paths$auroc_dir,
    OUT_CSV_NAME_NEW = basename(paths$auroc_csv),
    BASELINE_RDS = paths$baseline_rds,
    PATCHED_SCRIPT = paths$patched_eval_script
  )))
}

write_condition_eval <- function(paths, p_drop, seed) {
  if (!check_file(paths$auroc_csv, 1L)) {
    stop("Missing AUROC CSV for condition: ", paths$auroc_csv, call. = FALSE)
  }
  dt <- fread(paths$auroc_csv)
  dt <- dt[method == "reverse_atac_to_rna"]
  if (nrow(dt) == 0L) stop("No reverse_atac_to_rna rows in: ", paths$auroc_csv, call. = FALSE)
  out <- data.table(
    run_id = paths$run_id,
    method = "sctarf_reverse",
    dataset = "GRANJA_RNA",
    trait = dt$target,
    p_drop = as.numeric(p_drop),
    seed = as.integer(seed),
    atac_w = as.numeric(dt$atac_w),
    target = dt$target,
    n_positive = as.integer(dt$n_positive),
    n_total = as.integer(dt$n_total),
    auroc = as.numeric(dt$auroc),
    source_file = paths$auroc_csv,
    note = "scTARF reverse on degraded Granja RNA"
  )
  fwrite(out, paths$condition_eval_file)
  out
}

clean_sctarf_anchor <- function() {
  if (!check_file(CLEAN_AUROC_CSV, 1L)) {
    warning("Missing clean scTARF anchor: ", CLEAN_AUROC_CSV)
    return(data.table())
  }
  dt <- fread(CLEAN_AUROC_CSV)
  dt <- dt[method == "reverse_atac_to_rna" & target %in% ROBUST_TRAITS_TO_RUN]
  if (nrow(dt) == 0L) return(data.table())
  data.table(
    run_id = "clean_p0",
    method = "sctarf_reverse",
    dataset = "GRANJA_RNA",
    trait = dt$target,
    p_drop = 0,
    seed = 0L,
    atac_w = as.numeric(dt$atac_w),
    target = dt$target,
    n_positive = as.integer(dt$n_positive),
    n_total = as.integer(dt$n_total),
    auroc = as.numeric(dt$auroc),
    source_file = CLEAN_AUROC_CSV,
    note = "clean p_drop=0 anchor from main reverse run"
  )
}

summarize_sctarf <- function(all_eval) {
  summary_dt <- all_eval[is.finite(auroc), .(
    mean_auroc = mean(auroc, na.rm = TRUE),
    sd_auroc = if (.N > 1L) stats::sd(auroc, na.rm = TRUE) else NA_real_,
    n_runs = .N,
    n_positive = n_positive[1],
    n_total = n_total[1]
  ), by = .(method, dataset, trait, p_drop, atac_w)]
  summary_dt[, se_auroc := fifelse(n_runs > 1L, sd_auroc / sqrt(n_runs), NA_real_)]
  summary_dt[, lo95 := fifelse(n_runs > 1L, mean_auroc - 1.96 * se_auroc, mean_auroc)]
  summary_dt[, hi95 := fifelse(n_runs > 1L, mean_auroc + 1.96 * se_auroc, mean_auroc)]
  setorder(summary_dt, trait, atac_w, p_drop)
  summary_dt
}

write_combined_outputs <- function(sctarf_summary) {
  if (!check_file(PURE_SCDRS_SUMMARY, 1L)) {
    warning("Missing pure scDRS summary; combined plot was not written: ", PURE_SCDRS_SUMMARY)
    return(invisible(FALSE))
  }
  pure <- fread(PURE_SCDRS_SUMMARY)
  pure_plot <- data.table(
    method = "pure_scdrs",
    curve = "pure scDRS",
    trait = pure$trait,
    p_drop = pure$p_drop,
    atac_w = NA_real_,
    mean_auroc = pure$mean_auroc,
    lo95 = pure$lo95,
    hi95 = pure$hi95,
    n_runs = pure$n_runs,
    n_positive = pure$n_positive,
    n_total = pure$n_total
  )
  sctarf_plot <- data.table(
    method = "sctarf_reverse",
    curve = paste0("scTARF reverse w=", weight_label(sctarf_summary$atac_w)),
    trait = sctarf_summary$trait,
    p_drop = sctarf_summary$p_drop,
    atac_w = sctarf_summary$atac_w,
    mean_auroc = sctarf_summary$mean_auroc,
    lo95 = sctarf_summary$lo95,
    hi95 = sctarf_summary$hi95,
    n_runs = sctarf_summary$n_runs,
    n_positive = sctarf_summary$n_positive,
    n_total = sctarf_summary$n_total
  )
  combined <- rbindlist(list(pure_plot, sctarf_plot), use.names = TRUE, fill = TRUE)
  combined[, trait := factor(as.character(trait), levels = ROBUST_TRAITS_TO_RUN)]
  combined[, curve := factor(curve, levels = c(
    "pure scDRS",
    paste0("scTARF reverse w=", weight_label(ROBUST_ATAC_W_LIST))
  ))]
  setorder(combined, trait, curve, p_drop)
  fwrite(combined, COMBINED_SUMMARY)
  fwrite(combined, COMBINED_PLOT_DATA)

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("Package ggplot2 is not installed; combined CSVs were written but plots were skipped.")
    return(invisible(TRUE))
  }
  p <- ggplot2::ggplot(
    combined,
    ggplot2::aes(x = p_drop, y = mean_auroc, color = curve, linetype = curve, group = curve)
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lo95, ymax = hi95, fill = curve),
      alpha = 0.12,
      color = NA
    ) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::facet_wrap(~ trait, nrow = 1) +
    ggplot2::scale_x_continuous(
      breaks = seq(0, 0.9, by = 0.1),
      labels = function(x) paste0(round(100 * x), "%")
    ) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(
      title = "Reverse robustness under Granja RNA dropout",
      x = "RNA read dropout",
      y = "AUROC",
      color = NULL,
      fill = NULL,
      linetype = NULL
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      strip.background = ggplot2::element_rect(fill = "grey92", color = "grey70"),
      strip.text = ggplot2::element_text(face = "bold")
    )
  ggplot2::ggsave(COMBINED_PNG, p, width = 10.5, height = 4.6, dpi = 320)
  ggplot2::ggsave(COMBINED_PDF, p, width = 10.5, height = 4.6)
  invisible(TRUE)
}

cat(">>> scTARF reverse robustness under Granja RNA dropout\n\n")
cat("Traits      : ", paste(ROBUST_TRAITS_TO_RUN, collapse = ", "), "\n", sep = "")
cat("P_DROP_LIST : ", paste(ROBUST_P_DROP_LIST, collapse = ", "), "\n", sep = "")
cat("SEED_LIST   : ", paste(ROBUST_SEED_LIST, collapse = ", "), "\n", sep = "")
cat("Weights     : ", paste(ROBUST_ATAC_W_LIST, collapse = ", "), "\n", sep = "")
cat("Output      : ", SCTARF_OUT_DIR, "\n\n", sep = "")

required <- c(PYTHON_EXE, BASE_H5AD, BASE_META_RDS, DROP_HELPER, COUNT_HELPER)
missing <- required[!file.exists(required)]
if (length(missing) > 0L) {
  stop("Missing required file(s):\n  ", paste(missing, collapse = "\n  "), call. = FALSE)
}
if (any(!ROBUST_TRAITS_TO_RUN %in% names(TRAIT_IDS))) {
  stop("Unknown trait(s): ", paste(setdiff(ROBUST_TRAITS_TO_RUN, names(TRAIT_IDS)), collapse = ", "),
       call. = FALSE)
}

manifest_rows <- list()

for (p_drop in ROBUST_P_DROP_LIST) {
  for (seed in ROBUST_SEED_LIST) {
    paths <- condition_paths(p_drop, seed)
    dir.create(paths$run_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(paths$inputs_dir, recursive = TRUE, showWarnings = FALSE)

    cat("\n", strrep("=", 78), "\n", sep = "")
    cat("[condition] ", paths$run_id, "\n", sep = "")
    cat(strrep("=", 78), "\n", sep = "")

    missing_before <- condition_missing_items(paths)
    h5ad_file <- available_h5ad(paths)
    if (isTRUE(RESUME_SKIP_COMPLETED) && !FORCE_DROPOUT_H5AD && !FORCE_CONDITION_META &&
        !FORCE_REMAP && !FORCE_RETRANSFER && !FORCE_REBUILD_RNA_MAIN_Z &&
        length(missing_before) == 0L) {
      cat("  [skip] condition already complete.\n")
      manifest_rows[[paths$run_id]] <- manifest_row(
        paths, p_drop, seed, h5ad_file,
        status = "skipped_complete",
        missing_items = character()
      )
      fwrite(rbindlist(manifest_rows, use.names = TRUE, fill = TRUE), OUT_MANIFEST)
      next
    }
    if (length(missing_before) > 0L) {
      cat("  [resume] missing: ", paste(missing_before, collapse = ", "), "\n", sep = "")
    }

    h5ad_file <- ensure_dropout_h5ad(paths, p_drop, seed)
    meta_rds <- ensure_condition_meta(paths, h5ad_file)

    missing_pipeline <- setdiff(condition_missing_items(paths), "condition_eval")
    if (length(missing_pipeline) > 0L || !check_file(paths$auroc_csv, 1L) ||
        FORCE_REMAP || FORCE_RETRANSFER || FORCE_REBUILD_RNA_MAIN_Z) {
      run_condition_pipeline(paths, h5ad_file, meta_rds)
    } else {
      cat("  [pipeline] AUROC exists; only refreshing condition eval.\n")
    }

    write_condition_eval(paths, p_drop, seed)
    missing_after <- condition_missing_items(paths)
    manifest_rows[[paths$run_id]] <- manifest_row(
      paths, p_drop, seed, h5ad_file,
      status = if (length(missing_after) == 0L) "completed" else "incomplete",
      missing_items = missing_after
    )
    fwrite(rbindlist(manifest_rows, use.names = TRUE, fill = TRUE), OUT_MANIFEST)
  }
}

condition_eval_files <- list.files(
  RUNS_DIR,
  pattern = "^sctarf_reverse_condition_eval[.]csv$",
  recursive = TRUE,
  full.names = TRUE
)
stress_eval <- if (length(condition_eval_files) > 0L) {
  rbindlist(lapply(condition_eval_files, fread), use.names = TRUE, fill = TRUE)
} else {
  data.table()
}
all_eval <- rbindlist(list(clean_sctarf_anchor(), stress_eval), use.names = TRUE, fill = TRUE)
setorder(all_eval, trait, atac_w, p_drop, seed)
sctarf_summary <- summarize_sctarf(all_eval)

fwrite(all_eval, OUT_EVAL)
fwrite(sctarf_summary, OUT_SUMMARY)
fwrite(sctarf_summary, OUT_PLOT_DATA)
write_combined_outputs(sctarf_summary)

cat("\n=== scTARF reverse robustness finished ===\n")
cat("Wrote:\n")
cat("  ", OUT_EVAL, "\n", sep = "")
cat("  ", OUT_SUMMARY, "\n", sep = "")
cat("  ", OUT_PLOT_DATA, "\n", sep = "")
cat("  ", OUT_MANIFEST, "\n", sep = "")
cat("  ", COMBINED_SUMMARY, "\n", sep = "")
cat("  ", COMBINED_PLOT_DATA, "\n", sep = "")
cat("  ", COMBINED_PNG, "\n", sep = "")
cat("\nSummary:\n")
print(sctarf_summary[, .(
  trait,
  atac_w,
  p_drop,
  mean_auroc = round(mean_auroc, 6),
  lo95 = round(lo95, 6),
  hi95 = round(hi95, 6),
  n_runs
)])
