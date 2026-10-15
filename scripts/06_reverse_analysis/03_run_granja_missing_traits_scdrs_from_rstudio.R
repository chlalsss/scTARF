source("scripts/_paths.R")
# RStudio entry point: resume official scDRS for Granja missing traits.
#
# Usage in RStudio:
#   Run from the repository root with optional SCTARF_* path overrides.
#   source("reverse_atac_to_rna_experiment_20260608/scripts/run_granja_missing_traits_scdrs_from_rstudio.R")
#
# Behavior:
#   1. Runs only the official scDRS CLI.
#   2. Runs each trait separately so completed .score.gz/.full_score.gz files
#      are skipped on rerun.
#   3. Runs perform-downstream only if group-analysis files are missing.
#   4. Summarizes Granja missing-trait scores and refreshes the extended AUROC table.

if (!exists("ROOT")) ROOT <- Sys.getenv("SCTARF_WORK_ROOT", normalizePath(".", winslash = "/", mustWork = FALSE))
if (!exists("MODULE4_ROOT")) MODULE4_ROOT <- Sys.getenv("SCTARF_NONPAIRED_ROOT", file.path(ROOT, "outputs", "nonpaired"))
if (!exists("REV_ROOT")) REV_ROOT <- Sys.getenv("SCTARF_REVERSE_ROOT", file.path(ROOT, "outputs", "reverse"))
LOCAL_SCRIPT_DIR <- dirname(tryCatch(normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = FALSE),
                                      error = function(e) "scripts/06_reverse_analysis/03_run_granja_missing_traits_scdrs_from_rstudio.R"))

PYTHON_EXE <- sctarf_python()
SCDRS_CLI <- sctarf_scdrs()
H5AD <- file.path(
  MODULE4_ROOT,
  "01_data_prep/03_rna_scdrs_python/input/granja_rna_raw_counts.h5ad"
)
GS_MASTER <- file.path(
  MODULE4_ROOT,
  "01_data_prep_external_rna/03_cite_scdrs_python/cite_trait_screening.gs"
)
OUT_DIR <- file.path(
  MODULE4_ROOT,
  "01_data_prep/03_rna_scdrs_python/output_missing_traits"
)
GS_SPLIT_DIR <- file.path(
  MODULE4_ROOT,
  "01_data_prep/03_rna_scdrs_python/missing_trait_gs"
)
LOG_DIR <- file.path(REV_ROOT, "logs")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(GS_SPLIT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

N_CTRL <- 1000L

TRAITS <- c(
  EOS = "UKB_460K.blood_EOSINOPHIL_COUNT",
  MCH = "UKB_460K.blood_MEAN_CORPUSCULAR_HEMOGLOBIN",
  RDW = "UKB_460K.blood_RBC_DISTRIB_WIDTH",
  RBC = "UKB_460K.blood_RED_COUNT",
  WBC = "UKB_460K.blood_WHITE_COUNT"
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
    return(file.path(normalizePath(parent, winslash = "\\", mustWork = TRUE), basename(path)))
  }
  gsub("/", "\\\\", path, fixed = TRUE)
}

run_or_stop <- function(cmd, args, log_prefix) {
  stdout_log <- file.path(LOG_DIR, paste0(log_prefix, ".stdout.log"))
  stderr_log <- file.path(LOG_DIR, paste0(log_prefix, ".stderr.log"))
  cat("Running:\n  ", cmd, " ", paste(args, collapse = " "), "\n", sep = "")
  status <- system2(cmd, args = shQuote(args), stdout = stdout_log, stderr = stderr_log)
  if (!identical(status, 0L)) {
    cat("\nSTDOUT log: ", stdout_log, "\n", sep = "")
    cat("STDERR log: ", stderr_log, "\n", sep = "")
    stop("Command failed with status ", status, call. = FALSE)
  }
  cat("  OK. Logs:\n")
  cat("    ", stdout_log, "\n", sep = "")
  cat("    ", stderr_log, "\n", sep = "")
}

write_one_trait_gs <- function(trait_id) {
  lines <- readLines(GS_MASTER, warn = FALSE)
  header <- lines[1]
  hit <- grep(paste0("^", trait_id, "\t"), lines, value = TRUE)
  if (length(hit) != 1L) {
    stop("Expected exactly one gene-set row for ", trait_id,
         "; found ", length(hit), call. = FALSE)
  }
  out_gs <- file.path(GS_SPLIT_DIR, paste0(trait_id, ".gs"))
  writeLines(c(header, hit), out_gs, useBytes = TRUE)
  out_gs
}

cat(">>> Resume official scDRS for Granja missing traits\n\n")
for (path in c(PYTHON_EXE, SCDRS_CLI, H5AD, GS_MASTER)) {
  if (!file.exists(path)) stop("Missing required file: ", path, call. = FALSE)
}

Sys.setenv(NUMBA_CACHE_DIR = file.path(MODULE4_ROOT, ".numba_cache"))
dir.create(Sys.getenv("NUMBA_CACHE_DIR"), recursive = TRUE, showWarnings = FALSE)

for (short in names(TRAITS)) {
  trait_id <- TRAITS[[short]]
  score_file <- file.path(OUT_DIR, paste0(trait_id, ".score.gz"))
  full_file <- file.path(OUT_DIR, paste0(trait_id, ".full_score.gz"))

  cat("\n", strrep("=", 78), "\n", sep = "")
  cat("[compute-score] ", short, " / ", trait_id, "\n", sep = "")
  cat(strrep("=", 78), "\n", sep = "")

  if (check_file(score_file, 1000L) && check_file(full_file, 1000L)) {
    cat("Skip compute-score; existing files found:\n")
    cat("  ", score_file, " (", round(file.info(score_file)$size / 1e6, 2), " MB)\n", sep = "")
    cat("  ", full_file, " (", round(file.info(full_file)$size / 1e6, 2), " MB)\n", sep = "")
    next
  }

  one_trait_gs <- write_one_trait_gs(trait_id)
  run_or_stop(
    win_path(PYTHON_EXE),
    args = c(
      win_path(SCDRS_CLI),
      "compute-score",
      "--h5ad-file", win_path(H5AD),
      "--h5ad-species", "human",
      "--gs-file", win_path(one_trait_gs),
      "--gs-species", "human",
      "--out-folder", win_path(OUT_DIR),
      "--n-ctrl", as.character(N_CTRL)
    ),
    log_prefix = paste0("granja_missing_", short, "_compute_score")
  )
}

expected_group <- file.path(OUT_DIR, paste0(TRAITS, ".scdrs_group.cell_type"))
cat("\n", strrep("=", 78), "\n", sep = "")
cat("[perform-downstream] cell_type group analysis\n")
cat(strrep("=", 78), "\n", sep = "")

if (all(vapply(expected_group, check_file, logical(1), min_size = 1L))) {
  cat("Skip perform-downstream; all group-analysis files already exist.\n")
} else {
  missing_group <- expected_group[!vapply(expected_group, check_file, logical(1), min_size = 1L)]
  cat("Missing group-analysis files:\n")
  cat(paste0("  ", missing_group, collapse = "\n"), "\n", sep = "")

  full_files <- file.path(OUT_DIR, paste0(TRAITS, ".full_score.gz"))
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
      "--h5ad-file", win_path(H5AD),
      "--score-file", win_path(file.path(OUT_DIR, "@.full_score.gz"), mustWork = FALSE),
      "--out-folder", win_path(OUT_DIR),
      "--group-analysis", "cell_type"
    ),
    log_prefix = "granja_missing_perform_downstream"
  )
}

cat("\n>>> Summarizing new Granja scores and refreshing AUROC table\n")
source(file.path(LOCAL_SCRIPT_DIR, "05_summarize_granja_missing_trait_scdrs.R"))
source(file.path(LOCAL_SCRIPT_DIR, "02_step0_extended_traits_auroc.R"))

cat("\nAll done. Updated table:\n")
cat(file.path(REV_ROOT, "outputs/step0b_extended_rna_only_scdrs_baseline_auroc_wide.csv"), "\n")
