# module4 robustness: ATAC embedding stress sweep for MONO/RDW/RBC
# -----------------------------------------------------------------------------
# Purpose:
#   This is the module4 robustness test that can be run with the files currently
#   available in this folder. Module3's original robustness script sparsifies raw
#   ATAC peak counts and recomputes LSI. Module4 currently keeps the Granja ATAC
#   LSI/module3 input embedding, not the raw peak-by-cell count matrix, so this
#   script does an embedding-level ATAC stress test instead.
#
# What is stressed:
#   atac_z_block_v2.rds
#
# What is reused unchanged:
#   rna_z_block_v2_<trait>.rds
#   seed_v2_<trait>.rds
#   module3 KNN/RWR/TRS code via 03_pseudopair_pipeline/03_run_module3_v2_knn_rwr.R
#
# Stress rule:
#   stressed_atac = scale(atac_z + N(0, noise_sd))
#
# Recommended first run:
#   NOISE_SD_LIST <- c(0.00, 0.05, 0.10)
#   SEED_LIST <- c(1)
#
# Full sweep:
#   NOISE_SD_LIST <- c(0.00, 0.05, 0.10, 0.20, 0.30)
#   SEED_LIST <- c(1, 2, 3)
#
# Run in RStudio:
#   setwd("D:/RStudio/work/module4_traitaware_pseudopair")
#   source("05_robustness/02_run_top3_embedding_stress.R")
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

ROOT <- "D:/RStudio/work/module4_traitaware_pseudopair"
ROBUST_DIR <- file.path(ROOT, "05_robustness/output/top3_embedding_stress")
BASE_INPUT_DIR <- file.path(ROBUST_DIR, "base_module3_inputs")
RUNS_DIR <- file.path(ROBUST_DIR, "runs")
dir.create(RUNS_DIR, recursive = TRUE, showWarnings = FALSE)

if (!exists("ROBUST_TRAITS")) ROBUST_TRAITS <- c("MONO", "RDW", "RBC")
if (!exists("NOISE_SD_LIST")) NOISE_SD_LIST <- c(0.00, 0.05, 0.10, 0.20, 0.30)
if (!exists("SEED_LIST")) SEED_LIST <- c(1L, 2L, 3L)
if (!exists("ATAC_W_LIST")) ATAC_W_LIST <- c(0.7, 1.0)

stop_missing <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0L) {
    stop("Missing required file(s):\n  ", paste(missing, collapse = "\n  "),
         call. = FALSE)
  }
}

format_noise <- function(x) {
  gsub("\\.", "p", sprintf("%.2f", x))
}

scale_no_attrs <- function(x) {
  z <- scale(x)
  attributes(z)$`scaled:center` <- NULL
  attributes(z)$`scaled:scale` <- NULL
  z
}

prepare_run_dir <- function(base_dir, run_dir, noise_sd, seed) {
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

  base_atac_file <- file.path(base_dir, "atac_z_block_v2.rds")
  stop_missing(base_atac_file)
  base_atac <- readRDS(base_atac_file)

  atac_out <- file.path(run_dir, "atac_z_block_v2.rds")
  if (!file.exists(atac_out)) {
    if (noise_sd == 0) {
      stressed <- base_atac
    } else {
      set.seed(seed)
      stressed <- base_atac + matrix(
        stats::rnorm(length(base_atac), mean = 0, sd = noise_sd),
        nrow = nrow(base_atac),
        ncol = ncol(base_atac),
        dimnames = dimnames(base_atac)
      )
      stressed <- scale_no_attrs(stressed)
    }
    saveRDS(stressed, atac_out)
  }

  for (trait in ROBUST_TRAITS) {
    for (prefix in c("rna_z_block_v2_", "seed_v2_")) {
      from <- file.path(base_dir, paste0(prefix, trait, ".rds"))
      to <- file.path(run_dir, paste0(prefix, trait, ".rds"))
      stop_missing(from)
      if (!file.exists(to)) {
        ok <- file.copy(from, to, overwrite = TRUE)
        if (!ok) stop("Failed to copy: ", from, call. = FALSE)
      }
    }
  }

  data.table(
    noise_sd = noise_sd,
    stress_seed = seed,
    run_dir = run_dir,
    n_cells = nrow(base_atac),
    atac_dim = ncol(base_atac)
  )
}

cat(">>> Run top3 ATAC embedding stress robustness\n\n")
cat("Traits      : ", paste(ROBUST_TRAITS, collapse = ", "), "\n", sep = "")
cat("NOISE_SD    : ", paste(NOISE_SD_LIST, collapse = ", "), "\n", sep = "")
cat("SEED_LIST   : ", paste(SEED_LIST, collapse = ", "), "\n", sep = "")
cat("ATAC_W_LIST : ", paste(ATAC_W_LIST, collapse = ", "), "\n", sep = "")
cat("Base inputs : ", BASE_INPUT_DIR, "\n\n", sep = "")

required <- c(
  file.path(BASE_INPUT_DIR, "atac_z_block_v2.rds"),
  file.path(BASE_INPUT_DIR, paste0("rna_z_block_v2_", ROBUST_TRAITS, ".rds")),
  file.path(BASE_INPUT_DIR, paste0("seed_v2_", ROBUST_TRAITS, ".rds"))
)
stop_missing(required)

run_rows <- list()
for (noise_sd in NOISE_SD_LIST) {
  for (seed in SEED_LIST) {
    run_id <- paste0("noise", format_noise(noise_sd), "_seed", seed)
    run_dir <- file.path(RUNS_DIR, run_id)
    cat("\n", strrep("=", 72), "\n", sep = "")
    cat("Run: ", run_id, "\n", sep = "")
    cat(strrep("=", 72), "\n", sep = "")

    run_rows[[run_id]] <- prepare_run_dir(BASE_INPUT_DIR, run_dir, noise_sd, seed)

    MODULE3_COMPAT_OUT_DIR <- run_dir
    TRAITS_TO_RUN <- ROBUST_TRAITS
    source(file.path(ROOT, "03_pseudopair_pipeline/03_run_module3_v2_knn_rwr.R"))
  }
}

run_manifest <- rbindlist(run_rows, use.names = TRUE, fill = TRUE)
fwrite(run_manifest, file.path(ROBUST_DIR, "embedding_stress_run_manifest.csv"))

cat("\n=== ATAC embedding stress KNN/RWR/TRS finished ===\n")
cat("Outputs are in: ", RUNS_DIR, "\n", sep = "")
