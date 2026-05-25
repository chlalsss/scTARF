# module4 robustness: prepare one shared module3-compatible input directory
# -----------------------------------------------------------------------------
# Purpose:
#   The top 3 traits come from two completed module4 branches:
#
#     MONO  -> 03_pseudopair_pipeline/output/module3_v2_neighborhood_fdr_ge0p80
#     RDW   -> 03_pseudopair_pipeline/output/module3_v2_trait_screening_ge0p80
#     RBC   -> 03_pseudopair_pipeline/output/module3_v2_trait_screening_ge0p80
#
#   This script gathers their existing module3 inputs into one directory:
#
#     05_robustness/output/top3_embedding_stress/base_module3_inputs
#
#   It does not recompute mapping, transfer, PCA, KNN, RWR, or AUROC.
#
# Run in RStudio:
#   setwd("D:/RStudio/work/module4_traitaware_pseudopair")
#   source("05_robustness/01_prepare_top3_module3_inputs.R")
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

if (!exists("ROOT")) ROOT <- "D:/RStudio/work/module4_traitaware_pseudopair"
if (!exists("ROBUST_DIR")) {
  ROBUST_DIR <- file.path(ROOT, "05_robustness/output/top3_embedding_stress")
}
if (!exists("BASE_OUT_DIR")) {
  BASE_OUT_DIR <- file.path(ROBUST_DIR, "base_module3_inputs")
}
dir.create(BASE_OUT_DIR, recursive = TRUE, showWarnings = FALSE)

if (!exists("SELECTED_TRAITS")) SELECTED_TRAITS <- c("MONO", "RDW", "RBC")

if (!exists("MONO_DIR")) {
  MONO_DIR <- file.path(ROOT, "03_pseudopair_pipeline/output/module3_v2_neighborhood_fdr_ge0p80")
}
if (!exists("SCREEN_DIR")) {
  SCREEN_DIR <- file.path(ROOT, "03_pseudopair_pipeline/output/module3_v2_trait_screening_ge0p80")
}

if (!exists("TRAIT_SOURCE_DIR")) {
  TRAIT_SOURCE_DIR <- c(
    MONO = MONO_DIR,
    RDW = SCREEN_DIR,
    RBC = SCREEN_DIR
  )
}

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

cat(">>> Prepare top3 module3-compatible inputs for robustness\n\n")
cat("Selected traits: ", paste(SELECTED_TRAITS, collapse = ", "), "\n", sep = "")
cat("Output dir     : ", BASE_OUT_DIR, "\n\n", sep = "")

base_trait <- SELECTED_TRAITS[[1]]
base_atac_file <- file.path(TRAIT_SOURCE_DIR[[base_trait]], "atac_z_block_v2.rds")
stop_missing(base_atac_file)
base_atac <- readRDS(base_atac_file)

copy_required(base_atac_file, file.path(BASE_OUT_DIR, "atac_z_block_v2.rds"))

manifest <- list()
for (trait in SELECTED_TRAITS) {
  src <- TRAIT_SOURCE_DIR[[trait]]
  rna_file <- file.path(src, paste0("rna_z_block_v2_", trait, ".rds"))
  seed_file <- file.path(src, paste0("seed_v2_", trait, ".rds"))
  stop_missing(c(rna_file, seed_file))

  rna_z <- readRDS(rna_file)
  seed <- readRDS(seed_file)
  if (!identical(rownames(base_atac), rownames(rna_z))) {
    stop("ATAC and RNA rownames differ for ", trait, call. = FALSE)
  }
  if (!identical(rownames(base_atac), names(seed))) {
    stop("ATAC and seed names differ for ", trait, call. = FALSE)
  }

  copy_required(rna_file, file.path(BASE_OUT_DIR, paste0("rna_z_block_v2_", trait, ".rds")))
  copy_required(seed_file, file.path(BASE_OUT_DIR, paste0("seed_v2_", trait, ".rds")))

  manifest[[trait]] <- data.table(
    trait = trait,
    source_dir = src,
    n_cells = nrow(base_atac),
    atac_dim = ncol(base_atac),
    rna_dim = ncol(rna_z),
    n_seed_positive = sum(seed > 0),
    seed_min_positive = min(seed[seed > 0]),
    seed_median_positive = median(seed[seed > 0]),
    seed_max = max(seed)
  )
}

manifest_dt <- rbindlist(manifest, use.names = TRUE, fill = TRUE)
fwrite(manifest_dt, file.path(BASE_OUT_DIR, "top3_input_manifest.csv"))

cat("Prepared inputs:\n")
print(manifest_dt)
cat("\nOutputs are in: ", BASE_OUT_DIR, "\n", sep = "")
