# module4 trait screening: prepare module3 inputs and replace seeds
# -----------------------------------------------------------------------------
# This reuses:
#   03_pseudopair_pipeline/02_prepare_module3_v2_inputs.R
#
# Then replaces seed_v2_<trait>.rds with neighborhood-FDR mapping evidence:
#   candidate_group_weight_ge0p80
#   score = sig_group_source_weight for selected cells, 0 otherwise
#
# Run in RStudio:
#   setwd("D:/RStudio/work/module4_traitaware_pseudopair")
#   source("03_pseudopair_pipeline/07_prepare_trait_screening_module3_inputs.R")
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

ROOT <- "D:/RStudio/work/module4_traitaware_pseudopair"

if (!exists("MODULE3_TRANSFER_DIR")) {
  MODULE3_TRANSFER_DIR <- file.path(ROOT, "03_pseudopair_pipeline/output/transfer_trait_screening")
}
if (!exists("MODULE3_COMPAT_OUT_DIR")) {
  MODULE3_COMPAT_OUT_DIR <- file.path(ROOT, "03_pseudopair_pipeline/output/module3_v2_trait_screening_ge0p80")
}
if (!exists("TRAITS_TO_RUN")) TRAITS_TO_RUN <- c("MCH", "RDW", "RBC")

PY_HELPER <- file.path(ROOT, "03_pseudopair_pipeline/export_lsi_npy.py")
PYTHON_EXE <- "D:/python/python.exe"

if (!exists("SEED_CELL_FILE")) {
  SEED_CELL_FILE <- file.path(
    ROOT,
    "02_trait_mapping/output/trait_screening/diagnostics/trait_screening_atac_candidate_cell_level.csv.gz"
  )
}
if (!exists("SEED_SUMMARY_FILE")) {
  SEED_SUMMARY_FILE <- file.path(MODULE3_COMPAT_OUT_DIR, "seed_summary_neighborhood_fdr_ge0p80.csv")
}

cat(">>> Prepare module3 inputs for trait screening\n\n")
cat("Traits: ", paste(TRAITS_TO_RUN, collapse = ", "), "\n", sep = "")
cat("Transfer dir: ", MODULE3_TRANSFER_DIR, "\n", sep = "")
cat("Output dir  : ", MODULE3_COMPAT_OUT_DIR, "\n\n", sep = "")

source(file.path(ROOT, "03_pseudopair_pipeline/02_prepare_module3_v2_inputs.R"), local = TRUE)

cat("\n>>> Replace top1% pseudo-zRNA seeds with neighborhood-FDR ge0p80 seeds\n\n")
if (!file.exists(SEED_CELL_FILE)) {
  stop("Missing seed candidate cell-level file: ", SEED_CELL_FILE, call. = FALSE)
}
seed_dt <- fread(SEED_CELL_FILE)
required <- c(
  "trait", "atac_barcode", "sig_group_source_weight",
  "candidate_group_weight_ge0p80", "no_celltype_used_to_select"
)
if (!all(required %in% names(seed_dt))) {
  stop("Seed cell file missing required columns: ",
       paste(setdiff(required, names(seed_dt)), collapse = ", "),
       call. = FALSE)
}
if (!all(seed_dt$no_celltype_used_to_select %in% TRUE)) {
  stop("Seed candidate file indicates celltype may have been used; aborting.",
       call. = FALSE)
}

summary_rows <- list()
for (trait_name in TRAITS_TO_RUN) {
  rna_z_file <- file.path(MODULE3_COMPAT_OUT_DIR, paste0("rna_z_block_v2_", trait_name, ".rds"))
  if (!file.exists(rna_z_file)) stop("Missing rna_z file after module3 prep: ", rna_z_file, call. = FALSE)
  reference_cells <- rownames(readRDS(rna_z_file))
  dt <- seed_dt[trait == trait_name]
  dt <- dt[match(reference_cells, atac_barcode)]
  if (!identical(dt$atac_barcode, reference_cells)) {
    stop("Seed candidate cell order does not match module3 reference cells for ", trait_name,
         call. = FALSE)
  }

  candidate <- as.logical(dt$candidate_group_weight_ge0p80)
  candidate[is.na(candidate)] <- FALSE
  score <- dt$sig_group_source_weight
  score[!is.finite(score) | is.na(score)] <- 0
  score[!candidate] <- 0
  names(score) <- reference_cells

  if (sum(score > 0) == 0L) {
    stop("No positive ge0p80 seeds for ", trait_name, call. = FALSE)
  }
  saveRDS(score, file.path(MODULE3_COMPAT_OUT_DIR, paste0("seed_v2_", trait_name, ".rds")))

  summary_rows[[trait_name]] <- data.table(
    trait = trait_name,
    n_cells = length(score),
    n_seed_positive = sum(score > 0),
    seed_rule = "neighborhood_fdr_group_weight_ge0p80",
    score_min_positive = min(score[score > 0]),
    score_median_positive = median(score[score > 0]),
    score_max = max(score),
    no_celltype_used = TRUE
  )
  cat("  ", trait_name, ": n_seed_positive=", sum(score > 0),
      "; median=", signif(median(score[score > 0]), 4),
      "; max=", signif(max(score), 4), "\n", sep = "")
}

seed_summary <- rbindlist(summary_rows, use.names = TRUE, fill = TRUE)
fwrite(seed_summary, SEED_SUMMARY_FILE)

cat("\nWrote seed summary: ", SEED_SUMMARY_FILE, "\n", sep = "")
