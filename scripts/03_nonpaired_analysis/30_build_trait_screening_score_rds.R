# module4 trait screening: build CITE score RDS for transfer
# -----------------------------------------------------------------------------
# Purpose:
#   Existing transfer code expects one RDS containing:
#     cell, scdrs_input_total_counts, <TRAIT>_norm_score
#
#   The screening scDRS run wrote separate .score.gz files. This adapter merges
#   them with the existing CITE metadata so transfer can reuse 01_transfer_dual.R.
#
# Run from the repository root with optional SCTARF_* path overrides.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

source(file.path("scripts", "_paths.R"))
ROOT <- sctarf_nonpaired_root()
if (!exists("BASE_RDS")) {
  BASE_RDS <- file.path(ROOT, "01_data_prep_external_rna/output/cite_rna_scdrs_scores.rds")
}
if (!exists("SCORE_DIR")) {
  SCORE_DIR <- file.path(ROOT, "01_data_prep_external_rna/03_cite_scdrs_python/output_trait_screening")
}
if (!exists("OUT_RDS")) {
  OUT_RDS <- file.path(ROOT, "01_data_prep_external_rna/output/trait_screening/cite_rna_scdrs_scores_trait_screening.rds")
}
if (!exists("OUT_CSV")) {
  OUT_CSV <- file.path(ROOT, "01_data_prep_external_rna/output/trait_screening/cite_rna_scdrs_scores_trait_screening.csv.gz")
}
dir.create(dirname(OUT_RDS), recursive = TRUE, showWarnings = FALSE)

if (!exists("TRAIT_IDS")) {
  TRAIT_IDS <- c(
    MCH = "UKB_460K.blood_MEAN_CORPUSCULAR_HEMOGLOBIN",
    RDW = "UKB_460K.blood_RBC_DISTRIB_WIDTH",
    RBC = "UKB_460K.blood_RED_COUNT"
  )
}

score_col <- function(dt) {
  hit <- intersect(c("norm_score", "zscore", "score"), names(dt))[1]
  if (is.na(hit)) stop("No score column found in score file.", call. = FALSE)
  hit
}

cell_col <- function(dt) {
  hit <- intersect(c("cell", "cell_id", "barcode", "V1"), names(dt))[1]
  if (is.na(hit)) hit <- names(dt)[1]
  hit
}

cat(">>> Build trait-screening CITE score RDS\n\n")

if (!file.exists(BASE_RDS)) stop("Missing base CITE RDS: ", BASE_RDS, call. = FALSE)
base <- as.data.table(readRDS(BASE_RDS))
needed <- c("cell", "scdrs_input_total_counts")
if (!all(needed %in% names(base))) {
  stop("Base RDS must contain: ", paste(needed, collapse = ", "), call. = FALSE)
}
out <- copy(base)

for (short in names(TRAIT_IDS)) {
  trait_id <- TRAIT_IDS[[short]]
  f <- file.path(SCORE_DIR, paste0(trait_id, ".score.gz"))
  if (!file.exists(f)) stop("Missing screening score file: ", f, call. = FALSE)
  dt <- fread(f)
  sc <- score_col(dt)
  cc <- cell_col(dt)
  sdt <- data.table(cell = as.character(dt[[cc]]), score = as.numeric(dt[[sc]]))
  setnames(sdt, "score", paste0(short, "_norm_score"))
  out <- merge(out, sdt, by = "cell", all.x = TRUE, sort = FALSE)
  cat("  added ", short, "_norm_score from ", basename(f), "\n", sep = "")
}

score_cols <- paste0(names(TRAIT_IDS), "_norm_score")
missing_scores <- colSums(is.na(out[, ..score_cols]))
if (any(missing_scores > 0L)) {
  cat("\nNote: some CITE cells are absent from screening scDRS score files.\n")
  cat("These scores will be set to 0 so the cells can still be used for RNA profile transfer,\n")
  cat("but they do not contribute positive trait-score evidence.\n")
  print(missing_scores[missing_scores > 0L])
  for (col in score_cols) {
    set(out, i = which(is.na(out[[col]])), j = col, value = 0)
  }
}

saveRDS(as.data.frame(out), OUT_RDS)
fwrite(out, OUT_CSV)

cat("\nWrote:\n")
cat("  ", OUT_RDS, "\n", sep = "")
cat("  ", OUT_CSV, "\n", sep = "")
