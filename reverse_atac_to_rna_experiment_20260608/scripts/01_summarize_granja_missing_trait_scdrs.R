# Summarize Granja RNA missing-trait official scDRS outputs.
#
# This script does not compute scDRS. It reads output_missing_traits/*.score.gz
# produced by official scDRS and writes compact RDS/CSV score tables.

suppressPackageStartupMessages({
  library(data.table)
})

ROOT <- "D:/RStudio/work"
MODULE4_ROOT <- file.path(ROOT, "module4_traitaware_pseudopair")
SCDRS_OUT <- file.path(
  MODULE4_ROOT,
  "01_data_prep/03_rna_scdrs_python/output_missing_traits"
)
BASE_SCORE_RDS <- file.path(
  MODULE4_ROOT,
  "01_data_prep/output/granja_rna_scdrs_scores.rds"
)
META_RDS <- file.path(
  MODULE4_ROOT,
  "01_data_prep/output/granja_rna_cell_meta.rds"
)
SUMMARY_OUT <- file.path(MODULE4_ROOT, "01_data_prep/output")

TRAITS <- c(
  EOS = "UKB_460K.blood_EOSINOPHIL_COUNT",
  MCH = "UKB_460K.blood_MEAN_CORPUSCULAR_HEMOGLOBIN",
  RDW = "UKB_460K.blood_RBC_DISTRIB_WIDTH",
  RBC = "UKB_460K.blood_RED_COUNT",
  WBC = "UKB_460K.blood_WHITE_COUNT"
)

score_col <- function(dt) {
  candidates <- c("norm_score", "zscore", "score")
  hit <- intersect(candidates, names(dt))[1]
  if (is.na(hit)) stop("No score column found in: ", paste(names(dt), collapse = ", "))
  hit
}

aux_cols <- function(dt) {
  intersect(c("pval", "pval_ctrl", "pval_mc", "fdr"), names(dt))
}

cat(">>> Summarize Granja missing-trait official scDRS outputs\n\n")
stopifnot(file.exists(BASE_SCORE_RDS), file.exists(META_RDS))

meta <- readRDS(META_RDS)
if (is.null(rownames(meta))) {
  if (!("barcode" %in% colnames(meta))) {
    stop("Granja metadata lacks rownames and barcode column.")
  }
  rownames(meta) <- as.character(meta$barcode)
}

missing_scores <- data.table(
  barcode = rownames(meta),
  cell_type = as.character(meta$cell_type)
)

for (short in names(TRAITS)) {
  trait_id <- TRAITS[[short]]
  f <- file.path(SCDRS_OUT, paste0(trait_id, ".score.gz"))
  if (!file.exists(f)) {
    stop("Missing official scDRS score file: ", f)
  }
  cat("[trait] ", short, " <- ", basename(f), "\n", sep = "")
  dt <- fread(f)
  bc_col <- intersect(c("cell", "cell_id", "barcode", "V1"), names(dt))[1]
  if (is.na(bc_col)) bc_col <- names(dt)[1]
  sc_col <- score_col(dt)

  score_dt <- dt[, .(
    barcode = as.character(get(bc_col)),
    score = as.numeric(get(sc_col))
  )]
  setnames(score_dt, "score", paste0(short, "_norm_score"))
  missing_scores <- merge(
    missing_scores,
    score_dt,
    by = "barcode",
    all.x = TRUE,
    sort = FALSE
  )

  for (ac in aux_cols(dt)) {
    aux <- dt[, .(
      barcode = as.character(get(bc_col)),
      value = get(ac)
    )]
    setnames(aux, "value", paste0(short, "_", ac))
    missing_scores <- merge(
      missing_scores,
      aux,
      by = "barcode",
      all.x = TRUE,
      sort = FALSE
    )
  }
}

base_scores <- as.data.table(readRDS(BASE_SCORE_RDS))
base_scores[, barcode := as.character(barcode)]
score_cols_to_add <- grep("_norm_score$|_pval$|_fdr$|_pval_ctrl$|_pval_mc$", names(missing_scores), value = TRUE)
combined <- merge(
  base_scores,
  missing_scores[, c("barcode", score_cols_to_add), with = FALSE],
  by = "barcode",
  all.x = TRUE,
  sort = FALSE
)

missing_rds <- file.path(SUMMARY_OUT, "granja_rna_scdrs_scores_missing_traits.rds")
missing_csv <- file.path(SUMMARY_OUT, "granja_rna_scdrs_scores_missing_traits.csv.gz")
combined_rds <- file.path(SUMMARY_OUT, "granja_rna_scdrs_scores_with_missing_traits.rds")
combined_csv <- file.path(SUMMARY_OUT, "granja_rna_scdrs_scores_with_missing_traits.csv.gz")

saveRDS(missing_scores, missing_rds)
fwrite(missing_scores, missing_csv)
saveRDS(combined, combined_rds)
fwrite(combined, combined_csv)

cat("\nWrote:\n")
cat("  ", missing_rds, "\n", sep = "")
cat("  ", missing_csv, "\n", sep = "")
cat("  ", combined_rds, "\n", sep = "")
cat("  ", combined_csv, "\n", sep = "")
cat("\nGranja missing-trait scDRS summary OK.\n")
