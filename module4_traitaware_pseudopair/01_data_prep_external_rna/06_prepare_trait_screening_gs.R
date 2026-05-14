# module4 trait screening: prepare candidate blood-trait gene set for CITE scDRS
# -----------------------------------------------------------------------------
# Purpose:
#   Build a small .gs file for the next screening round from the existing
#   module3 74-trait MAGMA/scDRS gene-set file.
#
# This does not run scDRS. It only prepares:
#   01_data_prep_external_rna/03_cite_scdrs_python/cite_trait_screening.gs
#   01_data_prep_external_rna/output/trait_screening_available_traits.csv
#
# Run in RStudio:
#   setwd("D:/RStudio/work/module4_traitaware_pseudopair")
#   source("01_data_prep_external_rna/06_prepare_trait_screening_gs.R")
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

ROOT <- "D:/RStudio/work/module4_traitaware_pseudopair"
SOURCE_GS <- "D:/RStudio/work/module3_multiome/04_dual_modality/01_scdrs_python/magma_10kb_top1000_zscore.74_traits.rv1.gs"
CITE_SCDRS_DIR <- file.path(ROOT, "01_data_prep_external_rna/03_cite_scdrs_python")
OUT_GS <- file.path(CITE_SCDRS_DIR, "cite_trait_screening.gs")
OUT_DIR <- file.path(ROOT, "01_data_prep_external_rna/output")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

if (!exists("SCREENING_TRAIT_IDS")) {
  SCREENING_TRAIT_IDS <- c(
    "UKB_460K.blood_EOSINOPHIL_COUNT",
    "UKB_460K.blood_WHITE_COUNT",
    "UKB_460K.blood_RED_COUNT",
    "UKB_460K.blood_MEAN_CORPUSCULAR_HEMOGLOBIN",
    "UKB_460K.blood_RBC_DISTRIB_WIDTH"
  )
}

short_name <- function(trait_id) {
  x <- sub("^UKB_460K\\.blood_", "", trait_id)
  x <- gsub("_COUNT$", "", x)
  x <- gsub("MEAN_CORPUSCULAR_HEMOGLOBIN", "MCH", x)
  x <- gsub("RBC_DISTRIB_WIDTH", "RDW", x)
  x <- gsub("RED", "RBC", x)
  x
}

cat(">>> Prepare trait-screening .gs for CITE scDRS\n\n")
cat("Source .gs: ", SOURCE_GS, "\n", sep = "")
cat("Output .gs: ", OUT_GS, "\n\n", sep = "")

if (!file.exists(SOURCE_GS)) {
  stop("Missing source .gs file: ", SOURCE_GS, call. = FALSE)
}

gs <- fread(SOURCE_GS)
if (!all(c("TRAIT", "GENESET") %in% names(gs))) {
  stop("Source .gs must contain TRAIT and GENESET columns.", call. = FALSE)
}

available <- data.table(
  trait = SCREENING_TRAIT_IDS,
  short = vapply(SCREENING_TRAIT_IDS, short_name, character(1)),
  available_in_source_gs = SCREENING_TRAIT_IDS %in% gs$TRAIT
)
available[, n_genes := NA_integer_]
for (i in seq_len(nrow(available))) {
  if (available$available_in_source_gs[i]) {
    geneset <- gs[TRAIT == available$trait[i], GENESET][1]
    available$n_genes[i] <- length(strsplit(geneset, ",", fixed = TRUE)[[1]])
  }
}

missing <- available[available_in_source_gs == FALSE, trait]
if (length(missing) > 0L) {
  warning("Some requested screening traits are absent from source .gs and will be skipped:\n  ",
          paste(missing, collapse = "\n  "))
}

keep <- gs[TRAIT %in% available[available_in_source_gs == TRUE, trait]]
if (nrow(keep) == 0L) {
  stop("No requested screening traits are available in source .gs.", call. = FALSE)
}

setorder(keep, TRAIT)
fwrite(keep, OUT_GS, sep = "\t")
fwrite(available, file.path(OUT_DIR, "trait_screening_available_traits.csv"))

cat("Available screening traits:\n")
print(available)
cat("\nWrote:\n")
cat("  ", OUT_GS, "\n", sep = "")
cat("  ", file.path(OUT_DIR, "trait_screening_available_traits.csv"), "\n", sep = "")
