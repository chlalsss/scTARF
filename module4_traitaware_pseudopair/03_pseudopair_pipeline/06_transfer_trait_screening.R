# module4 trait screening: transfer CITE RNA information for MCH/RDW/RBC
# -----------------------------------------------------------------------------
# Thin wrapper around:
#   03_pseudopair_pipeline/01_transfer_dual.R
#
# Run in RStudio:
#   setwd("D:/RStudio/work/module4_traitaware_pseudopair")
#   source("03_pseudopair_pipeline/06_transfer_trait_screening.R")
# -----------------------------------------------------------------------------

ROOT <- "D:/RStudio/work/module4_traitaware_pseudopair"

if (!exists("TRANSFER_MAP_DIR")) {
  TRANSFER_MAP_DIR <- file.path(ROOT, "02_trait_mapping/output/trait_screening")
}
if (!exists("TRANSFER_OUT_DIR")) {
  TRANSFER_OUT_DIR <- file.path(ROOT, "03_pseudopair_pipeline/output/transfer_trait_screening")
}

if (!exists("CITE_MTX")) CITE_MTX <- file.path(TRANSFER_MAP_DIR, "cite_extract/cite_trait_counts.mtx")
if (!exists("CITE_GENES")) CITE_GENES <- file.path(TRANSFER_MAP_DIR, "cite_extract/cite_trait_counts_genes.txt")
if (!exists("CITE_CELLS")) CITE_CELLS <- file.path(TRANSFER_MAP_DIR, "cite_extract/cite_trait_counts_cells.txt")
if (!exists("CITE_SCORE_RDS")) {
  CITE_SCORE_RDS <- file.path(
    ROOT,
    "01_data_prep_external_rna/output/trait_screening/cite_rna_scdrs_scores_trait_screening.rds"
  )
}

if (!exists("TRAITS_ALL")) TRAITS_ALL <- c("MCH", "RDW", "RBC")
if (!exists("TRAITS_TO_RUN")) TRAITS_TO_RUN <- TRAITS_ALL

cat(">>> Run trait-screening transfer via existing transfer script\n\n")
cat("Traits: ", paste(TRAITS_TO_RUN, collapse = ", "), "\n", sep = "")
cat("Mapping dir : ", TRANSFER_MAP_DIR, "\n", sep = "")
cat("Transfer out: ", TRANSFER_OUT_DIR, "\n\n", sep = "")

source(file.path(ROOT, "03_pseudopair_pipeline/01_transfer_dual.R"))
