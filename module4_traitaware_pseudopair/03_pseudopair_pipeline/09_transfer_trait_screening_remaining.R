# module4 trait screening: transfer CITE RNA information for EOS/WBC
# -----------------------------------------------------------------------------
# Thin wrapper around:
#   03_pseudopair_pipeline/06_transfer_trait_screening.R
#
# Run in RStudio:
#   setwd("D:/RStudio/work/module4_traitaware_pseudopair")
#   source("03_pseudopair_pipeline/09_transfer_trait_screening_remaining.R")
# -----------------------------------------------------------------------------

ROOT <- "D:/RStudio/work/module4_traitaware_pseudopair"

TRANSFER_MAP_DIR <- file.path(ROOT, "02_trait_mapping/output/trait_screening_remaining")
TRANSFER_OUT_DIR <- file.path(ROOT, "03_pseudopair_pipeline/output/transfer_trait_screening_remaining")

CITE_MTX <- file.path(TRANSFER_MAP_DIR, "cite_extract/cite_trait_counts.mtx")
CITE_GENES <- file.path(TRANSFER_MAP_DIR, "cite_extract/cite_trait_counts_genes.txt")
CITE_CELLS <- file.path(TRANSFER_MAP_DIR, "cite_extract/cite_trait_counts_cells.txt")
CITE_SCORE_RDS <- file.path(
  ROOT,
  "01_data_prep_external_rna/output/trait_screening/cite_rna_scdrs_scores_trait_screening_remaining.rds"
)

TRAITS_ALL <- c("EOS", "WBC")
TRAITS_TO_RUN <- TRAITS_ALL

cat(">>> Run remaining trait-screening transfer via existing transfer script\n\n")
cat("Traits: ", paste(TRAITS_TO_RUN, collapse = ", "), "\n", sep = "")
cat("Mapping dir : ", TRANSFER_MAP_DIR, "\n", sep = "")
cat("Transfer out: ", TRANSFER_OUT_DIR, "\n\n", sep = "")

source(file.path(ROOT, "03_pseudopair_pipeline/06_transfer_trait_screening.R"))
