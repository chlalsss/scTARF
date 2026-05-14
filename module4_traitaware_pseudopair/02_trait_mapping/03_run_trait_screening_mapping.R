# module4 trait screening: thin wrapper for existing trait-aware mapping
# -----------------------------------------------------------------------------
# Purpose:
#   Reuse 02_trait_mapping/01_align_granja_cite_by_trait.R for the selected
#   screening traits. No mapping algorithm is reimplemented here.
#
# Current selected traits:
#   MCH, RDW, RBC
#
# Run in RStudio:
#   setwd("D:/RStudio/work/module4_traitaware_pseudopair")
#   source("02_trait_mapping/03_run_trait_screening_mapping.R")
# -----------------------------------------------------------------------------

ROOT <- "D:/RStudio/work/module4_traitaware_pseudopair"

GS_FILE <- file.path(
  ROOT,
  "01_data_prep_external_rna/03_cite_scdrs_python/cite_trait_screening.gs"
)
ALIGN_OUT_DIR <- file.path(ROOT, "02_trait_mapping/output/trait_screening")
GA_PATH <- "D:/SCAVENGE-reproducibility-main/blood_data/multiomics_output_v2/gene_activity_matrix.rds"
ATAC_META_PATH <- file.path(ROOT, "04_evaluation/inputs/granja_cell_meta.rds")
CITE_H5AD <- file.path(
  ROOT,
  "01_data_prep_external_rna/03_cite_scdrs_python/input/cite_rna_raw_counts_gex_only.h5ad"
)
CITE_SCORE_RDS <- file.path(ROOT, "01_data_prep_external_rna/output/cite_rna_scdrs_scores.rds")
PY_HELPER <- file.path(ROOT, "02_trait_mapping/extract_cite_trait_counts.py")
PYTHON_EXE <- "D:/python/python.exe"

TRAITS_ALL <- c(
  MCH = "UKB_460K.blood_MEAN_CORPUSCULAR_HEMOGLOBIN",
  RDW = "UKB_460K.blood_RBC_DISTRIB_WIDTH",
  RBC = "UKB_460K.blood_RED_COUNT"
)
TRAITS_TO_RUN <- c("MCH", "RDW", "RBC")

K_NEIGHBOR <- 5L
BATCH_SIZE <- 512L
USE_GWAS_WEIGHT <- TRUE
FORCE_REEXTRACT_CITE <- FALSE

cat(">>> Run trait-screening mapping via existing trait-aware mapper\n\n")
cat("Traits: ", paste(TRAITS_TO_RUN, collapse = ", "), "\n", sep = "")
cat("Gene set: ", GS_FILE, "\n", sep = "")
cat("Output : ", ALIGN_OUT_DIR, "\n\n", sep = "")

source(file.path(ROOT, "02_trait_mapping/01_align_granja_cite_by_trait.R"))
