# module4 trait screening: build CITE score RDS for EOS/WBC
# -----------------------------------------------------------------------------
# Purpose:
#   Reuse 10_build_trait_screening_score_rds.R for the two remaining screening
#   traits. The original trait IDs are long, but the pipeline short names are:
#     EOS = UKB_460K.blood_EOSINOPHIL_COUNT
#     WBC = UKB_460K.blood_WHITE_COUNT
#
# Run in RStudio:
#   Run from the repository root with optional SCTARF_* path overrides.
#   source("01_data_prep_external_rna/11_build_trait_screening_remaining_score_rds.R")
# -----------------------------------------------------------------------------

if (!exists("ROOT")) ROOT <- Sys.getenv("SCTARF_NONPAIRED_ROOT", file.path("outputs", "nonpaired"))
LOCAL_SCRIPT_DIR <- if (exists("SCTARF_SCRIPT_DIR")) {
  SCTARF_SCRIPT_DIR
} else {
  dirname(tryCatch(normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = FALSE),
                   error = function(e) "scripts/03_nonpaired_analysis/31_build_trait_screening_remaining_score_rds.R"))
}

TRAIT_IDS <- c(
  EOS = "UKB_460K.blood_EOSINOPHIL_COUNT",
  WBC = "UKB_460K.blood_WHITE_COUNT"
)

OUT_RDS <- file.path(
  ROOT,
  "01_data_prep_external_rna/output/trait_screening/cite_rna_scdrs_scores_trait_screening_remaining.rds"
)
OUT_CSV <- file.path(
  ROOT,
  "01_data_prep_external_rna/output/trait_screening/cite_rna_scdrs_scores_trait_screening_remaining.csv.gz"
)

cat(">>> Build CITE score RDS for remaining screening traits\n\n")
cat("Traits: ", paste(names(TRAIT_IDS), collapse = ", "), "\n", sep = "")
cat("Output: ", OUT_RDS, "\n\n", sep = "")

source(file.path(LOCAL_SCRIPT_DIR, "30_build_trait_screening_score_rds.R"))
