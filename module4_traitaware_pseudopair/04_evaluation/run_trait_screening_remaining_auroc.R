# module4 trait screening: AUROC for EOS/WBC
# -----------------------------------------------------------------------------
# Reuses:
#   04_evaluation/run_trait_screening_auroc.R
#
# Post hoc target proxies:
#   EOS is checked against a granulocyte-like ATAC proxy because this Granja ATAC
#   label set does not contain a clean eosinophil label.
#   WBC is checked against a broad non-erythroid leukocyte ATAC proxy.
#
# These celltype proxies are for evaluation only. They are not used in mapping,
# RNA neighborhood FDR, seed selection, KNN, RWR, or TRS.
#
# Run in RStudio after module3 finishes:
#   setwd("D:/RStudio/work/module4_traitaware_pseudopair")
#   source("04_evaluation/run_trait_screening_remaining_auroc.R")
# -----------------------------------------------------------------------------

ROOT <- "D:/RStudio/work/module4_traitaware_pseudopair"

V2_FILES_DIR_NEW <- file.path(
  ROOT,
  "03_pseudopair_pipeline/output/module3_v2_trait_screening_remaining_ge0p80"
)
OUT_DIR_NEW <- file.path(ROOT, "04_evaluation/trait_screening_remaining_ge0p80")
OUT_CSV_NAME_NEW <- "trait_screening_remaining_ge0p80_auroc.csv"
PATCHED_SCRIPT <- file.path(
  OUT_DIR_NEW,
  "eval_auroc.trait_screening_remaining_ge0p80_patched.R"
)

TARGET_REPLACEMENT <- c(
  'TARGETS <- list(',
  '  EOS = c("04_Early.Baso", "07_GMP", "08_GMP.Neut"),',
  '  WBC = c("04_Early.Baso",',
  '          "06_CLP.1", "07_GMP", "08_GMP.Neut", "09_pDC", "10_cDC",',
  '          "11_CD14.Mono.1", "12_CD14.Mono.2",',
  '          "15_CLP.2", "16_Pre.B", "17_B", "18_Plasma",',
  '          "19_CD8.N", "20_CD4.N1", "21_CD4.N2", "22_CD4.M",',
  '          "23_CD8.EM", "24_CD8.CM", "25_NK")',
  ')'
)
PRED_REPLACEMENT <- c(
  'PRED_TRAIT <- list(EOS = "EOS", WBC = "WBC")',
  ''
)

cat(">>> Run AUROC evaluation for remaining trait screening\n\n")
cat("TRS files dir: ", V2_FILES_DIR_NEW, "\n", sep = "")
cat("Output dir   : ", OUT_DIR_NEW, "\n\n", sep = "")

source(file.path(ROOT, "04_evaluation/run_trait_screening_auroc.R"))
