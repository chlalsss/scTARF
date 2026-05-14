# module4 trait screening: ATAC candidate diagnostics for EOS/WBC
# -----------------------------------------------------------------------------
# Purpose:
#   Reuse 02_trait_mapping/04_screen_trait_mapping_atac_candidates.R on the two
#   remaining mapping outputs.
#
# Important:
#   The leukocyte target below is only a post hoc proxy for checking whether
#   mapped RNA neighborhood signal lands on plausible ATAC lineages. It is not
#   used to select candidates, mapping neighbors, or seeds.
#
# Run in RStudio:
#   setwd("D:/RStudio/work/module4_traitaware_pseudopair")
#   source("02_trait_mapping/06_screen_trait_mapping_remaining_atac_candidates.R")
# -----------------------------------------------------------------------------

ROOT <- "D:/RStudio/work/module4_traitaware_pseudopair"

MAP_DIR <- file.path(ROOT, "02_trait_mapping/output/trait_screening_remaining")
OUT_DIR <- file.path(ROOT, "02_trait_mapping/output/trait_screening_remaining/diagnostics")
TRAITS_TO_SCREEN <- c("EOS", "WBC")

GROUP_TEST_TRAIT_MAP <- c(
  EOS = "EOSINOPHIL",
  WBC = "WHITE"
)

POSTHOC_TARGET_NAME <- "leukocyte_non_erythroid_proxy"
POSTHOC_TARGETS <- c(
  "04_Early.Baso",
  "06_CLP.1", "07_GMP", "08_GMP.Neut", "09_pDC", "10_cDC",
  "11_CD14.Mono.1", "12_CD14.Mono.2",
  "15_CLP.2", "16_Pre.B", "17_B", "18_Plasma",
  "19_CD8.N", "20_CD4.N1", "21_CD4.N2", "22_CD4.M",
  "23_CD8.EM", "24_CD8.CM", "25_NK"
)

GROUP_FDR_THRESHOLD <- 0.10
TOP_N_ATAC <- 350L

cat(">>> Screen remaining trait mappings into ATAC candidates\n\n")
cat("Traits: ", paste(TRAITS_TO_SCREEN, collapse = ", "), "\n", sep = "")
cat("Mapping dir: ", MAP_DIR, "\n", sep = "")
cat("Output dir : ", OUT_DIR, "\n\n", sep = "")

source(file.path(ROOT, "02_trait_mapping/04_screen_trait_mapping_atac_candidates.R"))
