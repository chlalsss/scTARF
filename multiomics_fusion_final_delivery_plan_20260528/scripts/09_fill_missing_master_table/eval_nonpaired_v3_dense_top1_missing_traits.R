# Evaluate AUROC for the non-paired dense top1 missing-traits run.
#
# Run this after run_nonpaired_v3_dense_top1_missing_traits.R has produced
# TRS_v2_<TRAIT>_w*.tsv.gz files. It reuses module4's existing evaluator.

ROOT <- "D:/RStudio/work/module4_traitaware_pseudopair"

if (!exists("TOP1_TRAITS_TO_EVAL")) {
  TOP1_TRAITS_TO_EVAL <- c("EOS", "LYM", "MCH", "PLT", "WBC")
}
if (!exists("DENSE_TOP1_TAG")) DENSE_TOP1_TAG <- "v3_dense_top1_missing_traits"
if (!exists("METHOD_LABEL_NEW")) METHOD_LABEL_NEW <- DENSE_TOP1_TAG

red_proxy <- c("02_Early.Eryth", "03_Late.Eryth", "04_Early.Baso")
lymph_proxy <- c(
  "16_Pre.B", "17_B", "18_Plasma",
  "19_CD8.N", "20_CD4.N1", "21_CD4.N2", "22_CD4.M",
  "23_CD8.EM", "24_CD8.CM", "25_NK"
)
granulocyte_proxy <- c("04_Early.Baso", "07_GMP", "08_GMP.Neut")
leukocyte_proxy <- c(
  "04_Early.Baso",
  "06_CLP.1", "07_GMP", "08_GMP.Neut", "09_pDC", "10_cDC",
  "11_CD14.Mono.1", "12_CD14.Mono.2",
  "15_CLP.2", "16_Pre.B", "17_B", "18_Plasma",
  "19_CD8.N", "20_CD4.N1", "21_CD4.N2", "22_CD4.M",
  "23_CD8.EM", "24_CD8.CM", "25_NK"
)
TARGETS_BY_TRAIT <- list(
  PLT = red_proxy,
  LYM = lymph_proxy,
  MCH = red_proxy,
  EOS = granulocyte_proxy,
  WBC = leukocyte_proxy
)

V2_FILES_DIR_NEW <- file.path(
  ROOT,
  "03_pseudopair_pipeline/output",
  DENSE_TOP1_TAG,
  "module3_neighborhood_fdr_ge0p80"
)
OUT_DIR_NEW <- file.path(ROOT, "04_evaluation", DENSE_TOP1_TAG, "neighborhood_fdr_ge0p80")
OUT_CSV_NAME_NEW <- paste0(DENSE_TOP1_TAG, "_neighborhood_fdr_ge0p80_auroc.csv")
PATCHED_SCRIPT <- file.path(OUT_DIR_NEW, paste0("eval_auroc.", DENSE_TOP1_TAG, "_patched.R"))

r_literal_chr <- function(x) paste(sprintf('"%s"', x), collapse = ", ")

make_target_replacement <- function(traits) {
  lines <- character(length(traits) + 2L)
  lines[1] <- "TARGETS <- list("
  for (i in seq_along(traits)) {
    trait <- traits[i]
    suffix <- if (i < length(traits)) "," else ""
    lines[i + 1L] <- paste0("  ", trait, " = c(",
                            r_literal_chr(TARGETS_BY_TRAIT[[trait]]), ")",
                            suffix)
  }
  lines[length(lines)] <- ")"
  lines
}

make_pred_replacement <- function(traits) {
  pairs <- paste0(traits, ' = "', traits, '"')
  c(paste0("PRED_TRAIT <- list(", paste(pairs, collapse = ", "), ")"), "")
}

unknown <- setdiff(TOP1_TRAITS_TO_EVAL, names(TARGETS_BY_TRAIT))
if (length(unknown) > 0L) {
  stop("Unknown TOP1_TRAITS_TO_EVAL: ", paste(unknown, collapse = ", "),
       call. = FALSE)
}

cat(">>> evaluate dense top1 missing-traits AUROC\n\n")
cat("Traits : ", paste(TOP1_TRAITS_TO_EVAL, collapse = ", "), "\n", sep = "")
cat("TRS dir: ", V2_FILES_DIR_NEW, "\n", sep = "")
cat("Out dir: ", OUT_DIR_NEW, "\n\n", sep = "")

source(file.path(ROOT, "04_evaluation/run_trait_screening_auroc.R"), local = list2env(list(
  V2_FILES_DIR_NEW = V2_FILES_DIR_NEW,
  OUT_DIR_NEW = OUT_DIR_NEW,
  OUT_CSV_NAME_NEW = OUT_CSV_NAME_NEW,
  PATCHED_SCRIPT = PATCHED_SCRIPT,
  METHOD_LABEL_NEW = METHOD_LABEL_NEW,
  TARGET_REPLACEMENT = make_target_replacement(TOP1_TRAITS_TO_EVAL),
  PRED_REPLACEMENT = make_pred_replacement(TOP1_TRAITS_TO_EVAL)
), parent = globalenv()))

cat("\nResult CSV: ", file.path(OUT_DIR_NEW, OUT_CSV_NAME_NEW), "\n", sep = "")
