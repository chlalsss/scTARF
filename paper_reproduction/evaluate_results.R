# Historical score reading and ROC conventions retained unchanged.
evaluate_paper <- function(data_root, replayed_root, output_file) {
bmmc <- file.path(data_root, "paired")
audit <- replayed_root
out <- output_file
suppressPackageStartupMessages({library(pROC); library(data.table)})
paired_targets <- list(MONO = c("CD14+ Mono", "CD16+ Mono"), PLT = "MK/E prog",
  LYM = c("B1 B", "Naive CD20+ B", "Transitional B", "Plasma cell"),
  LYM_TNK = c("NK", "ILC", "CD4+ T activated", "CD4+ T naive", "CD8+ T", "CD8+ T naive"),
  LYM_wide = c("NK", "ILC", "CD4+ T activated", "CD4+ T naive", "CD8+ T", "CD8+ T naive", "Lymph prog", "B1 B", "Naive CD20+ B", "Transitional B", "Plasma cell"))
red <- c("02_Early.Eryth", "03_Late.Eryth", "04_Early.Baso")
granja_targets <- list(MONO = c("11_CD14.Mono.1", "12_CD14.Mono.2"), RDW = red,
  LYM = c("16_Pre.B", "17_B", "18_Plasma", "19_CD8.N", "23_CD8.EM", "24_CD8.CM", "20_CD4.N1", "21_CD4.N2", "22_CD4.M", "25_NK"))
read_meta <- function(path) {
  x <- as.data.frame(readRDS(path))
  if (identical(rownames(x), as.character(seq_len(nrow(x)))) && "barcode" %in% names(x)) rownames(x) <- as.character(x$barcode)
  stopifnot(!anyDuplicated(rownames(x)), "cell_type" %in% names(x))
  x
}
records <- list()
add <- function(setting, target, method, weight, score, meta, positives, expected, direction) {
  stopifnot(!anyDuplicated(names(score)), nrow(expected) == 1L)
  cells <- intersect(names(score), rownames(meta))
  pred <- score[cells]
  labels <- as.integer(as.character(meta[cells, "cell_type"]) %in% positives)
  keep <- is.finite(pred)
  roc_args <- list(response = labels[keep], predictor = pred[keep], quiet = TRUE, levels = c(0L, 1L))
  if (direction != "auto") roc_args$direction <- direction
  roc <- suppressMessages(do.call(pROC::roc, roc_args))
  value <- as.numeric(pROC::auc(roc))
  records[[length(records) + 1L]] <<- data.frame(setting, target, method, weight,
    n_total = sum(keep), n_positive = sum(labels[keep]), auroc = value,
    roc_direction = roc$direction, reference_auroc = expected$auroc,
    abs_difference = abs(value - expected$auroc), matched = abs(value - expected$auroc) <= 1e-12)
}
read_score <- function(path, legacy_v1 = FALSE) {
  x <- fread(path)
  if (legacy_v1) {
    # Preserve the actual historical evaluator's column lookup and fallback.
    bc <- intersect(c("barcode", "V1", "rn"), names(x))[1]
    value <- intersect(c("TRS", "trs", "score"), names(x))[1]
    if (is.na(bc) || is.na(value)) { bc <- names(x)[1]; value <- names(x)[2] }
    return(setNames(as.numeric(x[[value]]), as.character(x[[bc]])))
  }
  bc <- intersect(c("barcode", "cell", "cell_id", "V1", "rn"), names(x))[1]
  if (is.na(bc)) bc <- names(x)[1]
  value <- intersect(c("TRS", "norm_score", "trs", "score"), names(x))[1]
  stopifnot(!is.na(value))
  setNames(as.numeric(x[[value]]), as.character(x[[bc]]))
}

meta <- read_meta(file.path(bmmc, "BMMC_cell_meta.rds"))
baseline <- as.data.frame(readRDS(file.path(bmmc, "TRS_matrix_baseline.rds")))
if ("cell" %in% names(baseline)) {
  rownames(baseline) <- as.character(baseline$cell)
  baseline$cell <- NULL
}
ref <- read.csv(file.path(bmmc, "eval_v2_auroc.csv"))
for (target in names(paired_targets)) {
  trait <- if (startsWith(target, "LYM")) "LYM" else target
  for (method in c("baseline_atac", "v1_paired", "v2")) {
    weights <- if (method == "v2") c(0.5, 0.6, 0.7, 0.8, 0.9, 1) else NA_real_
    for (w in weights) {
      score <- if (method == "baseline_atac") setNames(baseline[[trait]], rownames(baseline)) else
        if (method == "v1_paired") read_score(file.path(bmmc, paste0("TRS_dual_", trait, ".tsv.gz")), legacy_v1 = TRUE) else
        read_score(file.path(audit, "paired", sprintf("paired_%s_w%.1f.tsv", trait, w)))
      expected <- ref[ref$target == target & ref$method == method & (if (is.na(w)) is.na(ref$atac_w) else !is.na(ref$atac_w) & abs(ref$atac_w - w) < 1e-9), ]
      add("paired", target, method, w, score, meta, paired_targets[[target]], expected, "auto")
    }
  }
}
rna_ref <- read.csv(file.path(bmmc, "bmmc_paired_scdrs_only_auroc.csv"))
for (target in c("MONO", "PLT", "LYM_wide")) {
  trait_id <- c(MONO = "MONOCYTE_COUNT", PLT = "PLATELET_COUNT", LYM_wide = "LYMPHOCYTE_COUNT")[[target]]
  path <- file.path(bmmc, paste0("UKB_460K.blood_", trait_id, ".score.gz"))
  expected <- rna_ref[rna_ref$trait == target, ]
  expected$auroc <- expected$value
  add("paired_RNA_only", target, "scDRS-only", NA, read_score(path), meta, paired_targets[[target]], expected, "<")
}

m4 <- file.path(data_root, "nonpaired")
meta <- read_meta(file.path(m4, "04_evaluation/inputs/granja_cell_meta.rds"))
baseline <- readRDS(file.path(m4, "04_evaluation/inputs/granja_baseline_TRS_matrix.rds"))
ref_paths <- c(MONO = "v3_sparse_mapping_top5significant/mono_neighborhood_fdr_ge0p80/v3_sparse_mapping_top5significant_mono_neighborhood_fdr_ge0p80_auroc.csv",
  RDW = "v3_sparse_other_traits/neighborhood_fdr_ge0p80/v3_sparse_other_traits_neighborhood_fdr_ge0p80_auroc.csv",
  LYM = "v3_sparse_lym_granja_rna/neighborhood_fdr_ge0p80/v3_sparse_lym_granja_rna_neighborhood_fdr_ge0p80_auroc.csv")
for (trait in names(granja_targets)) {
  ref <- read.csv(file.path(m4, "04_evaluation", ref_paths[[trait]]))
  for (w in c(NA, 0.7, 1)) {
    score <- if (is.na(w)) setNames(baseline[[trait]], rownames(baseline)) else
      read_score(file.path(audit, "nonpaired", sprintf("nonpaired_%s_%s_w%.1f.tsv", trait, trait, w)))
    expected <- ref[ref$target == trait & (if (is.na(w)) is.na(ref$atac_w) else !is.na(ref$atac_w) & abs(ref$atac_w - w) < 1e-9), ]
    add("nonpaired", trait, if (is.na(w)) "baseline_atac" else "v2", w, score, meta, granja_targets[[trait]], expected, "<")
  }
}

reverse <- file.path(data_root, "reverse")
meta <- read_meta(file.path(m4, "01_data_prep/output/granja_rna_cell_meta.rds"))
baseline <- readRDS(file.path(reverse, "_inputs/granja_rna_scdrs_baseline_matrix.rds"))
ref <- read.csv(file.path(reverse, "auroc/reverse_granja_rna_main_plt_mch_rbc_auroc.csv"))
for (trait in c("PLT", "MCH", "RBC")) for (w in c(NA, 0.7, 1)) {
  score <- if (is.na(w)) setNames(baseline[[trait]], rownames(baseline)) else
    read_score(file.path(audit, "reverse", sprintf("reverse_%s_w%.1f.tsv", trait, w)))
  expected <- ref[ref$target == trait & (if (is.na(w)) is.na(ref$atac_w) else !is.na(ref$atac_w) & abs(ref$atac_w - w) < 1e-9), ]
  add("reverse", trait, if (is.na(w)) "baseline_rna_scdrs" else "reverse", w, score, meta, red, expected, "<")
}
result <- do.call(rbind, records)
write.csv(result, out, row.names = FALSE)
print(result[c("setting", "target", "method", "weight", "auroc", "abs_difference", "matched")], row.names = FALSE)
stopifnot(nrow(result) == 61L, all(result$matched))
cat("PASS: 61 AUROCs recomputed from real replayed scores and historical baseline scores\n")

}
