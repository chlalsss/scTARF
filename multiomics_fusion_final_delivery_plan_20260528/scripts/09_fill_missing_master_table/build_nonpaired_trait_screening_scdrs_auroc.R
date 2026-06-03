# Compute non-paired RNA-reference scDRS-only AUROC for missing traits.
#
# This script uses existing NeurIPS CITE RNA scDRS score RDS files. It does not
# run the non-paired mapping/KNN/RWR pipeline.

WORK_ROOT <- "D:/RStudio/work"
DELIVERY_ROOT <- file.path(WORK_ROOT, "multiomics_fusion_final_delivery_plan_20260528")
MODULE4_ROOT <- file.path(WORK_ROOT, "module4_traitaware_pseudopair")
OUT_DIR <- file.path(DELIVERY_ROOT, "final_outputs/09_missing_master_table")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

if (!exists("TRAITS_TO_RUN")) TRAITS_TO_RUN <- c("EOS", "MCH", "RBC", "RDW", "WBC")
if (!exists("OUT_CSV")) {
  OUT_CSV <- file.path(OUT_DIR, "nonpaired_trait_screening_scdrs_auroc.csv")
}

BASE_SCORE_RDS <- file.path(
  MODULE4_ROOT,
  "01_data_prep_external_rna/output/cite_rna_scdrs_scores.rds"
)
SCREEN_SCORE_RDS <- file.path(
  MODULE4_ROOT,
  "01_data_prep_external_rna/output/trait_screening/cite_rna_scdrs_scores_trait_screening.rds"
)
REMAINING_SCORE_RDS <- file.path(
  MODULE4_ROOT,
  "01_data_prep_external_rna/output/trait_screening/cite_rna_scdrs_scores_trait_screening_remaining.rds"
)

TRAIT_IDS <- c(
  EOS = "UKB_460K.blood_EOSINOPHIL_COUNT",
  MCH = "UKB_460K.blood_MEAN_CORPUSCULAR_HEMOGLOBIN",
  RBC = "UKB_460K.blood_RED_COUNT",
  RDW = "UKB_460K.blood_RBC_DISTRIB_WIDTH",
  WBC = "UKB_460K.blood_WHITE_COUNT"
)

erythroid_mk <- c("MK/E prog", "Proerythroblast", "Erythroblast",
                  "Normoblast", "Reticulocyte")
mono <- c("CD14+ Mono", "CD16+ Mono")
myeloid_other <- c("G/M prog", "cDC1", "cDC2", "pDC")
lymph <- c(
  "CD4+ T activated", "CD4+ T naive", "CD4+ T activated integrinB7+",
  "CD4+ T CD314+ CD45RA+", "CD8+ T naive",
  "CD8+ T naive CD127+ CD26- CD101-", "CD8+ T CD49f+",
  "CD8+ T CD57+ CD45RA+", "CD8+ T CD57+ CD45RO+",
  "CD8+ T CD69+ CD45RA+", "CD8+ T CD69+ CD45RO+",
  "CD8+ T TIGIT+ CD45RA+", "CD8+ T TIGIT+ CD45RO+",
  "T reg", "T prog cycling", "MAIT", "gdT CD158b+", "gdT TCRVD2+",
  "dnT", "NK", "NK CD158e1+", "ILC", "ILC1", "Lymph prog",
  "Naive CD20+ B IGKC+", "Naive CD20+ B IGKC-",
  "B1 B IGKC+", "B1 B IGKC-", "Transitional B",
  "Plasma cell IGKC+", "Plasma cell IGKC-",
  "Plasmablast IGKC+", "Plasmablast IGKC-"
)

TARGETS <- list(
  MCH = erythroid_mk,
  RBC = erythroid_mk,
  RDW = erythroid_mk,
  EOS = c(mono, myeloid_other, lymph),
  WBC = c(mono, myeloid_other, lymph)
)
TARGET_NAMES <- c(
  MCH = "erythroid_mk",
  RBC = "erythroid_mk",
  RDW = "erythroid_mk",
  EOS = "leukocyte_non_erythroid_proxy",
  WBC = "leukocyte_non_erythroid_proxy"
)

stop_missing <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0L) {
    stop("Missing required file(s):\n  ", paste(missing, collapse = "\n  "),
         call. = FALSE)
  }
}

merge_scores <- function(paths) {
  stop_missing(paths[1])
  out <- as.data.frame(readRDS(paths[1]), stringsAsFactors = FALSE)
  for (path in paths[-1]) {
    if (!file.exists(path)) next
    x <- as.data.frame(readRDS(path), stringsAsFactors = FALSE)
    score_cols <- grep("_norm_score$", names(x), value = TRUE)
    score_cols <- setdiff(score_cols, names(out))
    if (length(score_cols) == 0L) next
    out <- merge(out, x[, c("cell", score_cols), drop = FALSE],
                 by = "cell", all.x = TRUE, sort = FALSE)
  }
  out
}

compute_auc <- function(score, label) {
  keep <- is.finite(score) & !is.na(label)
  score <- score[keep]
  label <- as.integer(label[keep])
  if (length(score) < 20L || sum(label == 1L) < 10L || sum(label == 0L) < 10L) {
    return(NA_real_)
  }
  if (requireNamespace("pROC", quietly = TRUE)) {
    return(as.numeric(pROC::auc(pROC::roc(label, score, quiet = TRUE,
                                           direction = "<",
                                           levels = c(0L, 1L)))))
  }
  r <- rank(score, ties.method = "average")
  n_pos <- sum(label == 1L)
  n_neg <- sum(label == 0L)
  (sum(r[label == 1L]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

top_cell_types <- function(cell_type, score, n = 5L) {
  ord <- order(score, decreasing = TRUE)
  top <- cell_type[head(ord, max(1L, ceiling(0.01 * length(ord))))]
  tab <- sort(table(top), decreasing = TRUE)
  paste(names(head(tab, n)), collapse = ";")
}

cat(">>> Non-paired trait-screening scDRS-only AUROC\n\n")
scores <- merge_scores(c(BASE_SCORE_RDS, SCREEN_SCORE_RDS, REMAINING_SCORE_RDS))
if (!all(c("cell", "cell_type") %in% names(scores))) {
  stop("Merged CITE score table must contain cell and cell_type.", call. = FALSE)
}

rows <- list()
for (trait in TRAITS_TO_RUN) {
  score_name <- paste0(trait, "_norm_score")
  if (!score_name %in% names(scores)) {
    stop("Missing score column for ", trait, ": ", score_name, call. = FALSE)
  }
  score <- as.numeric(scores[[score_name]])
  label <- as.integer(scores$cell_type %in% TARGETS[[trait]])
  auc <- compute_auc(score, label)

  cat("[trait] ", trait, " target=", TARGET_NAMES[[trait]],
      " AUROC=", signif(auc, 6), " pos=", sum(label), "/", length(label), "\n",
      sep = "")

  rows[[trait]] <- data.frame(
    dataset = "RNA reference scDRS",
    pairing_setting = "rna_reference_only",
    trait = trait,
    target_definition = TARGET_NAMES[[trait]],
    method_group = "scdrs",
    method_version = "scdrs_only",
    method = "scDRS-only",
    rna_reference = "NeurIPS CITE RNA trait-screening",
    top_k = NA_real_,
    feature_mode = "rna_score_only",
    atac_w = NA_real_,
    metric = "AUROC",
    value = auc,
    ci_low = NA_real_,
    ci_high = NA_real_,
    n_positive = sum(label == 1L),
    n_total = length(label),
    source_file = paste(c(SCREEN_SCORE_RDS, REMAINING_SCORE_RDS), collapse = ";"),
    note = "Computed from CITE trait-screening scDRS scores; EOS/WBC use leukocyte_non_erythroid proxy because CITE has no eosinophil cell type",
    trait_id = TRAIT_IDS[[trait]],
    top5_cell_types = top_cell_types(scores$cell_type, score, 5L),
    stringsAsFactors = FALSE
  )
}

out <- do.call(rbind, rows)
utils::write.csv(out, OUT_CSV, row.names = FALSE, quote = TRUE)
cat("\nWrote: ", OUT_CSV, "\n", sep = "")
