# module4 / Phase 1.2c: summarize Granja RNA scDRS outputs
# ------------------------------------------------------------------
# Produces a compact cell-level score table and cell-type enrichment summary.

suppressPackageStartupMessages({
  library(data.table)
  library(pROC)
})

source(file.path("scripts", "_paths.R"))
ROOT <- sctarf_nonpaired_root()
SCDRS_DIR <- file.path(ROOT, "01_data_prep/03_rna_scdrs_python")
OUT_DIR <- file.path(SCDRS_DIR, "output")
META_FILE <- file.path(ROOT, "01_data_prep/output/granja_rna_cell_meta.rds")
SUMMARY_OUT <- file.path(ROOT, "01_data_prep/output")

TRAITS <- c(
  MONO = "UKB_460K.blood_MONOCYTE_COUNT",
  PLT  = "UKB_460K.blood_PLATELET_COUNT",
  LYM  = "UKB_460K.blood_LYMPHOCYTE_COUNT"
)

TARGETS <- list(
  MONO = c("11_CD14.Mono.1", "12_CD14.Mono.2"),
  PLT  = c("02_Early.Eryth", "03_Late.Eryth", "04_Early.Baso"),
  LYM  = c("16_Pre.B", "17_B", "18_Plasma",
           "19_CD8.N", "23_CD8.EM", "24_CD8.CM",
           "20_CD4.N1", "21_CD4.N2", "22_CD4.M",
           "25_NK")
)

score_col <- function(dt) {
  candidates <- c("norm_score", "zscore", "score")
  hit <- intersect(candidates, colnames(dt))[1]
  if (is.na(hit)) stop("No score column found in: ", paste(colnames(dt), collapse = ", "))
  hit
}

fdr_col <- function(dt) {
  candidates <- c("pval", "pval_ctrl", "pval_mc", "pval_nlog10", "fdr")
  intersect(candidates, colnames(dt))[1]
}

compute_auc <- function(score, label) {
  keep <- is.finite(score) & !is.na(label)
  score <- score[keep]
  label <- label[keep]
  if (length(score) < 20 || sum(label == 1L) < 10 || sum(label == 0L) < 10) {
    return(NA_real_)
  }
  as.numeric(pROC::auc(pROC::roc(label, score, quiet = TRUE,
                                 direction = "<", levels = c(0L, 1L))))
}

cat(">>> Phase 1.2c summarize Granja RNA scDRS\n\n")
stopifnot(file.exists(META_FILE))
meta <- readRDS(META_FILE)
if (is.null(rownames(meta))) rownames(meta) <- meta$barcode

cell_scores <- data.table(barcode = rownames(meta), cell_type = as.character(meta$cell_type))
summary_rows <- list()
celltype_rows <- list()

for (short in names(TRAITS)) {
  trait_id <- TRAITS[[short]]
  f <- file.path(OUT_DIR, paste0(trait_id, ".score.gz"))
  if (!file.exists(f)) stop("Missing scDRS score file: ", f)
  cat("[trait] ", short, " <- ", basename(f), "\n", sep = "")
  dt <- fread(f)

  bc_col <- intersect(c("cell", "cell_id", "barcode", "V1"), colnames(dt))[1]
  if (is.na(bc_col)) bc_col <- colnames(dt)[1]
  sc_col <- score_col(dt)
  fc <- fdr_col(dt)

  score <- dt[[sc_col]]
  names(score) <- as.character(dt[[bc_col]])
  common <- intersect(cell_scores$barcode, names(score))
  cell_scores[barcode %in% common, paste0(short, "_norm_score") := score[barcode]]

  if (!is.na(fc)) {
    aux <- dt[[fc]]
    names(aux) <- as.character(dt[[bc_col]])
    cell_scores[barcode %in% common, paste0(short, "_", fc) := aux[barcode]]
  }

  tmp <- data.table(
    barcode = common,
    cell_type = as.character(meta[common, "cell_type"]),
    score = as.numeric(score[common])
  )
  by_ct <- tmp[, .(
    n_cells = .N,
    mean_norm_score = mean(score, na.rm = TRUE),
    median_norm_score = median(score, na.rm = TRUE),
    q90_norm_score = as.numeric(quantile(score, 0.90, na.rm = TRUE))
  ), by = cell_type][order(-mean_norm_score)]
  by_ct[, trait := short]
  setcolorder(by_ct, c("trait", "cell_type", "n_cells",
                       "mean_norm_score", "median_norm_score", "q90_norm_score"))
  celltype_rows[[short]] <- by_ct

  label <- as.integer(tmp$cell_type %in% TARGETS[[short]])
  summary_rows[[short]] <- data.table(
    trait = short,
    trait_id = trait_id,
    n_cells = nrow(tmp),
    n_positive_target = sum(label),
    target_mean_norm_score = mean(tmp$score[label == 1L], na.rm = TRUE),
    non_target_mean_norm_score = mean(tmp$score[label == 0L], na.rm = TRUE),
    delta_target_minus_rest = mean(tmp$score[label == 1L], na.rm = TRUE) -
      mean(tmp$score[label == 0L], na.rm = TRUE),
    target_auroc = compute_auc(tmp$score, label),
    top5_cell_types = paste(head(by_ct$cell_type, 5), collapse = ";")
  )
}

summary_dt <- rbindlist(summary_rows)
celltype_dt <- rbindlist(celltype_rows)

fwrite(cell_scores, file.path(SUMMARY_OUT, "granja_rna_scdrs_scores.csv.gz"))
fwrite(summary_dt, file.path(SUMMARY_OUT, "granja_scdrs_target_enrichment.csv"))
fwrite(celltype_dt, file.path(SUMMARY_OUT, "granja_scdrs_celltype_enrichment.csv"))
saveRDS(cell_scores, file.path(SUMMARY_OUT, "granja_rna_scdrs_scores.rds"))

cat("\nTarget enrichment summary:\n")
print(summary_dt)
cat("\nWrote:\n")
cat("  ", file.path(SUMMARY_OUT, "granja_rna_scdrs_scores.csv.gz"), "\n", sep = "")
cat("  ", file.path(SUMMARY_OUT, "granja_scdrs_target_enrichment.csv"), "\n", sep = "")
cat("  ", file.path(SUMMARY_OUT, "granja_scdrs_celltype_enrichment.csv"), "\n", sep = "")
cat("  ", file.path(SUMMARY_OUT, "granja_rna_scdrs_scores.rds"), "\n", sep = "")

cat("\nPhase 1.2c OK.\n")
