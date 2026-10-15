# module4 / external RNA: summarize NeurIPS 2021 CITE RNA scDRS outputs
#
# Produces cell-level scores and cell-type enrichment summaries. This script
# assumes scDRS has already been run by 04_run_cite_scdrs.R or run_scdrs_cite.py.

suppressPackageStartupMessages({
  library(data.table)
  library(pROC)
})

source(file.path("scripts", "_paths.R"))
ROOT <- sctarf_nonpaired_root()
SCDRS_DIR <- file.path(ROOT, "01_data_prep_external_rna/03_cite_scdrs_python")
SCORE_DIR <- file.path(SCDRS_DIR, "output")
META_FILE <- file.path(SCDRS_DIR, "cite_scdrs_cell_meta.csv.gz")
SUMMARY_OUT <- file.path(ROOT, "01_data_prep_external_rna/output")
dir.create(SUMMARY_OUT, showWarnings = FALSE, recursive = TRUE)

TRAITS <- c(
  MONO = "UKB_460K.blood_MONOCYTE_COUNT",
  PLT  = "UKB_460K.blood_PLATELET_COUNT",
  LYM  = "UKB_460K.blood_LYMPHOCYTE_COUNT"
)

TARGETS <- list(
  MONO = c("CD14+ Mono", "CD16+ Mono"),
  LYM = c(
    "CD4+ T activated", "CD4+ T naive", "CD4+ T activated integrinB7+",
    "CD4+ T CD314+ CD45RA+", "CD8+ T naive", "CD8+ T naive CD127+ CD26- CD101-",
    "CD8+ T CD49f+", "CD8+ T CD57+ CD45RA+", "CD8+ T CD57+ CD45RO+",
    "CD8+ T CD69+ CD45RA+", "CD8+ T CD69+ CD45RO+",
    "CD8+ T TIGIT+ CD45RA+", "CD8+ T TIGIT+ CD45RO+",
    "T reg", "T prog cycling", "MAIT", "gdT CD158b+", "gdT TCRVD2+",
    "dnT", "NK", "NK CD158e1+", "ILC", "ILC1", "Lymph prog",
    "Naive CD20+ B IGKC+", "Naive CD20+ B IGKC-", "B1 B IGKC+",
    "B1 B IGKC-", "Transitional B", "Plasma cell IGKC+",
    "Plasma cell IGKC-", "Plasmablast IGKC+", "Plasmablast IGKC-"
  ),
  PLT = c("MK/E prog", "Proerythroblast", "Erythroblast", "Normoblast", "Reticulocyte")
)

PLT_DIAGNOSTIC_TARGETS <- list(
  PLT_MKE_strict = c("MK/E prog"),
  PLT_ery_mke_broad = c("MK/E prog", "Proerythroblast", "Erythroblast", "Normoblast", "Reticulocyte"),
  PLT_ery_mke_no_retic = c("MK/E prog", "Proerythroblast", "Erythroblast", "Normoblast"),
  PLT_myeloid_distractor = c("CD14+ Mono", "CD16+ Mono", "G/M prog", "cDC1", "cDC2", "pDC")
)

score_col <- function(dt) {
  candidates <- c("norm_score", "zscore", "score")
  hit <- intersect(candidates, colnames(dt))[1]
  if (is.na(hit)) stop("No score column found in: ", paste(colnames(dt), collapse = ", "))
  hit
}

aux_cols <- function(dt) {
  intersect(c("pval", "pval_ctrl", "pval_mc", "fdr"), colnames(dt))
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

summarize_target <- function(score, cell_type, target_types, label_name) {
  label <- as.integer(cell_type %in% target_types)
  data.table(
    target_def = label_name,
    n_cells = length(score),
    n_positive_target = sum(label, na.rm = TRUE),
    target_mean_norm_score = mean(score[label == 1L], na.rm = TRUE),
    non_target_mean_norm_score = mean(score[label == 0L], na.rm = TRUE),
    delta_target_minus_rest = mean(score[label == 1L], na.rm = TRUE) -
      mean(score[label == 0L], na.rm = TRUE),
    target_auroc = compute_auc(score, label)
  )
}

cat(">>> Summarize CITE RNA scDRS\n\n")
if (!file.exists(META_FILE)) {
  stop("Missing CITE scDRS metadata: ", META_FILE,
       "\nRun source(\"01_data_prep_external_rna/03_prepare_cite_scdrs_input.R\") first.",
       call. = FALSE)
}
meta <- fread(META_FILE)
if (!("cell" %in% colnames(meta)) || !("cell_type" %in% colnames(meta))) {
  stop("Metadata must contain cell and cell_type columns: ", META_FILE, call. = FALSE)
}

cell_scores <- copy(meta)
summary_rows <- list()
celltype_rows <- list()

for (short in names(TRAITS)) {
  trait_id <- TRAITS[[short]]
  f <- file.path(SCORE_DIR, paste0(trait_id, ".score.gz"))
  if (!file.exists(f)) stop("Missing scDRS score file: ", f)
  cat("[trait] ", short, " <- ", basename(f), "\n", sep = "")
  dt <- fread(f)

  bc_col <- intersect(c("cell", "cell_id", "barcode", "V1"), colnames(dt))[1]
  if (is.na(bc_col)) bc_col <- colnames(dt)[1]
  sc_col <- score_col(dt)
  score_dt <- dt[, .(cell = as.character(get(bc_col)), score = as.numeric(get(sc_col)))]
  setnames(score_dt, "score", paste0(short, "_norm_score"))
  cell_scores <- merge(cell_scores, score_dt, by = "cell", all.x = TRUE, sort = FALSE)

  for (ac in aux_cols(dt)) {
    aux <- dt[, .(cell = as.character(get(bc_col)), value = get(ac))]
    setnames(aux, "value", paste0(short, "_", ac))
    cell_scores <- merge(cell_scores, aux, by = "cell", all.x = TRUE, sort = FALSE)
  }

  tmp <- merge(meta[, .(cell, cell_type)], score_dt, by = "cell", all.x = FALSE, sort = FALSE)
  sc_name <- paste0(short, "_norm_score")
  by_ct <- tmp[, .(
    n_cells = .N,
    mean_norm_score = mean(get(sc_name), na.rm = TRUE),
    median_norm_score = median(get(sc_name), na.rm = TRUE),
    q90_norm_score = as.numeric(quantile(get(sc_name), 0.90, na.rm = TRUE))
  ), by = cell_type][order(-mean_norm_score)]
  by_ct[, trait := short]
  setcolorder(by_ct, c("trait", "cell_type", "n_cells",
                       "mean_norm_score", "median_norm_score", "q90_norm_score"))
  celltype_rows[[short]] <- by_ct

  s <- summarize_target(tmp[[sc_name]], tmp$cell_type, TARGETS[[short]], short)
  s[, `:=`(
    trait = short,
    trait_id = trait_id,
    top5_cell_types = paste(head(by_ct$cell_type, 5), collapse = ";")
  )]
  setcolorder(s, c("trait", "trait_id", "target_def"))
  summary_rows[[short]] <- s
}

summary_dt <- rbindlist(summary_rows, use.names = TRUE, fill = TRUE)
celltype_dt <- rbindlist(celltype_rows, use.names = TRUE, fill = TRUE)

plt_score_col <- "PLT_norm_score"
plt_dt <- cell_scores[is.finite(get(plt_score_col)), .(
  cell_type,
  score = get(plt_score_col)
)]
plt_diag <- rbindlist(lapply(names(PLT_DIAGNOSTIC_TARGETS), function(name) {
  summarize_target(plt_dt$score, plt_dt$cell_type, PLT_DIAGNOSTIC_TARGETS[[name]], name)
}))
plt_diag[, trait := "PLT"]
setcolorder(plt_diag, c("trait", "target_def"))

fwrite(cell_scores, file.path(SUMMARY_OUT, "cite_rna_scdrs_scores.csv.gz"))
fwrite(summary_dt, file.path(SUMMARY_OUT, "cite_scdrs_target_enrichment.csv"))
fwrite(celltype_dt, file.path(SUMMARY_OUT, "cite_scdrs_celltype_enrichment.csv"))
fwrite(plt_diag, file.path(SUMMARY_OUT, "cite_scdrs_plt_diagnostic_targets.csv"))
saveRDS(cell_scores, file.path(SUMMARY_OUT, "cite_rna_scdrs_scores.rds"))

cat("\nTarget enrichment summary:\n")
print(summary_dt)
cat("\nPLT diagnostic targets:\n")
print(plt_diag)

cat("\nWrote:\n")
cat("  ", file.path(SUMMARY_OUT, "cite_rna_scdrs_scores.csv.gz"), "\n", sep = "")
cat("  ", file.path(SUMMARY_OUT, "cite_scdrs_target_enrichment.csv"), "\n", sep = "")
cat("  ", file.path(SUMMARY_OUT, "cite_scdrs_celltype_enrichment.csv"), "\n", sep = "")
cat("  ", file.path(SUMMARY_OUT, "cite_scdrs_plt_diagnostic_targets.csv"), "\n", sep = "")
cat("  ", file.path(SUMMARY_OUT, "cite_rna_scdrs_scores.rds"), "\n", sep = "")

cat("\nCITE scDRS summary OK.\n")
