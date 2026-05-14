# module4 / Phase 1.2e: resolved PLT RNA seed score
# ------------------------------------------------------------------
# Diagnosis showed official PLT scDRS norm_score points to myeloid cells in
# Granja RNA. For downstream PLT seed transfer we keep the original scDRS score
# but add a PLT seed score based on megakaryocyte/platelet marker genes that are
# also present in the PLT GWAS gene set.

suppressPackageStartupMessages({
  library(Matrix)
  library(SummarizedExperiment)
  library(data.table)
  library(pROC)
})

ROOT <- "D:/RStudio/work/module4_traitaware_pseudopair"
OUT_DIR <- file.path(ROOT, "01_data_prep/output")
SCORE_FILE <- file.path(OUT_DIR, "granja_rna_scdrs_scores.rds")
RNA_RSE <- "D:/SCAVENGE-reproducibility-main/blood_data/scRNA-Healthy-Hematopoiesis-191120.rds"
GS_FILE <- file.path(ROOT, "01_data_prep/03_rna_scdrs_python/granja_3traits.gs")

PLT_TARGET <- c("02_Early.Eryth", "03_Late.Eryth", "04_Early.Baso")
MONO_TARGET <- c("11_CD14.Mono.1", "12_CD14.Mono.2", "13_CD16.Mono")
MYELOID_TARGET <- c("10_cDC", "11_CD14.Mono.1", "12_CD14.Mono.2", "13_CD16.Mono")

MK_MARKERS_CANDIDATE <- c(
  "ITGA2B", "ITGB3", "GP1BA", "GP1BB", "GP5", "GP9", "PF4", "PPBP",
  "MPL", "NFE2", "ZFPM1", "FLI1", "TAL1", "RUNX1", "PEAR1", "ARHGEF3",
  "SLFN14"
)

parse_gs <- function(path) {
  lines <- readLines(path, warn = FALSE)
  out <- list()
  for (line in lines[-1]) {
    parts <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(parts) < 2) next
    pairs <- strsplit(strsplit(parts[2], ",", fixed = TRUE)[[1]], ":", fixed = TRUE)
    out[[parts[1]]] <- data.table(
      gene = vapply(pairs, `[`, character(1), 1),
      weight = as.numeric(vapply(
        pairs,
        function(x) if (length(x) >= 2) x[2] else "1",
        character(1)
      ))
    )
  }
  out
}

zscore <- function(x) {
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s < 1e-8) return(rep(NA_real_, length(x)))
  as.numeric((x - mean(x, na.rm = TRUE)) / s)
}

auc_one <- function(score, target_cell_types, cell_type) {
  label <- as.integer(cell_type %in% target_cell_types)
  keep <- is.finite(score)
  if (sum(label[keep] == 1L) < 10 || sum(label[keep] == 0L) < 10) {
    return(NA_real_)
  }
  as.numeric(pROC::auc(pROC::roc(label[keep], score[keep], quiet = TRUE,
                                 direction = "<", levels = c(0L, 1L))))
}

summary_one <- function(name, score, cell_type, target) {
  label <- as.integer(cell_type %in% target)
  data.table(
    score_name = name,
    n_positive = sum(label),
    mean_target = mean(score[label == 1L], na.rm = TRUE),
    mean_rest = mean(score[label == 0L], na.rm = TRUE),
    delta_target_minus_rest = mean(score[label == 1L], na.rm = TRUE) -
      mean(score[label == 0L], na.rm = TRUE),
    auroc_current_plt_target = auc_one(score, PLT_TARGET, cell_type),
    auroc_mono_cdc_target = auc_one(score, MYELOID_TARGET, cell_type)
  )
}

cat(">>> Resolve PLT RNA seed score\n\n")
stopifnot(file.exists(SCORE_FILE), file.exists(RNA_RSE), file.exists(GS_FILE))

scores <- as.data.table(readRDS(SCORE_FILE))
rna <- readRDS(RNA_RSE)
counts <- as(SummarizedExperiment::assay(rna, "counts"), "CsparseMatrix")
stopifnot(identical(colnames(counts), scores$barcode))

gs <- parse_gs(GS_FILE)
plt_gs <- gs[["UKB_460K.blood_PLATELET_COUNT"]]
mk_markers <- intersect(MK_MARKERS_CANDIDATE, plt_gs$gene)
mk_markers <- intersect(mk_markers, rownames(counts))
if (length(mk_markers) < 5) {
  stop("Too few PLT GWAS-overlapping MK markers: ", paste(mk_markers, collapse = ", "))
}

cat("MK/platelet markers used for PLT rescue score:\n  ",
    paste(mk_markers, collapse = ", "), "\n\n", sep = "")

lib <- Matrix::colSums(counts)
lib[lib < 1] <- 1
mk_counts <- counts[mk_markers, , drop = FALSE]
mk_log <- log1p(Matrix::t(Matrix::t(mk_counts) / lib) * 1e4)
mk_raw <- Matrix::colMeans(mk_log)
mk_z <- zscore(as.numeric(mk_raw))

resolved <- copy(scores)
setnames(resolved, "PLT_norm_score", "PLT_scdrs_original_norm_score")
resolved[, PLT_mk_marker_raw_score := as.numeric(mk_raw)]
resolved[, PLT_mk_marker_zscore := mk_z]
resolved[, PLT_norm_score := PLT_mk_marker_zscore]
setcolorder(resolved, c(
  "barcode", "cell_type",
  "MONO_norm_score", "MONO_pval",
  "PLT_norm_score", "PLT_scdrs_original_norm_score", "PLT_mk_marker_raw_score",
  "PLT_mk_marker_zscore", "PLT_pval",
  "LYM_norm_score", "LYM_pval"
))

score_summary <- rbindlist(list(
  summary_one("PLT_scdrs_original_norm_score", scores$PLT_norm_score, scores$cell_type, PLT_TARGET),
  summary_one("PLT_mk_marker_zscore", mk_z, scores$cell_type, PLT_TARGET)
))

by_ct <- resolved[, .(
  n_cells = .N,
  mean_PLT_resolved = mean(PLT_norm_score, na.rm = TRUE),
  mean_PLT_original_scdrs = mean(PLT_scdrs_original_norm_score, na.rm = TRUE),
  mean_PLT_mk_marker_raw = mean(PLT_mk_marker_raw_score, na.rm = TRUE)
), by = cell_type][order(-mean_PLT_resolved)]

saveRDS(resolved, file.path(OUT_DIR, "granja_rna_trait_scores_resolved.rds"))
fwrite(resolved, file.path(OUT_DIR, "granja_rna_trait_scores_resolved.csv.gz"))
fwrite(score_summary, file.path(OUT_DIR, "plt_seed_score_resolution_summary.csv"))
fwrite(by_ct, file.path(OUT_DIR, "plt_seed_score_resolved_by_celltype.csv"))
writeLines(mk_markers, file.path(OUT_DIR, "plt_seed_score_marker_genes.txt"), useBytes = TRUE)

cat("Resolution summary:\n")
print(score_summary)
cat("\nTop cell types after PLT seed rescue:\n")
print(head(by_ct, 10))
cat("\nWrote:\n")
cat("  ", file.path(OUT_DIR, "granja_rna_trait_scores_resolved.rds"), "\n", sep = "")
cat("  ", file.path(OUT_DIR, "granja_rna_trait_scores_resolved.csv.gz"), "\n", sep = "")
cat("  ", file.path(OUT_DIR, "plt_seed_score_resolution_summary.csv"), "\n", sep = "")
cat("  ", file.path(OUT_DIR, "plt_seed_score_resolved_by_celltype.csv"), "\n", sep = "")
cat("  ", file.path(OUT_DIR, "plt_seed_score_marker_genes.txt"), "\n", sep = "")

cat("\nPhase 1.2e OK.\n")
