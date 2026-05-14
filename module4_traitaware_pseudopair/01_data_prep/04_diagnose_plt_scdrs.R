# module4 / Phase 1.2d: diagnose PLT RNA scDRS behavior
# ------------------------------------------------------------------
# PLT scDRS is unexpectedly high in myeloid cells. This script checks:
#   1. Which candidate cell-type target best matches PLT scDRS scores.
#   2. Whether platelet/megakaryocyte marker genes exist and where they express.
#   3. Whether the PLT GWAS gene set itself has a myeloid expression signature.

suppressPackageStartupMessages({
  library(Matrix)
  library(SummarizedExperiment)
  library(data.table)
  library(pROC)
})

ROOT <- "D:/RStudio/work/module4_traitaware_pseudopair"
OUT_DIR <- file.path(ROOT, "01_data_prep/output")
SCORE_FILE <- file.path(OUT_DIR, "granja_rna_scdrs_scores.csv.gz")
RNA_RSE <- "D:/SCAVENGE-reproducibility-main/blood_data/scRNA-Healthy-Hematopoiesis-191120.rds"
GS_FILE <- file.path(ROOT, "01_data_prep/03_rna_scdrs_python/granja_3traits.gs")

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

auc_one <- function(score, label) {
  keep <- is.finite(score) & !is.na(label)
  score <- score[keep]
  label <- label[keep]
  if (length(score) < 20 || sum(label == 1L) < 10 || sum(label == 0L) < 10) {
    return(NA_real_)
  }
  as.numeric(pROC::auc(pROC::roc(label, score, quiet = TRUE,
                                 direction = "<", levels = c(0L, 1L))))
}

program_score <- function(log_expr, genes) {
  g <- intersect(genes, rownames(log_expr))
  if (length(g) == 0) return(rep(NA_real_, ncol(log_expr)))
  Matrix::colMeans(log_expr[g, , drop = FALSE])
}

cat(">>> Diagnose PLT RNA scDRS\n\n")
stopifnot(file.exists(SCORE_FILE), file.exists(RNA_RSE), file.exists(GS_FILE))

scores <- fread(SCORE_FILE)
if (!"PLT_norm_score" %in% colnames(scores)) {
  stop("Missing PLT_norm_score in ", SCORE_FILE)
}

cat("[1/5] Candidate target definitions vs PLT score ...\n")
groups <- list(
  current_eryth_baso = c("02_Early.Eryth", "03_Late.Eryth", "04_Early.Baso"),
  eryth_only = c("02_Early.Eryth", "03_Late.Eryth"),
  baso_only = c("04_Early.Baso"),
  mono_cdc = c("10_cDC", "11_CD14.Mono.1", "12_CD14.Mono.2", "13_CD16.Mono"),
  mono_only = c("11_CD14.Mono.1", "12_CD14.Mono.2", "13_CD16.Mono"),
  myeloid_broad = c("07_GMP", "08_GMP.Neut", "09_pDC", "10_cDC",
                    "11_CD14.Mono.1", "12_CD14.Mono.2", "13_CD16.Mono", "14_Unk"),
  hspc_prog = c("01_HSC", "05_CMP.LMPP", "07_GMP"),
  lymphoid = c("16_Pre.B", "17_B", "18_Plasma", "19_CD8.N", "20_CD4.N1",
               "21_CD4.N2", "22_CD4.M", "23_CD8.EM", "24_CD8.CM", "25_NK")
)

target_diag <- rbindlist(lapply(names(groups), function(nm) {
  label <- as.integer(scores$cell_type %in% groups[[nm]])
  data.table(
    target_def = nm,
    n_positive = sum(label),
    mean_pos = mean(scores$PLT_norm_score[label == 1L], na.rm = TRUE),
    mean_rest = mean(scores$PLT_norm_score[label == 0L], na.rm = TRUE),
    delta_pos_minus_rest = mean(scores$PLT_norm_score[label == 1L], na.rm = TRUE) -
      mean(scores$PLT_norm_score[label == 0L], na.rm = TRUE),
    auroc = auc_one(scores$PLT_norm_score, label)
  )
}))[order(-auroc)]
print(target_diag)

cat("\n[2/5] Loading RNA counts for marker / gene-set expression diagnostics ...\n")
rna <- readRDS(RNA_RSE)
counts <- as(SummarizedExperiment::assay(rna, "counts"), "CsparseMatrix")
stopifnot(identical(colnames(counts), scores$barcode))
lib <- Matrix::colSums(counts)
lib[lib < 1] <- 1

marker_sets <- list(
  megakaryocyte_platelet = c("ITGA2B", "ITGB3", "GP1BA", "GP1BB", "GP5", "GP9",
                            "PF4", "PPBP", "MPL", "NFE2", "ZFPM1", "FLI1",
                            "TAL1", "RUNX1", "PEAR1", "ARHGEF3", "SLFN14"),
  erythroid = c("HBB", "HBA1", "HBA2", "ALAS2", "KLF1", "GYPA", "GYPB",
                "SLC4A1", "TFRC", "AHSP"),
  myeloid = c("LST1", "S100A8", "S100A9", "LYZ", "FCGR3A", "MS4A7",
              "CTSS", "ITGAX", "FCER1G", "SPI1", "CD33", "LGALS3",
              "QPCTL", "DNASE1L3")
)
marker_genes <- unique(unlist(marker_sets))
marker_genes_present <- intersect(marker_genes, rownames(counts))
marker_counts <- counts[marker_genes_present, , drop = FALSE]
marker_log <- log1p(Matrix::t(Matrix::t(marker_counts) / lib) * 1e4)

cat("[3/5] Marker program means by cell type ...\n")
program_dt <- data.table(barcode = scores$barcode, cell_type = scores$cell_type)
for (nm in names(marker_sets)) {
  program_dt[[paste0(nm, "_score")]] <- program_score(marker_log, marker_sets[[nm]])
}
program_by_ct <- program_dt[, lapply(.SD, mean, na.rm = TRUE),
                            by = cell_type,
                            .SDcols = patterns("_score$")]
for (nm in names(marker_sets)) {
  score_name <- paste0(nm, "_score")
  data.table::setorderv(program_by_ct, score_name, order = -1)
  cat("\nTop cell types for marker program: ", nm, "\n", sep = "")
  print(head(program_by_ct[, .(cell_type, score = get(score_name))], 8))
}

cat("\n[4/5] PLT GWAS gene-set expression signature ...\n")
gs <- parse_gs(GS_FILE)
plt_id <- "UKB_460K.blood_PLATELET_COUNT"
plt_gs <- gs[[plt_id]]
plt_gs <- plt_gs[gene %in% rownames(counts)]
setorder(plt_gs, -weight)
plt_counts <- counts[plt_gs$gene, , drop = FALSE]
plt_log <- log1p(Matrix::t(Matrix::t(plt_counts) / lib) * 1e4)
weighted_expr <- as.numeric(Matrix::crossprod(plt_gs$weight, plt_log)) / sum(plt_gs$weight)
plain_expr <- Matrix::colMeans(plt_log)
expr_dt <- data.table(
  barcode = scores$barcode,
  cell_type = scores$cell_type,
  PLT_norm_score = scores$PLT_norm_score,
  plt_gs_weighted_expr = weighted_expr,
  plt_gs_plain_expr = as.numeric(plain_expr)
)
expr_by_ct <- expr_dt[, .(
  n_cells = .N,
  mean_PLT_norm_score = mean(PLT_norm_score, na.rm = TRUE),
  mean_plt_gs_weighted_expr = mean(plt_gs_weighted_expr, na.rm = TRUE),
  mean_plt_gs_plain_expr = mean(plt_gs_plain_expr, na.rm = TRUE)
), by = cell_type][order(-mean_plt_gs_weighted_expr)]
print(head(expr_by_ct, 12))

cor_summary <- data.table(
  comparison = c("scDRS_norm_score_vs_weighted_expr", "scDRS_norm_score_vs_plain_expr"),
  spearman = c(
    suppressWarnings(cor(expr_dt$PLT_norm_score, expr_dt$plt_gs_weighted_expr,
                         method = "spearman", use = "complete.obs")),
    suppressWarnings(cor(expr_dt$PLT_norm_score, expr_dt$plt_gs_plain_expr,
                         method = "spearman", use = "complete.obs"))
  )
)
print(cor_summary)

cat("\n[5/5] PLT gene-set marker membership ...\n")
marker_membership <- rbindlist(lapply(names(marker_sets), function(nm) {
  g <- intersect(marker_sets[[nm]], plt_gs$gene)
  data.table(
    marker_program = nm,
    n_marker_genes = length(marker_sets[[nm]]),
    n_in_plt_gs = length(g),
    genes_in_plt_gs = paste(g, collapse = ";")
  )
}))
print(marker_membership)

top_plt_genes <- plt_gs[1:min(.N, 120)]

fwrite(target_diag, file.path(OUT_DIR, "plt_scdrs_candidate_targets.csv"))
fwrite(program_by_ct, file.path(OUT_DIR, "plt_marker_program_by_celltype.csv"))
fwrite(expr_by_ct, file.path(OUT_DIR, "plt_gwas_geneset_expr_by_celltype.csv"))
fwrite(cor_summary, file.path(OUT_DIR, "plt_scdrs_expr_correlation.csv"))
fwrite(marker_membership, file.path(OUT_DIR, "plt_gwas_marker_membership.csv"))
fwrite(top_plt_genes, file.path(OUT_DIR, "plt_gwas_top120_genes.csv"))

cat("\nWrote PLT diagnostics to ", OUT_DIR, "\n", sep = "")
cat("Phase 1.2d OK.\n")
