# module4 / 04_evaluation / 00: Granja 评估输入预处理
# ------------------------------------------------------------------
# 目的: 把 Granja 的纯 ATAC SCAVENGE 产物 (TRS matrix + cell_type meta)
#       标准化成跟 BMMC 同款格式, 供 eval_auroc.R 复用 (module3 接口).
#       不算 AUROC, 只做格式转换 + cell_type 列表打印.
#
# 输入 (来自 SCAVENGE-reproducibility-main):
#   blood_data/scavenge_out_small_all_traits/TRS_matrix_small_all_traits_fullcells.rds
#   blood_data/scavenge_out_small_all_traits/SCAVENGE_small_MONO.tsv.gz  (拿 cell_type)
#
# 输出 (04_evaluation/inputs/):
#   granja_baseline_TRS_matrix.rds   data.frame, rowname = cell barcode,
#                                    cols = MONO/PLT/LYM/... (15 traits)
#   granja_cell_meta.rds             data.frame, rowname = barcode,
#                                    cols 含 cell_type, 跟 BMMC_cell_meta.rds 同格式
#   granja_cell_types_summary.csv    每个 cell_type 的细胞数 (供后续定义 TARGETS)
# ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

BASE_DIR   <- "D:/SCAVENGE-reproducibility-main/blood_data/scavenge_out_small_all_traits"
OUT_DIR    <- "D:/RStudio/work/module4_traitaware_pseudopair/04_evaluation/inputs"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

REF_TRAIT  <- "MONO"   # 用哪个 SCAVENGE_small_*.tsv.gz 抽 cell_type meta
                       # (任意 trait 的 TSV 都含同一份 BioClassification)

cat(">>> Granja evaluation inputs preparation\n\n")

# ---- [1/4] load TRS matrix ----
cat("[1/4] Loading TRS_matrix_small_all_traits_fullcells.rds ...\n")
trs <- readRDS(file.path(BASE_DIR, "TRS_matrix_small_all_traits_fullcells.rds"))
cat("       class:", class(trs)[1], "  dim:",
    paste(dim(trs), collapse = " x "), "\n")
cat("       colnames (first 8):",
    paste(head(colnames(trs), 8), collapse = ", "), "\n")

# Normalize to data.frame with cell barcode as rownames (BMMC convention)
if (inherits(trs, "data.table")) {
  # has explicit cell_id-like column
  bc_col_cands <- c("cell", "cell_id", "barcode", "V1")
  bc_col       <- intersect(bc_col_cands, colnames(trs))[1]
  if (is.na(bc_col)) stop("data.table: no cell barcode column found")
  trs_df <- as.data.frame(trs)
  rownames(trs_df) <- as.character(trs_df[[bc_col]])
  trs_df[[bc_col]] <- NULL
  cat("       extracted barcode col -> rownames, dropped col '", bc_col,
      "'\n", sep = "")
} else if (inherits(trs, "matrix")) {
  if (is.null(rownames(trs))) stop("matrix has no rownames (cell barcodes)")
  trs_df <- as.data.frame(trs)
} else if (inherits(trs, "data.frame")) {
  bc_col_cands <- c("cell", "cell_id", "barcode")
  bc_col       <- intersect(bc_col_cands, colnames(trs))
  if (length(bc_col) > 0) {
    rownames(trs) <- as.character(trs[[bc_col[1]]])
    trs[[bc_col[1]]] <- NULL
    cat("       extracted barcode col '", bc_col[1], "' -> rownames\n",
        sep = "")
  } else if (is.null(rownames(trs)) ||
             identical(rownames(trs), as.character(seq_len(nrow(trs))))) {
    stop("data.frame has no barcode column and no meaningful rownames")
  }
  trs_df <- trs
} else {
  stop("unsupported TRS matrix class: ", class(trs)[1])
}
cat("       standardized -> data.frame, dim:",
    paste(dim(trs_df), collapse = " x "), "\n")
cat("       rowname head:", paste(head(rownames(trs_df), 2), collapse = ", "), "\n")
cat("       colnames     :", paste(colnames(trs_df), collapse = ", "), "\n")

# ---- [2/4] load cell_type from one SCAVENGE_small TSV ----
cat("\n[2/4] Loading cell_type from SCAVENGE_small_", REF_TRAIT,
    ".tsv.gz ...\n", sep = "")
meta_tsv <- fread(file.path(BASE_DIR,
                            paste0("SCAVENGE_small_", REF_TRAIT, ".tsv.gz")))
cat("       TSV colnames:", paste(colnames(meta_tsv), collapse = ", "), "\n")

bc_col_meta <- intersect(c("cell_id", "cell", "barcode", "V1"),
                         colnames(meta_tsv))[1]
ct_col      <- intersect(c("BioClassification", "cell_type", "celltype",
                            "ClusterName"),
                         colnames(meta_tsv))[1]
stopifnot(!is.na(bc_col_meta), !is.na(ct_col))
cat("       cell_id col -> '", bc_col_meta,
    "'; cell_type col -> '", ct_col, "'\n", sep = "")

# Build BMMC-style meta (data.frame, barcode = rownames, cell_type column)
meta_df <- data.frame(
  barcode   = as.character(meta_tsv[[bc_col_meta]]),
  cell_type = as.character(meta_tsv[[ct_col]]),
  stringsAsFactors = FALSE
)
# Also carry UMAP if present (for QC)
for (extra in c("UMAP1", "UMAP2", "UMAP_1", "UMAP_2")) {
  if (extra %in% colnames(meta_tsv)) meta_df[[extra]] <- meta_tsv[[extra]]
}
rownames(meta_df) <- meta_df$barcode
cat("       meta_df dim:", paste(dim(meta_df), collapse = " x "),
    "  cols:", paste(colnames(meta_df), collapse = ", "), "\n")

# ---- [3/4] cell_type summary ----
cat("\n[3/4] Unique cell_types in Granja (sorted by count desc):\n")
ct_summary <- as.data.frame(sort(table(meta_df$cell_type), decreasing = TRUE),
                            stringsAsFactors = FALSE)
colnames(ct_summary) <- c("cell_type", "n_cells")
ct_summary$cell_type <- as.character(ct_summary$cell_type)
ct_summary$pct <- round(100 * ct_summary$n_cells / sum(ct_summary$n_cells), 2)
print(ct_summary, row.names = FALSE)

# ---- [4/4] save ----
cat("\n[4/4] Writing outputs to:", OUT_DIR, "\n")
saveRDS(trs_df, file.path(OUT_DIR, "granja_baseline_TRS_matrix.rds"))
saveRDS(meta_df, file.path(OUT_DIR, "granja_cell_meta.rds"))
write.csv(ct_summary,
          file.path(OUT_DIR, "granja_cell_types_summary.csv"),
          row.names = FALSE)

# Sanity: TRS rownames vs meta rownames overlap
common <- intersect(rownames(trs_df), rownames(meta_df))
cat("\n=== Sanity check ===\n")
cat("  TRS n_cells       :", nrow(trs_df), "\n")
cat("  meta n_cells      :", nrow(meta_df), "\n")
cat("  common barcodes   :", length(common), "\n")
cat("  TRS-only barcodes :", nrow(trs_df) - length(common), "\n")
cat("  meta-only barcodes:", nrow(meta_df) - length(common), "\n")

cat("\n=== Granja inputs prepared, ready for eval_auroc.R ===\n")
cat("Outputs:\n",
    "  ", file.path(OUT_DIR, "granja_baseline_TRS_matrix.rds"), "\n",
    "  ", file.path(OUT_DIR, "granja_cell_meta.rds"), "\n",
    "  ", file.path(OUT_DIR, "granja_cell_types_summary.csv"), "\n", sep = "")
