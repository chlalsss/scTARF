# module4 / Phase 0.1: BMMC pairing shuffle
# ------------------------------------------------------------------
# 目的: 把 BMMC 真配对 RNA 矩阵的细胞顺序随机打乱, 制造"假非配对"用于
#       Phase 0.2 mapping 算法 sanity check. ATAC 矩阵保持不动.
#
# 输入:
#   D:/SCAVENGE-reproducibility-main/BMMC/01_data_qc/BMMC_rna_counts.rds
#   D:/SCAVENGE-reproducibility-main/BMMC/01_data_qc/BMMC_cell_meta.rds
#
# 输出 (00_sanity_bmmc/output/):
#   bmmc_rna_shuffled.rds         洗牌后 RNA counts (统一为 gene x cell sparse)
#   bmmc_meta_shuffled.rds        洗牌后 meta (行序同 rna_shuffled 列序)
#   true_pairing_groundtruth.csv  原 ATAC barcode -> 洗牌后 RNA 位置, 仅 0.3 评估用
# ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
})

in_qc_dir <- "D:/SCAVENGE-reproducibility-main/BMMC/01_data_qc"
out_dir   <- "D:/RStudio/work/module4_traitaware_pseudopair/00_sanity_bmmc/output"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- [1/6] load ----
cat("[1/6] Loading BMMC RNA counts ...\n")
rna_counts <- readRDS(file.path(in_qc_dir, "BMMC_rna_counts.rds"))
cat("       dim:", paste(dim(rna_counts), collapse = " x "),
    "  class:", class(rna_counts)[1], "\n")

cat("[2/6] Loading BMMC cell meta ...\n")
meta <- readRDS(file.path(in_qc_dir, "BMMC_cell_meta.rds"))
cat("       dim:", paste(dim(meta), collapse = " x "),
    "  class:", class(meta)[1], "\n")
cat("       columns:", paste(colnames(meta), collapse = ", "), "\n")

# ---- [3/6] detect cell_type column ----
celltype_candidates <- c("cell_type", "celltype", "CellType",
                         "BioClassification", "Annotation", "annotation",
                         "Cluster", "cluster", "celltype_l1", "celltype_l2")
hit_ct <- intersect(celltype_candidates, colnames(meta))
if (length(hit_ct) == 0) {
  stop("Cannot find cell_type column. Available: ",
       paste(colnames(meta), collapse = ", "))
}
celltype_col <- hit_ct[1]
cat("[3/6] cell_type column         :", celltype_col, "\n")

# ---- [4/6] detect barcode source + matrix orientation ----
# meta barcodes source: try common column names, then rownames as fallback
bc_cands <- c("barcode", "Barcode", "cell_id", "CellID", "cell", "cell_name")
hit_bc   <- intersect(bc_cands, colnames(meta))
if (length(hit_bc) > 0) {
  meta_barcodes <- as.character(meta[[hit_bc[1]]])
  bc_src <- paste0("meta$", hit_bc[1])
} else if (!is.null(rownames(meta))) {
  meta_barcodes <- rownames(meta)
  bc_src <- "rownames(meta)"
} else {
  stop("Cannot find barcode in meta. Columns: ",
       paste(colnames(meta), collapse = ", "))
}
cat("       barcode source         :", bc_src,
    " (head:", paste(head(meta_barcodes, 3), collapse=","), ")\n")

# Detect which axis of rna_counts holds cells.
# Cells axis = the one that matches meta barcodes.
rn <- rownames(rna_counts); if (is.null(rn)) rn <- character(0)
cn <- colnames(rna_counts); if (is.null(cn)) cn <- character(0)
r_axis_match <- length(intersect(rn, meta_barcodes))
c_axis_match <- length(intersect(cn, meta_barcodes))

cat("       barcode matches rownames(rna):", r_axis_match,
    "; colnames(rna):", c_axis_match, "\n")

if (c_axis_match >= r_axis_match && c_axis_match > 0.9 * length(meta_barcodes)) {
  orientation <- "gene_x_cell"
} else if (r_axis_match > 0.9 * length(meta_barcodes)) {
  orientation <- "cell_x_gene"
} else {
  stop("Neither axis of rna_counts matches meta barcodes (need >90%).")
}
cat("[4/6] orientation              :", orientation, "\n")

# Standardize to gene x cell internally (cells are columns)
if (orientation == "cell_x_gene") {
  cat("       transposing rna_counts to gene x cell ...\n")
  rna_counts <- Matrix::t(rna_counts)
}
cat("       after standardize, dim :", paste(dim(rna_counts), collapse = " x "), "\n")
cat("       sample gene names      :", paste(head(rownames(rna_counts), 3), collapse=","), "\n")
cat("       sample cell barcodes   :", paste(head(colnames(rna_counts), 3), collapse=","), "\n")

# ---- [5/6] align meta to rna columns ----
rna_cells <- colnames(rna_counts)
common    <- intersect(rna_cells, meta_barcodes)
cat("[5/6] aligning meta to rna columns\n")
cat("       common barcodes        :", length(common), "\n")
cat("       rna cells              :", length(rna_cells), "\n")
cat("       meta rows              :", length(meta_barcodes), "\n")
stopifnot(length(common) > 0.95 * length(rna_cells))

rna_counts <- rna_counts[, common, drop = FALSE]
meta       <- as.data.frame(meta)
meta       <- meta[match(common, meta_barcodes), , drop = FALSE]
rownames(meta) <- common

# ---- [6/6] shuffle ----
cat("[6/6] Shuffling cell order with set.seed(2026) ...\n")
set.seed(2026)
n_cells       <- ncol(rna_counts)
perm          <- sample.int(n_cells)            # new_position -> original_position
orig_barcodes <- colnames(rna_counts)
new_barcodes  <- orig_barcodes[perm]

rna_shuffled  <- rna_counts[, perm, drop = FALSE]
meta_shuffled <- meta[perm, , drop = FALSE]
rownames(meta_shuffled) <- new_barcodes

stopifnot(identical(colnames(rna_shuffled),  new_barcodes))
stopifnot(identical(rownames(meta_shuffled), new_barcodes))

# ---- ground truth: ATAC orig_pos -> RNA shuffled_pos ----
# In BMMC paired: ATAC col i  barcode == RNA col i barcode.
# After shuffling RNA, RNA cell with that barcode is at column inv[i].
inv <- order(perm)                              # original_position -> new_position
gt <- data.frame(
  atac_orig_position        = seq_len(n_cells),
  atac_barcode              = orig_barcodes,
  true_rna_pos_in_shuffled  = inv,
  true_rna_barcode_shuffled = new_barcodes[inv],
  stringsAsFactors          = FALSE
)
stopifnot(identical(gt$true_rna_barcode_shuffled, gt$atac_barcode))

# ---- save ----
saveRDS(rna_shuffled,  file.path(out_dir, "bmmc_rna_shuffled.rds"))
saveRDS(meta_shuffled, file.path(out_dir, "bmmc_meta_shuffled.rds"))
write.csv(gt, file.path(out_dir, "true_pairing_groundtruth.csv"),
          row.names = FALSE)

# ---- summary ----
cat("\n=== Phase 0.1 shuffle DONE ===\n")
cat("rna_shuffled         : ", paste(dim(rna_shuffled),  collapse = " x "),
    " (gene x cell)\n", sep = "")
cat("meta_shuffled        : ", paste(dim(meta_shuffled), collapse = " x "), "\n", sep = "")
cat("ground truth         : ", nrow(gt), " rows\n", sep = "")
cat("celltype column used : ", celltype_col, "\n", sep = "")
cat("seed                 :  2026\n")

cat("\n--- cell_type frequency after shuffle (top 10) ---\n")
print(head(sort(table(meta_shuffled[[celltype_col]]), decreasing = TRUE), 10))

identity_rate <- mean(perm == seq_len(n_cells))
cat("\nshuffle check: fraction at original position =",
    round(identity_rate, 6),
    "  (theoretical 1/N =", round(1/n_cells, 6), ")\n")
