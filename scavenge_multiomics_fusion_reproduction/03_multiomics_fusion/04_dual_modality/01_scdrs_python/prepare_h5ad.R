# prepare_h5ad.R
# Step 2 of SCDRS pipeline: BMMC RNA counts (.rds, dgCMatrix) -> .h5ad
#
# 直接在 RStudio 里 source(). 必须用全新 R session (Session > Restart R).
# 通过 reticulate 直连系统 Python (D:/python, 已装 anndata/scipy/pandas/scdrs),
# 绕开 zellkonverter + basilisk (basilisk 在 Windows 国内会去 GitHub 拉 pyenv-win, 拉不下来).

# --- 0. 锁定 Python 解释器 (必须在 library(reticulate) 之前) ---
Sys.setenv(RETICULATE_PYTHON = "D:/python/python.exe")

suppressPackageStartupMessages({
  library(Matrix)
  library(reticulate)
})

cat("[0/6] python sanity check ...\n")
cat("  python:  ", reticulate::py_config()$python, "\n")
ad <- import("anndata")
sp <- import("scipy.sparse", convert = FALSE)  # keep Python object so $T$tocsr() works
np <- import("numpy",        convert = FALSE)  # keep Python object for sp.csc_matrix args
cat("  anndata: ", ad$`__version__`, "\n\n")

# --- 1. 路径 + 读 R 数据 ---
in_counts <- "D:/SCAVENGE-reproducibility-main/BMMC/01_data_qc/BMMC_rna_counts.rds"
in_genes  <- "D:/SCAVENGE-reproducibility-main/BMMC/01_data_qc/BMMC_rna_features.rds"
in_meta   <- "D:/SCAVENGE-reproducibility-main/BMMC/01_data_qc/BMMC_cell_meta.rds"
out_h5ad  <- "D:/SCAVENGE-reproducibility-main/BMMC/04_dual_modality/BMMC_rna.h5ad"

stopifnot(file.exists(in_counts), file.exists(in_genes), file.exists(in_meta))

cat("[1/6] reading inputs ...\n")
counts <- readRDS(in_counts)
genes  <- readRDS(in_genes)
meta   <- readRDS(in_meta)
cat("  counts:    ", class(counts)[1], paste(dim(counts), collapse = " x "), "\n")
cat("  genes:     ", class(genes)[1],  length(genes), "\n")
cat("  meta:      ", class(meta)[1],   paste(dim(meta), collapse = " x "), "\n")
cat("  meta cols: ", paste(colnames(meta), collapse = ", "), "\n\n")

# --- 2. orientation: gene x cell ---
cat("[2/6] orientation check ...\n")
n_cells <- nrow(meta)
if (nrow(counts) == n_cells && ncol(counts) != n_cells) {
  cat("  cell x gene -> transpose to gene x cell\n")
  counts <- t(counts)
} else if (ncol(counts) == n_cells) {
  cat("  gene x cell (already)\n")
} else {
  stop("counts dim ", paste(dim(counts), collapse = "x"),
       " doesn't match #cells ", n_cells)
}
counts <- as(counts, "CsparseMatrix")
cat("  final (gene x cell): ", paste(dim(counts), collapse = " x "), "\n")
nnz <- length(counts@x)
cat(sprintf("  nnz: %d  (sparsity %.2f%%)\n", nnz,
            100 * nnz / (as.numeric(nrow(counts)) * as.numeric(ncol(counts)))))
cat("\n")

# --- 3. align row/col names ---
cat("[3/6] aligning row/col names ...\n")
if (is.null(rownames(counts))) {
  stopifnot(length(genes) == nrow(counts))
  rownames(counts) <- as.character(genes)
}
if (is.null(colnames(counts))) {
  if (!is.null(rownames(meta))) {
    colnames(counts) <- rownames(meta)
  } else {
    bc_col <- intersect(c("barcode", "cell_id", "cellID", "CellID"),
                        colnames(meta))[1]
    if (is.na(bc_col)) stop("no cell barcode source in meta")
    colnames(counts) <- as.character(meta[[bc_col]])
  }
}
if (!is.null(rownames(meta)) && all(colnames(counts) %in% rownames(meta))) {
  meta <- meta[colnames(counts), , drop = FALSE]
}
cat("  rownames example: ", paste(head(rownames(counts), 3), collapse = ", "), "\n")
cat("  colnames example: ", paste(head(colnames(counts), 3), collapse = ", "), "\n\n")

# --- 4. dgCMatrix -> scipy.sparse.csc, then .T to get cell x gene CSR ---
cat("[4/6] transferring sparse matrix to Python ...\n")
t0 <- Sys.time()
mat_py <- sp$csc_matrix(
  tuple(
    np$asarray(counts@x,  dtype = "float32"),
    np$asarray(counts@i,  dtype = "int32"),
    np$asarray(counts@p,  dtype = "int32")
  ),
  shape = tuple(as.integer(nrow(counts)), as.integer(ncol(counts)))
)
mat_py <- mat_py$T$tocsr()  # cell x gene
cat("  shape (cell x gene): ",
    paste(py_to_r(mat_py$shape), collapse = " x "), "\n")
cat("  elapsed: ",
    round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), " s\n\n")

# --- 5. ship obs/var to Python ---
cat("[5/6] shipping obs/var to Python ...\n")
py$mat_py       <- mat_py
py$barcodes     <- as.character(colnames(counts))
py$gene_symbols <- as.character(rownames(counts))
py$meta_df      <- as.data.frame(meta, stringsAsFactors = FALSE)
py$out_h5ad     <- out_h5ad
cat("  shipped\n\n")

# --- 6. build AnnData & write h5ad (Python side) ---
cat("[6/6] constructing AnnData & writing h5ad ...\n")
t0 <- Sys.time()
py_run_string('
import os
import anndata as ad
import pandas as pd

obs = pd.DataFrame(meta_df).copy()
obs.index = pd.Index(list(barcodes), name="barcode")
var = pd.DataFrame(index=pd.Index(list(gene_symbols), name="gene_symbol"))

adata = ad.AnnData(X=mat_py, obs=obs, var=var)
adata.write_h5ad(out_h5ad, compression="gzip")

_size = os.path.getsize(out_h5ad) / 1e6
print(f"  AnnData shape: {adata.shape}")
print(f"  obs cols:      {list(adata.obs.columns)}")
print(f"  file size:     {_size:.1f} MB")
')
cat("  elapsed: ",
    round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min\n\n")

cat("prepare_h5ad.R OK.\n")
