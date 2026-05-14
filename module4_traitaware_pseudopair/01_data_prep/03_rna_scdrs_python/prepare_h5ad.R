# module4 / Phase 1.2b: Granja RNA counts -> h5ad for scDRS
# ------------------------------------------------------------------
# Uses reticulate + D:/python to write a raw-count AnnData object. scDRS will
# do its own size-factor normalization and log1p transform.

Sys.setenv(RETICULATE_PYTHON = "D:/python/python.exe")
Sys.setenv(NUMBA_CACHE_DIR = "D:/RStudio/work/module4_traitaware_pseudopair/.numba_cache")

suppressPackageStartupMessages({
  library(Matrix)
  library(SummarizedExperiment)
  library(reticulate)
})

ROOT <- "D:/RStudio/work/module4_traitaware_pseudopair"
RNA_RSE <- "D:/SCAVENGE-reproducibility-main/blood_data/scRNA-Healthy-Hematopoiesis-191120.rds"
OUT_DIR <- file.path(ROOT, "01_data_prep/03_rna_scdrs_python/input")
OUT_H5AD <- file.path(OUT_DIR, "granja_rna_raw_counts.h5ad")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(Sys.getenv("NUMBA_CACHE_DIR"), recursive = TRUE, showWarnings = FALSE)

cat(">>> Phase 1.2b Granja RNA h5ad preparation\n\n")

cat("[1/6] Python sanity check ...\n")
cfg <- reticulate::py_config()
cat("  python: ", cfg$python, "\n", sep = "")
ad <- import("anndata")
sp <- import("scipy.sparse", convert = FALSE)
np <- import("numpy", convert = FALSE)
cat("  anndata imported\n\n")

cat("[2/6] Reading RNA RangedSummarizedExperiment ...\n")
rna <- readRDS(RNA_RSE)
if (!inherits(rna, "SummarizedExperiment")) {
  stop("RNA object is not a SummarizedExperiment: ", paste(class(rna), collapse = ", "))
}
counts <- SummarizedExperiment::assay(rna, "counts")
counts <- as(counts, "CsparseMatrix")
cat("  counts: ", paste(dim(counts), collapse = " x "),
    " (gene x RNA cell), nnz = ", length(counts@x), "\n", sep = "")

cat("[3/6] Building obs / var metadata ...\n")
obs <- as.data.frame(SummarizedExperiment::colData(rna), stringsAsFactors = FALSE)
obs$barcode <- colnames(rna)
if ("BioClassification" %in% colnames(obs)) {
  obs$cell_type <- as.character(obs$BioClassification)
} else {
  obs$cell_type <- NA_character_
}
rownames(obs) <- obs$barcode
var <- as.data.frame(SummarizedExperiment::rowData(rna), stringsAsFactors = FALSE)
var$gene_symbol <- rownames(rna)
rownames(var) <- rownames(rna)
cat("  obs: ", paste(dim(obs), collapse = " x "), "\n", sep = "")
cat("  var: ", paste(dim(var), collapse = " x "), "\n", sep = "")

cat("[4/6] Transferring sparse matrix to Python ...\n")
t0 <- Sys.time()
mat_py <- sp$csc_matrix(
  tuple(
    np$asarray(counts@x, dtype = "float32"),
    np$asarray(counts@i, dtype = "int32"),
    np$asarray(counts@p, dtype = "int32")
  ),
  shape = tuple(as.integer(nrow(counts)), as.integer(ncol(counts)))
)
mat_py <- mat_py$T$tocsr()
cat("  shape (cell x gene): ",
    paste(reticulate::py_to_r(mat_py$shape), collapse = " x "), "\n", sep = "")
cat("  elapsed: ",
    round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), " s\n", sep = "")

cat("[5/6] Constructing AnnData ...\n")
py$mat_py <- mat_py
py$obs_df <- obs
py$var_df <- var
py$barcodes <- as.character(colnames(counts))
py$genes <- as.character(rownames(counts))
py$out_h5ad <- OUT_H5AD

cat("[6/6] Writing h5ad ...\n")
t0 <- Sys.time()
py_run_string('
import os
import anndata as ad
import pandas as pd

obs = pd.DataFrame(obs_df).copy()
obs.index = pd.Index(list(barcodes), name="barcode")
var = pd.DataFrame(var_df).copy()
var.index = pd.Index(list(genes), name="gene_symbol")

adata = ad.AnnData(X=mat_py, obs=obs, var=var)
adata.write_h5ad(out_h5ad, compression="gzip")

print(f"  AnnData shape: {adata.shape}")
print(f"  obs cols:      {list(adata.obs.columns)}")
print(f"  file size:     {os.path.getsize(out_h5ad) / 1e6:.1f} MB")
')
cat("  elapsed: ",
    round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min\n", sep = "")

cat("\nPhase 1.2b OK: ", OUT_H5AD, "\n", sep = "")
