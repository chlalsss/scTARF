# module1_method_B_export_scSHEFT_input.R
# 模块 1 · 方法 B：为官方 scSHEFT 代码导出输入中间件
#
# 官方 scSHEFT 需要：
#   1. raw ATAC peak h5ad
#   2. ATAC Gene Activity h5ad
#   3. RNA gene expression h5ad
#   4. ATAC LSI npy
#   5. ATAC structural PCA npy
#
# 本脚本只负责从现有 RDS 产物导出共享基因层面的 RNA / Gene Activity
# MatrixMarket + 元数据。h5ad / npy 在 scSHEFT Python 环境中生成。
#
# 运行方式：
#   source("D:/RStudio/work/module1_method_B_export_scSHEFT_input.R")

required_pkgs <- c("Seurat", "Matrix", "data.table", "jsonlite")
missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("缺少包: ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(jsonlite)
})

set.seed(9527)
script_start <- Sys.time()

data_dir <- "D:/SCAVENGE-reproducibility-main/blood_data"
in_dir   <- file.path(data_dir, "multiomics_output_v2")
out_dir  <- file.path(data_dir, "multiomics_module1_pairing", "method_B_scSHEFT", "scsheft_input")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("========== scSHEFT 输入导出 ==========\n")
cat("  输入:", in_dir, "\n")
cat("  输出:", out_dir, "\n")

cat("\n========== 加载 RDS ==========\n")
t0 <- Sys.time()
ga_mat     <- readRDS(file.path(in_dir, "gene_activity_matrix.rds"))
seurat_rna <- readRDS(file.path(in_dir, "seurat_rna_normalized.rds"))
atac_meta  <- readRDS(file.path(in_dir, "atac_metadata.rds"))
rna_meta   <- as.data.frame(seurat_rna@meta.data)
rna_counts <- GetAssayData(seurat_rna, assay = "RNA", layer = "counts")

cat("  Gene Activity:", nrow(ga_mat), "genes x", ncol(ga_mat), "ATAC cells\n")
cat("  RNA counts   :", nrow(rna_counts), "genes x", ncol(rna_counts), "RNA cells\n")
cat(sprintf("  [计时] 加载: %.1f 秒\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

cat("\n========== 共享基因对齐 ==========\n")
shared_genes <- intersect(rownames(rna_counts), rownames(ga_mat))
shared_genes <- sort(shared_genes)
if (length(shared_genes) < 100) {
  stop("共享基因过少: ", length(shared_genes))
}
rna_counts <- rna_counts[shared_genes, , drop = FALSE]
ga_mat     <- ga_mat[shared_genes, , drop = FALSE]
cat("  共享基因数:", length(shared_genes), "\n")

cat("\n========== 写 MatrixMarket ==========\n")
t0 <- Sys.time()
Matrix::writeMM(rna_counts, file.path(out_dir, "rna_counts_shared_genes.mtx"))
Matrix::writeMM(ga_mat,     file.path(out_dir, "gene_activity_counts_shared_genes.mtx"))
cat(sprintf("  [计时] 写矩阵: %.1f 秒\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

cat("\n========== 写基因与元数据 ==========\n")
gene_df <- data.frame(gene = shared_genes, stringsAsFactors = FALSE)
fwrite(gene_df, file.path(out_dir, "shared_genes.csv"))

rna_cells <- data.frame(
  cell = rownames(rna_meta),
  cell_type = as.character(rna_meta$BioClassification),
  stringsAsFactors = FALSE
)
atac_cells <- data.frame(
  cell = rownames(atac_meta),
  cell_type = as.character(atac_meta$BioClassification),
  stringsAsFactors = FALSE
)
fwrite(rna_cells,  file.path(out_dir, "rna_cells.csv"))
fwrite(atac_cells, file.path(out_dir, "atac_cells.csv"))

manifest <- list(
  method = "scSHEFT",
  source = "module1_method_B_export_scSHEFT_input.R",
  input_dir = in_dir,
  output_dir = out_dir,
  n_shared_genes = length(shared_genes),
  n_rna_cells = ncol(rna_counts),
  n_atac_cells = ncol(ga_mat),
  files = list(
    rna_counts_shared_genes_mtx = file.path(out_dir, "rna_counts_shared_genes.mtx"),
    gene_activity_counts_shared_genes_mtx = file.path(out_dir, "gene_activity_counts_shared_genes.mtx"),
    shared_genes_csv = file.path(out_dir, "shared_genes.csv"),
    rna_cells_csv = file.path(out_dir, "rna_cells.csv"),
    atac_cells_csv = file.path(out_dir, "atac_cells.csv")
  ),
  note = paste(
    "Use the Python scSHEFT environment to convert these MatrixMarket files",
    "and glue_input peak files into h5ad/npy inputs required by the official repo."
  ),
  total_seconds = as.numeric(difftime(Sys.time(), script_start, units = "secs"))
)
jsonlite::write_json(manifest, file.path(out_dir, "export_manifest.json"),
                     pretty = TRUE, auto_unbox = TRUE)

cat("\n======================================================\n")
cat("  scSHEFT 输入中间件导出完成\n")
cat("======================================================\n")
cat("  共享基因:", length(shared_genes), "\n")
cat("  RNA cells:", ncol(rna_counts), "\n")
cat("  ATAC cells:", ncol(ga_mat), "\n")
cat("  输出目录:", out_dir, "\n")
cat("======================================================\n")
