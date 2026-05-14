# module4 / Phase 0.1.5: BMMC Gene Activity (hg38, Granja-同款几何)
# ------------------------------------------------------------------
# 目的: 给 BMMC ATAC 算 gene × cell Gene Activity 矩阵, 用作 Phase 0.2 mapping
#       ATAC 端在 G_trait 子集上的表征. 算法与 Granja
#       (D:/RStudio/work/multiomics_phase2_shared_space.R 行 105-163) 完全一致,
#       只把 EnsDb 从 v75 (hg19) 换成 v86 (hg38).
#
# 算法: 自实现 binary peak-gene overlap, 不用 Signac::GeneActivity
#   1. EnsDb v86 取 protein_coding genes, seqlevels -> UCSC, dedup by gene_name
#   2. strand-aware 2kb upstream promoter 延伸
#   3. findOverlaps(gene_regions, peak_gr) -> sparse gene×peak 0/1 indicator
#   4. gene_activity = gene_peak_map %*% atac_counts  (gene × cell)
#   5. drop rowSums == 0 的基因
#
# 输入:
#   D:/SCAVENGE-reproducibility-main/BMMC/03_scavenge_baseline/SE_atac_gcbiased.rds
#     SummarizedExperiment: rowRanges = peak GRanges, assay = peak × cell counts
#
# 输出:
#   00_sanity_bmmc/output/BMMC_gene_activity.rds  (gene × cell dgCMatrix)
#   00_sanity_bmmc/output/BMMC_gene_activity_summary.csv
# ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Matrix)
  library(SummarizedExperiment)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(IRanges)
  library(S4Vectors)
  library(ensembldb)
  library(EnsDb.Hsapiens.v86)   # hg38
})

SE_PATH <- "D:/SCAVENGE-reproducibility-main/BMMC/03_scavenge_baseline/SE_atac_gcbiased.rds"
OUT_DIR <- "D:/RStudio/work/module4_traitaware_pseudopair/00_sanity_bmmc/output"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

t0 <- Sys.time()

# ---- [1/6] load BMMC ATAC SE ----
cat("[1/6] Loading BMMC SE_atac_gcbiased.rds ...\n")
SE <- readRDS(SE_PATH)
peak_gr     <- rowRanges(SE)
atac_counts <- assay(SE)
cat("       peak_gr length     :", length(peak_gr), "\n")
cat("       atac_counts dim    :", paste(dim(atac_counts), collapse = " x "), "\n")
cat("       atac_counts class  :", class(atac_counts)[1], "\n")
cat("       peak_gr seqlevels (first 5):",
    paste(head(seqlevels(peak_gr), 5), collapse = ", "), "\n")

# Sanity: peaks should be rows of atac_counts (peak × cell)
stopifnot(length(peak_gr) == nrow(atac_counts))

# Ensure UCSC-style seqlevels on peaks (chr1, chr2, ...)
if (!grepl("^chr", seqlevels(peak_gr)[1])) {
  seqlevelsStyle(peak_gr) <- "UCSC"
  cat("       converted peak seqlevels to UCSC\n")
}

# ---- [2/6] extract hg38 protein-coding genes ----
cat("[2/6] Extracting hg38 protein_coding genes from EnsDb.Hsapiens.v86 ...\n")
gene_coords <- ensembldb::genes(
  EnsDb.Hsapiens.v86,
  filter  = ~ gene_biotype == "protein_coding",
  columns = c("gene_name", "gene_biotype")
)
cat("       raw protein_coding genes:", length(gene_coords), "\n")

# UCSC style + standard chroms
seqlevelsStyle(gene_coords) <- "UCSC"
gene_coords <- keepStandardChromosomes(gene_coords, pruning.mode = "coarse")
cat("       after keepStandardChromosomes:", length(gene_coords), "\n")

# Dedup by gene_name (keep first occurrence, same as Granja)
gene_coords <- gene_coords[!duplicated(gene_coords$gene_name)]
names(gene_coords) <- gene_coords$gene_name
cat("       after dedup by gene_name:", length(gene_coords), "\n")

# ---- [3/6] strand-aware 2kb upstream promoter extension ----
cat("[3/6] Extending 2kb upstream (strand-aware) ...\n")
gene_regions <- gene_coords
sp <- as.character(strand(gene_regions))
strand_pos  <- sp == "+"
strand_neg  <- sp == "-"
strand_star <- sp == "*"

start(gene_regions)[strand_pos | strand_star] <-
  pmax(start(gene_regions)[strand_pos | strand_star] - 2000L, 1L)
end(gene_regions)[strand_neg] <-
  end(gene_regions)[strand_neg] + 2000L

cat("       gene_regions: ", length(gene_regions),
    " (pos:", sum(strand_pos),
    ", neg:", sum(strand_neg),
    ", star:", sum(strand_star), ")\n", sep = "")

# ---- [4/6] peak-gene overlap ----
cat("[4/6] Finding peak-gene overlaps ...\n")
ov <- GenomicRanges::findOverlaps(gene_regions, peak_gr)
cat("       overlaps           :", length(ov), "\n")
cat("       unique genes hit   :", length(unique(queryHits(ov))), "\n")
cat("       unique peaks hit   :", length(unique(subjectHits(ov))), "\n")

# ---- [5/6] sparse gene×peak indicator, then matmul ----
cat("[5/6] Building gene×peak indicator and computing GA ...\n")
n_genes <- length(gene_regions)
n_peaks <- nrow(atac_counts)
gene_peak_map <- Matrix::sparseMatrix(
  i    = queryHits(ov),
  j    = subjectHits(ov),
  x    = 1,
  dims = c(n_genes, n_peaks)
)
cat("       gene_peak_map dim  :", paste(dim(gene_peak_map), collapse = " x "), "\n")
cat("       indicator nnz      :", length(ov), "\n")

# Force atac_counts to sparse if not already
if (!inherits(atac_counts, "sparseMatrix")) {
  cat("       converting atac_counts to dgCMatrix\n")
  atac_counts <- as(atac_counts, "CsparseMatrix")
}

t_mm <- Sys.time()
gene_activity <- gene_peak_map %*% atac_counts   # gene × cell
cat("       matmul elapsed     :",
    round(as.numeric(difftime(Sys.time(), t_mm, units = "secs")), 1), "s\n")

rownames(gene_activity) <- names(gene_regions)
# colnames inherited from atac_counts (cell barcodes); verify
if (is.null(colnames(gene_activity))) {
  colnames(gene_activity) <- colnames(atac_counts)
}

# Drop all-zero genes
gene_keep     <- Matrix::rowSums(gene_activity) > 0
n_kept        <- sum(gene_keep)
gene_activity <- gene_activity[gene_keep, , drop = FALSE]
cat("       kept genes (rowSum > 0):", n_kept, " / ", n_genes, "\n", sep = "")
cat("       final GA dim       :", paste(dim(gene_activity), collapse = " x "),
    " (gene × cell)\n")

# ---- [6/6] save + summary ----
cat("[6/6] Saving outputs ...\n")
saveRDS(gene_activity,
        file.path(OUT_DIR, "BMMC_gene_activity.rds"))

# quick summary
ga_cell_total <- Matrix::colSums(gene_activity)
ga_gene_total <- Matrix::rowSums(gene_activity)
ga_nnz_per_cell <- Matrix::colSums(gene_activity > 0)

summary_df <- data.frame(
  metric = c("n_genes", "n_cells", "matrix_nnz",
             "median_GA_per_cell", "median_genes_active_per_cell",
             "median_GA_per_gene", "mean_GA_per_cell"),
  value  = c(nrow(gene_activity), ncol(gene_activity),
             length(gene_activity@x),
             median(ga_cell_total),
             median(ga_nnz_per_cell),
             median(ga_gene_total),
             mean(ga_cell_total))
)
write.csv(summary_df,
          file.path(OUT_DIR, "BMMC_gene_activity_summary.csv"),
          row.names = FALSE)

cat("\n=== Phase 0.1.5 BMMC Gene Activity DONE ===\n")
print(summary_df)
cat("\nTotal elapsed: ",
    round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "s\n")
cat("Outputs in: ", OUT_DIR, "\n")

# Sanity print: are key marker genes present?
markers <- c("CD14", "FCGR3A", "GATA1", "GP1BA", "CD8A", "CD4", "MS4A1", "GNLY")
present <- markers[markers %in% rownames(gene_activity)]
absent  <- setdiff(markers, present)
cat("\nMarker gene presence check (",
    length(present), "/", length(markers), " present)\n", sep = "")
if (length(present) > 0) cat("  present: ", paste(present, collapse = ", "), "\n")
if (length(absent)  > 0) cat("  absent : ", paste(absent,  collapse = ", "), "\n")
