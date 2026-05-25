# 03_scavenge_baseline/01_build_atac_se.R
#
# Pack BMMC ATAC counts (cells x peaks dgCMatrix) + peak coordinates into a
# RangedSummarizedExperiment with rowRanges = GRanges(peaks). This is the
# input format expected by chromVAR (addGCBias / getBackgroundPeaks / etc) and
# by SCAVENGE downstream.
#
# Inputs (from 01_data_qc/01_extract_matrices_from_h5ad.R):
#   BMMC_atac_counts.rds     dgCMatrix [cells x peaks]
#   BMMC_atac_features.rds   data.table peak_id, chr, start, end
#   BMMC_cell_meta.rds       data.table per-cell metadata (cell_type, donor, batch, ...)
#
# Output:
#   <out_dir>/03_scavenge_baseline/SE_atac.rds   RangedSummarizedExperiment

source("D:/RStudio/work/scavenge_multiomics_fusion_reproduction/reproduction_outputs/patched_scripts/module3_multiome/_bootstrap.R")
required_pkgs(c("Matrix", "data.table", "SummarizedExperiment", "GenomicRanges"))

paths <- m3_paths()
in_dir  <- file.path(paths$out_dir, "01_data_qc")
out_dir <- file.path(paths$out_dir, "03_scavenge_baseline")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

atac     <- readRDS(file.path(in_dir, "BMMC_atac_counts.rds"))
features <- readRDS(file.path(in_dir, "BMMC_atac_features.rds"))
meta     <- readRDS(file.path(in_dir, "BMMC_cell_meta.rds"))
message("Loaded: atac ", nrow(atac), " x ", ncol(atac),
        ", features ", nrow(features),
        ", meta ", nrow(meta))

# Re-parse peak coords from peak_id directly using the shared parser.
# This sidesteps any earlier parse failures cached in BMMC_atac_features.rds.
reparsed <- parse_peak_names(features$peak_id)
features$chr   <- reparsed$chr
features$start <- reparsed$start
features$end   <- reparsed$end

ok <- !is.na(features$chr) & !is.na(features$start) & !is.na(features$end)
if (!any(ok)) stop("All peaks unparseable. Sample peak_id: ",
                   paste(head(features$peak_id, 5), collapse = " | "))
if (any(!ok)) {
  message("Dropping ", sum(!ok), " peaks with unparseable coords (alt contigs)")
  features <- features[ok]
  atac <- atac[, ok, drop = FALSE]
}
message("Peaks kept: ", nrow(features))

# chromVAR expects peaks-as-rows, cells-as-columns. Our matrix is cells x peaks.
counts <- Matrix::t(atac)
rownames(counts) <- features$peak_id
colnames(counts) <- meta$barcode

peak_gr <- GenomicRanges::GRanges(
  seqnames = features$chr,
  ranges   = IRanges::IRanges(start = features$start + 1L, end = features$end),
  peak_id  = features$peak_id
)

# Use barcode as rowname so colData lookup is consistent
meta_df <- as.data.frame(meta)
rownames(meta_df) <- meta_df$barcode

SE <- SummarizedExperiment::SummarizedExperiment(
  assays    = list(counts = counts),
  rowRanges = peak_gr,
  colData   = S4Vectors::DataFrame(meta_df)
)

out_path <- file.path(out_dir, "SE_atac.rds")
saveRDS(SE, out_path)
message("Wrote: ", out_path,
        " | dim: ", nrow(SE), " peaks x ", ncol(SE), " cells")
