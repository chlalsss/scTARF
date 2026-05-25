# Build weighted GWAS variant -> ATAC peak -> gene links.
#
# This script creates the gene-mapping layer needed before trait-level gene risk
# aggregation. It keeps distal regulatory peaks by adding a nearest-TSS link
# within a fixed cis window with distance-decay weighting.

source("D:/RStudio/work/module2_gwas_gene_risk/_bootstrap.R")

required_pkgs(c(
  "data.table", "SummarizedExperiment", "GenomicRanges", "IRanges", "jsonlite"
))
suppressPackageStartupMessages({
  library(data.table)
  library(SummarizedExperiment)
  library(GenomicRanges)
  library(IRanges)
  library(jsonlite)
})

paths <- gwas_gene_paths()
ensure_gwas_gene_dirs(paths)

cis_window <- env_integer("MODULE2_CIS_WINDOW_BP", 500000)
decay_scale <- env_integer("MODULE2_DISTANCE_DECAY_BP", 50000)
promoter_upstream <- env_integer("MODULE2_PROMOTER_UPSTREAM_BP", 2000)
if (cis_window <= 0 || cis_window > 1000000) {
  stop("MODULE2_CIS_WINDOW_BP must be in (0, 1000000]", call. = FALSE)
}
if (decay_scale <= 0) {
  stop("MODULE2_DISTANCE_DECAY_BP must be positive", call. = FALSE)
}

out_dir <- file.path(paths$out_dir, "gwas_peak_gene_map")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("========== Module 2 GWAS gene-risk: build peak-gene map ==========\n")
cat("  ATAC SE       :", paths$atac_se, "\n")
cat("  GWAS          :", paths$gwas, "\n")
cat("  RNA genes     :", paths$rna_genes, "\n")
cat("  cis window bp :", cis_window, "\n")
cat("  decay bp      :", decay_scale, "\n")
cat("  promoter bp   :", promoter_upstream, "\n")
cat("  output dir    :", out_dir, "\n")

stopifnot(file.exists(paths$atac_se), file.exists(paths$gwas), file.exists(paths$rna_genes))

cat("\n========== Load inputs ==========\n")
SE_atac <- readRDS(paths$atac_se)
peak_gr <- SummarizedExperiment::rowRanges(SE_atac)
gwas <- data.table::fread(paths$gwas)
genes <- data.table::fread(paths$rna_genes)

required_gene_cols <- c("gene_id", "chrom", "chromStart", "chromEnd", "strand")
if (!all(required_gene_cols %in% colnames(genes))) {
  stop("RNA gene table must contain: ", paste(required_gene_cols, collapse = ", "), call. = FALSE)
}
if (!all(c("var", "trait", "PP") %in% colnames(gwas))) {
  stop("GWAS table must contain var, trait, and PP columns", call. = FALSE)
}

peak_id <- names(peak_gr)
if (is.null(peak_id) || any(!nzchar(peak_id))) {
  peak_id <- paste0(as.character(GenomeInfoDb::seqnames(peak_gr)), "_", start(peak_gr), "_", end(peak_gr))
}
peak_dt <- data.table(
  peak_idx = seq_along(peak_gr),
  peak_id = peak_id,
  peak_chr = as.character(GenomeInfoDb::seqnames(peak_gr)),
  peak_start = start(peak_gr),
  peak_end = end(peak_gr),
  peak_center = floor((start(peak_gr) + end(peak_gr)) / 2)
)

genes <- genes[
  !is.na(gene_id) & nzchar(gene_id) &
    !is.na(chrom) & nzchar(chrom) &
    is.finite(chromStart) & is.finite(chromEnd) &
    chromEnd > chromStart
]
genes[, `:=`(
  gene_start = as.integer(pmax(chromStart + 1L, 1L)),
  gene_end = as.integer(chromEnd),
  strand = ifelse(strand %in% c("+", "-"), strand, "+")
)]
genes[, tss := fifelse(strand == "-", gene_end, gene_start)]
genes[, body_promoter_start := fifelse(
  strand == "+",
  pmax(gene_start - promoter_upstream, 1L),
  gene_start
)]
genes[, body_promoter_end := fifelse(
  strand == "+",
  gene_end,
  gene_end + promoter_upstream
)]

gene_body_promoter_gr <- GRanges(
  seqnames = genes$chrom,
  ranges = IRanges(start = genes$body_promoter_start, end = genes$body_promoter_end),
  strand = "*"
)
gene_tss_gr <- GRanges(
  seqnames = genes$chrom,
  ranges = IRanges(start = genes$tss, end = genes$tss),
  strand = "*"
)
peak_center_gr <- GRanges(
  seqnames = peak_dt$peak_chr,
  ranges = IRanges(start = peak_dt$peak_center, end = peak_dt$peak_center),
  strand = "*"
)

cat("  peaks:", length(peak_gr), "\n")
cat("  genes:", nrow(genes), "\n")
cat("  GWAS rows:", nrow(gwas), "\n")

cat("\n========== Build peak -> gene links ==========\n")
body_hits <- GenomicRanges::findOverlaps(peak_gr, gene_body_promoter_gr, ignore.strand = TRUE)
body_links <- data.table(
  peak_idx = S4Vectors::queryHits(body_hits),
  gene_idx = S4Vectors::subjectHits(body_hits),
  link_type = "body_promoter",
  link_weight = 1.0
)

if (nrow(body_links) > 0) {
  body_links <- merge(body_links, peak_dt[, .(peak_idx, peak_center)], by = "peak_idx", all.x = TRUE)
  body_links[, distance_to_tss := abs(peak_center - genes$tss[gene_idx])]
  body_links[, peak_center := NULL]
} else {
  body_links[, distance_to_tss := numeric()]
}

nearest_idx <- GenomicRanges::nearest(peak_center_gr, gene_tss_gr, ignore.strand = TRUE)
has_nearest <- !is.na(nearest_idx)
nearest_distance <- rep(NA_integer_, length(nearest_idx))
nearest_distance[has_nearest] <- GenomicRanges::distance(
  peak_center_gr[has_nearest],
  gene_tss_gr[nearest_idx[has_nearest]]
)

nearest_links <- data.table(
  peak_idx = which(has_nearest),
  gene_idx = nearest_idx[has_nearest],
  distance_to_tss = as.integer(nearest_distance[has_nearest])
)
nearest_links <- nearest_links[distance_to_tss <= cis_window]
nearest_links[, `:=`(
  link_type = "nearest_tss_decay",
  link_weight = exp(-distance_to_tss / decay_scale)
)]

peak_gene <- rbindlist(list(body_links, nearest_links), use.names = TRUE, fill = TRUE)
if (nrow(peak_gene) == 0) {
  stop("No peak-to-gene links generated", call. = FALSE)
}

peak_gene[, priority := fifelse(link_type == "body_promoter", 1L, 2L)]
setorder(peak_gene, peak_idx, gene_idx, priority, -link_weight)
peak_gene <- peak_gene[, .SD[1], by = .(peak_idx, gene_idx)]
peak_gene[, priority := NULL]

peak_gene <- merge(
  peak_gene,
  peak_dt[, .(peak_idx, peak_id, peak_chr, peak_start, peak_end, peak_center)],
  by = "peak_idx",
  all.x = TRUE
)
gene_annot <- genes[, .(
  gene_idx = .I,
  gene_id,
  gene_chr = chrom,
  gene_start,
  gene_end,
  strand,
  tss
)]
peak_gene <- merge(peak_gene, gene_annot, by = "gene_idx", all.x = TRUE)
setcolorder(peak_gene, c(
  "peak_id", "peak_chr", "peak_start", "peak_end", "peak_center",
  "gene_id", "gene_chr", "gene_start", "gene_end", "strand", "tss",
  "link_type", "distance_to_tss", "link_weight"
))
peak_gene <- peak_gene[order(peak_chr, peak_start, gene_id)]

peak_gene_out <- file.path(out_dir, "peak_gene_weighted_links.csv.gz")
data.table::fwrite(peak_gene, peak_gene_out)

peak_gene_summary <- peak_gene[, .(
  n_links = .N,
  n_peaks = uniqueN(peak_id),
  n_genes = uniqueN(gene_id),
  median_distance_to_tss = median(distance_to_tss, na.rm = TRUE),
  mean_link_weight = mean(link_weight, na.rm = TRUE)
), by = link_type][order(link_type)]
peak_gene_summary_out <- file.path(out_dir, "peak_gene_link_summary.csv")
data.table::fwrite(peak_gene_summary, peak_gene_summary_out)
print(peak_gene_summary)

cat("\n========== Build GWAS variant -> peak -> gene links ==========\n")
if (!all(c("chr", "start", "end") %in% colnames(gwas))) {
  parsed <- parse_var_to_bed(gwas$var)
  bad_n <- sum(!parsed$ok)
  if (bad_n > 0) {
    stop("Failed to parse ", bad_n, " GWAS variants from var column", call. = FALSE)
  }
  gwas[, `:=`(chr = parsed$chr, start = parsed$start, end = parsed$end)]
  gwas[end <= start, end := start + 1L]
}
gwas[, gwas_idx := .I]

variant_gr <- GRanges(
  seqnames = gwas$chr,
  ranges = IRanges(start = as.integer(gwas$start + 1L), end = as.integer(gwas$end)),
  strand = "*"
)
variant_peak_hits <- GenomicRanges::findOverlaps(variant_gr, peak_gr, ignore.strand = TRUE)
variant_peak <- data.table(
  gwas_idx = S4Vectors::queryHits(variant_peak_hits),
  peak_idx = S4Vectors::subjectHits(variant_peak_hits)
)
if (nrow(variant_peak) == 0) {
  stop("No GWAS variants overlap ATAC peaks in SE rowRanges", call. = FALSE)
}

variant_peak <- merge(
  variant_peak,
  gwas[, .(gwas_idx, var, trait, PP = as.numeric(PP), ancestry, chr, start, end)],
  by = "gwas_idx",
  all.x = TRUE
)
variant_peak <- merge(
  variant_peak,
  peak_dt[, .(peak_idx, peak_id, peak_chr, peak_start, peak_end)],
  by = "peak_idx",
  all.x = TRUE
)
variant_peak_gene <- merge(
  variant_peak,
  peak_gene[, .(peak_idx, peak_id, gene_id, link_type, distance_to_tss, link_weight)],
  by = c("peak_idx", "peak_id"),
  allow.cartesian = TRUE
)
variant_peak_gene[, weighted_PP := PP * link_weight]
setcolorder(variant_peak_gene, c(
  "trait", "var", "PP", "weighted_PP", "ancestry",
  "chr", "start", "end", "peak_id", "peak_chr", "peak_start", "peak_end",
  "gene_id", "link_type", "distance_to_tss", "link_weight"
))
variant_peak_gene <- variant_peak_gene[order(trait, var, peak_id, gene_id)]

variant_peak_gene_out <- file.path(out_dir, "gwas_variant_peak_gene_links.csv.gz")
data.table::fwrite(variant_peak_gene, variant_peak_gene_out)

trait_summary <- variant_peak_gene[, .(
  n_gwas_rows_linked = uniqueN(gwas_idx),
  n_variants = uniqueN(var),
  n_peaks = uniqueN(peak_id),
  n_genes = uniqueN(gene_id),
  n_links = .N,
  sum_PP = sum(PP, na.rm = TRUE),
  sum_weighted_PP = sum(weighted_PP, na.rm = TRUE),
  median_distance_to_tss = median(distance_to_tss, na.rm = TRUE)
), by = trait][order(trait)]
trait_summary_out <- file.path(out_dir, "gwas_peak_gene_trait_summary.csv")
data.table::fwrite(trait_summary, trait_summary_out)

link_type_summary <- variant_peak_gene[, .(
  n_links = .N,
  n_variants = uniqueN(var),
  n_peaks = uniqueN(peak_id),
  n_genes = uniqueN(gene_id),
  sum_weighted_PP = sum(weighted_PP, na.rm = TRUE)
), by = .(trait, link_type)][order(trait, link_type)]
link_type_summary_out <- file.path(out_dir, "gwas_peak_gene_link_type_summary.csv")
data.table::fwrite(link_type_summary, link_type_summary_out)

manifest <- list(
  script = "02_build_gwas_peak_gene_map.R",
  atac_se = paths$atac_se,
  gwas = paths$gwas,
  rna_genes = paths$rna_genes,
  cis_window_bp = cis_window,
  distance_decay_bp = decay_scale,
  promoter_upstream_bp = promoter_upstream,
  peak_gene_links = peak_gene_out,
  peak_gene_summary = peak_gene_summary_out,
  gwas_variant_peak_gene_links = variant_peak_gene_out,
  trait_summary = trait_summary_out,
  link_type_summary = link_type_summary_out,
  created_at = as.character(Sys.time())
)
jsonlite::write_json(
  manifest,
  file.path(out_dir, "gwas_peak_gene_map_manifest.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

cat("\n========== Done ==========\n")
cat("  peak-gene links :", peak_gene_out, "\n")
cat("  GWAS-gene links :", variant_peak_gene_out, "\n")
cat("  trait summary   :", trait_summary_out, "\n")
