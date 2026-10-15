# 00_precheck/03_variant_peak_overlap.R
#
# Compute fine-mapped variant <-> BMMC peak overlap, weighted by PP, per trait.
# Gating: per trait, sum(PP of variants overlapping any peak) / sum(PP of all
# variants) >= M3_MIN_PP_OVERLAP (0.30).
#
# Input:
#   paths$gwas_tsv  CS.trans_all.PP01.strong_atac_overlap.tsv (existing project file)
#   paths$peaks_bed BMMC peak bed (3-col: chr start end)
#                   Caller responsibility: drop this file in place before running.
#                   See 00_precheck/README.md for sourcing options.
#
# Output:
#   <out_dir>/00_precheck/BMMC_finemap_overlap_summary.csv
#     columns: trait, n_variants, n_overlap, pp_sum_total, pp_sum_overlap,
#              pp_overlap_frac, pass

source(file.path("scripts", "_bootstrap.R"))
required_pkgs(c("data.table", "GenomicRanges"))

paths  <- m3_paths()
gating <- m3_gating()

if (!file.exists(paths$gwas_tsv)) {
  stop("Missing GWAS TSV: ", paths$gwas_tsv,
       ". Set M3_GWAS_TSV or copy CS.trans_all.PP01.strong_atac_overlap.tsv into place.")
}
if (!file.exists(paths$peaks_bed)) {
  stop("Missing BMMC peak bed: ", paths$peaks_bed,
       ". See 00_precheck/README.md for how to obtain it.")
}

# --- 1. load GWAS variants -------------------------------------------------

gwas <- data.table::fread(paths$gwas_tsv)
message("GWAS rows: ", nrow(gwas), " | cols: ", paste(colnames(gwas), collapse = ", "))

# Older copies of this TSV use "var" instead of "variant" for the variant ID
# column. Normalize so downstream code can assume "variant".
if ("var" %in% colnames(gwas) && !"variant" %in% colnames(gwas)) {
  data.table::setnames(gwas, "var", "variant")
}

# expected columns from CS.trans_all.PP01.strong_atac_overlap.tsv:
# variant, trait, PP (or pip), ...
needed <- c("variant", "trait")
miss <- setdiff(needed, colnames(gwas))
if (length(miss) > 0) {
  stop("GWAS file missing expected columns: ", paste(miss, collapse = ", "))
}
pp_col <- intersect(c("PP", "pip", "PIP", "posterior"), colnames(gwas))[1]
if (is.na(pp_col)) {
  stop("Could not find a posterior-probability column (looked for PP/pip/PIP/posterior).")
}
message("Using PP column: ", pp_col)

# Trait subset: only the three blood traits we plan to study downstream.
keep_traits <- c("MONO", "PLT", "LYM")
gwas <- gwas[trait %in% keep_traits]
message("After trait filter (MONO/PLT/LYM): ", nrow(gwas))

bed_dt <- parse_var_to_bed(gwas$variant)
bad <- !bed_dt$ok
if (any(bad)) {
  message("Dropping ", sum(bad), " unparseable variants")
  gwas   <- gwas[!bad]
  bed_dt <- bed_dt[!bad]
}

variant_gr <- GenomicRanges::GRanges(
  seqnames = bed_dt$chr,
  ranges   = IRanges::IRanges(start = bed_dt$start + 1L, end = bed_dt$end),
  variant  = gwas$variant,
  trait    = gwas$trait,
  pp       = gwas[[pp_col]]
)

# --- 2. load BMMC peaks ----------------------------------------------------

peak_dt <- data.table::fread(
  paths$peaks_bed,
  header = FALSE,
  col.names = c("chr", "start", "end")[seq_len(min(3, ncol(data.table::fread(paths$peaks_bed, header = FALSE, nrows = 1))))]
)
# Force exactly 3 cols
peak_dt <- peak_dt[, 1:3]
data.table::setnames(peak_dt, c("chr", "start", "end"))
peak_dt <- peak_dt[!is.na(start) & !is.na(end)]

if (!any(grepl("^chr", peak_dt$chr))) {
  peak_dt[, chr := paste0("chr", chr)]
}

peak_gr <- GenomicRanges::GRanges(
  seqnames = peak_dt$chr,
  ranges   = IRanges::IRanges(start = peak_dt$start + 1L, end = peak_dt$end)
)
message("BMMC peaks: ", length(peak_gr))

# --- 3. overlap ------------------------------------------------------------

hits <- GenomicRanges::findOverlaps(variant_gr, peak_gr)
overlap_idx <- unique(S4Vectors::queryHits(hits))
variant_gr$overlap <- seq_along(variant_gr) %in% overlap_idx

# --- 4. per-trait summary --------------------------------------------------

dt <- data.table::data.table(
  trait    = variant_gr$trait,
  pp       = as.numeric(variant_gr$pp),
  overlap  = variant_gr$overlap
)

summary_dt <- dt[, .(
  n_variants      = .N,
  n_overlap       = sum(overlap),
  pp_sum_total    = sum(pp, na.rm = TRUE),
  pp_sum_overlap  = sum(pp[overlap], na.rm = TRUE)
), by = trait]
summary_dt[, pp_overlap_frac := pp_sum_overlap / pp_sum_total]
summary_dt[, threshold := gating$min_pp_sum_overlap]
summary_dt[, pass := pp_sum_overlap >= gating$min_pp_sum_overlap]

data.table::fwrite(summary_dt, paths$overlap_csv)
message("Wrote ", paths$overlap_csv)
print(summary_dt)

# --- 5. verdict ------------------------------------------------------------

if (all(summary_dt$pass)) {
  message("\nPRECHECK 2 PASS — all traits meet PP overlap threshold.")
} else {
  failing <- summary_dt[pass == FALSE, trait]
  message(
    "\nPRECHECK 2 BLOCKED — failing trait(s): ", paste(failing, collapse = ", "),
    ". Options: (a) LD-expand variant intervals (+/- 1kb) and rerun; ",
    "(b) drop the trait from downstream; (c) switch dataset."
  )
}
