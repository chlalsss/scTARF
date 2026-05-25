# 02: Subtract RNA-incompatible cells from the original ATAC z-score seed.
#
# Inputs:
#   - cell_trait_rna_target_fraction.csv     # from 01
#   - seed_signal_by_cell.csv.gz             # original ATAC z-score + seed flag
#                                            # (from module2_gwas_gene_risk)
#
# Outputs:
#   - seed_with_rna_fraction.csv             # one row per (trait, original seed cell)
#   - seed_filter_summary.csv                # one row per (trait, threshold)
#   - filtered_seed_cells_<trait>_t<NN>.csv  # filtered seed list per (trait, threshold)
#   - seed_filter_manifest.json
#
# Parameters:
#   MODULE2_TRAITS                : traits to process (default "MONO,PLT,LYM")
#   MODULE2_SEED_FILTER_THRESHOLDS: thresholds to scan (default
#                                   "0.30,0.40,0.50,0.60,0.70,0.80")
#   MODULE2_SEED_FILTER_PRIMARY   : primary threshold for downstream
#                                   diagnostics (default 0.50)

source("D:/RStudio/work/module2_seed_rna_filter/_bootstrap.R")

required_pkgs(c("data.table", "jsonlite"))
suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

paths <- seed_filter_paths()
ensure_seed_filter_dirs(paths)

traits <- parse_csv_env("MODULE2_TRAITS", "MONO,PLT,LYM")
threshold_strs <- parse_csv_env("MODULE2_SEED_FILTER_THRESHOLDS", "0.30,0.40,0.50,0.60,0.70,0.80")
thresholds <- as.numeric(threshold_strs)
if (any(!is.finite(thresholds)) || any(thresholds < 0) || any(thresholds > 1)) {
  stop("MODULE2_SEED_FILTER_THRESHOLDS must contain values in [0, 1]", call. = FALSE)
}
primary_threshold <- env_numeric("MODULE2_SEED_FILTER_PRIMARY", 0.50)

cat("========== Module 2 seed RNA filter: 02 filter seed ==========\n")
cat("  traits     :", paste(traits, collapse = ", "), "\n")
cat("  thresholds :", paste(thresholds, collapse = ", "), "\n")
cat("  primary    :", primary_threshold, "\n")
cat("  out_dir    :", paths$out_dir, "\n")

# Load --------------------------------------------------------------------

trait_target_path <- file.path(paths$out_dir, "cell_trait_rna_target_fraction.csv")
if (!file.exists(trait_target_path)) {
  stop("Missing 01 output: ", trait_target_path,
       ". Run 01_build_rna_inferred_type.R first.", call. = FALSE)
}
if (!file.exists(paths$seed_signal_by_cell)) {
  stop("Missing seed_signal_by_cell: ", paths$seed_signal_by_cell, call. = FALSE)
}

cat("\n[1/4] Load inputs\n")
trait_fraction_dt <- fread(trait_target_path)
seed_signal_dt <- fread(paths$seed_signal_by_cell)
cat("  trait_fraction rows:", nrow(trait_fraction_dt), "\n")
cat("  seed_signal rows   :", nrow(seed_signal_dt), "\n")

required_signal_cols <- c("trait", "cell", "z_score", "seed", "BioClassification", "target_type")
missing_cols <- setdiff(required_signal_cols, names(seed_signal_dt))
if (length(missing_cols) > 0) {
  stop("seed_signal_by_cell missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

# Restrict to requested traits
trait_fraction_dt <- trait_fraction_dt[trait %in% traits]
seed_signal_dt <- seed_signal_dt[trait %in% traits]
seed_signal_dt[, seed := as.logical(seed)]
seed_signal_dt[, target_type := as.logical(target_type)]

# Join: every seed cell gets an rna target_fraction for its own trait ------

cat("\n[2/4] Join seed signal + RNA target fraction\n")
joined_dt <- merge(
  seed_signal_dt,
  trait_fraction_dt[, .(atac_cell, trait, n_neighbors, n_target_neighbors, target_fraction)],
  by.x = c("cell", "trait"),
  by.y = c("atac_cell", "trait"),
  all.x = TRUE
)
n_seed_no_neighbors <- joined_dt[seed == TRUE & is.na(target_fraction), .N]
if (n_seed_no_neighbors > 0) {
  warning("Found ", n_seed_no_neighbors, " seed cells without RNA neighbors; treating target_fraction as 0.")
  joined_dt[is.na(target_fraction), target_fraction := 0]
}

setnames(joined_dt, "target_type", "target_type_atac")
setnames(joined_dt, "target_fraction", "rna_target_fraction")

# Restrict export to seed cells only (for size; full join not needed)
seed_with_rna_dt <- joined_dt[seed == TRUE, .(
  trait, cell, z_score, BioClassification,
  target_type_atac,
  n_neighbors, n_target_neighbors, rna_target_fraction
)]
seed_out <- file.path(paths$out_dir, "seed_with_rna_fraction.csv")
fwrite(seed_with_rna_dt, seed_out)
cat("  wrote:", seed_out, "(", nrow(seed_with_rna_dt), "rows)\n")

# Threshold sweep summary -------------------------------------------------

cat("\n[3/4] Threshold sweep\n")
summary_list <- list()
filtered_cells_list <- list()
for (tr in traits) {
  for (th in thresholds) {
    sub <- seed_with_rna_dt[trait == tr]
    n_orig <- nrow(sub)
    if (n_orig == 0) next
    keep_idx <- sub$rna_target_fraction >= th
    n_filt <- sum(keep_idx)
    purity_orig <- mean(sub$target_type_atac, na.rm = TRUE)
    purity_filt <- if (n_filt > 0) mean(sub$target_type_atac[keep_idx], na.rm = TRUE) else NA_real_
    median_z_orig <- median(sub$z_score, na.rm = TRUE)
    median_z_filt <- if (n_filt > 0) median(sub$z_score[keep_idx], na.rm = TRUE) else NA_real_

    summary_list[[paste(tr, th, sep = "_")]] <- data.table(
      trait = tr,
      threshold = th,
      n_seed_original = n_orig,
      n_seed_filtered = n_filt,
      kept_pct = if (n_orig > 0) n_filt / n_orig else NA_real_,
      target_purity_original = purity_orig,
      target_purity_filtered = purity_filt,
      delta_purity_pp = if (is.finite(purity_filt)) (purity_filt - purity_orig) * 100 else NA_real_,
      median_z_original = median_z_orig,
      median_z_filtered = median_z_filt,
      target_definition = target_definition(tr)
    )

    filtered_cells_dt <- sub[keep_idx, .(trait, cell, z_score, BioClassification,
                                         target_type_atac, rna_target_fraction)]
    th_label <- gsub("\\.", "p", sprintf("%.2f", th))
    filtered_out <- file.path(
      paths$out_dir,
      sprintf("filtered_seed_cells_%s_t%s.csv", tr, th_label)
    )
    fwrite(filtered_cells_dt, filtered_out)
    filtered_cells_list[[paste(tr, th, sep = "_")]] <- list(
      trait = tr, threshold = th, file = filtered_out, n = nrow(filtered_cells_dt)
    )
  }
}
summary_dt <- rbindlist(summary_list, use.names = TRUE)
setorder(summary_dt, trait, threshold)
summary_out <- file.path(paths$out_dir, "seed_filter_summary.csv")
fwrite(summary_dt, summary_out)
cat("  wrote:", summary_out, "\n")
cat("\n  Threshold sweep summary:\n")
print(summary_dt)

# Manifest ----------------------------------------------------------------

cat("\n[4/4] Write manifest\n")
manifest <- list(
  script = "02_filter_seed_by_rna.R",
  traits = traits,
  thresholds = thresholds,
  primary_threshold = primary_threshold,
  inputs = list(
    cell_trait_target_fraction = trait_target_path,
    seed_signal_by_cell = paths$seed_signal_by_cell
  ),
  outputs = list(
    seed_with_rna_fraction = seed_out,
    threshold_sweep_summary = summary_out,
    filtered_seed_files = filtered_cells_list
  ),
  notes = c(
    "Filter rule: keep seed iff rna_target_fraction >= threshold.",
    "rna_target_fraction = fraction of top-K RNA neighbors with cell_type matching trait target rule.",
    "Independent verification only; full TRS run is downstream of seed-purity acceptance."
  ),
  created_at = as.character(Sys.time())
)
manifest_out <- file.path(paths$out_dir, "seed_filter_manifest.json")
write_json(manifest, manifest_out, pretty = TRUE, auto_unbox = TRUE)

cat("\n========== Done ==========\n")
cat("  primary threshold for 03 diagnostics:", primary_threshold, "\n")
cat("  manifest:", manifest_out, "\n")
