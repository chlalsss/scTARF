# 00_precheck/02_celltype_coverage.R
#
# Compute trait x cell_type x donor coverage for BMMC multiome metadata.
# Apply gating thresholds from m3_gating().
#
# Input:
#   <out_dir>/00_precheck/BMMC_meta_raw.csv  (from 01_fetch_metadata.R)
#
# Outputs:
#   <out_dir>/00_precheck/BMMC_celltype_coverage.csv
#     columns: trait, cell_type, donor, n_cells
#   <out_dir>/00_precheck/BMMC_celltype_coverage_summary.csv
#     columns: trait, n_cells_total, n_donors_with_>=min, pass
#   <out_dir>/00_precheck/figures/celltype_coverage.png
#
# Gating: per trait, total target cells >= M3_MIN_CELLS_TOTAL (200) AND each of
# 4 donors >= M3_MIN_CELLS_DONOR (50).

source(file.path(dirname(sys.frame(1)$ofile %||% "."), "..", "_bootstrap.R"))
required_pkgs(c("data.table", "ggplot2"))

paths  <- m3_paths()
gating <- m3_gating()

if (!file.exists(paths$meta_raw)) {
  stop("Missing ", paths$meta_raw, " — run 01_fetch_metadata.R first.")
}
meta <- data.table::fread(paths$meta_raw)
message("Loaded metadata: ", nrow(meta), " cells")

# --- 1. verify cell_type vocabulary matches is_bmmc_target() ---------------

obs_types <- sort(unique(as.character(meta$cell_type)))
target_types <- bmmc_target_types()
declared_types <- unique(unlist(target_types, use.names = FALSE))

unseen_declared <- setdiff(declared_types, obs_types)
if (length(unseen_declared) > 0) {
  message(
    "Warning: declared target types not present in metadata -> ",
    paste(unseen_declared, collapse = ", "),
    ". Update bmmc_target_types() in _config.R after inspecting unique cell_type values:"
  )
  print(obs_types)
}

# --- 2. expand to (trait, cell_type, donor, n_cells) -----------------------

per_donor <- meta[, .N, by = .(cell_type, donor)]

rows <- list()
for (.tt in names(target_types)) {
  hit <- per_donor[cell_type %in% target_types[[.tt]]]
  if (nrow(hit) == 0) next
  hit[, trait := .tt]
  rows[[.tt]] <- hit
}
coverage <- data.table::rbindlist(rows, fill = TRUE)
data.table::setcolorder(coverage, c("trait", "cell_type", "donor", "N"))
data.table::setnames(coverage, "N", "n_cells")
data.table::fwrite(coverage, paths$coverage_csv)
message("Wrote ", paths$coverage_csv, " (", nrow(coverage), " rows)")

# --- 3. summary + gating ---------------------------------------------------

donors_all <- sort(unique(meta$donor))
summary_rows <- list()
for (.tt in names(target_types)) {
  trait_rows <- coverage[trait == .tt]
  per_donor_n <- trait_rows[, .(n = sum(n_cells)), by = donor]
  total_n <- sum(per_donor_n$n)
  donors_pass <- sum(per_donor_n$n >= gating$min_cells_per_target_donor)
  donors_total <- length(donors_all)
  pass <- (total_n >= gating$min_cells_per_target_total) &&
    (donors_pass == donors_total)

  summary_rows[[.tt]] <- data.table::data.table(
    trait                  = .tt,
    n_cells_total          = total_n,
    n_donors_total         = donors_total,
    n_donors_with_min      = donors_pass,
    threshold_total        = gating$min_cells_per_target_total,
    threshold_per_donor    = gating$min_cells_per_target_donor,
    pass                   = pass
  )
}
summary_dt <- data.table::rbindlist(summary_rows)
summary_path <- sub("\\.csv$", "_summary.csv", paths$coverage_csv)
data.table::fwrite(summary_dt, summary_path)
message("Wrote ", summary_path)
print(summary_dt)

# --- 4. plot ---------------------------------------------------------------

plot_dt <- coverage[, .(n_cells = sum(n_cells)), by = .(trait, donor)]
p <- ggplot2::ggplot(
  plot_dt,
  ggplot2::aes(x = trait, y = n_cells, fill = donor)
) +
  ggplot2::geom_col(position = ggplot2::position_stack()) +
  ggplot2::geom_hline(
    yintercept = gating$min_cells_per_target_total,
    linetype = "dashed", color = "red"
  ) +
  ggplot2::labs(
    title = "BMMC target cell type coverage by trait",
    subtitle = sprintf(
      "dashed = total threshold (%d); each donor must individually have >= %d",
      gating$min_cells_per_target_total,
      gating$min_cells_per_target_donor
    ),
    x = "Trait", y = "Cells"
  ) +
  ggplot2::theme_bw()

fig_path <- file.path(paths$fig_dir, "celltype_coverage.png")
ggplot2::ggsave(fig_path, p, width = 6, height = 4, dpi = 150)
message("Wrote ", fig_path)

# --- 5. final verdict ------------------------------------------------------

if (all(summary_dt$pass)) {
  message("\nPRECHECK 1 PASS — all traits meet coverage thresholds.")
} else {
  failing <- summary_dt[pass == FALSE, trait]
  message(
    "\nPRECHECK 1 BLOCKED — failing trait(s): ", paste(failing, collapse = ", "),
    ". Either expand bmmc_target_types() in _config.R or change trait/dataset."
  )
}
