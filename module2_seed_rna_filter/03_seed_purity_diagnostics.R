# 03: Independent seed-purity diagnostics. No full TRS run.
#
# Two acceptance signals on filtered vs original seed:
#   (a) target_type purity on BioClassification (i.e. "is the seed in the
#       trait target cell type?"); expect filtered >= original + 5pp.
#   (b) hop1 target_fraction on the original LSI mutualKNN graph; expect
#       filtered >= original + 1pp.
#
# Inputs:
#   - seed_with_rna_fraction.csv   # from 02
#   - seed_filter_summary.csv      # from 02 (used to know thresholds scanned)
#   - mutualknn_k30.rds            # original LSI mutualKNN graph
#   - seed_signal_by_cell.csv.gz   # for full target_type vector across all cells
#
# Outputs:
#   - seed_purity_summary.csv
#   - seed_purity_diagnostics_manifest.json
#   - figures/Fig_seed_purity_threshold_sweep.png
#
# Parameters:
#   MODULE2_TRAITS                : traits to process
#   MODULE2_SEED_FILTER_THRESHOLDS: thresholds to evaluate

source("D:/RStudio/work/module2_seed_rna_filter/_bootstrap.R")

required_pkgs(c("data.table", "Matrix", "ggplot2", "jsonlite"))
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(ggplot2)
  library(jsonlite)
})

paths <- seed_filter_paths()
ensure_seed_filter_dirs(paths)

traits <- parse_csv_env("MODULE2_TRAITS", "MONO,PLT,LYM")
threshold_strs <- parse_csv_env("MODULE2_SEED_FILTER_THRESHOLDS", "0.30,0.40,0.50,0.60,0.70,0.80")
thresholds <- as.numeric(threshold_strs)

cat("========== Module 2 seed RNA filter: 03 purity diagnostics ==========\n")
cat("  traits     :", paste(traits, collapse = ", "), "\n")
cat("  thresholds :", paste(thresholds, collapse = ", "), "\n")

seed_with_rna_path <- file.path(paths$out_dir, "seed_with_rna_fraction.csv")
if (!file.exists(seed_with_rna_path)) {
  stop("Missing 02 output: ", seed_with_rna_path, ". Run 02_filter_seed_by_rna.R first.", call. = FALSE)
}
if (!file.exists(paths$original_knn)) {
  stop("Missing original KNN: ", paths$original_knn, call. = FALSE)
}
if (!file.exists(paths$seed_signal_by_cell)) {
  stop("Missing seed_signal_by_cell: ", paths$seed_signal_by_cell, call. = FALSE)
}

# Load --------------------------------------------------------------------

cat("\n[1/4] Load inputs\n")
seed_with_rna_dt <- fread(seed_with_rna_path)
seed_signal_dt <- fread(paths$seed_signal_by_cell)
seed_signal_dt[, target_type := as.logical(target_type)]
original_knn <- readRDS(paths$original_knn)
adj <- as_binary_graph(original_knn)
cells_in_graph <- rownames(adj)
cat("  graph cells:", length(cells_in_graph), "\n")

# Hop1 helper -------------------------------------------------------------

hop1_target_fraction <- function(seed_cells, target_vec_named) {
  seed_in_graph <- intersect(seed_cells, cells_in_graph)
  if (length(seed_in_graph) == 0) return(list(n_seed = 0, n_hop1 = 0, frac = NA_real_))
  seed_idx <- match(seed_in_graph, cells_in_graph)
  hop1 <- Matrix::rowSums(adj[, seed_idx, drop = FALSE]) > 0
  hop1_logical <- hop1
  hop1_logical[seed_idx] <- FALSE
  n_hop1 <- sum(hop1_logical)
  if (n_hop1 == 0) return(list(n_seed = length(seed_in_graph), n_hop1 = 0, frac = NA_real_))
  hop1_cells <- cells_in_graph[hop1_logical]
  tg <- target_vec_named[hop1_cells]
  list(
    n_seed = length(seed_in_graph),
    n_hop1 = n_hop1,
    frac = mean(tg, na.rm = TRUE)
  )
}

# Per-trait diagnostics ---------------------------------------------------

cat("\n[2/4] Compute purity + hop1 enrichment for original and each threshold\n")
out_rows <- list()
for (tr in traits) {
  cat("  trait:", tr, "\n")
  signal_tr <- seed_signal_dt[trait == tr]
  if (nrow(signal_tr) == 0) {
    warning("seed_signal_by_cell has no rows for trait ", tr)
    next
  }
  target_vec_named <- setNames(as.logical(signal_tr$target_type), signal_tr$cell)
  background_target <- mean(target_vec_named, na.rm = TRUE)

  seed_tr <- seed_with_rna_dt[trait == tr]
  if (nrow(seed_tr) == 0) next

  # Original (no filter)
  orig_cells <- seed_tr$cell
  orig_purity <- mean(seed_tr$target_type_atac, na.rm = TRUE)
  orig_hop1 <- hop1_target_fraction(orig_cells, target_vec_named)
  out_rows[[paste(tr, "original", sep = "_")]] <- data.table(
    trait = tr, threshold = NA_real_, label = "original",
    n_seed = orig_hop1$n_seed,
    target_purity = orig_purity,
    delta_purity_pp = 0,
    hop1_n = orig_hop1$n_hop1,
    hop1_target_fraction = orig_hop1$frac,
    delta_hop1_pp = 0,
    background_target = background_target,
    target_definition = target_definition(tr)
  )

  # Each threshold
  for (th in thresholds) {
    keep <- seed_tr$rna_target_fraction >= th
    if (sum(keep) == 0) {
      out_rows[[paste(tr, th, sep = "_")]] <- data.table(
        trait = tr, threshold = th, label = sprintf("filt_%.2f", th),
        n_seed = 0, target_purity = NA_real_, delta_purity_pp = NA_real_,
        hop1_n = 0, hop1_target_fraction = NA_real_, delta_hop1_pp = NA_real_,
        background_target = background_target,
        target_definition = target_definition(tr)
      )
      next
    }
    filt_cells <- seed_tr$cell[keep]
    filt_purity <- mean(seed_tr$target_type_atac[keep], na.rm = TRUE)
    filt_hop1 <- hop1_target_fraction(filt_cells, target_vec_named)
    out_rows[[paste(tr, th, sep = "_")]] <- data.table(
      trait = tr, threshold = th, label = sprintf("filt_%.2f", th),
      n_seed = filt_hop1$n_seed,
      target_purity = filt_purity,
      delta_purity_pp = (filt_purity - orig_purity) * 100,
      hop1_n = filt_hop1$n_hop1,
      hop1_target_fraction = filt_hop1$frac,
      delta_hop1_pp = if (is.finite(filt_hop1$frac) && is.finite(orig_hop1$frac)) (filt_hop1$frac - orig_hop1$frac) * 100 else NA_real_,
      background_target = background_target,
      target_definition = target_definition(tr)
    )
  }
}
purity_dt <- rbindlist(out_rows, use.names = TRUE, fill = TRUE)
setorder(purity_dt, trait, threshold)
purity_out <- file.path(paths$out_dir, "seed_purity_summary.csv")
fwrite(purity_dt, purity_out)
cat("\n  Purity summary:\n")
print(purity_dt)
cat("  wrote:", purity_out, "\n")

# Pre-registered acceptance check -----------------------------------------

acceptance_threshold_purity_pp <- 5.0
acceptance_threshold_hop1_pp <- 1.0
purity_dt[, accept_purity := delta_purity_pp >= acceptance_threshold_purity_pp]
purity_dt[, accept_hop1 := delta_hop1_pp >= acceptance_threshold_hop1_pp]
cat("\n  Acceptance flags (>=", acceptance_threshold_purity_pp, "pp purity AND >=", acceptance_threshold_hop1_pp, "pp hop1):\n")
print(purity_dt[label != "original", .(trait, threshold, n_seed, delta_purity_pp, delta_hop1_pp, accept_purity, accept_hop1)])

# Plot --------------------------------------------------------------------

cat("\n[3/4] Plot threshold sweep\n")
plot_dt <- purity_dt[label != "original"]
plot_dt[, metric_purity := target_purity * 100]
plot_dt[, metric_hop1 := hop1_target_fraction * 100]
plot_long <- melt(
  plot_dt,
  id.vars = c("trait", "threshold", "n_seed"),
  measure.vars = c("metric_purity", "metric_hop1"),
  variable.name = "metric",
  value.name = "value_pct"
)
plot_long[, metric := factor(
  metric,
  levels = c("metric_purity", "metric_hop1"),
  labels = c("Seed target purity (%)", "Hop1 target fraction (%)")
)]

# overlay original baselines as horizontal lines
baseline_dt <- purity_dt[label == "original", .(trait, baseline_purity = target_purity * 100, baseline_hop1 = hop1_target_fraction * 100)]
baseline_long <- melt(
  baseline_dt,
  id.vars = "trait",
  measure.vars = c("baseline_purity", "baseline_hop1"),
  variable.name = "metric",
  value.name = "baseline_pct"
)
baseline_long[, metric := factor(
  metric,
  levels = c("baseline_purity", "baseline_hop1"),
  labels = c("Seed target purity (%)", "Hop1 target fraction (%)")
)]

p <- ggplot(plot_long, aes(x = threshold, y = value_pct, color = trait)) +
  geom_line() +
  geom_point() +
  geom_hline(data = baseline_long, aes(yintercept = baseline_pct, color = trait), linetype = "dashed") +
  facet_wrap(~ metric, scales = "free_y") +
  labs(
    title = "Seed RNA-target-fraction filter sweep",
    subtitle = "Solid = filtered seeds; dashed = original (unfiltered) seeds",
    x = "RNA target fraction threshold",
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "top")

fig_out <- file.path(paths$fig_dir, "Fig_seed_purity_threshold_sweep.png")
ggsave(fig_out, p, width = 10, height = 5, dpi = 300)
cat("  wrote:", fig_out, "\n")

# Manifest ----------------------------------------------------------------

cat("\n[4/4] Write manifest\n")
manifest <- list(
  script = "03_seed_purity_diagnostics.R",
  traits = traits,
  thresholds = thresholds,
  acceptance_thresholds = list(
    purity_delta_pp = acceptance_threshold_purity_pp,
    hop1_delta_pp = acceptance_threshold_hop1_pp
  ),
  inputs = list(
    seed_with_rna_fraction = seed_with_rna_path,
    seed_signal_by_cell = paths$seed_signal_by_cell,
    original_knn = paths$original_knn
  ),
  outputs = list(
    purity_summary = purity_out,
    figure = fig_out
  ),
  notes = c(
    "Hop1 is computed on the original LSI mutualKNN graph used by the published SCAVENGE pipeline.",
    "Acceptance: choose the smallest threshold whose delta_purity_pp >= 5 AND delta_hop1_pp >= 1 for at least 1 trait.",
    "If no threshold passes for any trait, direction 1 is negative for this dataset; do not run the full TRS."
  ),
  created_at = as.character(Sys.time())
)
manifest_out <- file.path(paths$out_dir, "seed_purity_diagnostics_manifest.json")
write_json(manifest, manifest_out, pretty = TRUE, auto_unbox = TRUE)

cat("\n========== Done ==========\n")
cat("  manifest:", manifest_out, "\n")
