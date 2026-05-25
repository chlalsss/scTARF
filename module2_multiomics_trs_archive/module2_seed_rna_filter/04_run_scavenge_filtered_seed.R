# 04: Run full SCAVENGE TRS using RNA-filtered seeds + original LSI mutualKNN.
#
# Difference from the original SCAVENGE pipeline:
#   - Seed cells = filtered_seed_cells_<trait>_t<NN>.csv (from 02), instead of
#     z_score top 5%.
#   - KNN graph, gamma, scale_factor recipe are unchanged.
#
# Inputs:
#   filtered_seed_cells_<trait>_t<NN>.csv  # from 02
#   seed_signal_by_cell.csv.gz             # has z_score per cell per trait
#   mutualknn_k30.rds                      # original LSI mutualKNN
#   TRS_matrix_small_all_traits_fullcells.rds (optional, for original baseline)
#
# Outputs (under <out_dir>/trs_filtered_t<NN>/):
#   SCAVENGE_filtered_<trait>.tsv.gz       # per-cell z_score, seed flag, TRS
#   trs_top5_enrichment_summary.csv        # filtered vs original top 5% target fraction
#   isolated_cells_removed_before_random_walk.csv
#   trait_summary.csv
#   run_manifest.json
#
# Parameters:
#   MODULE2_TRAITS                : traits (default "MONO,PLT,LYM")
#   MODULE2_SEED_FILTER_PRIMARY   : threshold to use (default 0.50)
#   MODULE2_RW_GAMMA              : RWR restart probability (default 0.05)
#   MODULE2_TOP_PERCENT           : top-x percent for target enrichment (default 0.05)

source("D:/RStudio/work/module2_seed_rna_filter/_bootstrap.R")

required_pkgs(c(
  "SCAVENGE", "data.table", "Matrix", "jsonlite"
))
suppressPackageStartupMessages({
  library(SCAVENGE)
  library(data.table)
  library(Matrix)
  library(jsonlite)
})

paths <- seed_filter_paths()
ensure_seed_filter_dirs(paths)

traits <- parse_csv_env("MODULE2_TRAITS", "MONO,PLT,LYM")
threshold <- env_numeric("MODULE2_SEED_FILTER_PRIMARY", 0.50)
gamma <- env_numeric("MODULE2_RW_GAMMA", 0.05)
top_pct <- env_numeric("MODULE2_TOP_PERCENT", 0.05)
cap_q <- env_numeric("MODULE2_CAP_QUANTILE", 0.95)
scale_pct <- env_numeric("MODULE2_SCALE_PERCENT", 0.01)

th_label <- gsub("\\.", "p", sprintf("%.2f", threshold))
out_dir <- file.path(paths$out_dir, sprintf("trs_filtered_t%s", th_label))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("========== Module 2 seed RNA filter: 04 run SCAVENGE on filtered seed ==========\n")
cat("  traits     :", paste(traits, collapse = ", "), "\n")
cat("  threshold  :", threshold, "(label", th_label, ")\n")
cat("  gamma      :", gamma, "\n")
cat("  cap quant  :", cap_q, "\n")
cat("  scale pct  :", scale_pct, "\n")
cat("  top pct    :", top_pct, "\n")
cat("  out_dir    :", out_dir, "\n")

# Load --------------------------------------------------------------------

cat("\n[1/5] Load inputs\n")
if (!file.exists(paths$seed_signal_by_cell)) {
  stop("Missing seed_signal_by_cell: ", paths$seed_signal_by_cell, call. = FALSE)
}
if (!file.exists(paths$original_knn)) {
  stop("Missing original LSI KNN: ", paths$original_knn, call. = FALSE)
}

seed_signal_dt <- fread(paths$seed_signal_by_cell)
seed_signal_dt[, target_type := as.logical(target_type)]
seed_signal_dt[, seed := as.logical(seed)]
original_knn <- readRDS(paths$original_knn)
original_knn <- as_binary_graph(original_knn)
cells_in_graph <- rownames(original_knn)
cat("  graph cells:", length(cells_in_graph), "\n")

# Optional: original SCAVENGE TRS for baseline comparison
original_trs_path <- file.path(paths$original_dir, "TRS_matrix_small_all_traits_fullcells.rds")
original_trs <- if (file.exists(original_trs_path)) readRDS(original_trs_path) else NULL
if (is.null(original_trs)) {
  warning("Missing original TRS at ", original_trs_path, "; original baseline top5 will be recomputed from seed_signal.")
}

# Remove isolated cells once
graph_degree <- Matrix::rowSums(original_knn)
graph_cells <- cells_in_graph[graph_degree > 0]
isolated_cells <- cells_in_graph[graph_degree == 0]
if (length(isolated_cells) > 0) {
  cat("  isolated cells in original graph:", length(isolated_cells), "\n")
  fwrite(
    data.table(cell = isolated_cells),
    file.path(out_dir, "isolated_cells_removed_before_random_walk.csv")
  )
}
graph <- original_knn[graph_cells, graph_cells]

# Run ---------------------------------------------------------------------

cat("\n[2/5] Run RWR for each trait\n")
trait_summary_list <- list()
top_summary_list <- list()
for (tr in traits) {
  cat("\n[Trait]", tr, "\n")
  signal_tr <- seed_signal_dt[trait == tr]
  if (nrow(signal_tr) == 0) {
    warning("seed_signal_by_cell has no rows for trait ", tr); next
  }
  z_score <- setNames(as.numeric(signal_tr$z_score), signal_tr$cell)
  target_vec_named <- setNames(as.logical(signal_tr$target_type), signal_tr$cell)

  # filtered seed
  filtered_path <- file.path(
    paths$out_dir,
    sprintf("filtered_seed_cells_%s_t%s.csv", tr, th_label)
  )
  if (!file.exists(filtered_path)) {
    warning("Missing filtered seed file: ", filtered_path, "; skipping ", tr)
    next
  }
  filtered_dt <- fread(filtered_path)
  seed_total <- as.character(filtered_dt$cell)
  seed_in_graph <- intersect(seed_total, graph_cells)
  seed_dropped_isolated <- setdiff(seed_total, seed_in_graph)
  if (length(seed_in_graph) == 0) {
    warning("No filtered seed cells in graph for trait ", tr); next
  }

  # Use the same scale_factor recipe as the original pipeline:
  #   top scale_pct of z_score (e.g. top 1%) mean.
  z_for_scale <- z_score[!is.na(z_score)]
  cap_n <- max(1L, ceiling(length(z_for_scale) * scale_pct))
  z_top <- sort(z_for_scale, decreasing = TRUE)[seq_len(cap_n)]
  scale_factor <- mean(z_top, na.rm = TRUE)

  cat("  filtered seed total:", length(seed_total),
      "  in graph:", length(seed_in_graph),
      "  isolated dropped:", length(seed_dropped_isolated), "\n")
  cat("  scale_factor (mean top", scale_pct, "z):", round(scale_factor, 4), "\n")

  np_score <- randomWalk_sparse(
    intM = graph,
    queryCells = seed_in_graph,
    gamma = gamma
  )

  omit_idx <- np_score == 0
  np_use <- np_score[!omit_idx]
  keep_cells <- names(np_use)
  TRS <- capOutlierQuantile(x = np_use, q_ceiling = cap_q) |> max_min_scale()
  TRS <- TRS * scale_factor

  out_df <- data.table(
    cell = keep_cells,
    z_score = as.numeric(z_score[keep_cells]),
    seed = keep_cells %in% seed_in_graph,
    target_type = target_vec_named[keep_cells],
    np_score = as.numeric(np_use),
    TRS = as.numeric(TRS)
  )
  out_file <- file.path(out_dir, sprintf("SCAVENGE_filtered_%s.tsv.gz", tr))
  fwrite(out_df, out_file, sep = "\t")
  cat("  saved:", out_file, "(", nrow(out_df), "rows)\n")

  # top-x% target enrichment, filtered
  n_top <- max(1L, ceiling(nrow(out_df) * top_pct))
  setorder(out_df, -TRS)
  top_filtered <- out_df[seq_len(n_top)]
  top_filtered_target_fraction <- mean(top_filtered$target_type, na.rm = TRUE)

  # original baseline: prefer original_trs if available, else recompute from
  # seed_signal (z_score top 5% as seed).
  if (!is.null(original_trs) && tr %in% colnames(original_trs)) {
    orig_dt <- data.table(
      cell = as.character(original_trs[["cell"]]),
      TRS = as.numeric(original_trs[[tr]])
    )
    orig_dt <- merge(
      orig_dt,
      data.table(cell = signal_tr$cell, target_type = target_vec_named[signal_tr$cell]),
      by = "cell"
    )
    orig_dt <- orig_dt[!is.na(TRS)]
    setorder(orig_dt, -TRS)
    n_top_orig <- max(1L, ceiling(nrow(orig_dt) * top_pct))
    top_orig <- orig_dt[seq_len(n_top_orig)]
    top_original_target_fraction <- mean(top_orig$target_type, na.rm = TRUE)
    n_orig_cells <- nrow(orig_dt)
    baseline_source <- "TRS_matrix_small_all_traits_fullcells.rds"
  } else {
    # Recompute baseline TRS using z-score top 5% as seed on the same graph.
    orig_seed <- signal_tr[seed == TRUE, cell]
    orig_seed_in_graph <- intersect(orig_seed, graph_cells)
    if (length(orig_seed_in_graph) == 0) {
      warning("Cannot recompute baseline: no original seed in graph for ", tr)
      top_original_target_fraction <- NA_real_
      n_orig_cells <- 0
      baseline_source <- "missing"
    } else {
      np_orig <- randomWalk_sparse(
        intM = graph,
        queryCells = orig_seed_in_graph,
        gamma = gamma
      )
      omit_orig <- np_orig == 0
      np_orig_use <- np_orig[!omit_orig]
      keep_orig <- names(np_orig_use)
      TRS_orig <- capOutlierQuantile(x = np_orig_use, q_ceiling = cap_q) |> max_min_scale()
      TRS_orig <- TRS_orig * scale_factor
      orig_dt_recomp <- data.table(
        cell = keep_orig,
        TRS = as.numeric(TRS_orig),
        target_type = target_vec_named[keep_orig]
      )
      setorder(orig_dt_recomp, -TRS)
      n_top_orig <- max(1L, ceiling(nrow(orig_dt_recomp) * top_pct))
      top_orig <- orig_dt_recomp[seq_len(n_top_orig)]
      top_original_target_fraction <- mean(top_orig$target_type, na.rm = TRUE)
      n_orig_cells <- nrow(orig_dt_recomp)
      baseline_source <- "recomputed_z_top_seed_same_graph"
    }
  }

  background_target <- mean(target_vec_named, na.rm = TRUE)

  top_summary_list[[tr]] <- data.table(
    trait = tr,
    threshold = threshold,
    n_filtered_seed_total = length(seed_total),
    n_filtered_seed_in_graph = length(seed_in_graph),
    n_filtered_cells_in_TRS = nrow(out_df),
    n_top_filtered = n_top,
    top_filtered_target_fraction = top_filtered_target_fraction,
    n_orig_cells_in_TRS = n_orig_cells,
    top_original_target_fraction = top_original_target_fraction,
    delta_top_pp = if (is.finite(top_original_target_fraction))
      (top_filtered_target_fraction - top_original_target_fraction) * 100 else NA_real_,
    background_target = background_target,
    baseline_source = baseline_source,
    target_definition = target_definition(tr)
  )

  trait_summary_list[[tr]] <- data.table(
    trait = tr,
    n_filtered_seed_total = length(seed_total),
    n_filtered_seed_in_graph = length(seed_in_graph),
    n_seed_isolated_dropped = length(seed_dropped_isolated),
    n_graph_cells = length(graph_cells),
    n_omit = sum(omit_idx),
    pct_omit = mean(omit_idx),
    scale_factor = scale_factor,
    mean_np_score = mean(np_use),
    mean_TRS = mean(TRS),
    median_TRS = median(TRS),
    output_file = out_file
  )
  cat("  top", top_pct, "filtered target frac:",
      round(top_filtered_target_fraction * 100, 2), "%",
      " original:", round(top_original_target_fraction * 100, 2), "%",
      " delta_pp:", round((top_filtered_target_fraction - top_original_target_fraction) * 100, 2),
      "(baseline:", baseline_source, ")\n")
}

# Summary tables ----------------------------------------------------------

cat("\n[3/5] Write summary tables\n")
trait_summary_dt <- rbindlist(trait_summary_list, use.names = TRUE, fill = TRUE)
top_summary_dt <- rbindlist(top_summary_list, use.names = TRUE, fill = TRUE)
trait_summary_path <- file.path(out_dir, "trait_summary.csv")
top_summary_path <- file.path(out_dir, "trs_top5_enrichment_summary.csv")
fwrite(trait_summary_dt, trait_summary_path)
fwrite(top_summary_dt, top_summary_path)
cat("\n  Top-", top_pct, " enrichment summary:\n", sep = "")
print(top_summary_dt[, .(trait, threshold, n_filtered_seed_in_graph,
                         top_filtered_target_fraction, top_original_target_fraction,
                         delta_top_pp, baseline_source)])

# Acceptance flags --------------------------------------------------------

cat("\n[4/5] Acceptance check (delta_top_pp >= 1 pp on at least 1 trait)\n")
top_summary_dt[, accept := delta_top_pp >= 1]
print(top_summary_dt[, .(trait, threshold, top_original_target_fraction,
                         top_filtered_target_fraction, delta_top_pp, accept)])

# Manifest ----------------------------------------------------------------

cat("\n[5/5] Write manifest\n")
manifest <- list(
  script = "04_run_scavenge_filtered_seed.R",
  threshold = threshold,
  traits = traits,
  parameters = list(
    gamma = gamma,
    cap_quantile = cap_q,
    scale_percent = scale_pct,
    top_percent = top_pct
  ),
  inputs = list(
    seed_signal_by_cell = paths$seed_signal_by_cell,
    original_knn = paths$original_knn,
    original_trs = if (file.exists(original_trs_path)) original_trs_path else NA
  ),
  outputs = list(
    out_dir = out_dir,
    trait_summary = trait_summary_path,
    top_enrichment_summary = top_summary_path
  ),
  notes = c(
    "RWR is run on the original LSI mutualKNN graph; only the seed set differs from baseline.",
    "scale_factor follows the original recipe: mean of top-1% z_score values across all cells.",
    "Cap quantile = 0.95 (matches module2_wnn_lite_fusion/03 and the published blood_small workflow).",
    "Acceptance: delta_top_pp >= 1 pp on at least 1 trait."
  ),
  created_at = as.character(Sys.time())
)
manifest_out <- file.path(out_dir, "run_manifest.json")
write_json(manifest, manifest_out, pretty = TRUE, auto_unbox = TRUE)

cat("\n========== Done ==========\n")
cat("  manifest:", manifest_out, "\n")
