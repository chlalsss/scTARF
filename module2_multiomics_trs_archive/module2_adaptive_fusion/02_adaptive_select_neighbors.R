# Convert top-k distances to per-ATAC relative similarity and select neighbors.
#
# Optional cell-type consistency filter (MODULE2_TYPE_CONSISTENCY, default "true"):
#   after similarity cutoff and rank-1 fallback, drop selected RNA neighbors whose
#   cell_type differs from the ATAC cell's cell_type. ATAC cells with zero same-type
#   neighbors keep their nearest RNA neighbor regardless of type so RNA_ref still
#   has at least one source.

source("D:/RStudio/work/module2_adaptive_fusion/_bootstrap.R")

required_pkgs(c("data.table", "jsonlite"))
suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

paths <- module2_paths()
ensure_module2_dirs(paths)

k <- as.integer(Sys.getenv("MODULE2_TOPK", "1"))
cutoffs <- get_cutoffs()
type_consistency <- tolower(Sys.getenv("MODULE2_TYPE_CONSISTENCY", "true")) %in% c("true", "1", "yes", "on")
topk_csv <- file.path(paths$out_dir, sprintf("atac_top%d_rna_neighbors.csv", k))
if (!file.exists(topk_csv)) {
  stop("Missing top-k neighbors. Run 01_build_top100_neighbors.R first: ", topk_csv)
}

cat("========== Module 2: adaptive RNA neighbor selection ==========\n")
cat("  input             :", topk_csv, "\n")
cat("  cutoffs           :", paste(cutoffs, collapse = ", "), "\n")
cat("  type consistency  :", type_consistency, "\n")

dt <- data.table::fread(topk_csv)
needed <- c("atac_cell", "rank", "rna_cell", "distance")
if (!all(needed %in% colnames(dt))) {
  stop("Neighbor table must contain: ", paste(needed, collapse = ", "))
}
if (max(dt$rank) == 1L) {
  out_label <- "top1"
  selected <- dt[, .(atac_cell, rank = 1L, rna_cell, distance,
                     similarity = 1.0, weight = 1.0)]
  sel_out <- file.path(paths$out_dir,
                       sprintf("adaptive_neighbors_cutoff_%s.csv", out_label))
  by_cell_out <- file.path(paths$out_dir,
                           sprintf("adaptive_selection_by_cell_cutoff_%s.csv", out_label))
  fwrite(selected, sel_out)
  fwrite(data.table(atac_cell = selected$atac_cell, n_selected = 1L,
                    mean_similarity = 1.0, scale = NA_real_),
         by_cell_out)
  jsonlite::write_json(
    list(script = "02_adaptive_select_neighbors.R", mode = "top1",
         input = topk_csv, output_csv = sel_out,
         created_at = as.character(Sys.time())),
    file.path(paths$out_dir, "adaptive_selection_manifest.json"),
    pretty = TRUE, auto_unbox = TRUE
  )
  cat("top1 mode: skipping cutoff / type filter\n")
  cat("========== Done ==========\n")
} else {
setorder(dt, atac_cell, rank)

dt[, d_min := min(distance), by = atac_cell]
dt[, delta := distance - d_min]
scale_dt <- dt[, .(scale = safe_scale_delta(delta)), by = atac_cell]
dt <- merge(dt, scale_dt, by = "atac_cell", all.x = TRUE, sort = FALSE)
dt[, similarity := exp(-delta / scale)]

rna_lab <- NULL
atac_lab <- NULL
if (file.exists(paths$rna_cells) && file.exists(paths$atac_cells)) {
  rna_lab <- fread(paths$rna_cells)
  atac_lab <- fread(paths$atac_cells)
  if (all(c("cell", "cell_type") %in% colnames(rna_lab)) &&
      all(c("cell", "cell_type") %in% colnames(atac_lab))) {
    setnames(rna_lab, c("cell", "cell_type"), c("rna_cell", "rna_type"))
    setnames(atac_lab, c("cell", "cell_type"), c("atac_cell", "atac_type"))
  } else {
    rna_lab <- NULL
    atac_lab <- NULL
  }
}

if (type_consistency && (is.null(rna_lab) || is.null(atac_lab))) {
  warning("Type consistency requested but cell-type labels not available; filter disabled.")
  type_consistency <- FALSE
}

selected_cols <- c("atac_cell", "rank", "rna_cell", "distance", "d_min", "delta", "scale", "similarity")

apply_type_filter <- function(selected, dt, rna_lab, atac_lab) {
  before_n <- nrow(selected)
  all_atac <- unique(selected$atac_cell)

  sel_lab <- merge(selected, rna_lab, by = "rna_cell", all.x = TRUE, sort = FALSE)
  sel_lab <- merge(sel_lab, atac_lab, by = "atac_cell", all.x = TRUE, sort = FALSE)
  sel_lab[, type_match := !is.na(rna_type) & !is.na(atac_type) & (rna_type == atac_type)]

  same_type <- sel_lab[type_match == TRUE, ..selected_cols]
  missing_typed <- setdiff(all_atac, unique(same_type$atac_cell))
  if (length(missing_typed) > 0) {
    fallback <- dt[atac_cell %in% missing_typed & rank == 1, ..selected_cols]
    same_type <- rbind(same_type, fallback, use.names = TRUE)
  }
  setorder(same_type, atac_cell, rank)

  list(
    selected = same_type,
    n_before = before_n,
    n_after = nrow(same_type),
    n_cells_no_same_type = length(missing_typed)
  )
}

all_summary <- list()
for (cutoff in cutoffs) {
  label <- cutoff_label(cutoff)
  cat("  cutoff", cutoff, "\n")

  selected <- dt[similarity >= cutoff]

  # Keep at least the nearest RNA cell for every ATAC cell.
  missing_cells <- setdiff(unique(dt$atac_cell), unique(selected$atac_cell))
  if (length(missing_cells) > 0) {
    selected <- rbind(
      selected,
      dt[atac_cell %in% missing_cells & rank == 1],
      use.names = TRUE
    )
  }
  setorder(selected, atac_cell, rank)

  filter_info <- list(
    enabled = type_consistency,
    n_before_filter = nrow(selected),
    n_after_filter = NA_integer_,
    n_cells_no_same_type = NA_integer_,
    n_cross_type_dropped = NA_integer_
  )
  if (type_consistency) {
    res <- apply_type_filter(selected, dt, rna_lab, atac_lab)
    selected <- res$selected
    filter_info$n_after_filter <- res$n_after
    filter_info$n_cells_no_same_type <- res$n_cells_no_same_type
    filter_info$n_cross_type_dropped <- res$n_before - res$n_after + res$n_cells_no_same_type
    cat(sprintf(
      "    type filter: %d -> %d rows; cross-type dropped %d; %d cells used nearest fallback\n",
      res$n_before, res$n_after,
      filter_info$n_cross_type_dropped,
      res$n_cells_no_same_type
    ))
  }

  selected[, weight := similarity / sum(similarity), by = atac_cell]

  by_cell <- selected[, .(
    n_selected = .N,
    max_selected_rank = max(rank),
    mean_similarity = mean(similarity),
    min_similarity = min(similarity),
    mean_distance = mean(distance),
    best_distance = min(distance),
    scale = scale[1]
  ), by = atac_cell]

  if (!is.null(rna_lab) && !is.null(atac_lab)) {
    selected_labeled <- merge(selected, rna_lab, by = "rna_cell", all.x = TRUE)
    selected_labeled <- merge(selected_labeled, atac_lab, by = "atac_cell", all.x = TRUE)
    label_stats <- selected_labeled[, .(
      selected_same_type_frac = mean(rna_type == atac_type, na.rm = TRUE),
      selected_weighted_same_type_frac = sum(weight * as.numeric(rna_type == atac_type), na.rm = TRUE)
    ), by = atac_cell]
    by_cell <- merge(by_cell, label_stats, by = "atac_cell", all.x = TRUE)
  }

  selected_out <- file.path(paths$out_dir, sprintf("adaptive_neighbors_cutoff_%s.csv", label))
  by_cell_out <- file.path(paths$out_dir, sprintf("adaptive_selection_by_cell_cutoff_%s.csv", label))
  fwrite(selected[, .(atac_cell, rank, rna_cell, distance, similarity, weight)], selected_out)
  fwrite(by_cell, by_cell_out)

  summary_row <- data.table(
    cutoff = cutoff,
    type_consistency = type_consistency,
    n_atac = nrow(by_cell),
    n_before_type_filter = filter_info$n_before_filter,
    n_after_type_filter = filter_info$n_after_filter,
    n_cross_type_dropped = filter_info$n_cross_type_dropped,
    n_cells_no_same_type = filter_info$n_cells_no_same_type,
    mean_n_selected = mean(by_cell$n_selected),
    median_n_selected = median(by_cell$n_selected),
    q10_n_selected = as.numeric(quantile(by_cell$n_selected, 0.10)),
    q90_n_selected = as.numeric(quantile(by_cell$n_selected, 0.90)),
    pct_single_neighbor = mean(by_cell$n_selected <= 1),
    pct_ge_80_neighbors = mean(by_cell$n_selected >= 80),
    mean_selected_same_type_frac = if ("selected_same_type_frac" %in% colnames(by_cell)) {
      mean(by_cell$selected_same_type_frac, na.rm = TRUE)
    } else NA_real_,
    mean_selected_weighted_same_type_frac = if ("selected_weighted_same_type_frac" %in% colnames(by_cell)) {
      mean(by_cell$selected_weighted_same_type_frac, na.rm = TRUE)
    } else NA_real_,
    selected_file = selected_out,
    by_cell_file = by_cell_out
  )
  all_summary[[label]] <- summary_row
}

summary_dt <- rbindlist(all_summary, use.names = TRUE)
summary_out <- file.path(paths$out_dir, "adaptive_selection_summary.csv")
fwrite(summary_dt, summary_out)
cat("  saved summary:", summary_out, "\n")

jsonlite::write_json(
  list(
    script = "02_adaptive_select_neighbors.R",
    input = topk_csv,
    cutoffs = cutoffs,
    type_consistency = type_consistency,
    similarity = "exp(-(distance - nearest_distance) / per_cell_median_delta)",
    type_filter = if (type_consistency) {
      "drop selected RNA neighbors whose cell_type != ATAC cell_type; cells with zero same-type neighbors fall back to rank-1 nearest"
    } else {
      "disabled"
    },
    output_summary = summary_out,
    created_at = as.character(Sys.time())
  ),
  file.path(paths$out_dir, "adaptive_selection_manifest.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

cat("========== Done ==========\n")
}
