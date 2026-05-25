# Common builder for WNN-lite gated mutualKNN graphs.
#
# Entry scripts must define:
#   gated_version
#   gated_label
#   gated_weights_file

if (!exists("wnn_gated_code_dir", inherits = TRUE)) {
  source("D:/RStudio/work/module2_wnn_lite_gated/_bootstrap.R")
}

required_pkgs(c("data.table", "Matrix", "RANN", "ggplot2", "jsonlite"))
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(RANN)
  library(ggplot2)
  library(jsonlite)
})

if (!exists("gated_version", inherits = TRUE)) {
  stop("Entry script must define gated_version", call. = FALSE)
}
if (!exists("gated_label", inherits = TRUE)) {
  stop("Entry script must define gated_label", call. = FALSE)
}
if (!exists("gated_weights_file", inherits = TRUE)) {
  stop("Entry script must define gated_weights_file", call. = FALSE)
}

paths <- gated_paths(gated_version)
ensure_gated_dirs(paths)

cutoff <- get_single_cutoff()
k <- env_integer("WNN_LITE_KNN_K", env_integer("MODULE2_KNN_K", 30))
candidate_k <- env_integer("WNN_LITE_CANDIDATE_K", k)
l2_norm <- env_flag("WNN_LITE_L2_NORM", FALSE)

if (k < 1) stop("WNN_LITE_KNN_K must be positive", call. = FALSE)
if (candidate_k < k) {
  stop("WNN_LITE_CANDIDATE_K must be >= WNN_LITE_KNN_K", call. = FALSE)
}

rna_ref_file <- rna_reference_path(paths, cutoff)
if (!file.exists(gated_weights_file)) {
  stop("Missing gated weights: ", gated_weights_file, call. = FALSE)
}
if (!file.exists(rna_ref_file)) {
  stop("Missing RNA reference latent: ", rna_ref_file, call. = FALSE)
}

cat("========== WNN-lite gated: build mutualKNN ==========\n")
cat("  version     :", gated_version, "\n")
cat("  k           :", k, "\n")
cat("  candidate k :", candidate_k, "\n")
cat("  l2_norm     :", l2_norm, "\n")
cat("  weights     :", gated_weights_file, "\n")
cat("  RNA ref     :", rna_ref_file, "\n")
cat("  output dir  :", paths$out_dir, "\n")

loaded <- load_shared_embedding(paths$shared_embedding)
emb <- as.data.table(loaded$emb)
latent_cols <- loaded$latent_cols
atac_dt <- emb[modality == "atac"]
if (nrow(atac_dt) == 0) stop("No ATAC rows in shared embedding", call. = FALSE)
if (anyDuplicated(atac_dt$cell)) stop("Duplicated ATAC cell barcodes", call. = FALSE)

rna_ref <- as.data.table(readRDS(rna_ref_file))
rna_ref_cols <- grep("^rna_ref_", colnames(rna_ref), value = TRUE)
if (!"atac_cell" %in% colnames(rna_ref) || length(rna_ref_cols) == 0) {
  stop("RNA reference latent must contain atac_cell and rna_ref_* columns", call. = FALSE)
}
if (anyDuplicated(rna_ref$atac_cell)) {
  stop("Duplicated atac_cell values in RNA reference latent", call. = FALSE)
}

weights <- fread(gated_weights_file)
needed_weight_cols <- c("atac_cell", "rna_weight", "atac_weight")
if (!all(needed_weight_cols %in% colnames(weights))) {
  stop("Weight table must contain: ", paste(needed_weight_cols, collapse = ", "), call. = FALSE)
}
if (anyDuplicated(weights$atac_cell)) {
  stop("Duplicated atac_cell values in weight table", call. = FALSE)
}
if (!"gated" %in% names(weights)) {
  weights[, gated := FALSE]
}
weights[, gated := as.logical(gated)]

common <- Reduce(intersect, list(atac_dt$cell, rna_ref$atac_cell, weights$atac_cell))
if (length(common) == 0) {
  stop("No common ATAC cells across shared embedding, RNA reference, and weights", call. = FALSE)
}
if (length(common) < nrow(atac_dt)) {
  warning("Using common ATAC cells only: ", length(common), " / ", nrow(atac_dt))
}

atac_dt <- atac_dt[match(common, cell)]
rna_ref <- rna_ref[match(common, atac_cell)]
weights <- weights[match(common, atac_cell)]
if (!all(atac_dt$cell == rna_ref$atac_cell) || !all(atac_dt$cell == weights$atac_cell)) {
  stop("Internal cell-order alignment failed", call. = FALSE)
}
if (length(latent_cols) != length(rna_ref_cols)) {
  warning("ATAC latent dims (", length(latent_cols), ") != RNA ref dims (", length(rna_ref_cols), ")")
}

gated_weight_check <- weights[gated == TRUE, .(
  n_gated = .N,
  min_rna_weight = min(rna_weight),
  max_rna_weight = max(rna_weight),
  min_atac_weight = min(atac_weight),
  max_atac_weight = max(atac_weight)
)]
if (nrow(gated_weight_check) == 1 && gated_weight_check$n_gated > 0) {
  if (gated_weight_check$max_rna_weight != 0 || gated_weight_check$min_atac_weight != 1) {
    stop("Gated cells are not forced to rna_weight=0 and atac_weight=1", call. = FALSE)
  }
}

atac_mat <- as.matrix(atac_dt[, ..latent_cols])
rna_mat <- as.matrix(rna_ref[, ..rna_ref_cols])
storage.mode(atac_mat) <- "double"
storage.mode(rna_mat) <- "double"
atac_mat[!is.finite(atac_mat)] <- 0
rna_mat[!is.finite(rna_mat)] <- 0

if (l2_norm) {
  atac_mat <- l2_normalize_rows(atac_mat)
  rna_mat <- l2_normalize_rows(rna_mat)
}
rownames(atac_mat) <- common
rownames(rna_mat) <- common

n <- nrow(atac_mat)
cat("  cells       :", n, "\n")
cat("  gated cells :", sum(weights$gated), "\n")
cat("  ATAC dims   :", ncol(atac_mat), "\n")
cat("  RNA ref dims:", ncol(rna_mat), "\n")

t0 <- Sys.time()
nn_atac_raw <- RANN::nn2(data = atac_mat, query = atac_mat, k = candidate_k + 1)
nn_rna_raw <- RANN::nn2(data = rna_mat, query = rna_mat, k = candidate_k + 1)
nn_elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("  RANN elapsed: %.1f sec\n", nn_elapsed))

nn_atac <- strip_self_neighbors(nn_atac_raw$nn.idx, nn_atac_raw$nn.dists, candidate_k)
nn_rna <- strip_self_neighbors(nn_rna_raw$nn.idx, nn_rna_raw$nn.dists, candidate_k)

atac_norm2 <- rowSums(atac_mat^2)
rna_norm2 <- rowSums(rna_mat^2)

out_to_idx <- matrix(NA_integer_, nrow = n, ncol = k)
out_rank <- matrix(rep(seq_len(k), times = n), nrow = n, ncol = k, byrow = TRUE)
out_score <- matrix(NA_real_, nrow = n, ncol = k)
out_atac_sim <- matrix(NA_real_, nrow = n, ncol = k)
out_rna_sim <- matrix(NA_real_, nrow = n, ncol = k)
out_atac_component <- matrix(NA_real_, nrow = n, ncol = k)
out_rna_component <- matrix(NA_real_, nrow = n, ncol = k)
out_from_atac <- matrix(FALSE, nrow = n, ncol = k)
out_from_rna <- matrix(FALSE, nrow = n, ncol = k)

view_overlap <- numeric(n)
candidate_count <- integer(n)
selected_atac_fraction <- numeric(n)
selected_rna_fraction <- numeric(n)
selected_both_fraction <- numeric(n)
selected_atac_component_fraction <- numeric(n)

cat("  scoring gated WNN-lite candidates...\n")
pb <- utils::txtProgressBar(min = 0, max = n, style = 3)
for (i in seq_len(n)) {
  atac_candidates <- nn_atac$idx[i, ]
  rna_candidates <- nn_rna$idx[i, ]
  candidates <- sort(unique(c(atac_candidates, rna_candidates)))
  candidates <- candidates[!is.na(candidates) & candidates != i]
  if (length(candidates) < k) {
    stop("Only ", length(candidates), " candidates for row ", i, "; need ", k, call. = FALSE)
  }

  from_atac <- candidates %in% atac_candidates
  from_rna <- candidates %in% rna_candidates
  atac_dist <- candidate_distances(atac_mat, atac_norm2, i, candidates)
  rna_dist <- candidate_distances(rna_mat, rna_norm2, i, candidates)
  atac_sim <- distance_to_similarity(atac_dist)
  rna_sim <- distance_to_similarity(rna_dist)
  atac_component <- weights$atac_weight[i] * atac_sim
  rna_component <- weights$rna_weight[i] * rna_sim
  score <- atac_component + rna_component

  ord <- order(-score, atac_dist, rna_dist, candidates)
  keep <- ord[seq_len(k)]

  out_to_idx[i, ] <- candidates[keep]
  out_score[i, ] <- score[keep]
  out_atac_sim[i, ] <- atac_sim[keep]
  out_rna_sim[i, ] <- rna_sim[keep]
  out_atac_component[i, ] <- atac_component[keep]
  out_rna_component[i, ] <- rna_component[keep]
  out_from_atac[i, ] <- from_atac[keep]
  out_from_rna[i, ] <- from_rna[keep]

  view_overlap[i] <- length(intersect(atac_candidates, rna_candidates)) /
    length(union(atac_candidates, rna_candidates))
  candidate_count[i] <- length(candidates)
  selected_atac_fraction[i] <- mean(from_atac[keep])
  selected_rna_fraction[i] <- mean(from_rna[keep])
  selected_both_fraction[i] <- mean(from_atac[keep] & from_rna[keep])
  selected_atac_component_fraction[i] <- mean(
    atac_component[keep] / pmax(score[keep], 1e-12)
  )

  if (i %% 200 == 0 || i == n) utils::setTxtProgressBar(pb, i)
}
close(pb)
score_elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs")) - nn_elapsed
cat(sprintf("  scoring elapsed: %.1f sec\n", score_elapsed))

edge_dt <- data.table(
  from_idx = rep(seq_len(n), each = k),
  to_idx = as.vector(t(out_to_idx)),
  rank = as.vector(t(out_rank)),
  wnn_score = as.vector(t(out_score)),
  atac_similarity = as.vector(t(out_atac_sim)),
  rna_similarity = as.vector(t(out_rna_sim)),
  atac_component = as.vector(t(out_atac_component)),
  rna_component = as.vector(t(out_rna_component)),
  candidate_from_atac = as.vector(t(out_from_atac)),
  candidate_from_rna = as.vector(t(out_from_rna))
)
edge_dt[, from_cell := common[from_idx]]
edge_dt[, to_cell := common[to_idx]]
setcolorder(edge_dt, c(
  "from_cell", "to_cell", "rank", "wnn_score",
  "atac_similarity", "rna_similarity",
  "atac_component", "rna_component",
  "candidate_from_atac", "candidate_from_rna",
  "from_idx", "to_idx"
))

directed_out <- file.path(paths$out_dir, sprintf("wnn_lite_%s_knn_directed_k%d.csv", gated_label, k))
fwrite(edge_dt[, !c("from_idx", "to_idx")], directed_out)

edge_dt[, edge_key := paste(pmin(from_idx, to_idx), pmax(from_idx, to_idx), sep = "_")]
edge_counts <- edge_dt[, .N, by = edge_key]
mutual_keys <- edge_counts[N >= 2, edge_key]
mutual_dt <- edge_dt[edge_key %in% mutual_keys]
mutual_pairs <- mutual_dt[, .(
  i = min(from_idx),
  j = max(to_idx),
  mean_rank = mean(rank),
  min_rank = min(rank),
  mean_wnn_score = mean(wnn_score),
  mean_atac_component = mean(atac_component),
  mean_rna_component = mean(rna_component),
  both_view_fraction = mean(candidate_from_atac & candidate_from_rna)
), by = edge_key]

adj <- Matrix::sparseMatrix(
  i = c(mutual_pairs$i, mutual_pairs$j),
  j = c(mutual_pairs$j, mutual_pairs$i),
  x = 1,
  dims = c(n, n),
  dimnames = list(common, common)
)

degree <- Matrix::rowSums(adj)
degree_dt <- data.table(
  atac_cell = common,
  degree = as.numeric(degree),
  gated = weights$gated,
  rna_weight = weights$rna_weight,
  atac_weight = weights$atac_weight,
  view_jaccard_atac_vs_rna = view_overlap,
  n_candidates = candidate_count,
  selected_from_atac_frac = selected_atac_fraction,
  selected_from_rna_frac = selected_rna_fraction,
  selected_from_both_frac = selected_both_fraction,
  selected_atac_component_frac = selected_atac_component_fraction,
  selected_rna_component_frac = 1 - selected_atac_component_fraction
)

cell_type_map <- load_atac_cell_types(paths)
same_type_summary <- NULL
if (!is.null(cell_type_map)) {
  degree_dt[, cell_type := cell_type_map[atac_cell]]
  directed_type <- edge_dt[, .(
    from_cell,
    to_cell,
    from_type = cell_type_map[from_cell],
    to_type = cell_type_map[to_cell]
  )]
  directed_type[, same_type := !is.na(from_type) & !is.na(to_type) & from_type == to_type]
  same_type_by_cell <- directed_type[, .(
    directed_same_type_frac = mean(same_type, na.rm = TRUE)
  ), by = from_cell]
  degree_dt[, directed_same_type_frac := same_type_by_cell$directed_same_type_frac[
    match(atac_cell, same_type_by_cell$from_cell)
  ]]

  mutual_edge_dt <- data.table(
    from_cell = c(common[mutual_pairs$i], common[mutual_pairs$j]),
    to_cell = c(common[mutual_pairs$j], common[mutual_pairs$i])
  )
  mutual_edge_dt[, from_type := cell_type_map[from_cell]]
  mutual_edge_dt[, to_type := cell_type_map[to_cell]]
  mutual_edge_dt[, same_type := !is.na(from_type) & !is.na(to_type) & from_type == to_type]
  mutual_same_by_cell <- mutual_edge_dt[, .(
    mutual_same_type_frac = mean(same_type, na.rm = TRUE)
  ), by = from_cell]
  degree_dt[, mutual_same_type_frac := mutual_same_by_cell$mutual_same_type_frac[
    match(atac_cell, mutual_same_by_cell$from_cell)
  ]]

  same_type_summary <- degree_dt[
    !is.na(cell_type),
    .(
      n_cells = .N,
      n_gated = sum(gated),
      pct_gated = mean(gated),
      mean_degree = mean(degree),
      mean_directed_same_type_frac = mean(directed_same_type_frac, na.rm = TRUE),
      mean_mutual_same_type_frac = mean(mutual_same_type_frac, na.rm = TRUE)
    ),
    by = cell_type
  ]
  same_type_out <- file.path(paths$out_dir, sprintf("wnn_lite_%s_same_type_by_celltype_k%d.csv", gated_label, k))
  fwrite(same_type_summary, same_type_out)
} else {
  same_type_out <- NA_character_
}

diagnostics <- data.table(
  version = gated_version,
  cutoff = cutoff,
  k = k,
  candidate_k = candidate_k,
  l2_norm = l2_norm,
  n_cells = n,
  n_gated = sum(weights$gated),
  pct_gated = mean(weights$gated),
  n_directed_edges = nrow(edge_dt),
  n_mutual_pairs = nrow(mutual_pairs),
  n_undirected_edges_in_sparse_matrix = length(adj@x) / 2,
  min_degree = min(degree),
  q10_degree = as.numeric(stats::quantile(degree, 0.10)),
  q25_degree = as.numeric(stats::quantile(degree, 0.25)),
  median_degree = as.numeric(stats::median(degree)),
  mean_degree = mean(degree),
  q75_degree = as.numeric(stats::quantile(degree, 0.75)),
  q90_degree = as.numeric(stats::quantile(degree, 0.90)),
  max_degree = max(degree),
  n_isolated = sum(degree == 0),
  pct_isolated = mean(degree == 0),
  mean_rna_weight = mean(weights$rna_weight),
  median_rna_weight = stats::median(weights$rna_weight),
  mean_rna_weight_gated = if (any(weights$gated)) mean(weights[gated == TRUE]$rna_weight) else NA_real_,
  mean_rna_weight_non_gated = if (any(!weights$gated)) mean(weights[gated == FALSE]$rna_weight) else NA_real_,
  mean_view_jaccard_atac_vs_rna = mean(view_overlap),
  median_view_jaccard_atac_vs_rna = stats::median(view_overlap),
  mean_selected_from_atac_frac = mean(selected_atac_fraction),
  mean_selected_from_rna_frac = mean(selected_rna_fraction),
  mean_selected_from_both_frac = mean(selected_both_fraction),
  mean_selected_atac_component_frac = mean(selected_atac_component_fraction),
  mean_selected_rna_component_frac = 1 - mean(selected_atac_component_fraction),
  mean_selected_rna_component_frac_gated = if (any(weights$gated)) {
    mean((1 - selected_atac_component_fraction)[weights$gated])
  } else {
    NA_real_
  },
  mean_selected_rna_component_frac_non_gated = if (any(!weights$gated)) {
    mean((1 - selected_atac_component_fraction)[!weights$gated])
  } else {
    NA_real_
  },
  rann_seconds = nn_elapsed,
  scoring_seconds = score_elapsed
)

adj_out <- file.path(paths$out_dir, sprintf("wnn_lite_%s_mutualknn_k%d.rds", gated_label, k))
diag_out <- file.path(paths$out_dir, sprintf("wnn_lite_%s_knn_diagnostics_k%d.csv", gated_label, k))
degree_out <- file.path(paths$out_dir, sprintf("wnn_lite_%s_knn_degree_by_cell_k%d.csv", gated_label, k))
mutual_out <- file.path(paths$out_dir, sprintf("wnn_lite_%s_mutualknn_edges_k%d.csv", gated_label, k))

saveRDS(adj, adj_out, compress = FALSE)
fwrite(diagnostics, diag_out)
fwrite(degree_dt, degree_out)
fwrite(mutual_pairs[, .(
  from_cell = common[i],
  to_cell = common[j],
  mean_rank,
  min_rank,
  mean_wnn_score,
  mean_atac_component,
  mean_rna_component,
  both_view_fraction
)], mutual_out)

p_degree <- ggplot(degree_dt, aes(x = degree, fill = gated)) +
  geom_histogram(binwidth = 1, color = "white", linewidth = 0.15, alpha = 0.85, position = "identity") +
  geom_vline(xintercept = diagnostics$median_degree, color = "#D62728", linewidth = 0.8) +
  scale_fill_manual(values = c("FALSE" = "#4C78A8", "TRUE" = "#F58518")) +
  labs(
    title = sprintf("WNN-lite %s mutualKNN degree (k=%d)", gated_label, k),
    subtitle = sprintf(
      "n_isolated=%d (%.3f%%), median=%s, gated=%.2f%%",
      diagnostics$n_isolated,
      diagnostics$pct_isolated * 100,
      round(diagnostics$median_degree, 2),
      diagnostics$pct_gated * 100
    ),
    x = "Mutual KNN degree",
    y = "ATAC cell count",
    fill = "Gated"
  ) +
  theme_bw(base_size = 12)
fig_degree <- file.path(paths$fig_dir, sprintf("Fig_wnn_lite_%s_degree_k%d.png", gated_label, k))
ggsave(fig_degree, p_degree, width = 7, height = 4.8, dpi = 300)

p_contrib <- ggplot(degree_dt, aes(x = selected_rna_component_frac, fill = gated)) +
  geom_histogram(bins = 70, color = "white", linewidth = 0.15, alpha = 0.85, position = "identity") +
  scale_fill_manual(values = c("FALSE" = "#54A24B", "TRUE" = "#F58518")) +
  labs(
    title = sprintf("RNA component fraction after %s gating", gated_label),
    x = "RNA component / WNN-lite score",
    y = "ATAC cell count",
    fill = "Gated"
  ) +
  theme_bw(base_size = 12)
fig_contrib <- file.path(paths$fig_dir, sprintf("Fig_wnn_lite_%s_component_fraction_k%d.png", gated_label, k))
ggsave(fig_contrib, p_contrib, width = 6.5, height = 4.6, dpi = 300)

manifest <- list(
  script = "_build_gated_knn_common.R",
  entry_version = gated_version,
  gated_label = gated_label,
  cutoff = cutoff,
  k = k,
  candidate_k = candidate_k,
  l2_norm = l2_norm,
  shared_embedding = paths$shared_embedding,
  rna_reference = rna_ref_file,
  weights = gated_weights_file,
  graph_construction = "ATAC-view and RNA-reference-view KNN candidate union; candidates rescored by gated cell-specific weighted similarities; final top-k converted to mutualKNN",
  gating_behavior = "Cells marked gated keep candidate union but score RNA component with rna_weight=0 and ATAC component with atac_weight=1",
  similarity = "exp(-(distance - per_cell_min_distance) / per_cell_median_delta)",
  output_rds = adj_out,
  diagnostics_csv = diag_out,
  degree_by_cell_csv = degree_out,
  directed_knn_csv = directed_out,
  mutual_edges_csv = mutual_out,
  same_type_by_celltype_csv = same_type_out,
  degree_figure = fig_degree,
  contribution_figure = fig_contrib,
  created_at = as.character(Sys.time())
)
jsonlite::write_json(
  manifest,
  file.path(paths$out_dir, sprintf("wnn_lite_%s_knn_manifest_k%d.json", gated_label, k)),
  pretty = TRUE,
  auto_unbox = TRUE
)

cat("\n[Diagnostics]\n")
print(diagnostics)
cat("\n[Outputs]\n")
cat("  graph      :", adj_out, "\n")
cat("  diagnostics:", diag_out, "\n")
cat("  degree     :", degree_out, "\n")
cat("  same type  :", same_type_out, "\n")
cat("========== Done ==========\n")

if (diagnostics$pct_isolated > 0.01) {
  warning("More than 1% isolated nodes. Inspect degree CSV before running TRS.")
}
