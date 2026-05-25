# Build fused ATAC-ATAC mutual KNN graph from module 2 fused embedding.

source("D:/RStudio/work/module2_adaptive_fusion/_bootstrap.R")

required_pkgs(c("data.table", "Matrix", "RANN", "ggplot2", "jsonlite"))
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(RANN)
  library(ggplot2)
  library(jsonlite)
})

paths <- module2_paths()
ensure_module2_dirs(paths)

label <- get_active_label()
cutoff <- if (label == "top1" || startsWith(label, "top1disease")) NA_real_ else get_single_cutoff()
npcs_raw <- Sys.getenv("MODULE2_FUSED_NPCS", "auto")
npcs <- if (tolower(npcs_raw) == "auto") {
  read_auto_value(
    file.path(paths$out_dir, sprintf("fused_embedding_auto_npcs_cutoff_%s.txt", label)),
    "fused PCA dimension"
  )
} else {
  as.integer(npcs_raw)
}
if (!is.finite(npcs) || npcs < 1) stop("MODULE2_FUSED_NPCS must be 'auto' or a positive integer")

k_raw <- Sys.getenv("MODULE2_KNN_K", "auto")
k_auto <- tolower(k_raw) == "auto"
k <- if (k_auto) NA_integer_ else as.integer(k_raw)
if (!k_auto && (!is.finite(k) || k < 1)) {
  stop("MODULE2_KNN_K must be 'auto' or a positive integer")
}

fused_rds <- file.path(paths$out_dir, sprintf("fused_embedding_cutoff_%s_pca%d.rds", label, npcs))
if (!file.exists(fused_rds)) {
  stop("Missing fused embedding. Run 04_build_fused_embedding.R first: ", fused_rds)
}

cat("========== Module 2: build fused mutualKNN ==========\n")
cat("  fused embedding:", fused_rds, "\n")
cat("  cutoff         :", cutoff, "\n")
cat("  PCA dims       :", npcs, "\n")
cat("  KNN k          :", if (k_auto) "auto" else k, "\n")

fused <- as.data.table(readRDS(fused_rds))
pc_cols <- grep("^fusedPC_", colnames(fused), value = TRUE)
if (!"atac_cell" %in% colnames(fused) || length(pc_cols) == 0) {
  stop("Fused embedding must contain atac_cell and fusedPC_* columns")
}
if (anyDuplicated(fused$atac_cell)) {
  stop("Duplicated atac_cell values in fused embedding")
}

mat <- as.matrix(fused[, ..pc_cols])
storage.mode(mat) <- "double"
rownames(mat) <- fused$atac_cell

cat("  cells:", nrow(mat), "\n")
cat("  dims :", ncol(mat), "\n")

auto_k_summary <- NULL
auto_k_summary_csv <- NA_character_
if (k_auto) {
  target_iso <- as.numeric(Sys.getenv("MODULE2_KNN_TARGET_ISO", "0.01"))
  if (!is.finite(target_iso) || target_iso <= 0 || target_iso >= 1) {
    stop("MODULE2_KNN_TARGET_ISO must be between 0 and 1")
  }
  k_search <- parse_int_sequence_env("MODULE2_KNN_K_SEARCH", 10:50)
  auto_k <- auto_select_knn_k(mat, target_iso = target_iso, k_search = k_search)
  k <- auto_k$k
  auto_k_summary <- as.data.table(auto_k$search)
  auto_k_summary[, selected := k == auto_k$k]
  auto_k_summary_csv <- file.path(
    paths$out_dir,
    sprintf("fused_knn_auto_k_search_cutoff_%s_pca%d.csv", label, npcs)
  )
  fwrite(auto_k_summary, auto_k_summary_csv)
  cat("  auto KNN k    :", k, "\n")
}
k_auto_file <- file.path(paths$out_dir, sprintf("fused_knn_auto_k_cutoff_%s_pca%d.txt", label, npcs))
write_auto_value(k_auto_file, k)

t0 <- Sys.time()
nn <- RANN::nn2(data = mat, query = mat, k = k + 1)
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("  RANN elapsed: %.1f sec\n", elapsed))

# Remove self neighbor if present. Keep exactly k non-self neighbors per cell.
n <- nrow(mat)
neighbor_idx <- matrix(NA_integer_, nrow = n, ncol = k)
neighbor_dist <- matrix(NA_real_, nrow = n, ncol = k)
for (i in seq_len(n)) {
  idx_i <- nn$nn.idx[i, ]
  dist_i <- nn$nn.dists[i, ]
  keep <- idx_i != i
  idx_i <- idx_i[keep]
  dist_i <- dist_i[keep]
  if (length(idx_i) < k) {
    stop("Could not find enough non-self neighbors for row ", i)
  }
  neighbor_idx[i, ] <- idx_i[seq_len(k)]
  neighbor_dist[i, ] <- dist_i[seq_len(k)]
}

edge_dt <- data.table(
  from_idx = rep(seq_len(n), each = k),
  to_idx = as.vector(t(neighbor_idx)),
  rank = rep(seq_len(k), times = n),
  distance = as.vector(t(neighbor_dist))
)
edge_dt[, from_cell := rownames(mat)[from_idx]]
edge_dt[, to_cell := rownames(mat)[to_idx]]
edge_dt[, sim_scale := stats::median(distance), by = from_idx]
edge_dt[!is.finite(sim_scale) | sim_scale <= 1e-12, sim_scale := 1]
edge_dt[, similarity := exp(-distance / sim_scale)]

directed_out <- file.path(paths$out_dir, sprintf("fused_knn_directed_cutoff_%s_pca%d_k%d.csv", label, npcs, k))
fwrite(edge_dt[, .(from_cell, to_cell, rank, distance)], directed_out)

# Mutual edges: keep i -> j only if j -> i exists. Store symmetric adjacency.
edge_dt[, edge_key := paste(pmin(from_idx, to_idx), pmax(from_idx, to_idx), sep = "_")]
edge_counts <- edge_dt[, .N, by = edge_key]
mutual_keys <- edge_counts[N >= 2, edge_key]
mutual_dt <- edge_dt[edge_key %in% mutual_keys]

# One undirected row per mutual pair, distance averaged across both directions.
mutual_pairs <- mutual_dt[, .(
  i = min(from_idx),
  j = max(to_idx),
  mean_distance = mean(distance),
  mean_similarity = mean(similarity),
  min_rank = min(rank)
), by = edge_key]

adj <- Matrix::sparseMatrix(
  i = c(mutual_pairs$i, mutual_pairs$j),
  j = c(mutual_pairs$j, mutual_pairs$i),
  x = c(mutual_pairs$mean_similarity, mutual_pairs$mean_similarity),
  dims = c(n, n),
  dimnames = list(rownames(mat), rownames(mat))
)

degree <- Matrix::rowSums(adj)
edge_w <- mutual_pairs$mean_similarity
edge_w_p <- edge_w / sum(edge_w)
edge_w_entropy <- -sum(edge_w_p * log(edge_w_p + 1e-30))
edge_w_rel_entropy <- edge_w_entropy / log(length(edge_w))
diag_dt <- data.table(
  atac_cell = rownames(mat),
  degree = as.numeric(degree)
)

diagnostics <- data.table(
  cutoff = cutoff,
  npcs = npcs,
  k = k,
  n_cells = n,
  n_directed_edges = nrow(edge_dt),
  n_mutual_pairs = nrow(mutual_pairs),
  n_undirected_edges_in_sparse_matrix = length(adj@x) / 2,
  edge_weight_entropy = edge_w_entropy,
  edge_weight_relative_entropy = edge_w_rel_entropy,
  k_mode = if (k_auto) "auto_min_isolated" else "manual",
  k_requested = k_raw,
  k_auto_summary_csv = auto_k_summary_csv,
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
  rann_seconds = elapsed
)

adj_out <- file.path(paths$out_dir, sprintf("fused_mutualknn_cutoff_%s_pca%d_k%d.rds", label, npcs, k))
diag_out <- file.path(paths$out_dir, sprintf("fused_knn_diagnostics_cutoff_%s_pca%d_k%d.csv", label, npcs, k))
degree_out <- file.path(paths$out_dir, sprintf("fused_knn_degree_by_cell_cutoff_%s_pca%d_k%d.csv", label, npcs, k))
mutual_out <- file.path(paths$out_dir, sprintf("fused_mutualknn_edges_cutoff_%s_pca%d_k%d.csv", label, npcs, k))

saveRDS(adj, adj_out, compress = FALSE)
fwrite(diagnostics, diag_out)
fwrite(diag_dt, degree_out)
fwrite(mutual_pairs[, .(
  from_cell = rownames(mat)[i],
  to_cell = rownames(mat)[j],
  mean_distance,
  mean_similarity,
  min_rank
)], mutual_out)

fig <- ggplot(diag_dt, aes(x = degree)) +
  geom_histogram(binwidth = 1, fill = "#4C78A8", color = "white", linewidth = 0.15) +
  geom_vline(xintercept = diagnostics$median_degree, color = "#D62728", linewidth = 0.8) +
  labs(
    title = sprintf("Fused mutualKNN degree distribution (k=%d)", k),
    subtitle = sprintf(
      "n_isolated=%d (%.3f%%), median=%s, q10=%s",
      diagnostics$n_isolated,
      diagnostics$pct_isolated * 100,
      round(diagnostics$median_degree, 2),
      round(diagnostics$q10_degree, 2)
    ),
    x = "Mutual KNN degree",
    y = "ATAC cell count"
  ) +
  theme_bw(base_size = 12)

fig_out <- file.path(paths$fig_dir, sprintf("Fig_fused_knn_degree_cutoff_%s_pca%d_k%d.png", label, npcs, k))
ggsave(fig_out, fig, width = 7, height = 4.8, dpi = 300)

manifest <- list(
  script = "06_build_fused_knn.R",
  fused_embedding = fused_rds,
  cutoff = cutoff,
  npcs = npcs,
  k = k,
  k_mode = if (k_auto) "auto_min_isolated" else "manual",
  k_requested = k_raw,
  k_auto_file = k_auto_file,
  k_auto_summary_csv = auto_k_summary_csv,
  graph = "mutual KNN from Euclidean nearest neighbors on fused PCA embedding",
  output_rds = adj_out,
  diagnostics_csv = diag_out,
  degree_by_cell_csv = degree_out,
  directed_knn_csv = directed_out,
  mutual_edges_csv = mutual_out,
  degree_figure = fig_out,
  created_at = as.character(Sys.time())
)
jsonlite::write_json(
  manifest,
  file.path(paths$out_dir, sprintf("fused_knn_manifest_cutoff_%s_pca%d_k%d.json", label, npcs, k)),
  pretty = TRUE,
  auto_unbox = TRUE
)

cat("  saved graph      :", adj_out, "\n")
cat("  saved diagnostics:", diag_out, "\n")
cat("  saved figure     :", fig_out, "\n")
cat("\n  degree summary\n")
print(diagnostics)
if (diagnostics$pct_isolated > 0.01) {
  warning("More than 1% isolated nodes. Consider increasing k or using simpleKNN/fallback edges before TRS.")
}
cat("========== Done ==========\n")

