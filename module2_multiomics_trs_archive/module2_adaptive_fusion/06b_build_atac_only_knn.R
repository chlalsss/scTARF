# Build mutualKNN graph from scSHEFT ATAC latent only (no PCA, no RNA reference).
#
# Purpose: baseline that uses only training-time multi-omics fusion (the scSHEFT
# shared encoder), with no inference-time fusion. Used to compare against:
#   - original LSI mutualKNN (no fusion at all)
#   - fused mutualKNN from 06 (training-time + inference-time fusion)

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

k_raw <- Sys.getenv("MODULE2_KNN_K", "auto")
k_auto <- tolower(k_raw) == "auto"
k <- if (k_auto) NA_integer_ else as.integer(k_raw)
if (!k_auto && (!is.finite(k) || k < 1)) {
  stop("MODULE2_KNN_K must be 'auto' or a positive integer")
}

cat("========== Module 2: build scSHEFT ATAC-only mutualKNN ==========\n")
cat("  shared embedding:", paths$shared_embedding, "\n")
cat("  KNN k           :", if (k_auto) "auto" else k, "\n")

stopifnot(file.exists(paths$shared_embedding))

loaded <- load_shared_embedding(paths$shared_embedding)
emb <- as.data.table(loaded$emb)
latent_cols <- loaded$latent_cols
atac_dt <- emb[modality == "atac"]
if (nrow(atac_dt) == 0) stop("No ATAC rows in shared embedding")
if (anyDuplicated(atac_dt$cell)) stop("Duplicated ATAC cell barcodes")

mat <- as.matrix(atac_dt[, ..latent_cols])
storage.mode(mat) <- "double"
rownames(mat) <- atac_dt$cell
cat("  cells:", nrow(mat), "\n")
cat("  dims :", ncol(mat), "\n")

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
  auto_k_summary_csv <- file.path(paths$out_dir, "atac_only_knn_auto_k_search.csv")
  fwrite(auto_k_summary, auto_k_summary_csv)
  write_auto_value(file.path(paths$out_dir, "atac_only_knn_auto_k.txt"), k)
  cat("  auto KNN k      :", k, "\n")
}

t0 <- Sys.time()
nn <- RANN::nn2(data = mat, query = mat, k = k + 1)
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("  RANN elapsed: %.1f sec\n", elapsed))

n <- nrow(mat)
neighbor_idx <- matrix(NA_integer_, nrow = n, ncol = k)
neighbor_dist <- matrix(NA_real_, nrow = n, ncol = k)
for (i in seq_len(n)) {
  idx_i <- nn$nn.idx[i, ]
  dist_i <- nn$nn.dists[i, ]
  keep <- idx_i != i
  idx_i <- idx_i[keep]
  dist_i <- dist_i[keep]
  if (length(idx_i) < k) stop("Could not find enough non-self neighbors for row ", i)
  neighbor_idx[i, ] <- idx_i[seq_len(k)]
  neighbor_dist[i, ] <- dist_i[seq_len(k)]
}

edge_dt <- data.table(
  from_idx = rep(seq_len(n), each = k),
  to_idx = as.vector(t(neighbor_idx)),
  rank = rep(seq_len(k), times = n),
  distance = as.vector(t(neighbor_dist))
)
edge_dt[, edge_key := paste(pmin(from_idx, to_idx), pmax(from_idx, to_idx), sep = "_")]
edge_counts <- edge_dt[, .N, by = edge_key]
mutual_keys <- edge_counts[N >= 2, edge_key]
mutual_dt <- edge_dt[edge_key %in% mutual_keys]
mutual_pairs <- mutual_dt[, .(
  i = min(from_idx),
  j = max(to_idx),
  mean_distance = mean(distance),
  min_rank = min(rank)
), by = edge_key]

adj <- Matrix::sparseMatrix(
  i = c(mutual_pairs$i, mutual_pairs$j),
  j = c(mutual_pairs$j, mutual_pairs$i),
  x = 1,
  dims = c(n, n),
  dimnames = list(rownames(mat), rownames(mat))
)

degree <- Matrix::rowSums(adj)
diag_dt <- data.table(atac_cell = rownames(mat), degree = as.numeric(degree))

diagnostics <- data.table(
  variant = "atac_only_scSHEFT_latent",
  k = k,
  k_mode = if (k_auto) "auto_min_isolated" else "manual",
  k_requested = k_raw,
  k_auto_summary_csv = auto_k_summary_csv,
  n_cells = n,
  n_mutual_pairs = nrow(mutual_pairs),
  min_degree = min(degree),
  q10_degree = as.numeric(stats::quantile(degree, 0.10)),
  median_degree = as.numeric(stats::median(degree)),
  mean_degree = mean(degree),
  q90_degree = as.numeric(stats::quantile(degree, 0.90)),
  max_degree = max(degree),
  n_isolated = sum(degree == 0),
  pct_isolated = mean(degree == 0),
  rann_seconds = elapsed
)

adj_out <- file.path(paths$out_dir, sprintf("atac_only_mutualknn_k%d.rds", k))
diag_out <- file.path(paths$out_dir, sprintf("atac_only_knn_diagnostics_k%d.csv", k))
degree_out <- file.path(paths$out_dir, sprintf("atac_only_knn_degree_by_cell_k%d.csv", k))

saveRDS(adj, adj_out, compress = FALSE)
fwrite(diagnostics, diag_out)
fwrite(diag_dt, degree_out)

fig <- ggplot(diag_dt, aes(x = degree)) +
  geom_histogram(binwidth = 1, fill = "#54A24B", color = "white", linewidth = 0.15) +
  geom_vline(xintercept = diagnostics$median_degree, color = "#D62728", linewidth = 0.8) +
  labs(
    title = sprintf("scSHEFT ATAC-only mutualKNN degree (k=%d)", k),
    subtitle = sprintf(
      "n_isolated=%d (%.3f%%), median=%s, q10=%s",
      diagnostics$n_isolated,
      diagnostics$pct_isolated * 100,
      round(diagnostics$median_degree, 2),
      round(diagnostics$q10_degree, 2)
    ),
    x = "Degree", y = "ATAC cell count"
  ) +
  theme_bw(base_size = 12)
fig_out <- file.path(paths$fig_dir, sprintf("Fig_atac_only_knn_degree_k%d.png", k))
ggsave(fig_out, fig, width = 7, height = 4.8, dpi = 300)

manifest <- list(
  script = "06b_build_atac_only_knn.R",
  variant = "atac_only_scSHEFT_latent",
  shared_embedding = paths$shared_embedding,
  k = k,
  k_mode = if (k_auto) "auto_min_isolated" else "manual",
  k_requested = k_raw,
  k_auto_summary_csv = auto_k_summary_csv,
  graph_construction = "RANN::nn2 on 128-dim scSHEFT ATAC latent + mutualKNN",
  output_rds = adj_out,
  diagnostics_csv = diag_out,
  degree_by_cell_csv = degree_out,
  degree_figure = fig_out,
  created_at = as.character(Sys.time())
)
jsonlite::write_json(
  manifest,
  file.path(paths$out_dir, sprintf("atac_only_knn_manifest_k%d.json", k)),
  pretty = TRUE,
  auto_unbox = TRUE
)

cat("  saved graph      :", adj_out, "\n")
cat("  saved diagnostics:", diag_out, "\n")
cat("  saved figure     :", fig_out, "\n")
cat("\n  degree summary\n")
print(diagnostics)
if (diagnostics$pct_isolated > 0.01) {
  warning("More than 1% isolated nodes. Consider increasing k.")
}
cat("========== Done ==========\n")
