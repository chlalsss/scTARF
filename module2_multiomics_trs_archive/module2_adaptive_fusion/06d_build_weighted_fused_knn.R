# Build mutualKNN graph from weighted-concat fused embedding.
#
# Difference from 04 + 06:
#   - No scale() — preserves the natural norm ratio between ATAC and RNA_ref.
#   - Configurable RNA_ref weight via MODULE2_FUSION_WEIGHT (default 0.30).
#   - PCA optional: MODULE2_FUSION_PCA = "none" (default, 256-dim KNN) or an integer.

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

cutoff <- get_single_cutoff()
label <- cutoff_label(cutoff)
k <- as.integer(Sys.getenv("MODULE2_KNN_K", "30"))
if (!is.finite(k) || k < 1) stop("MODULE2_KNN_K must be a positive integer")

w <- as.numeric(Sys.getenv("MODULE2_FUSION_WEIGHT", "0.30"))
if (!is.finite(w) || w < 0) stop("MODULE2_FUSION_WEIGHT must be non-negative numeric")
w_label <- gsub("\\.", "p", sprintf("%.2f", w))

pca_str <- Sys.getenv("MODULE2_FUSION_PCA", "none")
do_pca <- !(pca_str %in% c("", "none", "NONE", "None"))
npcs <- if (do_pca) as.integer(pca_str) else NA_integer_
pca_label <- if (do_pca) sprintf("pca%d", npcs) else "raw256"

cat("========== Module 2: build weighted fused mutualKNN ==========\n")
cat("  cutoff (rna_ref):", cutoff, "\n")
cat("  RNA weight w    :", w, "  (label:", w_label, ")\n")
cat("  PCA mode        :", if (do_pca) sprintf("PCA %d", npcs) else "none", "  (label:", pca_label, ")\n")
cat("  KNN k           :", k, "\n")

rna_ref_rds <- file.path(paths$out_dir, sprintf("rna_reference_latent_cutoff_%s.rds", label))
stopifnot(file.exists(paths$shared_embedding), file.exists(rna_ref_rds))

loaded <- load_shared_embedding(paths$shared_embedding)
emb <- as.data.table(loaded$emb)
latent_cols <- loaded$latent_cols
atac_dt <- emb[modality == "atac"]
rna_ref <- as.data.table(readRDS(rna_ref_rds))
rna_ref_cols <- grep("^rna_ref_", colnames(rna_ref), value = TRUE)

if (length(latent_cols) != length(rna_ref_cols)) {
  stop("Latent dim mismatch: ATAC ", length(latent_cols), " vs RNA_ref ", length(rna_ref_cols))
}

common <- intersect(atac_dt$cell, rna_ref$atac_cell)
if (length(common) == 0) stop("No common ATAC cells between shared embedding and RNA reference")
atac_dt <- atac_dt[match(common, cell)]
rna_ref <- rna_ref[match(common, atac_cell)]

atac_mat <- as.matrix(atac_dt[, ..latent_cols])
ref_mat <- as.matrix(rna_ref[, ..rna_ref_cols])
storage.mode(atac_mat) <- "double"
storage.mode(ref_mat) <- "double"

fused_mat <- cbind(atac_mat, w * ref_mat)
fused_mat[!is.finite(fused_mat)] <- 0
rownames(fused_mat) <- common
cat("  cells :", nrow(fused_mat), "\n")
cat("  pre-PCA dims :", ncol(fused_mat), "\n")

if (do_pca) {
  if (requireNamespace("irlba", quietly = TRUE)) {
    pca <- irlba::prcomp_irlba(fused_mat, n = npcs, center = FALSE, scale. = FALSE)
  } else {
    pca <- stats::prcomp(fused_mat, rank. = npcs, center = FALSE, scale. = FALSE)
  }
  pc_mat <- pca$x[, seq_len(min(npcs, ncol(pca$x))), drop = FALSE]
  rownames(pc_mat) <- common
  fused_mat <- pc_mat
  cat("  post-PCA dims:", ncol(fused_mat), "\n")
}

t0 <- Sys.time()
nn <- RANN::nn2(data = fused_mat, query = fused_mat, k = k + 1)
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("  RANN elapsed: %.1f sec\n", elapsed))

n <- nrow(fused_mat)
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
  dimnames = list(common, common)
)

degree <- Matrix::rowSums(adj)
diag_dt <- data.table(atac_cell = common, degree = as.numeric(degree))

diagnostics <- data.table(
  variant = sprintf("weighted_w%s_%s", w_label, pca_label),
  rna_weight = w,
  pca_dims = if (do_pca) npcs else NA_integer_,
  k = k,
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

base <- sprintf("weighted_w%s_%s", w_label, pca_label)
adj_out <- file.path(paths$out_dir, sprintf("%s_mutualknn_k%d.rds", base, k))
diag_out <- file.path(paths$out_dir, sprintf("%s_knn_diagnostics_k%d.csv", base, k))
degree_out <- file.path(paths$out_dir, sprintf("%s_knn_degree_by_cell_k%d.csv", base, k))

saveRDS(adj, adj_out, compress = FALSE)
fwrite(diagnostics, diag_out)
fwrite(diag_dt, degree_out)

fig <- ggplot(diag_dt, aes(x = degree)) +
  geom_histogram(binwidth = 1, fill = "#F58518", color = "white", linewidth = 0.15) +
  geom_vline(xintercept = diagnostics$median_degree, color = "#D62728", linewidth = 0.8) +
  labs(
    title = sprintf("Weighted fused mutualKNN degree (w=%.2f, %s, k=%d)",
                    w, pca_label, k),
    subtitle = sprintf(
      "n_isolated=%d (%.3f%%), median=%s, q10=%s",
      diagnostics$n_isolated, diagnostics$pct_isolated * 100,
      round(diagnostics$median_degree, 2), round(diagnostics$q10_degree, 2)
    ),
    x = "Degree", y = "ATAC cell count"
  ) +
  theme_bw(base_size = 12)
fig_out <- file.path(paths$fig_dir, sprintf("Fig_%s_knn_degree_k%d.png", base, k))
ggsave(fig_out, fig, width = 7, height = 4.8, dpi = 300)

manifest <- list(
  script = "06d_build_weighted_fused_knn.R",
  variant = diagnostics$variant,
  rna_weight = w,
  pca = if (do_pca) npcs else "none",
  shared_embedding = paths$shared_embedding,
  rna_reference = rna_ref_rds,
  k = k,
  fusion_formula = "cbind(atac_latent, w * rna_ref_latent), no scale()",
  output_rds = adj_out,
  diagnostics_csv = diag_out,
  degree_by_cell_csv = degree_out,
  degree_figure = fig_out,
  created_at = as.character(Sys.time())
)
jsonlite::write_json(
  manifest,
  file.path(paths$out_dir, sprintf("%s_knn_manifest_k%d.json", base, k)),
  pretty = TRUE,
  auto_unbox = TRUE
)

cat("  saved graph      :", adj_out, "\n")
cat("  saved diagnostics:", diag_out, "\n")
cat("\n  degree summary\n")
print(diagnostics)
if (diagnostics$pct_isolated > 0.01) {
  warning("More than 1% isolated nodes. Consider increasing k or w.")
}
cat("========== Done ==========\n")
