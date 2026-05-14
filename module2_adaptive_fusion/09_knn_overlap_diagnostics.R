# Diagnose why fused KNN did not improve TRS over original.
#
# Reads existing artifacts only; writes to a new diagnostics subdirectory.
#   1. Edge Jaccard between fused mutualKNN and original mutualKNN.
#   2. Per-cell same-celltype neighbor fraction (overall + by celltype).
#   3. RNA reference residual ||RNA_ref - ATAC_latent|| distribution.

source("D:/RStudio/work/module2_adaptive_fusion/_bootstrap.R")

required_pkgs(c("data.table", "Matrix", "ggplot2", "jsonlite", "SummarizedExperiment"))
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(ggplot2)
  library(jsonlite)
  library(SummarizedExperiment)
})

paths <- module2_paths()
ensure_module2_dirs(paths)

cutoff <- get_single_cutoff()
label <- cutoff_label(cutoff)
npcs <- as.integer(Sys.getenv("MODULE2_FUSED_NPCS", "50"))
k <- as.integer(Sys.getenv("MODULE2_KNN_K", "30"))

data_dir <- paths$data_dir
fused_knn_path <- file.path(
  paths$out_dir,
  sprintf("fused_mutualknn_cutoff_%s_pca%d_k%d.rds", label, npcs, k)
)
original_knn_path <- Sys.getenv(
  "MODULE2_ORIGINAL_KNN",
  file.path(data_dir, "scavenge_out_small_all_traits", "mutualknn_k30.rds")
)
se_path <- Sys.getenv(
  "MODULE2_ATAC_SE",
  file.path(data_dir, "scATAC-Healthy-Hematopoiesis-191120.rds")
)
rna_ref_path <- file.path(
  paths$out_dir,
  sprintf("rna_reference_latent_cutoff_%s.rds", label)
)

out_dir <- file.path(
  paths$out_dir,
  sprintf("diagnostics_knn_cutoff_%s_pca%d_k%d", label, npcs, k)
)
fig_dir <- file.path(out_dir, "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

cat("========== Module 2: KNN diagnostics ==========\n")
cat("  fused KNN   :", fused_knn_path, "\n")
cat("  original KNN:", original_knn_path, "\n")
cat("  output dir  :", out_dir, "\n")

stopifnot(
  file.exists(fused_knn_path),
  file.exists(original_knn_path),
  file.exists(se_path),
  file.exists(rna_ref_path),
  file.exists(paths$shared_embedding)
)

# =========================
# Load + align cell sets
# =========================
fused <- readRDS(fused_knn_path)
original <- readRDS(original_knn_path)
SE <- readRDS(se_path)
loaded <- load_shared_embedding(paths$shared_embedding)
emb <- as.data.table(loaded$emb)
latent_cols <- loaded$latent_cols
rna_ref <- as.data.table(readRDS(rna_ref_path))
rna_ref_cols <- grep("^rna_ref_", colnames(rna_ref), value = TRUE)

if (nrow(fused) != ncol(fused)) stop("Fused KNN must be square")
if (nrow(original) != ncol(original)) stop("Original KNN must be square")

fused_cells <- rownames(fused)
original_cells <- rownames(original)
common_cells <- intersect(fused_cells, original_cells)
cat("  cells fused   :", length(fused_cells), "\n")
cat("  cells original:", length(original_cells), "\n")
cat("  cells common  :", length(common_cells), "\n")
if (length(common_cells) == 0) stop("No common cells between fused and original KNN")

fused_aligned <- fused[common_cells, common_cells]
original_aligned <- original[common_cells, common_cells]

meta <- as.data.frame(colData(SE))
meta$cell <- colnames(SE)
celltype_map <- setNames(as.character(meta$BioClassification), meta$cell)
celltype_vec <- celltype_map[common_cells]
n_unmapped <- sum(is.na(celltype_vec))
if (n_unmapped > 0) {
  warning("Some cells missing BioClassification: ", n_unmapped)
  celltype_vec[is.na(celltype_vec)] <- "Unknown"
}

# =========================
# 1) Edge Jaccard
# =========================
extract_undirected_edges <- function(adj) {
  if (!inherits(adj, "Matrix")) adj <- as(adj, "CsparseMatrix")
  upper <- triu(adj, k = 1)
  upper <- as(upper, "TsparseMatrix")
  data.table(i = upper@i + 1L, j = upper@j + 1L)
}

f_edges <- extract_undirected_edges(fused_aligned)
o_edges <- extract_undirected_edges(original_aligned)
f_keys <- paste(common_cells[f_edges$i], common_cells[f_edges$j], sep = "|")
o_keys <- paste(common_cells[o_edges$i], common_cells[o_edges$j], sep = "|")

n_f <- length(f_keys)
n_o <- length(o_keys)
n_inter <- length(intersect(f_keys, o_keys))
n_union <- n_f + n_o - n_inter
jaccard <- if (n_union > 0) n_inter / n_union else NA_real_

edge_summary <- data.table(
  n_edges_fused = n_f,
  n_edges_original = n_o,
  n_edges_common = n_inter,
  n_edges_union = n_union,
  n_fused_only = n_f - n_inter,
  n_original_only = n_o - n_inter,
  jaccard = jaccard,
  pct_fused_overlap_with_original = if (n_f > 0) n_inter / n_f else NA_real_,
  pct_original_overlap_with_fused = if (n_o > 0) n_inter / n_o else NA_real_
)
fwrite(edge_summary, file.path(out_dir, "knn_edge_overlap.csv"))
cat("\n[Edge overlap]\n")
print(edge_summary)

# =========================
# 2) Same-celltype neighbor fraction
# =========================
neighbor_same_frac_per_cell <- function(adj, types) {
  if (!inherits(adj, "Matrix")) adj <- as(adj, "CsparseMatrix")
  adj_t <- as(adj, "TsparseMatrix")
  i_vec <- adj_t@i + 1L
  j_vec <- adj_t@j + 1L
  keep <- i_vec != j_vec
  i_vec <- i_vec[keep]
  j_vec <- j_vec[keep]
  type_i <- types[i_vec]
  type_j <- types[j_vec]
  same <- type_i == type_j
  n <- length(types)
  same_per_cell <- tabulate(i_vec[same], nbins = n)
  total_per_cell <- tabulate(i_vec, nbins = n)
  data.table(
    cell_idx = seq_len(n),
    type = types,
    n_same = same_per_cell,
    n_total = total_per_cell,
    same_frac = ifelse(total_per_cell > 0, same_per_cell / total_per_cell, NA_real_)
  )
}

f_per_cell <- neighbor_same_frac_per_cell(fused_aligned, celltype_vec)
o_per_cell <- neighbor_same_frac_per_cell(original_aligned, celltype_vec)
f_per_cell[, method := "fused"]
o_per_cell[, method := "original"]
per_cell_dt <- rbind(f_per_cell, o_per_cell)
fwrite(per_cell_dt, file.path(out_dir, "knn_same_type_per_cell.csv"))

overall <- per_cell_dt[
  n_total > 0,
  .(
    mean_same_type_frac = mean(same_frac, na.rm = TRUE),
    median_same_type_frac = median(same_frac, na.rm = TRUE),
    mean_cross_type_frac = 1 - mean(same_frac, na.rm = TRUE),
    n_cells_with_edges = .N
  ),
  by = method
]
fwrite(overall, file.path(out_dir, "knn_same_type_overall.csv"))
cat("\n[Same-type fraction overall]\n")
print(overall)

per_type <- per_cell_dt[
  n_total > 0,
  .(
    n_cells = .N,
    mean_same_type_frac = mean(same_frac, na.rm = TRUE),
    median_same_type_frac = median(same_frac, na.rm = TRUE)
  ),
  by = .(type, method)
]
per_type_wide <- dcast(
  per_type,
  type ~ method,
  value.var = c("n_cells", "mean_same_type_frac", "median_same_type_frac")
)
fwrite(per_type_wide, file.path(out_dir, "knn_same_type_by_celltype.csv"))

# =========================
# 3) RNA_ref residual distribution
# =========================
atac_emb <- emb[modality == "atac"]
common_ref_cells <- intersect(atac_emb$cell, rna_ref$atac_cell)
atac_emb <- atac_emb[match(common_ref_cells, cell)]
rna_ref <- rna_ref[match(common_ref_cells, atac_cell)]

if (length(latent_cols) != length(rna_ref_cols)) {
  stop("Latent dim mismatch: ATAC ", length(latent_cols), " vs RNA_ref ", length(rna_ref_cols))
}

atac_mat <- as.matrix(atac_emb[, ..latent_cols])
ref_mat <- as.matrix(rna_ref[, ..rna_ref_cols])
storage.mode(atac_mat) <- "double"
storage.mode(ref_mat) <- "double"

resid <- ref_mat - atac_mat
resid_norm <- sqrt(rowSums(resid^2))
atac_norm <- sqrt(rowSums(atac_mat^2))
ref_norm <- sqrt(rowSums(ref_mat^2))
rel_norm <- resid_norm / pmax(atac_norm, 1e-12)
cos_sim <- rowSums(atac_mat * ref_mat) / pmax(atac_norm * ref_norm, 1e-12)

resid_dt <- data.table(
  atac_cell = common_ref_cells,
  resid_norm = resid_norm,
  atac_norm = atac_norm,
  ref_norm = ref_norm,
  rel_norm = rel_norm,
  cos_sim = cos_sim
)
fwrite(resid_dt, file.path(out_dir, "rna_ref_residual_per_cell.csv"))

quantile_row <- function(x, name) {
  q <- as.numeric(quantile(x, c(0.05, 0.25, 0.50, 0.75, 0.95), na.rm = TRUE))
  data.table(
    metric = name,
    mean = mean(x, na.rm = TRUE),
    p05 = q[1], p25 = q[2], p50 = q[3], p75 = q[4], p95 = q[5]
  )
}
resid_summary <- rbind(
  quantile_row(resid_norm, "resid_norm"),
  quantile_row(atac_norm, "atac_norm"),
  quantile_row(ref_norm, "ref_norm"),
  quantile_row(rel_norm, "rel_norm"),
  quantile_row(cos_sim, "cos_sim")
)
fwrite(resid_summary, file.path(out_dir, "rna_ref_residual_summary.csv"))
cat("\n[Residual distribution]\n")
print(resid_summary)

# =========================
# 4) Plots
# =========================
overlap_long <- data.table(
  set = factor(
    c("fused only", "common", "original only"),
    levels = c("fused only", "common", "original only")
  ),
  count = c(n_f - n_inter, n_inter, n_o - n_inter)
)
p_edges <- ggplot(overlap_long, aes(x = set, y = count, fill = set)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c(
    "fused only" = "#4C78A8",
    "common" = "#54A24B",
    "original only" = "#8C8C8C"
  )) +
  labs(
    title = sprintf("Edge overlap (Jaccard = %.3f)", jaccard),
    x = NULL,
    y = "Number of undirected edges"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")
ggsave(file.path(fig_dir, "Fig_edge_overlap.png"),
       p_edges, width = 6, height = 4.5, dpi = 300)

p_same <- ggplot(per_cell_dt[n_total > 0], aes(x = method, y = same_frac, fill = method)) +
  geom_violin(linewidth = 0.2, alpha = 0.85) +
  geom_boxplot(width = 0.18, outlier.shape = NA, fill = "white", linewidth = 0.2) +
  scale_fill_manual(values = c(original = "#8C8C8C", fused = "#4C78A8")) +
  labs(
    title = "Per-cell same-celltype neighbor fraction",
    x = NULL,
    y = "Fraction of neighbors sharing BioClassification"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")
ggsave(file.path(fig_dir, "Fig_same_type_fraction.png"),
       p_same, width = 5, height = 5, dpi = 300)

per_type_plot <- dcast(per_type, type ~ method, value.var = "mean_same_type_frac")
setorder(per_type_plot, -original)
type_order <- per_type_plot$type
per_type_long <- melt(
  per_type_plot,
  id.vars = "type",
  measure.vars = c("original", "fused"),
  variable.name = "method",
  value.name = "same_frac"
)
per_type_long[, type := factor(type, levels = type_order)]
p_per_type <- ggplot(per_type_long, aes(x = type, y = same_frac, fill = method)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  coord_flip() +
  scale_fill_manual(values = c(original = "#8C8C8C", fused = "#4C78A8")) +
  labs(
    title = "Mean same-celltype neighbor fraction by cell type",
    x = NULL,
    y = "Mean same-celltype fraction"
  ) +
  theme_bw(base_size = 11)
ggsave(file.path(fig_dir, "Fig_same_type_by_celltype.png"),
       p_per_type, width = 8, height = 7, dpi = 300)

p_resid <- ggplot(resid_dt, aes(x = resid_norm)) +
  geom_histogram(bins = 80, fill = "#4C78A8", color = "white", linewidth = 0.15) +
  geom_vline(
    xintercept = median(resid_dt$resid_norm),
    color = "#D62728", linewidth = 0.7
  ) +
  labs(
    title = "RNA_ref - ATAC_latent residual norm",
    subtitle = sprintf(
      "median = %.3f, p25 = %.3f, p75 = %.3f",
      median(resid_dt$resid_norm),
      quantile(resid_dt$resid_norm, 0.25),
      quantile(resid_dt$resid_norm, 0.75)
    ),
    x = "||RNA_ref - ATAC_latent||_2",
    y = "Number of ATAC cells"
  ) +
  theme_bw(base_size = 12)
ggsave(file.path(fig_dir, "Fig_residual_norm.png"),
       p_resid, width = 6, height = 4.5, dpi = 300)

p_rel <- ggplot(resid_dt, aes(x = rel_norm)) +
  geom_histogram(bins = 80, fill = "#F58518", color = "white", linewidth = 0.15) +
  geom_vline(
    xintercept = median(resid_dt$rel_norm),
    color = "#D62728", linewidth = 0.7
  ) +
  labs(
    title = "Relative residual ||RNA_ref - ATAC|| / ||ATAC||",
    subtitle = sprintf(
      "median = %.3f, p25 = %.3f, p75 = %.3f",
      median(resid_dt$rel_norm),
      quantile(resid_dt$rel_norm, 0.25),
      quantile(resid_dt$rel_norm, 0.75)
    ),
    x = "Relative residual",
    y = "Number of ATAC cells"
  ) +
  theme_bw(base_size = 12)
ggsave(file.path(fig_dir, "Fig_residual_relative.png"),
       p_rel, width = 6, height = 4.5, dpi = 300)

# =========================
# 5) Manifest
# =========================
manifest <- list(
  script = "09_knn_overlap_diagnostics.R",
  cutoff = cutoff,
  npcs = npcs,
  k = k,
  fused_knn = fused_knn_path,
  original_knn = original_knn_path,
  shared_embedding = paths$shared_embedding,
  rna_reference = rna_ref_path,
  output_dir = out_dir,
  jaccard = jaccard,
  n_common_cells = length(common_cells),
  created_at = as.character(Sys.time())
)
jsonlite::write_json(
  manifest,
  file.path(out_dir, "diagnostics_manifest.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

cat("\n========== Done ==========\n")
cat("  edge overlap:", file.path(out_dir, "knn_edge_overlap.csv"), "\n")
cat("  same-type   :", file.path(out_dir, "knn_same_type_overall.csv"), "\n")
cat("  residual    :", file.path(out_dir, "rna_ref_residual_summary.csv"), "\n")
cat("  figures     :", fig_dir, "\n")
