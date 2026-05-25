# 04_dual_modality/03_atac_rna_weighted_fusion.R
#
# §2.4 ATAC-major / RNA-minor 非对称加权融合表征. ATAC LSI 30 维直接复用
# baseline 产物; RNA 端只取 G_trait 子集 log1p 表达, 单独 PCA 到 K_rna 维
# (auto: 累计方差 >= 80% 或 kneedle elbow 取较小, bracket [10, 30]).
# 两块各自 z-score 后用 sqrt(weight) 非对称 cbind, 默认 0.7/0.3.
#
# 输出三个独立 block (atac_z + rna_z_<trait>) 而不是直接出最终 fused 矩阵,
# 因为 §2.5 会按 isolated 细胞比例反复调权重 (0.7/0.3 起步, 上限 0.5/0.5);
# 把 cbind 留在 §2.5 做, §2.4 只负责把可重用的 block 算好.
#
# Inputs:
#   <out>/03_scavenge_baseline/atac_lsi.rds      n_cells x 30 (cells x dims)
#   <out>/01_data_qc/BMMC_rna_counts.rds         sparse counts
#   <out>/01_data_qc/BMMC_rna_features.rds
#   <out>/04_dual_modality/Gtrait_<trait>.rds    per trait gene list
#
# Outputs:
#   <out>/04_dual_modality/atac_z_block.rds      n_cells x 30 z-scored LSI
#   <out>/04_dual_modality/rna_z_block_<tr>.rds  n_cells x K_rna z-scored PC scores
#   <out>/04_dual_modality/fusion_default_<tr>.rds  reference 0.7/0.3 fusion (n_cells x (30+K_rna))
#   <out>/04_dual_modality/fusion_summary.csv
#   <out>/04_dual_modality/figures/rna_pca_var_explained.png
#
# Red line: NO cell_type / cell_meta in this script.

source("D:/RStudio/work/scavenge_multiomics_fusion_reproduction/reproduction_outputs/patched_scripts/module3_multiome/_bootstrap.R")
required_pkgs(c("data.table", "Matrix", "matrixStats", "irlba", "ggplot2"))

paths   <- m3_paths()
out_dir <- file.path(paths$out_dir, "04_dual_modality")
fig_dir <- file.path(out_dir, "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# --- locked params --------------------------------------------------------
K_BRACKET   <- c(10L, 30L)
CUMVAR_TGT  <- 0.80
ATAC_W_DEF  <- 0.7
RNA_W_DEF   <- 0.3
TRAITS      <- c("MONO", "PLT", "LYM")

# --- 1. ATAC LSI block (z-scored, shared across traits) -------------------

atac_lsi <- readRDS(file.path(paths$out_dir, "03_scavenge_baseline",
                              "atac_lsi.rds"))
message("ATAC LSI: ", nrow(atac_lsi), " x ", ncol(atac_lsi))
if (is.null(rownames(atac_lsi))) {
  stop("atac_lsi has no rownames; cannot align with RNA matrix.")
}

# --- 2. RNA counts (cells x genes) + log1p --------------------------------

rna_feat <- readRDS(file.path(paths$out_dir, "01_data_qc",
                              "BMMC_rna_features.rds"))
counts   <- readRDS(file.path(paths$out_dir, "01_data_qc",
                              "BMMC_rna_counts.rds"))
n_genes_feat <- nrow(rna_feat)
if (nrow(counts) == n_genes_feat && ncol(counts) != n_genes_feat) {
  counts <- Matrix::t(counts)
}
stopifnot(ncol(counts) == n_genes_feat)
if (is.null(colnames(counts)) || all(colnames(counts) == "")) {
  colnames(counts) <- rna_feat$gene_symbol
} else if (!identical(colnames(counts), rna_feat$gene_symbol) &&
           all(colnames(counts) %in% rna_feat$gene_id)) {
  map <- setNames(rna_feat$gene_symbol, rna_feat$gene_id)
  colnames(counts) <- unname(map[colnames(counts)])
}
message("RNA counts (cells x genes): ", nrow(counts), " x ", ncol(counts))

# --- align ATAC and RNA on barcode order ---------------------------------
common_cells <- intersect(rownames(atac_lsi), rownames(counts))
if (length(common_cells) < 0.99 * nrow(atac_lsi)) {
  warning("Only ", length(common_cells), " / ", nrow(atac_lsi),
          " cells match between ATAC LSI and RNA counts.")
}
atac_lsi <- atac_lsi[common_cells, , drop = FALSE]
counts   <- counts[common_cells, , drop = FALSE]
n_cells  <- length(common_cells)
message("Aligned cell count: ", n_cells)

# log1p normalize (locked, same recipe as §2.3)
lib_size    <- Matrix::rowSums(counts)
scale_med   <- stats::median(lib_size[lib_size > 0])
norm_factor <- ifelse(lib_size > 0, scale_med / lib_size, 0)
expr <- Matrix::Diagonal(x = norm_factor) %*% counts
expr <- methods::as(expr, "CsparseMatrix")
expr@x <- log1p(expr@x)
gene_symbols <- colnames(expr)
message("expr (log1p): ", nrow(expr), " x ", ncol(expr),
        " | nnz=", length(expr@x))

# --- 3. z-score ATAC LSI per column --------------------------------------

# matrix already centered-ish from LSI; z-score makes the column scale comparable
atac_mu <- matrixStats::colMeans2(atac_lsi)
atac_sd <- matrixStats::colSds(atac_lsi)
atac_sd[atac_sd == 0] <- 1
atac_z <- sweep(atac_lsi, 2L, atac_mu, "-")
atac_z <- sweep(atac_z,   2L, atac_sd, "/")
saveRDS(atac_z, file.path(out_dir, "atac_z_block.rds"))
message("atac_z_block: ", nrow(atac_z), " x ", ncol(atac_z),
        " | per-col sd in [", round(min(matrixStats::colSds(atac_z)), 3),
        ", ", round(max(matrixStats::colSds(atac_z)), 3), "]")

# --- 4. per-trait RNA PCA + K_rna auto + z-score --------------------------

# Kneedle on a decreasing variance curve (same algorithm as §2.2 G_trait):
# argmax(1 - x - y) where x is normalized rank, y is normalized magnitude.
kneedle_pick <- function(v, bracket) {
  n <- length(v)
  if (n < bracket[1]) return(n)
  s <- sort(v, decreasing = TRUE)
  x <- (seq_len(n) - 1) / (n - 1)
  s_min <- min(s); s_max <- max(s)
  y <- if (s_max > s_min) (s - s_min) / (s_max - s_min) else rep(1, n)
  knee <- which.max(1 - x - y)
  max(bracket[1], min(bracket[2], knee))
}

summary_rows <- list()
var_curve    <- list()

for (tr in TRAITS) {
  message("\n[Trait] ", tr)
  pick <- readRDS(file.path(out_dir, paste0("Gtrait_", tr, ".rds")))
  d_idx <- match(pick$gene_symbol, gene_symbols)
  d_idx <- d_idx[!is.na(d_idx)]
  if (length(d_idx) == 0) {
    warning("Trait ", tr, ": no overlap with RNA features; skip"); next
  }

  rna_sub <- expr[, d_idx, drop = FALSE]
  # densify for PCA (rna_sub has ~100 cols, n_cells rows; ~55 MB)
  rna_dense <- as.matrix(rna_sub)
  rna_mu    <- matrixStats::colMeans2(rna_dense)
  rna_centered <- sweep(rna_dense, 2L, rna_mu, "-")

  k_max <- min(K_BRACKET[2], ncol(rna_centered) - 1L)
  pca <- irlba::irlba(rna_centered, nv = k_max, work = k_max + 8L)
  scores <- pca$u %*% diag(pca$d)        # n_cells x k_max PC scores
  total_var <- sum(rna_centered^2) / (nrow(rna_centered) - 1L)
  pc_var    <- (pca$d^2) / (nrow(rna_centered) - 1L)
  cumvar    <- cumsum(pc_var) / total_var

  K_cumvar_idx <- as.integer(which(cumvar >= CUMVAR_TGT)[1])
  K_kneedle    <- kneedle_pick(pc_var, K_BRACKET)
  if (is.na(K_cumvar_idx)) {
    # cumvar never reached target within k_max; fall back to upper bracket.
    # kneedle alone is unsafe: for cell-type-structured PCA on ~100 disease
    # genes PC1 dominates so kneedle reports knee=1 (clamped to 10), losing
    # sub-cell-type structure that downstream KNN needs.
    K_pre               <- K_BRACKET[2]
    K_cumvar_summary    <- NA_integer_
    sel_rule            <- "cumvar_unreached_use_upper"
  } else {
    K_pre               <- min(K_cumvar_idx, K_kneedle)
    K_cumvar_summary    <- K_cumvar_idx
    sel_rule            <- "min_cumvar_kneedle"
  }
  K_rna <- max(K_BRACKET[1], min(K_BRACKET[2], K_pre))

  scores_K <- scores[, seq_len(K_rna), drop = FALSE]
  rna_mu_K <- matrixStats::colMeans2(scores_K)
  rna_sd_K <- matrixStats::colSds(scores_K)
  rna_sd_K[rna_sd_K == 0] <- 1
  rna_z <- sweep(scores_K, 2L, rna_mu_K, "-")
  rna_z <- sweep(rna_z,    2L, rna_sd_K, "/")
  rownames(rna_z) <- common_cells

  saveRDS(rna_z, file.path(out_dir, paste0("rna_z_block_", tr, ".rds")))

  # default reference fusion (sqrt(w) → variance scales as w)
  fused_def <- cbind(
    atac_z * sqrt(ATAC_W_DEF),
    rna_z  * sqrt(RNA_W_DEF)
  )
  rownames(fused_def) <- common_cells
  saveRDS(fused_def, file.path(out_dir, paste0("fusion_default_", tr, ".rds")))

  message(sprintf(
    "  K_cumvar=%s (cumvar>=%.2f), K_kneedle=%d, rule=%s, K_rna picked=%d | cumvar@K=%.3f",
    if (is.na(K_cumvar_summary)) "NA" else as.character(K_cumvar_summary),
    CUMVAR_TGT, K_kneedle, sel_rule, K_rna, cumvar[K_rna]
  ))

  summary_rows[[tr]] <- data.table::data.table(
    trait        = tr,
    n_disease    = length(d_idx),
    K_cumvar     = K_cumvar_summary,
    K_kneedle    = K_kneedle,
    sel_rule     = sel_rule,
    K_rna_picked = K_rna,
    cumvar_at_K  = cumvar[K_rna],
    fused_dim    = ncol(fused_def),
    atac_w       = ATAC_W_DEF,
    rna_w        = RNA_W_DEF
  )

  var_curve[[tr]] <- data.table::data.table(
    trait      = tr,
    pc_idx     = seq_along(pc_var),
    var_share  = pc_var / total_var,
    cumvar     = cumvar,
    knee_at    = K_rna
  )
}

summary_dt <- data.table::rbindlist(summary_rows)
data.table::fwrite(summary_dt, file.path(out_dir, "fusion_summary.csv"))
print(summary_dt)

curve_dt <- data.table::rbindlist(var_curve)
p <- ggplot2::ggplot(curve_dt,
                     ggplot2::aes(x = pc_idx, y = cumvar)) +
  ggplot2::geom_line() +
  ggplot2::geom_point(size = 0.8) +
  ggplot2::geom_vline(
    data = unique(curve_dt[, .(trait, knee_at)]),
    ggplot2::aes(xintercept = knee_at), color = "red", linetype = 2
  ) +
  ggplot2::geom_hline(yintercept = CUMVAR_TGT,
                      color = "grey40", linetype = 3) +
  ggplot2::facet_wrap(~ trait) +
  ggplot2::labs(
    title = "RNA G_trait subset PCA — cumulative variance explained",
    x = "PC index", y = "Cumulative variance share"
  ) +
  ggplot2::theme_bw()
ggplot2::ggsave(
  file.path(fig_dir, "rna_pca_var_explained.png"),
  p, width = 10, height = 3.2, dpi = 150
)

message("\nFusion blocks built. Outputs in: ", out_dir)
