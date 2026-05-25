# 04_dual_modality/02_scdrs_montecarlo_zrna.R
#
# §2.3 SCDRS Monte Carlo cell-level z_rna scalar. For each trait, output one
# scalar per cell: z_rna(c) = (s_disease(c) - mean_k s_control_k(c)) / sd_k
# s_control_k(c). The MC dispersion of weighted-average expression on
# (mean,var)-matched control gene sets is what makes this robust to per-cell
# sparsity.
#
# Inputs:
#   <out>/01_data_qc/BMMC_rna_counts.rds       sparse counts (orientation auto-detected)
#   <out>/01_data_qc/BMMC_rna_features.rds     gene_id / gene_symbol
#   <out>/04_dual_modality/Gtrait_<trait>.rds  gene_symbol + w_g per trait
#
# Outputs:
#   <out>/04_dual_modality/zrna_<trait>.rds            named numeric vector (cell_barcode -> z)
#   <out>/04_dual_modality/zrna_summary.csv            per-trait stats
#   <out>/04_dual_modality/figures/zrna_distribution.png  histogram per trait
#
# Red line (§2.10): does NOT load BMMC_cell_meta.rds or any cell_type field.
# Cells are identified only by barcode (counts matrix rownames/colnames).

source("D:/RStudio/work/scavenge_multiomics_fusion_reproduction/reproduction_outputs/patched_scripts/module3_multiome/_bootstrap.R")
required_pkgs(c("data.table", "Matrix", "matrixStats", "ggplot2"))

paths   <- m3_paths()
out_dir <- file.path(paths$out_dir, "04_dual_modality")
fig_dir <- file.path(out_dir, "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# --- locked params --------------------------------------------------------
N_MC          <- 200L                  # SCDRS default
CONTROL_FLOOR <- 100L                  # min non-disease genes per (mean,var) bin
N_BINS_GRID   <- c(20L, 15L, 12L, 10L, 8L, 6L, 5L, 4L)
TRAITS        <- c("MONO", "PLT", "LYM")

# --- 1. load + orient counts ---------------------------------------------

rna_feat <- readRDS(file.path(paths$out_dir, "01_data_qc",
                              "BMMC_rna_features.rds"))
counts <- readRDS(file.path(paths$out_dir, "01_data_qc",
                            "BMMC_rna_counts.rds"))
message("Raw counts: ", nrow(counts), " x ", ncol(counts), " (rows x cols)")
n_genes_feat <- nrow(rna_feat)

# Target orientation: rows = cells, cols = genes (so cell-level ops are colMeans-free).
if (nrow(counts) == n_genes_feat && ncol(counts) != n_genes_feat) {
  message("Detected genes x cells; transposing to cells x genes.")
  counts <- Matrix::t(counts)
} else if (ncol(counts) == n_genes_feat && nrow(counts) != n_genes_feat) {
  message("Detected cells x genes; keeping as is.")
} else {
  stop("Cannot infer counts orientation: dim=",
       paste(dim(counts), collapse = "x"),
       ", n_genes_feat=", n_genes_feat)
}
n_cells <- nrow(counts); n_genes <- ncol(counts)
stopifnot(n_genes == n_genes_feat)

# Align column names with rna_feat$gene_symbol when possible.
if (is.null(colnames(counts)) || all(colnames(counts) == "")) {
  colnames(counts) <- rna_feat$gene_symbol
} else if (!identical(colnames(counts), rna_feat$gene_symbol)) {
  if (all(colnames(counts) %in% rna_feat$gene_id)) {
    map <- setNames(rna_feat$gene_symbol, rna_feat$gene_id)
    colnames(counts) <- unname(map[colnames(counts)])
    message("Remapped counts colnames from gene_id -> gene_symbol.")
  }
}
gene_symbols <- colnames(counts)
cell_barcodes <- rownames(counts)
message("Aligned: ", n_cells, " cells x ", n_genes, " genes")

# --- 2. normalize + log1p -------------------------------------------------

lib_size <- Matrix::rowSums(counts)
if (any(lib_size == 0)) {
  warning(sum(lib_size == 0), " cells have zero RNA counts; their z_rna = NA")
}
scale_med <- stats::median(lib_size[lib_size > 0])
norm_factor <- ifelse(lib_size > 0, scale_med / lib_size, 0)
# Diagonal scaling is sparse-friendly: row-wise multiply preserves zeros.
expr <- Matrix::Diagonal(x = norm_factor) %*% counts
expr <- methods::as(expr, "CsparseMatrix")
expr@x <- log1p(expr@x)
message("log1p-normalized expr: ", nrow(expr), " x ", ncol(expr),
        " | nnz=", length(expr@x))

# --- 3. gene mean / var on log1p data ------------------------------------

# Sparse-friendly: mean = sum/n, var = (sumSq - n*mean^2) / (n-1)
n_c <- nrow(expr)
gene_sum   <- Matrix::colSums(expr)
gene_sumSq <- Matrix::colSums(expr * expr)
gene_mean  <- gene_sum / n_c
gene_var   <- (gene_sumSq - n_c * gene_mean^2) / (n_c - 1L)
gene_var[gene_var < 0] <- 0   # numerical guard

expressed_mask <- gene_mean > 0   # genes that ever fire
message("Genes expressed (mean > 0): ", sum(expressed_mask), " / ", n_genes)

# --- 4. 2D rank-based bins, auto-adjust until min pool >= CONTROL_FLOOR --

mean_rank <- rank(gene_mean, ties.method = "average") / n_genes
var_rank  <- rank(gene_var,  ties.method = "average") / n_genes

make_bin_id <- function(N_bins) {
  bm <- pmin(N_bins, pmax(1L, ceiling(mean_rank * N_bins)))
  bv <- pmin(N_bins, pmax(1L, ceiling(var_rank  * N_bins)))
  (bm - 1L) * N_bins + bv
}

# Disease symbol union across traits, used to mask out non-disease genes
all_disease <- unique(unlist(lapply(TRAITS, function(tr) {
  pick <- readRDS(file.path(out_dir, paste0("Gtrait_", tr, ".rds")))
  pick$gene_symbol
})))

# auto-pick N_bins: smallest N_bins (= largest bins) is permissive; we want
# the largest N_bins that still leaves CONTROL_FLOOR non-disease genes per
# bin that contains a disease gene.
chosen_N_bins <- NA_integer_
for (N_bins in N_BINS_GRID) {
  bin_id <- make_bin_id(N_bins)
  is_disease <- gene_symbols %in% all_disease
  # only check bins that hold disease genes AND filter to expressed genes
  disease_bins <- unique(bin_id[is_disease & expressed_mask])
  pool_sizes <- vapply(disease_bins, function(b) {
    sum(bin_id == b & !is_disease & expressed_mask)
  }, integer(1))
  message(sprintf(
    "  N_bins=%d -> %d disease-bins | min pool=%d | median pool=%d",
    N_bins, length(disease_bins), min(pool_sizes),
    as.integer(stats::median(pool_sizes))
  ))
  if (length(pool_sizes) > 0 && min(pool_sizes) >= CONTROL_FLOOR) {
    chosen_N_bins <- N_bins
    break
  }
}
if (is.na(chosen_N_bins)) {
  chosen_N_bins <- tail(N_BINS_GRID, 1L)
  warning("Even N_bins=", chosen_N_bins,
          " cannot satisfy CONTROL_FLOOR; using it anyway.")
}
N_BINS <- chosen_N_bins
bin_id <- make_bin_id(N_BINS)
message("Chosen N_bins = ", N_BINS)

# --- 5. SCDRS Monte Carlo per trait --------------------------------------

# Helper: efficient rowSds on a dense matrix (cells x mc).
row_sds <- function(M) matrixStats::rowSds(M)

zrna_summary <- list()
hist_rows    <- list()

for (tr in TRAITS) {
  message("\n[Trait] ", tr)
  pick <- readRDS(file.path(out_dir, paste0("Gtrait_", tr, ".rds")))
  disease_sym <- pick$gene_symbol
  w_g         <- pick$w_g
  # match to expr columns
  d_idx <- match(disease_sym, gene_symbols)
  ok    <- !is.na(d_idx) & expressed_mask[d_idx]
  if (sum(!ok) > 0) {
    message("  dropping ", sum(!ok), " disease genes (NA or zero-mean)")
  }
  d_idx <- d_idx[ok]; w_g <- w_g[ok]
  n_d   <- length(d_idx)
  if (n_d == 0) {
    warning("Trait ", tr, ": no usable disease genes; skip"); next
  }

  # control pools (non-disease + expressed in same bin)
  is_disease_tr <- gene_symbols %in% disease_sym
  pool_list <- lapply(d_idx, function(i) {
    which(bin_id == bin_id[i] & !is_disease_tr & expressed_mask)
  })
  pool_sizes <- vapply(pool_list, length, integer(1))
  message(sprintf(
    "  n_disease=%d | pool size: min=%d median=%d max=%d",
    n_d, min(pool_sizes),
    as.integer(stats::median(pool_sizes)), max(pool_sizes)
  ))
  if (min(pool_sizes) == 0L) {
    bad <- which(pool_sizes == 0)
    warning("Trait ", tr, ": ", length(bad),
            " disease genes have empty control pool; dropping them.")
    keep_d <- pool_sizes > 0
    d_idx  <- d_idx[keep_d]; w_g <- w_g[keep_d]
    pool_list <- pool_list[keep_d]; n_d <- length(d_idx)
    if (n_d == 0) { warning("All disease genes dropped; skip"); next }
  }

  # MC sampling: n_mc x n_d matrix of column indices into expr
  set.seed(202604L + which(TRAITS == tr))
  ctrl_mat <- vapply(seq_len(n_d), function(j) {
    sample(pool_list[[j]], N_MC, replace = TRUE)
  }, integer(N_MC))   # N_MC x n_d

  # s_disease vector (length n_cells)
  w_norm <- w_g / sum(w_g)   # absorb sum into weights, makes everything matrix-vec
  expr_d <- expr[, d_idx, drop = FALSE]
  s_disease <- as.numeric(expr_d %*% w_norm)

  # s_control matrix (n_cells x N_MC). For iter k, the k-th column gathers
  # the k-th sampled control per disease gene.
  s_ctrl <- matrix(0.0, nrow = n_cells, ncol = N_MC)
  # vectorize per iteration: extract n_d cols of expr, multiply by w_norm
  for (k in seq_len(N_MC)) {
    cols_k <- ctrl_mat[k, ]
    expr_k <- expr[, cols_k, drop = FALSE]
    s_ctrl[, k] <- as.numeric(expr_k %*% w_norm)
  }

  m_k  <- rowMeans(s_ctrl)
  sd_k <- row_sds(s_ctrl)
  z    <- (s_disease - m_k) / sd_k
  z[!is.finite(z)] <- NA_real_

  # save
  z_named <- z; names(z_named) <- cell_barcodes
  saveRDS(z_named, file.path(out_dir, paste0("zrna_", tr, ".rds")))
  message(sprintf(
    "  z_rna: mean=%.3f sd=%.3f q95=%.3f q99=%.3f n_NA=%d  saved",
    mean(z, na.rm = TRUE), stats::sd(z, na.rm = TRUE),
    as.numeric(stats::quantile(z, 0.95, na.rm = TRUE)),
    as.numeric(stats::quantile(z, 0.99, na.rm = TRUE)),
    sum(is.na(z))
  ))

  zrna_summary[[tr]] <- data.table::data.table(
    trait        = tr,
    n_disease    = n_d,
    n_cells      = n_cells,
    n_NA         = sum(is.na(z)),
    z_mean       = mean(z, na.rm = TRUE),
    z_sd         = stats::sd(z, na.rm = TRUE),
    z_q50        = as.numeric(stats::quantile(z, 0.50, na.rm = TRUE)),
    z_q95        = as.numeric(stats::quantile(z, 0.95, na.rm = TRUE)),
    z_q99        = as.numeric(stats::quantile(z, 0.99, na.rm = TRUE)),
    z_max        = max(z, na.rm = TRUE),
    pool_min     = min(pool_sizes),
    pool_median  = as.integer(stats::median(pool_sizes))
  )

  hist_rows[[tr]] <- data.table::data.table(trait = tr, z = z[is.finite(z)])
}

summary_dt <- data.table::rbindlist(zrna_summary)
data.table::fwrite(summary_dt, file.path(out_dir, "zrna_summary.csv"))
print(summary_dt)

hist_dt <- data.table::rbindlist(hist_rows)
p <- ggplot2::ggplot(hist_dt, ggplot2::aes(x = z)) +
  ggplot2::geom_histogram(bins = 80, fill = "#2c7fb8", alpha = 0.85) +
  ggplot2::geom_vline(xintercept = 0, color = "grey40", linetype = 2) +
  ggplot2::facet_wrap(~ trait, scales = "free_y") +
  ggplot2::labs(
    title = "SCDRS Monte Carlo z_rna distribution per trait",
    x = "z_rna", y = "Cells"
  ) +
  ggplot2::theme_bw()
ggplot2::ggsave(
  file.path(fig_dir, "zrna_distribution.png"),
  p, width = 10, height = 3.2, dpi = 150
)

message("\nz_rna built. Outputs in: ", out_dir)
