# module4 / Phase 0.2: BMMC trait-aware mapping (sanity check)
# ------------------------------------------------------------------
# 目的: 在 G_trait 子集空间, 给每个 ATAC cell 找 top-5 cosine 最近的 RNA cell
#       (洗牌后), 输出 mapping. 用 BMMC 已知 ground truth 在 Phase 0.3 评估
#       cell-type 一致率.
#
# 输入:
#   00_sanity_bmmc/output/BMMC_gene_activity.rds   (gene × cell, Phase 0.1.5)
#   00_sanity_bmmc/output/bmmc_rna_shuffled.rds    (gene × cell, Phase 0.1)
#   D:/SCAVENGE-reproducibility-main/BMMC/01_data_qc/BMMC_cell_meta.rds
#     (ATAC barcode 用; ATAC 端列序来自 GA, 来自 SE_atac_gcbiased.rds)
#   D:/RStudio/work/module3_multiome/04_dual_modality/01_scdrs_python/bmmc_3traits.gs
#
# 输出:
#   00_sanity_bmmc/output/bmmc_mapping_<trait>.rds  (list 含 5 个 field)
#   00_sanity_bmmc/output/bmmc_mapping_summary.csv  (距离 / G_trait 使用情况)
#
# 算法 (README §8.2 + WORKFLOW §0.2):
#   ATAC 表征: GA[G_trait_∩, ] -> 行 z-score (gene 维度, across cells)
#              -> 行 * w_g -> 转置 cell × gene -> 行 L2 normalize (每 cell 单位)
#   RNA 同上但用 log1p(library-normalized counts) 代 GA
#   cosine top-5 via Euclidean (L2 normalize 后等价), softmax weight per-cell
# ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
})
# 注: 不用 RANN::nn2. 在 |G_trait| ~700 维上 kd-tree 退化成接近暴力扫描
# (curse of dimensionality), 且 RANN 单线程, 实测在 69k×678 + k=5 上跑不完.
# 改用 batched 稠密 GEMM + partial sort, 利用 BLAS 多线程, ~5 min/trait.

# ---- paths ----
GA_PATH    <- "D:/RStudio/work/module4_traitaware_pseudopair/00_sanity_bmmc/output/BMMC_gene_activity.rds"
RNA_PATH   <- "D:/RStudio/work/module4_traitaware_pseudopair/00_sanity_bmmc/output/bmmc_rna_shuffled.rds"
GS_FILE    <- "D:/RStudio/work/module3_multiome/04_dual_modality/01_scdrs_python/bmmc_3traits.gs"
OUT_DIR    <- "D:/RStudio/work/module4_traitaware_pseudopair/00_sanity_bmmc/output"

# ---- experiment switch ----
# USE_WG = TRUE  : 原方案, 行向乘 w_g (MAGMA z-score)
# USE_WG = FALSE : 对照实验 (A0), 不乘 w_g, gene 等权
USE_WG         <- FALSE
OUTPUT_SUFFIX  <- if (USE_WG) "" else "_noweight"
cat(">>> Phase 0.2 run mode: USE_WG = ", USE_WG,
    "  ->  output suffix = '", OUTPUT_SUFFIX, "'\n", sep = "")

TRAITS_ALL <- c(
  MONO = "UKB_460K.blood_MONOCYTE_COUNT",
  PLT  = "UKB_460K.blood_PLATELET_COUNT",
  LYM  = "UKB_460K.blood_LYMPHOCYTE_COUNT"
)
# 选要跑哪些 trait. 用于跨 session 续跑 (避免重复跑已完成的).
# 默认全跑; 单 trait 续跑示例: TRAITS_TO_RUN <- c("LYM")
if (!exists("TRAITS_TO_RUN")) TRAITS_TO_RUN <- names(TRAITS_ALL)
TRAITS <- TRAITS_ALL[TRAITS_TO_RUN]
cat(">>> will run traits: ", paste(names(TRAITS), collapse = ", "), "\n", sep = "")

K          <- 5L
EPS        <- 1e-8

# ------------------------------------------------------------------
# helpers
# ------------------------------------------------------------------

parse_gs_file <- function(gs_path) {
  ln <- readLines(gs_path)
  out <- list()
  for (i in 2:length(ln)) {
    parts <- strsplit(ln[i], "\t", fixed = TRUE)[[1]]
    if (length(parts) < 2) next
    trait <- parts[1]
    pairs <- strsplit(strsplit(parts[2], ",", fixed = TRUE)[[1]],
                       ":", fixed = TRUE)
    out[[trait]] <- list(
      genes   = vapply(pairs, `[`, character(1), 1),
      weights = as.numeric(vapply(pairs,
                                  function(x) if (length(x) >= 2) x[2] else "1",
                                  character(1)))
    )
  }
  out
}

# Build cell × |G| feature block:
#   1. gene-level z-score (each gene zero-mean unit-var across cells)
#   2. gene-level w_g multiplication (broadcast across cells)
#   3. transpose -> cell × gene
#   4. cell-level L2 normalize (each row unit-L2)
# Modality-agnostic: same function for GA (ATAC) and log1p_norm (RNA).
build_trait_block <- function(mat_g_by_c, w_g, use_wg = TRUE) {
  # mat_g_by_c : numeric / dgCMatrix, gene × cell, rows aligned to w_g order
  # w_g        : numeric vector, length == nrow(mat_g_by_c)
  # use_wg     : if FALSE, skip step 2 (gene equally weighted)

  # 1. gene-level z-score across cells (dense, |G| rows manageable)
  mu <- Matrix::rowMeans(mat_g_by_c)
  # rowSds via rowMeans (using formula Var = E[X^2] - E[X]^2; needs squared mat)
  m2 <- Matrix::rowMeans(mat_g_by_c * mat_g_by_c)
  sg <- sqrt(pmax(m2 - mu * mu, 0))
  sg[sg < EPS] <- EPS
  # broadcast (mat - mu) / sg row-wise
  Z <- as.matrix(mat_g_by_c)                # dense, |G|×N (|G|~1k OK)
  Z <- (Z - mu) / sg

  # 2. (optional) gene-level w_g weighting (row-wise multiply)
  if (isTRUE(use_wg)) {
    Z <- Z * w_g
  }

  # 3. transpose -> cell × gene
  X <- t(Z)                                  # N × |G|

  # 4. cell-level L2 normalize
  norms <- sqrt(rowSums(X * X))
  norms[norms < EPS] <- EPS
  X <- X / norms

  X
}

# Batched top-k cosine search via dense GEMM (BLAS multi-threaded).
# Both X_query and X_ref must be row-L2-normalized (then cosine = inner product,
# Euclidean d² = 2 - 2·cos).
# Returns same shape as RANN::nn2 result: list(nn.idx, nn.dists).
batched_top_k_cosine <- function(X_query, X_ref, k, batch = 1024L,
                                  verbose = TRUE) {
  n_q <- nrow(X_query); n_r <- nrow(X_ref); d <- ncol(X_query)
  stopifnot(ncol(X_ref) == d)
  X_ref_t  <- t(X_ref)
  idx_all  <- matrix(0L, n_q, k)
  dist_all <- matrix(0,  n_q, k)
  n_batch  <- as.integer(ceiling(n_q / batch))
  t0 <- Sys.time()
  for (bi in seq_len(n_batch)) {
    b_start <- (bi - 1L) * batch + 1L
    b_end   <- min(bi * batch, n_q)
    Xq      <- X_query[b_start:b_end, , drop = FALSE]
    sim     <- Xq %*% X_ref_t                       # nb × n_r, cosine sim
    sim     <- as.matrix(sim)                       # ensure dense for max.col
    nb      <- nrow(sim)

    # top-k via repeated max.col + mask-to-(-Inf). max.col is vectorized C.
    idx_b  <- matrix(0L, nb, k)
    dist_b <- matrix(0,  nb, k)
    row_ix <- seq_len(nb)
    for (i in seq_len(k)) {
      mx <- max.col(sim, ties.method = "first")      # nb integer indices
      idx_b[, i]  <- mx
      vals        <- sim[cbind(row_ix, mx)]
      dist_b[, i] <- sqrt(pmax(2 - 2 * vals, 0))
      sim[cbind(row_ix, mx)] <- -Inf                 # mask for next iter
    }
    idx_all[b_start:b_end, ]  <- idx_b
    dist_all[b_start:b_end, ] <- dist_b

    if (verbose && (bi %% 5L == 0L || bi == n_batch)) {
      el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      cat(sprintf("       batch %d/%d done (cells %d/%d, %.1fs elapsed)\n",
                  bi, n_batch, b_end, n_q, el))
      flush.console()
    }
    # free this batch's sim matrix before next allocation
    rm(sim); invisible(gc(verbose = FALSE))
  }
  list(nn.idx = idx_all, nn.dists = dist_all)
}

# softmax with per-cell tau = median(d_top_k)
softmax_per_cell <- function(d_mat) {
  # d_mat : N × k Euclidean distances
  tau <- apply(d_mat, 1, function(x) {
    m <- stats::median(x)
    if (m < EPS) EPS else m
  })
  S <- exp(-d_mat / tau)
  S <- S / rowSums(S)
  list(weight = S, tau = tau)
}

# ------------------------------------------------------------------
# load inputs
# ------------------------------------------------------------------

t_total <- Sys.time()
cat("[1/4] Loading inputs ...\n")

ga          <- readRDS(GA_PATH)            # gene × ATAC_cell, dgCMatrix
rna_shuf    <- readRDS(RNA_PATH)           # gene × RNA_cell (shuffled), dgCMatrix

cat("       ga       : ", paste(dim(ga),       collapse = " x "),
    "  (gene x ATAC_cell)\n", sep = "")
cat("       rna_shuf : ", paste(dim(rna_shuf), collapse = " x "),
    "  (gene x RNA_cell)\n", sep = "")

# Check ATAC and RNA cell counts match (BMMC paired, should be aligned via 0.1)
stopifnot(ncol(ga) == ncol(rna_shuf))
atac_cells <- colnames(ga)
rna_cells  <- colnames(rna_shuf)
cat("       atac_cells head: ", paste(head(atac_cells, 2), collapse = ", "), "\n")
cat("       rna_cells head : ", paste(head(rna_cells,  2), collapse = ", "), "\n")

# RNA library size for log1p normalization (compute once, reuse across traits)
rna_lib    <- Matrix::colSums(rna_shuf)
rna_lib[rna_lib < 1] <- 1
cat("       rna_lib median : ", median(rna_lib), "\n")

# Parse .gs
gs <- parse_gs_file(GS_FILE)
cat("       parsed .gs traits: ", paste(names(gs), collapse = ", "), "\n")

# ------------------------------------------------------------------
# per-trait mapping
# ------------------------------------------------------------------

summary_rows <- list()

for (tname in names(TRAITS)) {
  trait_id <- TRAITS[[tname]]
  cat("\n[2/4] === Trait ", tname, " (", trait_id, ") ===\n", sep = "")

  g_full <- gs[[trait_id]]$genes
  w_full <- gs[[trait_id]]$weights
  stopifnot(length(g_full) == length(w_full))
  cat("       G_trait full size : ", length(g_full), "\n")

  # intersect with ATAC GA gene set & RNA gene set
  g_atac <- intersect(g_full, rownames(ga))
  g_rna  <- intersect(g_full, rownames(rna_shuf))
  g_use  <- intersect(g_atac, g_rna)        # final shared G_trait
  cat("       G ∩ GA            : ", length(g_atac), "\n")
  cat("       G ∩ RNA           : ", length(g_rna),  "\n")
  cat("       G shared (use)    : ", length(g_use),  "\n")
  stopifnot(length(g_use) >= 50)

  # align w_g to g_use order
  idx_in_full <- match(g_use, g_full)
  w_g <- w_full[idx_in_full]

  # ---- ATAC side: GA on g_use ----
  cat("       building ATAC block (use_wg = ", USE_WG, ") ...\n", sep = "")
  ga_sub <- ga[g_use, , drop = FALSE]       # |G_use| × N (dgCMatrix)
  X_atac <- build_trait_block(ga_sub, w_g, use_wg = USE_WG)

  # ---- RNA side: log1p(library-normalized) on g_use ----
  cat("       building RNA block (use_wg = ", USE_WG, ") ...\n", sep = "")
  rna_sub <- rna_shuf[g_use, , drop = FALSE]
  # library size normalize to 1e4, log1p
  rna_sub_norm <- Matrix::t(Matrix::t(rna_sub) / rna_lib) * 1e4
  rna_sub_log  <- log1p(rna_sub_norm)
  X_rna <- build_trait_block(rna_sub_log, w_g, use_wg = USE_WG)

  cat("       X_atac dim: ", paste(dim(X_atac), collapse = " x "), "\n")
  cat("       X_rna  dim: ", paste(dim(X_rna),  collapse = " x "), "\n")

  # ---- cosine top-K via batched GEMM (BLAS multi-threaded) ----
  cat("       running batched GEMM top-", K,
      " (batch=1024 ATAC cells/iter) ...\n", sep = "")
  t0 <- Sys.time()
  nn <- batched_top_k_cosine(X_query = X_atac, X_ref = X_rna,
                              k = K, batch = 1024L, verbose = TRUE)
  dt_nn <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat("       knn elapsed: ", round(dt_nn, 1), "s\n")

  nbr_idx  <- nn$nn.idx          # N_atac × K, integer
  nbr_dist <- nn$nn.dists        # N_atac × K, Euclidean

  # softmax weights (per-cell tau)
  sm <- softmax_per_cell(nbr_dist)

  # neighbor barcode lookup
  nbr_bc <- matrix(rna_cells[as.integer(nbr_idx)],
                   nrow = nrow(nbr_idx), ncol = K)

  # ---- save ----
  mapping <- list(
    trait                 = tname,
    trait_id              = trait_id,
    k                     = K,
    g_trait_used          = g_use,
    w_g                   = w_g,
    n_g_trait_used        = length(g_use),
    atac_barcodes         = atac_cells,
    rna_neighbor_idx      = nbr_idx,        # N_atac × K
    rna_neighbor_barcode  = nbr_bc,         # N_atac × K
    neighbor_distance     = nbr_dist,       # N_atac × K, Euclidean
    neighbor_weight       = sm$weight,      # N_atac × K, softmax normalized
    tau_per_cell          = sm$tau          # N_atac
  )

  out_file <- file.path(OUT_DIR,
                        paste0("bmmc_mapping_", tname, OUTPUT_SUFFIX, ".rds"))
  saveRDS(mapping, out_file)
  cat("       saved ", basename(out_file), "\n", sep = "")

  # ---- per-trait stats for summary ----
  d_med   <- median(nbr_dist)
  d_q05   <- as.numeric(quantile(nbr_dist, 0.05))
  d_q95   <- as.numeric(quantile(nbr_dist, 0.95))
  d_top1  <- nbr_dist[, 1]
  iso_n   <- sum(d_top1 > quantile(d_top1, 0.95) * 1.1)  # heuristic
  wsum    <- rowSums(sm$weight)
  cat(sprintf(
    "       distance median: %.3f  [q05=%.3f, q95=%.3f]; tau median: %.3f\n",
    d_med, d_q05, d_q95, median(sm$tau)
  ))
  cat(sprintf(
    "       softmax weight rowsum range: [%.6f, %.6f]\n",
    min(wsum), max(wsum)
  ))

  summary_rows[[tname]] <- data.frame(
    trait               = tname,
    n_g_trait_full      = length(g_full),
    n_g_trait_shared    = length(g_use),
    n_atac_cells        = nrow(X_atac),
    n_rna_cells         = nrow(X_rna),
    k                   = K,
    dist_median         = d_med,
    dist_q05            = d_q05,
    dist_q95            = d_q95,
    tau_median          = median(sm$tau),
    nn2_seconds         = dt_nn,
    stringsAsFactors    = FALSE
  )

  # free large objects before the next trait iteration
  rm(X_atac, X_rna, ga_sub, rna_sub, rna_sub_norm, rna_sub_log,
     nn, nbr_idx, nbr_dist, nbr_bc, sm, mapping)
  invisible(gc(verbose = FALSE))
  cat("       memory freed; gc done\n")
  flush.console()
}

# ------------------------------------------------------------------
# summary CSV
# ------------------------------------------------------------------
cat("\n[3/4] Writing summary CSV ...\n")
summary_df <- do.call(rbind, summary_rows)
summary_df$use_wg <- USE_WG
write.csv(summary_df,
          file.path(OUT_DIR,
                    paste0("bmmc_mapping_summary", OUTPUT_SUFFIX, ".csv")),
          row.names = FALSE)

cat("\n=== Phase 0.2 mapping DONE ===\n")
print(summary_df)
cat("\nTotal elapsed: ",
    round(as.numeric(difftime(Sys.time(), t_total, units = "secs")), 1),
    "s\n", sep = "")
cat("Outputs in: ", OUT_DIR, "\n", sep = "")
