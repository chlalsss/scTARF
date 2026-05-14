# 04_dual_modality/04_build_weighted_knn.R
#
# §2.5 带权 mutual KNN 图. 起步 k=30, ATAC/RNA 权重 0.7/0.3 (§2.4 默认).
# 每细胞出边相似度 s_ij = exp(-d_ij / median_d_per_cell_i); mutual filter
# (i→j 和 j→i 都有才保留, 相似度取均值); 行归一化让每细胞出边和=1
# (RWR 转移概率直接对应边权).
#
# Auto-tune (§2.5 step 2): isolated 比例 > 1% 时
#   先回 §2.4 调权重: atac_w -= 0.05, rna_w += 0.05, 上限 0.5/0.5
#   权重已到上限仍 isolated > 1% → k ← ceiling(k × 1.5), 上限 K_HARD_CAP
#
# Inputs:
#   <out>/04_dual_modality/atac_z_block.rds       n_cells x 30 z-scored ATAC
#   <out>/04_dual_modality/rna_z_block_<tr>.rds   n_cells x K_rna z-scored RNA per trait
#
# Outputs:
#   <out>/04_dual_modality/knn_graph_<tr>.rds     row-normalized sparse n x n graph
#   <out>/04_dual_modality/knn_summary.csv        per-trait final params + iso_frac
#
# Red line: NO cell_type / cell_meta in this script.

source(file.path(dirname(sys.frame(1)$ofile %||% "."), "..", "_bootstrap.R"))
required_pkgs(c("data.table", "Matrix", "matrixStats", "RANN"))

paths   <- m3_paths()
out_dir <- file.path(paths$out_dir, "04_dual_modality")

# --- locked / starting params ---------------------------------------------
TRAITS         <- c("MONO", "PLT", "LYM")
K_INIT         <- 30L
K_HARD_CAP     <- 90L
ATAC_W_START   <- 0.7
RNA_W_START    <- 0.3
W_STEP         <- 0.05
W_FLOOR_ATAC   <- 0.5    # ATAC must stay >= 0.5 (主辅 principle, §2.4)
ISO_TOL        <- 0.01   # 1% isolated allowed without retune
# §2.5 doc says "weights first then k" but BMMC empirics: weight tuning
# barely moves iso_frac (~0.1pp per step) while k bump from 30→45 cuts
# isolated by half. Cap weight iterations at MAX_W_ITERS — if iso hasn't
# crossed the tol after that, jump straight to k.
MAX_W_ITERS    <- 1L

# --- helpers --------------------------------------------------------------

# Build weighted directed KNN, then mutual-filter into a symmetric weighted
# sparse graph (BEFORE row-normalization).
build_mutual_knn_weighted <- function(fused, k) {
  n <- nrow(fused)
  # +1 because nn2 returns the cell itself as neighbor 0 (distance 0)
  nn <- RANN::nn2(fused, k = k + 1L)
  nn_idx   <- nn$nn.idx[, 2:(k + 1L), drop = FALSE]
  nn_dists <- nn$nn.dists[, 2:(k + 1L), drop = FALSE]

  # per-cell median distance over the k neighbors → similarity scale
  med_d <- matrixStats::rowMedians(nn_dists)
  med_d[med_d <= 0] <- 1

  sim <- exp(-nn_dists / med_d)

  # directed sparse: i has k out-edges; row-major flatten
  i_vec <- rep(seq_len(n), each = k)
  j_vec <- as.vector(t(nn_idx))
  s_vec <- as.vector(t(sim))

  # drop any (i, i) self-loops that might have leaked through
  keep <- i_vec != j_vec
  i_vec <- i_vec[keep]; j_vec <- j_vec[keep]; s_vec <- s_vec[keep]

  M_dir <- Matrix::sparseMatrix(
    i = i_vec, j = j_vec, x = s_vec, dims = c(n, n)
  )

  # mutual filter via triplet intersection. Forces double arithmetic — for
  # n > sqrt(2^31) ≈ 46341 the integer linear index (j-1)*n+i overflows and
  # silently produces NAs (which then break the %in% check below).
  T_dir <- Matrix::summary(M_dir)
  n_d   <- as.numeric(n)
  lin     <- (as.numeric(T_dir$j) - 1) * n_d + as.numeric(T_dir$i)
  lin_rev <- (as.numeric(T_dir$i) - 1) * n_d + as.numeric(T_dir$j)
  mut_mask <- lin_rev %in% lin
  T_mut    <- T_dir[mut_mask, , drop = FALSE]

  # average i→j and j→i similarities
  match_idx   <- match(lin_rev[mut_mask], lin)
  T_mut$x_avg <- (T_mut$x + T_dir$x[match_idx]) / 2

  M_mut <- Matrix::sparseMatrix(
    i = T_mut$i, j = T_mut$j, x = T_mut$x_avg,
    dims = c(n, n), symmetric = FALSE
  )
  M_mut
}

row_normalize <- function(M) {
  deg <- Matrix::rowSums(M)
  deg_safe <- pmax(deg, 1e-12)
  D_inv <- Matrix::Diagonal(x = 1 / deg_safe)
  D_inv %*% M
}

# --- per-trait auto-tune loop --------------------------------------------

atac_z <- readRDS(file.path(out_dir, "atac_z_block.rds"))
n_cells <- nrow(atac_z)
cell_barcodes <- rownames(atac_z)
message("ATAC block: ", n_cells, " x ", ncol(atac_z))

summary_rows <- list()

for (tr in TRAITS) {
  message("\n[Trait] ", tr)
  rna_z <- readRDS(file.path(out_dir, paste0("rna_z_block_", tr, ".rds")))
  stopifnot(identical(rownames(rna_z), cell_barcodes))
  message("  RNA block: ", nrow(rna_z), " x ", ncol(rna_z))

  k       <- K_INIT
  atac_w  <- ATAC_W_START
  rna_w   <- RNA_W_START
  iter    <- 0L
  history <- list()

  M_mut   <- NULL
  iso_frac <- NA_real_

  repeat {
    iter <- iter + 1L
    fused <- cbind(atac_z * sqrt(atac_w), rna_z * sqrt(rna_w))

    t0 <- Sys.time()
    M_mut <- build_mutual_knn_weighted(fused, k)
    deg <- Matrix::rowSums(M_mut)
    iso_n    <- sum(deg == 0)
    iso_frac <- iso_n / n_cells
    dt <- round(as.numeric(Sys.time() - t0, units = "secs"))

    message(sprintf(
      "  iter %d: k=%d, atac_w=%.2f, rna_w=%.2f | nnz=%d, isolated=%d (%.2f%%) | %ds",
      iter, k, atac_w, rna_w, length(M_mut@x), iso_n, 100 * iso_frac, dt
    ))
    history[[iter]] <- data.table::data.table(
      iter = iter, k = k, atac_w = atac_w, rna_w = rna_w,
      nnz = length(M_mut@x), isolated = iso_n, iso_frac = iso_frac
    )

    if (iso_frac <= ISO_TOL) break

    weight_iters_done <- sum(vapply(history,
                                    function(h) h$atac_w < ATAC_W_START - 1e-9,
                                    logical(1)))

    # tune order: weights first (down to ATAC=0.5, capped at MAX_W_ITERS),
    # then k. The cap is a pragmatic short-circuit: BMMC weight tuning here
    # only moves iso by ~0.1pp per step while k bump moves it by half a pp.
    if (weight_iters_done < MAX_W_ITERS &&
        atac_w - W_STEP >= W_FLOOR_ATAC - 1e-9) {
      atac_w <- atac_w - W_STEP
      rna_w  <- rna_w  + W_STEP
      next
    }
    if (k < K_HARD_CAP) {
      k <- min(K_HARD_CAP, as.integer(ceiling(k * 1.5)))
      next
    }
    message("  reached weight cap AND k cap (", K_HARD_CAP,
            "); accepting iso_frac=", round(iso_frac, 4))
    break
  }

  # row-normalize + assign barcode dimnames
  M_norm <- row_normalize(M_mut)
  rownames(M_norm) <- cell_barcodes
  colnames(M_norm) <- cell_barcodes
  saveRDS(M_norm, file.path(out_dir, paste0("knn_graph_", tr, ".rds")))

  hist_dt <- data.table::rbindlist(history)
  data.table::fwrite(hist_dt, file.path(out_dir,
                                        paste0("knn_tune_history_", tr, ".csv")))

  deg <- Matrix::rowSums(M_mut)
  summary_rows[[tr]] <- data.table::data.table(
    trait        = tr,
    n_cells      = n_cells,
    k_final      = k,
    atac_w_final = atac_w,
    rna_w_final  = rna_w,
    iter_total   = iter,
    iso_frac     = iso_frac,
    nnz          = length(M_mut@x),
    deg_mean     = mean(deg),
    deg_median   = stats::median(deg)
  )
}

summary_dt <- data.table::rbindlist(summary_rows)
data.table::fwrite(summary_dt, file.path(out_dir, "knn_summary.csv"))
print(summary_dt)

message("\nWeighted mutual KNN built. Outputs in: ", out_dir)
