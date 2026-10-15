# Reusable scTARF core definitions; calculation parameters are preserved.
# Use examples/run_core_demo.R to load graph and propagation in separate scopes.
# Run paper_reproduction/run_paper.R for the supplied manuscript data.

suppressPackageStartupMessages({
  library(Matrix)
  library(BiocNeighbors)
})

# K_INIT 是初始近邻数；如果 mutual KNN 后孤立细胞太多，可以把 k 放大。
# ISO_TOL 控制允许的孤立细胞比例；MAX_ITER 控制最多放大几轮。
K_INIT   <- 30
K_MAX    <- 90
ISO_TOL  <- 0.02
MAX_ITER <- 1

# emb 是已经融合好的细胞嵌入矩阵，行是细胞，列是特征。
# 返回值 adj 是行归一化的 mutual KNN 邻接矩阵。
build_knn <- function(emb, k_init = K_INIT, k_max = K_MAX,
                      iso_tol = ISO_TOL, max_iter = MAX_ITER) {
  n <- nrow(emb)
  k <- k_init
  iter <- 0
  repeat {
    cat("    KNN k=", k, " ...\n", sep = "")

    # 使用 Annoy 近似最近邻以提高速度；nn$index 是每个细胞的 k 个邻居。
    nn <- BiocNeighbors::findKNN(emb, k = k, BNPARAM = AnnoyParam())

    # 将距离转换为边权重。每个细胞用自身邻居距离中位数做尺度归一化，
    # 可减少不同细胞局部密度差异对边权重的影响。
    median_d <- apply(nn$distance, 1, median)
    median_d[median_d == 0] <- 1
    edge_w <- exp(-nn$distance / median_d)

    # 把邻居列表转换成稀疏邻接矩阵 A。
    # A[i, j] 表示细胞 i 指向细胞 j 的近邻边权重。
    i_idx <- rep(seq_len(n), each = k)
    j_idx <- as.integer(t(nn$index))
    w_idx <- as.numeric(t(edge_w))
    A <- sparseMatrix(i = i_idx, j = j_idx, x = w_idx, dims = c(n, n))

    # mutual filter: 只有 i 选 j 且 j 也选 i 时才保留边。
    # 保留后的边权使用两个方向权重的平均值，使邻接矩阵对称。
    At <- t(A)
    mutual_mask <- (A != 0) * (At != 0)   # binary sparse
    A_sym <- (A + At) / 2
    A_mut <- A_sym * mutual_mask          # zeros out non-mutual entries
    A_mut <- drop0(A_mut)                 # remove explicit zeros

    # deg 是每个细胞保留下来的 mutual 邻居数。
    # deg 为 0 的细胞是孤立细胞，后续随机游走无法从图结构获得传播信息。
    deg <- diff(A_mut@p)
    iso_frac <- mean(deg == 0)
    cat("    iso_frac=", round(iso_frac, 4), "\n", sep = "")

    # 如果孤立细胞比例已经足够低，或者达到最大迭代/最大 k，则停止调参。
    if (iso_frac <= iso_tol || iter >= max_iter || k >= k_max) break
    k <- min(ceiling(k * 1.5), k_max)
    iter <- iter + 1
  }

  # 行归一化后，每一行可以视为从当前细胞跳转到邻居细胞的概率分布。
  # 孤立细胞的行和为 0，这里临时设为 1，避免除零。
  rs <- rowSums(A_mut)
  rs[rs == 0] <- 1
  A_norm <- Diagonal(n, 1 / rs) %*% A_mut
  rownames(A_norm) <- rownames(emb)
  colnames(A_norm) <- rownames(emb)
  list(adj = A_norm, k_final = k, iso_frac = iso_frac)
}
