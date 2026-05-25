# 03_build_weighted_knn.R
# §2.C: weighted mutual KNN (per trait, per atac_w).
# 顶部 ATAC_W_LIST 控制要算哪些权重. Step 1 = c(0.7). Sweep 时改成 c(0.5,0.6,0.7,0.8,0.9).
# Skip-if-exists: 同名 .rds 已在跳过.
#
# 脚本功能：
# 读取 ATAC 表征和每个 trait 对应的 RNA 表征，将两种模态按 atac_w 加权拼接。
# 在融合后的细胞嵌入空间中寻找最近邻，并只保留互为最近邻的边。
# 输出的是行归一化后的稀疏邻接矩阵，供后续随机游走直接使用。

ATAC_W_LIST <- c(1.0)  # 对照: 纯 ATAC (RNA 项归零), 跟 atac_w=0.7 比看融合贡献

suppressPackageStartupMessages({
  library(Matrix)
  library(BiocNeighbors)
})

OUT_DIR <- "D:/SCAVENGE-reproducibility-main/BMMC/04_dual_modality"
TRAITS  <- c("MONO", "PLT", "LYM")

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

atac_z_full <- readRDS(file.path(OUT_DIR, "atac_z_block_v2.rds"))
cat("atac_z_full rows: ", nrow(atac_z_full), "\n")

knn_log <- data.frame()
for (tname in TRAITS) {
  # 每个 trait 的 RNA 表征可能不同，因此按 trait 分别读取 RNA block。
  # rownames(rna_z) 定义本轮建图使用的共同细胞顺序。
  rna_z <- readRDS(file.path(OUT_DIR, paste0("rna_z_block_v2_", tname, ".rds")))
  common <- rownames(rna_z)  # SCDRS-aligned 68303
  atac_z <- atac_z_full[common, , drop = FALSE]
  for (w in ATAC_W_LIST) {
    out_file <- file.path(OUT_DIR, sprintf("knn_v2_%s_w%.1f.rds", tname, w))
    if (file.exists(out_file)) {
      cat("[knn] skip exists: ", basename(out_file), "\n")
      next
    }
    cat("\n=== ", tname, "  atac_w=", w, "  cells=", length(common), " ===\n", sep = "")

    # 使用 sqrt(w) 和 sqrt(1-w) 缩放两个模态后再拼接。
    # 这样欧氏距离中的平方项近似对应 w * ATAC距离 + (1-w) * RNA距离。
    fusion <- cbind(sqrt(w) * atac_z, sqrt(1 - w) * rna_z)
    t0 <- Sys.time()
    res <- build_knn(fusion)
    cat("  elapsed: ",
        round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2),
        " min\n", sep = "")
    saveRDS(res$adj, out_file)
    cat("  saved ", basename(out_file),
        "  k_final=", res$k_final,
        " iso=", round(res$iso_frac, 4), "\n", sep = "")
    knn_log <- rbind(knn_log, data.frame(
      trait = tname, atac_w = w, k_final = res$k_final,
      iso_frac = round(res$iso_frac, 4)
    ))
  }
}

sum_file <- file.path(OUT_DIR, "knn_summary_v2.csv")
if (file.exists(sum_file)) {
  # 如果之前已经跑过其他权重或 trait，保留旧记录并追加本轮记录。
  old <- read.csv(sum_file)
  knn_log <- unique(rbind(old, knn_log))
}
write.csv(knn_log, sum_file, row.names = FALSE)
cat("\n")
print(knn_log)
cat("\n03_build_weighted_knn.R OK.\n")
