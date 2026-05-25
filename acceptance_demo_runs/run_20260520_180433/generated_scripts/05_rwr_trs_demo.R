# 05_rwr_trs.R
# §2.E: Random Walk with Restart + TRS post-process (per trait, per atac_w).
# 顶部 ATAC_W_LIST 跟 03_build_weighted_knn.R 保持一致.
# Skip-if-exists: 同名 TRS_v2_<trait>_w<w>.tsv.gz 在跳过.
#
# 脚本功能：
# 读取 seed 向量和对应的 KNN 图，在图上执行带重启随机游走。
# 随机游走结果表示 seed 信号传播到每个细胞后的稳态分数。
# 后处理会截断极端高值并缩放到统一范围，得到最终 TRS。

ATAC_W_LIST <- c(0.7)

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
})

OUT_DIR <- "D:/RStudio/work/acceptance_demo_runs/run_20260520_180433/outputs/algorithm_core_mono"
TRAITS  <- c("MONO")

# GAMMA 是重启概率：每轮有 gamma 的权重回到 seed 分布。
# TOL 是收敛阈值，MAX_ITER 是最大迭代次数。
# CAP_Q 用于截断 TRS 的高分极端值，避免少数异常点压缩整体尺度。
GAMMA   <- 0.05
TOL     <- 1e-6
MAX_ITER  <- 200
CAP_Q     <- 0.95
SCALE_FAC <- 1.0

# A_norm 是行归一化的邻接矩阵；seed_vec 是 0/1 seed 向量。
# s 是归一化后的 seed 概率分布，x 是当前迭代的传播分数。
rwr <- function(A_norm, seed_vec, gamma = GAMMA, tol = TOL, max_iter = MAX_ITER) {
  s <- seed_vec / sum(seed_vec)
  x <- s
  for (iter in seq_len(max_iter)) {
    # 带重启随机游走更新公式：
    # x_new = gamma * seed分布 + (1 - gamma) * 图传播后的上一轮分数。
    x_new <- gamma * s + (1 - gamma) * as.numeric(A_norm %*% x)
    delta <- sum(abs(x_new - x))
    x <- x_new
    if (delta < tol) {
      cat("    rwr converged iter=", iter, " delta=",
          format(delta, scientific = TRUE, digits = 3), "\n", sep = "")
      return(x)
    }
  }
  cat("    rwr maxed iter=", max_iter, " delta=",
      format(delta, scientific = TRUE, digits = 3), "\n", sep = "")
  x
}

# 对原始 RWR 分数做后处理：
# 先按 cap_q 分位数截断高值，再做 min-max 归一化。
# 这样不同 trait/权重下输出的 TRS 尺度更容易比较。
post_process <- function(trs, cap_q = CAP_Q, scale_factor = SCALE_FAC) {
  cap <- as.numeric(quantile(trs, cap_q))
  trs_c <- pmin(trs, cap)
  rng <- range(trs_c)
  if (rng[2] - rng[1] > 0)
    (trs_c - rng[1]) / (rng[2] - rng[1]) * scale_factor
  else
    rep(0, length(trs))
}

trs_log <- data.frame()
for (tname in TRAITS) {
  # 每个 trait 使用自己的 seed 向量。
  seed <- readRDS(file.path(OUT_DIR, paste0("seed_v2_", tname, ".rds")))
  for (w in ATAC_W_LIST) {
    out_file <- file.path(OUT_DIR, sprintf("TRS_v2_%s_w%.1f.tsv.gz", tname, w))
    if (file.exists(out_file)) {
      cat("[rwr] skip exists: ", basename(out_file), "\n")
      next
    }
    cat("\n=== ", tname, "  atac_w=", w, " ===\n", sep = "")
    knn_file <- file.path(OUT_DIR, sprintf("knn_v2_%s_w%.1f.rds", tname, w))
    if (!file.exists(knn_file))
      stop("missing knn file (run 03 first): ", knn_file)

    # KNN 文件中保存的是行归一化邻接矩阵，行列名都是细胞 barcode。
    A_norm <- readRDS(knn_file)
    cat("  loaded knn (", paste(dim(A_norm), collapse = " x "), ")\n", sep = "")

    cat("  RWR (gamma=", GAMMA, ") ...\n", sep = "")
    t0 <- Sys.time()
    trs_raw <- rwr(A_norm, seed)
    cat("  elapsed: ",
        round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1),
        "s\n", sep = "")

    trs_p <- post_process(trs_raw)
    names(trs_p) <- rownames(A_norm)

    # 输出每个细胞一行：barcode + TRS。
    # 文件压缩为 gzip，便于保存完整细胞级分数。
    out <- data.frame(barcode = names(trs_p), TRS = trs_p)
    fwrite(out, out_file, sep = "\t", compress = "gzip")
    cat("  saved ", basename(out_file), "\n", sep = "")

    trs_log <- rbind(trs_log, data.frame(
      trait = tname, atac_w = w, n_cell = length(trs_p),
      mean_TRS = round(mean(trs_p), 4),
      max_TRS  = round(max(trs_p),  4)
    ))
  }
}

sum_file <- file.path(OUT_DIR, "trs_summary_v2.csv")
if (file.exists(sum_file)) {
  # 保留之前 trait/权重的摘要记录，并追加本轮记录。
  old <- read.csv(sum_file)
  trs_log <- unique(rbind(old, trs_log))
}
write.csv(trs_log, sum_file, row.names = FALSE)
cat("\n")
print(trs_log)
cat("\n05_rwr_trs.R OK.\n")
