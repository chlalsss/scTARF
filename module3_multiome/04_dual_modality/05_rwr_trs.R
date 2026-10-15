# Reusable scTARF core definitions; calculation parameters are preserved.
# Use examples/run_core_demo.R to load graph and propagation in separate scopes.
# Run paper_reproduction/run_paper.R for the supplied manuscript data.

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
})

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
