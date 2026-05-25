# 05b_rwr_trs_dualseed.R
# §2.E' RWR + TRS using DUAL-seed (ATAC ∪ SCDRS top 1% union seed).
# 跟 05_rwr_trs.R 完全一致, 只换 seed 文件来源 + 输出文件名.
# 复用阶段 2 已存的 knn_MultiOmicsFusion_<trait>_w<W>.rds, 不重算 KNN.
# Skip-if-exists 续传安全.

ATAC_W_LIST <- c(0.7, 1.0)  # 两个权重都跑, 跟 05_rwr_trs.R 主流水线对齐

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
})

OUT_DIR <- "D:/RStudio/work/scavenge_multiomics_fusion_reproduction/reproduction_outputs/BMMC/04_dual_modality"
TRAITS  <- c("MONO", "PLT", "LYM")
GAMMA   <- 0.05
TOL     <- 1e-6
MAX_ITER  <- 200
CAP_Q     <- 0.95
SCALE_FAC <- 1.0

rwr <- function(A_norm, seed_vec, gamma = GAMMA, tol = TOL, max_iter = MAX_ITER) {
  s <- seed_vec / sum(seed_vec)
  x <- s
  for (iter in seq_len(max_iter)) {
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
  seed_file <- file.path(OUT_DIR, paste0("seed_dual_", tname, ".rds"))
  if (!file.exists(seed_file)) {
    stop("missing dual seed file (run 04b_seed_dual.R first): ", seed_file)
  }
  seed <- readRDS(seed_file)
  cat("\n=== ", tname, " (n_seed = ", sum(seed), ") ===\n", sep = "")

  for (w in ATAC_W_LIST) {
    out_file <- file.path(OUT_DIR,
      sprintf("TRS_dualseed_%s_w%.1f.tsv.gz", tname, w))
    if (file.exists(out_file)) {
      cat("[rwr] skip exists: ", basename(out_file), "\n", sep = "")
      next
    }
    knn_file <- file.path(OUT_DIR, sprintf("knn_MultiOmicsFusion_%s_w%.1f.rds", tname, w))
    if (!file.exists(knn_file))
      stop("missing knn file (阶段 2 应当已存): ", knn_file)
    A_norm <- readRDS(knn_file)
    cat("  loaded knn (", paste(dim(A_norm), collapse = " x "),
        "), atac_w=", w, "\n", sep = "")

    cat("  RWR (gamma=", GAMMA, ") with dual seed ...\n", sep = "")
    t0 <- Sys.time()
    trs_raw <- rwr(A_norm, seed)
    cat("  elapsed: ",
        round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1),
        "s\n", sep = "")

    trs_p <- post_process(trs_raw)
    names(trs_p) <- rownames(A_norm)

    out <- data.frame(barcode = names(trs_p), TRS = trs_p)
    fwrite(out, out_file, sep = "\t", compress = "gzip")
    cat("  saved ", basename(out_file), "\n", sep = "")

    trs_log <- rbind(trs_log, data.frame(
      trait = tname, atac_w = w, n_seed = sum(seed), n_cell = length(trs_p),
      mean_TRS = round(mean(trs_p), 4),
      max_TRS  = round(max(trs_p),  4)
    ))
  }
}

sum_file <- file.path(OUT_DIR, "trs_dualseed_summary.csv")
if (file.exists(sum_file)) {
  old <- read.csv(sum_file)
  trs_log <- unique(rbind(old, trs_log))
}
write.csv(trs_log, sum_file, row.names = FALSE)
cat("\n")
print(trs_log)
cat("\n05b_rwr_trs_dualseed.R OK.\n")
