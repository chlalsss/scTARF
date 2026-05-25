# 05_rwr_trs.R
# §2.E: Random Walk with Restart + TRS post-process (per trait, per atac_w).
# 顶部 ATAC_W_LIST 跟 03_build_weighted_knn.R 保持一致.
# Skip-if-exists: 同名 TRS_v2_<trait>_w<w>.tsv.gz 在跳过.

ATAC_W_LIST <- c(0.7)

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
})

OUT_DIR <- "D:/RStudio/work/acceptance_demo_runs/run_20260520_152900/outputs/algorithm_core_mono"
TRAITS  <- c("MONO")
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
  old <- read.csv(sum_file)
  trs_log <- unique(rbind(old, trs_log))
}
write.csv(trs_log, sum_file, row.names = FALSE)
cat("\n")
print(trs_log)
cat("\n05_rwr_trs.R OK.\n")
