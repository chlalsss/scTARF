# 03_build_weighted_knn.R
# §2.C: weighted mutual KNN (per trait, per atac_w).
# 顶部 ATAC_W_LIST 控制要算哪些权重. Step 1 = c(0.7). Sweep 时改成 c(0.5,0.6,0.7,0.8,0.9).
# Skip-if-exists: 同名 .rds 已在跳过.

ATAC_W_LIST <- c(0.7)

suppressPackageStartupMessages({
  library(Matrix)
  library(BiocNeighbors)
})

OUT_DIR <- "D:/RStudio/work/acceptance_demo_runs/run_20260519_192610/outputs/algorithm_core_mono"
TRAITS  <- c("MONO")
K_INIT   <- 30
K_MAX    <- 90
ISO_TOL  <- 0.02
MAX_ITER <- 1

build_knn <- function(emb, k_init = K_INIT, k_max = K_MAX,
                      iso_tol = ISO_TOL, max_iter = MAX_ITER) {
  n <- nrow(emb)
  k <- k_init
  iter <- 0
  repeat {
    cat("    KNN k=", k, " ...\n", sep = "")
    nn <- BiocNeighbors::findKNN(emb, k = k, BNPARAM = AnnoyParam())
    median_d <- apply(nn$distance, 1, median)
    median_d[median_d == 0] <- 1
    edge_w <- exp(-nn$distance / median_d)

    i_idx <- rep(seq_len(n), each = k)
    j_idx <- as.integer(t(nn$index))
    w_idx <- as.numeric(t(edge_w))
    A <- sparseMatrix(i = i_idx, j = j_idx, x = w_idx, dims = c(n, n))

    # mutual filter, keep symmetric weight (A + At)/2 only at mutual positions
    At <- t(A)
    mutual_mask <- (A != 0) * (At != 0)   # binary sparse
    A_sym <- (A + At) / 2
    A_mut <- A_sym * mutual_mask          # zeros out non-mutual entries
    A_mut <- drop0(A_mut)                 # remove explicit zeros

    deg <- diff(A_mut@p)
    iso_frac <- mean(deg == 0)
    cat("    iso_frac=", round(iso_frac, 4), "\n", sep = "")

    if (iso_frac <= iso_tol || iter >= max_iter || k >= k_max) break
    k <- min(ceiling(k * 1.5), k_max)
    iter <- iter + 1
  }

  # row-normalize
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
  old <- read.csv(sum_file)
  knn_log <- unique(rbind(old, knn_log))
}
write.csv(knn_log, sum_file, row.names = FALSE)
cat("\n")
print(knn_log)
cat("\n03_build_weighted_knn.R OK.\n")
