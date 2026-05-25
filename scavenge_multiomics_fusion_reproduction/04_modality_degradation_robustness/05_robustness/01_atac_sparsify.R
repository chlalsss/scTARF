# 01_atac_sparsify.R
# §3.1: ATAC binomial downsample sweep + 完整下游 (LSI 重算 + KNN + RWR + AUROC).
#
# Step 1 sanity = P_DROP_LIST=c(0.0), SEED_LIST=c(1). 验证 AUROC 等价 MultiOmicsFusion.
# Step 2 sweep  = P_DROP_LIST=c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5), SEED_LIST=c(1,2,3).
#
# Skip-if-exists: 同名 atac_z / knn / TRS 已在则跳过, 续传安全.
# 预期 sanity 数字 (atac_w=0.7): MONO 0.976 / PLT 0.866 / LYM_TNK 0.886.

# ============================================================================
# 参数 (Step 1 sanity 配置, 切到 sweep 时改这里即可)
# ============================================================================
P_DROP_LIST <- c(0.6, 0.7)
SEED_LIST   <- c(1, 2, 3)
ATAC_W_LIST <- c(0.7, 1.0)
TRAITS      <- c("MONO", "PLT", "LYM")

suppressPackageStartupMessages({
  library(Matrix)
  library(SummarizedExperiment)
  library(SCAVENGE)
  library(BiocNeighbors)
  library(data.table)
  library(pROC)
})

# ============================================================================
# 路径
# ============================================================================
SE_FILE   <- "D:/SCAVENGE-reproducibility-main/BMMC/03_scavenge_baseline/SE_atac.rds"
DUAL_DIR  <- "D:/SCAVENGE-reproducibility-main/BMMC/04_dual_modality"
META_FILE <- "D:/SCAVENGE-reproducibility-main/BMMC/01_data_qc/BMMC_cell_meta.rds"
OUT_DIR   <- "D:/SCAVENGE-reproducibility-main/BMMC/05_robustness"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# 评估靶 (跟 04_dual_modality/06_evaluate_auroc.R 一致)
# ============================================================================
TARGETS <- list(
  MONO     = c("CD14+ Mono", "CD16+ Mono"),
  PLT      = c("MK/E prog"),
  LYM      = c("B1 B", "Naive CD20+ B", "Transitional B", "Plasma cell"),
  LYM_TNK  = c("NK", "ILC", "CD4+ T activated", "CD4+ T naive",
               "CD8+ T", "CD8+ T naive"),
  LYM_wide = c("NK", "ILC", "CD4+ T activated", "CD4+ T naive",
               "CD8+ T", "CD8+ T naive", "Lymph prog",
               "B1 B", "Naive CD20+ B", "Transitional B", "Plasma cell")
)
PRED_TRAIT <- list(MONO = "MONO", PLT = "PLT", LYM = "LYM",
                   LYM_TNK = "LYM", LYM_wide = "LYM")

# ============================================================================
# KNN / RWR 参数 (跟阶段 2 MultiOmicsFusion 完全一致)
# ============================================================================
K_INIT       <- 30
K_MAX        <- 90
ISO_TOL      <- 0.02
MAX_ITER_KNN <- 1
GAMMA        <- 0.05
TOL          <- 1e-6
MAX_ITER_RWR <- 200
CAP_Q        <- 0.95

# ============================================================================
# helpers
# ============================================================================

# binomial read-level downsample (peak x cell counts)
sparsify_binomial <- function(counts, p_keep, seed) {
  if (p_keep >= 1.0) return(counts)
  set.seed(seed)
  out <- counts
  out@x <- as.numeric(rbinom(length(counts@x), size = as.integer(counts@x),
                             prob = p_keep))
  drop0(out)
}

# weighted mutual KNN (cf. 04_dual_modality/03_build_weighted_knn.R)
build_knn <- function(emb, k_init = K_INIT, k_max = K_MAX,
                      iso_tol = ISO_TOL, max_iter = MAX_ITER_KNN) {
  n <- nrow(emb)
  k <- k_init
  iter <- 0
  repeat {
    nn <- BiocNeighbors::findKNN(emb, k = k, BNPARAM = AnnoyParam())
    median_d <- apply(nn$distance, 1, median)
    median_d[median_d == 0] <- 1
    edge_w <- exp(-nn$distance / median_d)

    i_idx <- rep(seq_len(n), each = k)
    j_idx <- as.integer(t(nn$index))
    w_idx <- as.numeric(t(edge_w))
    A <- sparseMatrix(i = i_idx, j = j_idx, x = w_idx, dims = c(n, n))

    At <- t(A)
    mutual_mask <- (A != 0) * (At != 0)
    A_sym <- (A + At) / 2
    A_mut <- drop0(A_sym * mutual_mask)

    deg <- diff(A_mut@p)
    iso_frac <- mean(deg == 0)
    if (iso_frac <= iso_tol || iter >= max_iter || k >= k_max) break
    k <- min(ceiling(k * 1.5), k_max)
    iter <- iter + 1
  }
  rs <- rowSums(A_mut)
  rs[rs == 0] <- 1
  A_norm <- Diagonal(n, 1 / rs) %*% A_mut
  rownames(A_norm) <- rownames(emb)
  colnames(A_norm) <- rownames(emb)
  list(adj = A_norm, k_final = k, iso_frac = iso_frac)
}

# RWR (cf. 04_dual_modality/05_rwr_trs.R)
rwr <- function(A_norm, seed_vec, gamma = GAMMA, tol = TOL,
                max_iter = MAX_ITER_RWR) {
  s <- seed_vec / sum(seed_vec)
  x <- s
  for (iter in seq_len(max_iter)) {
    x_new <- gamma * s + (1 - gamma) * as.numeric(A_norm %*% x)
    delta <- sum(abs(x_new - x))
    x <- x_new
    if (delta < tol) return(x)
  }
  x
}

post_process <- function(trs, cap_q = CAP_Q) {
  cap <- as.numeric(quantile(trs, cap_q))
  trs_c <- pmin(trs, cap)
  rng <- range(trs_c)
  if (rng[2] - rng[1] > 0)
    (trs_c - rng[1]) / (rng[2] - rng[1])
  else
    rep(0, length(trs))
}

compute_auroc <- function(predictor, labels) {
  if (sum(labels == 1) < 10 || sum(labels == 0) < 10) return(NA_real_)
  roc_obj <- suppressMessages(pROC::roc(labels, predictor, quiet = TRUE))
  as.numeric(pROC::auc(roc_obj))
}

# ============================================================================
# 0. meta + SE_atac (一次性加载)
# ============================================================================
cat("[0] loading meta + SE_atac (~613MB on disk) ...\n")
meta <- readRDS(META_FILE)
if (inherits(meta, "data.table")) meta <- as.data.frame(meta)
if (is.null(rownames(meta)) ||
    identical(rownames(meta), as.character(seq_len(nrow(meta))))) {
  if ("barcode" %in% colnames(meta)) rownames(meta) <- meta$barcode
}

SE <- readRDS(SE_FILE)
counts_full <- SummarizedExperiment::assay(SE, "counts")
all_cells   <- colnames(counts_full)
cat("  SE peak x cell:", nrow(counts_full), "x", ncol(counts_full),
    "; nnz=", length(counts_full@x), "\n")

# ============================================================================
# 1. RNA z + seed (atac_w 无关, sparsify 不影响, sweep 全程 reuse)
# ============================================================================
cat("[1] loading rna_z + seed for", paste(TRAITS, collapse = ","), "...\n")
rna_z_list <- list()
seed_list  <- list()
for (tn in TRAITS) {
  rna_z_list[[tn]] <- readRDS(file.path(DUAL_DIR,
                                         paste0("rna_z_block_MultiOmicsFusion_", tn, ".rds")))
  seed_list[[tn]]  <- readRDS(file.path(DUAL_DIR,
                                         paste0("seed_MultiOmicsFusion_", tn, ".rds")))
}
common_cells <- rownames(rna_z_list[[1]])  # 68303 SCDRS-aligned
cat("  common (SCDRS-aligned) cells:", length(common_cells), "\n")
stopifnot(all(common_cells %in% all_cells))

# ============================================================================
# 2. main loop
# ============================================================================
eval_results <- data.frame()

for (p_drop in P_DROP_LIST) {
  p_keep <- 1 - p_drop
  for (sd in SEED_LIST) {
    cat("\n", strrep("=", 70), "\n", sep = "")
    cat("[main] p_drop=", p_drop, " seed=", sd, "\n", sep = "")
    cat(strrep("=", 70), "\n", sep = "")

    # 2a. sparsify + LSI
    lsi_file <- file.path(OUT_DIR, sprintf("atac_lsi_p%.2f_seed%d.rds",
                                            p_drop, sd))
    z_file   <- file.path(OUT_DIR, sprintf("atac_z_block_p%.2f_seed%d.rds",
                                            p_drop, sd))

    if (file.exists(z_file)) {
      cat("  [2a] skip exists: ", basename(z_file), "\n", sep = "")
      atac_z_full <- readRDS(z_file)
    } else {
      cat("  [2a] sparsify_binomial p_keep=", p_keep, " seed=", sd, " ...\n",
          sep = "")
      t0 <- Sys.time()
      counts_sp <- sparsify_binomial(counts_full, p_keep, sd)
      cat("    nnz: ", length(counts_full@x), " -> ", length(counts_sp@x),
          " (", sprintf("%.1f",
                        length(counts_sp@x) / length(counts_full@x) * 100),
          "%)\n", sep = "")

      cat("  [2a] tfidf + do_lsi (dims=30) ...\n")
      tfidf_mat <- SCAVENGE::tfidf(bmat = counts_sp,
                                    mat_binary = TRUE, TF = TRUE, log_TF = TRUE)
      lsi_mat <- SCAVENGE::do_lsi(mat = tfidf_mat, dims = 30)

      # orient cell x 30
      if (ncol(lsi_mat) != 30 && nrow(lsi_mat) == 30) lsi_mat <- t(lsi_mat)
      stopifnot(ncol(lsi_mat) == 30)

      # cell rownames 强制对齐 SE 列顺序
      if (is.null(rownames(lsi_mat)) || nrow(lsi_mat) != length(all_cells)) {
        if (nrow(lsi_mat) == length(all_cells)) {
          rownames(lsi_mat) <- all_cells
        } else {
          stop("LSI cell dim mismatch: ", nrow(lsi_mat),
               " vs SE ", length(all_cells))
        }
      }

      saveRDS(lsi_mat, lsi_file)
      atac_z_full <- scale(lsi_mat)
      attributes(atac_z_full)$`scaled:center` <- NULL
      attributes(atac_z_full)$`scaled:scale`  <- NULL
      saveRDS(atac_z_full, z_file)
      cat("    elapsed: ",
          round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2),
          " min\n", sep = "")
      cat("    saved ", basename(z_file),
          " (", paste(dim(atac_z_full), collapse = " x "), ")\n", sep = "")
    }

    # 2b. per (trait, w): KNN + RWR + TRS + AUROC
    for (tn in TRAITS) {
      rna_z    <- rna_z_list[[tn]]
      seed_vec <- seed_list[[tn]]
      atac_z   <- atac_z_full[common_cells, , drop = FALSE]
      stopifnot(identical(rownames(atac_z), rownames(rna_z)))

      for (w in ATAC_W_LIST) {
        knn_file <- file.path(OUT_DIR,
          sprintf("knn_p%.2f_seed%d_%s_w%.1f.rds", p_drop, sd, tn, w))
        trs_file <- file.path(OUT_DIR,
          sprintf("TRS_p%.2f_seed%d_%s_w%.1f.tsv.gz", p_drop, sd, tn, w))

        if (file.exists(trs_file)) {
          cat("    [2b] skip exists: ", basename(trs_file), "\n", sep = "")
          MultiOmicsFusion <- fread(trs_file)
          trs_p <- MultiOmicsFusion$TRS
          names(trs_p) <- MultiOmicsFusion$barcode
        } else {
          cat("\n  [2b] trait=", tn, " atac_w=", w, " ...\n", sep = "")
          fusion <- cbind(sqrt(w) * atac_z, sqrt(1 - w) * rna_z)

          if (file.exists(knn_file)) {
            cat("    knn cached, loading: ", basename(knn_file), "\n", sep = "")
            A_norm <- readRDS(knn_file)
          } else {
            t1 <- Sys.time()
            res <- build_knn(fusion)
            cat("    KNN k_final=", res$k_final, " iso=",
                round(res$iso_frac, 4), " (",
                round(as.numeric(difftime(Sys.time(), t1, units = "mins")), 2),
                " min)\n", sep = "")
            A_norm <- res$adj
            saveRDS(A_norm, knn_file)
          }

          t2 <- Sys.time()
          trs_raw <- rwr(A_norm, seed_vec)
          trs_p <- post_process(trs_raw)
          names(trs_p) <- rownames(A_norm)
          fwrite(data.frame(barcode = names(trs_p), TRS = trs_p),
                 trs_file, sep = "\t", compress = "gzip")
          cat("    RWR+TRS done (",
              round(as.numeric(difftime(Sys.time(), t2, units = "secs")), 1),
              "s)\n", sep = "")
        }

        # 2c. AUROC per applicable target
        for (tgt_name in names(TARGETS)) {
          if (PRED_TRAIT[[tgt_name]] != tn) next
          ct_pos <- TARGETS[[tgt_name]]
          common_eval <- intersect(names(trs_p), rownames(meta))
          lbl <- as.integer(as.character(meta[common_eval, "cell_type"]) %in% ct_pos)
          a <- compute_auroc(trs_p[common_eval], lbl)
          cat(sprintf("    AUROC[%s -> %-9s w=%.1f] = %.4f  (pos=%d)\n",
                      tn, tgt_name, w, a, sum(lbl)))
          eval_results <- rbind(eval_results, data.frame(
            p_drop = p_drop, seed = sd, trait = tn, target = tgt_name,
            atac_w = w, auroc = a
          ))
        }
      }
    }
  }
}

# ============================================================================
# 3. eval csv (append)
# ============================================================================
eval_file <- file.path(OUT_DIR, "atac_sparsify_eval.csv")
if (file.exists(eval_file)) {
  old <- read.csv(eval_file)
  eval_results <- unique(rbind(old, eval_results))
}
write.csv(eval_results, eval_file, row.names = FALSE)
cat("\nsaved ", basename(eval_file), " (", nrow(eval_results), " rows)\n",
    sep = "")

# ============================================================================
# 4. sanity headline (p_drop=0)
# ============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("SANITY HEADLINE  expected MultiOmicsFusion atac_w=0.7: MONO 0.976 / PLT 0.866 / LYM_TNK 0.886\n")
cat(strrep("=", 70), "\n", sep = "")
sanity <- eval_results[eval_results$p_drop == 0, ]
sanity_keys <- list(c("MONO", "MONO"), c("PLT", "PLT"), c("LYM", "LYM_TNK"))
exp_MultiOmicsFusion_w07 <- list("MONO/MONO" = 0.976, "PLT/PLT" = 0.866, "LYM/LYM_TNK" = 0.886)
exp_MultiOmicsFusion_w10 <- list("MONO/MONO" = 0.973, "PLT/PLT" = 0.864, "LYM/LYM_TNK" = 0.873)
for (w in ATAC_W_LIST) {
  cat("\n  atac_w =", w, ":\n")
  exp_tbl <- if (abs(w - 0.7) < 1e-6) exp_MultiOmicsFusion_w07 else
             if (abs(w - 1.0) < 1e-6) exp_MultiOmicsFusion_w10 else NULL
  for (kk in sanity_keys) {
    sub <- sanity[sanity$trait == kk[1] & sanity$target == kk[2] &
                  sanity$atac_w == w, ]
    if (nrow(sub) == 0) next
    key <- paste0(kk[1], "/", kk[2])
    obs <- mean(sub$auroc)
    delta_msg <- if (!is.null(exp_tbl) && key %in% names(exp_tbl)) {
      sprintf(" (delta=%+.4f vs expected %.3f)",
              obs - exp_tbl[[key]], exp_tbl[[key]])
    } else ""
    cat(sprintf("    %-12s = %.4f%s\n", key, obs, delta_msg))
  }
}

cat("\n01_atac_sparsify.R OK.\n")
