# 02_fusion_autotune.R
# §2.B of v2: ATAC LSI + RNA PCA fusion, auto-tune atac_w by silhouette of SCDRS top 1%
# 直接在 RStudio source(). 用全新 session.

suppressPackageStartupMessages({
  library(Matrix)
  library(irlba)        # truncated PCA
  library(bluster)      # approxSilhouette (KNN-based, O(n*k))
  library(data.table)   # fread for .score.gz
})

# ---- paths ----
ATAC_LSI   <- "D:/SCAVENGE-reproducibility-main/BMMC/03_scavenge_baseline/atac_lsi.rds"
RNA_COUNTS <- "D:/SCAVENGE-reproducibility-main/BMMC/01_data_qc/BMMC_rna_counts.rds"
RNA_FEATS  <- "D:/SCAVENGE-reproducibility-main/BMMC/01_data_qc/BMMC_rna_features.rds"
META       <- "D:/SCAVENGE-reproducibility-main/BMMC/01_data_qc/BMMC_cell_meta.rds"
GS_FILE    <- "D:/RStudio/work/module3_multiome/04_dual_modality/01_scdrs_python/bmmc_3traits.gs"
SCORE_DIR  <- "D:/RStudio/work/module3_multiome/04_dual_modality/01_scdrs_python/output"
OUT_DIR    <- "D:/SCAVENGE-reproducibility-main/BMMC/04_dual_modality"

TRAITS <- c(
  MONO = "UKB_460K.blood_MONOCYTE_COUNT",
  PLT  = "UKB_460K.blood_PLATELET_COUNT",
  LYM  = "UKB_460K.blood_LYMPHOCYTE_COUNT"
)

ATAC_W_GRID   <- c(0.5, 0.6, 0.7, 0.8, 0.9)
RNA_PCA_MAX_K <- 30
RNA_PCA_MIN_K <- 10
LABEL_TOP_PCT <- 0.01   # SCDRS norm_score top 1% as A class

# ---- 1. ATAC LSI -> z-scored atac_z_block (cell x 30) ----
cat("[1] loading ATAC LSI ...\n")
atac_lsi <- readRDS(ATAC_LSI)
cat("  class: ", class(atac_lsi)[1],
    " dim: ", paste(dim(atac_lsi), collapse = " x "), "\n")
if (ncol(atac_lsi) != 30) {
  if (nrow(atac_lsi) == 30) {
    atac_lsi <- t(atac_lsi)
    cat("  transposed -> cell x 30\n")
  } else {
    stop("expected 30-dim LSI; got dim ",
         paste(dim(atac_lsi), collapse = " x "))
  }
}
atac_z_block <- scale(atac_lsi)
attributes(atac_z_block)$`scaled:center` <- NULL
attributes(atac_z_block)$`scaled:scale`  <- NULL
saveRDS(atac_z_block, file.path(OUT_DIR, "atac_z_block_v2.rds"))
cat("  saved atac_z_block_v2.rds (", paste(dim(atac_z_block), collapse = " x "), ")\n\n")

# ---- 2. load RNA counts + .gs gene sets ----
cat("[2] loading RNA counts + meta + gs ...\n")
rna_counts <- readRDS(RNA_COUNTS)
rna_feats  <- readRDS(RNA_FEATS)
meta       <- readRDS(META)
# orient to gene x cell
if (nrow(rna_counts) == nrow(meta) && ncol(rna_counts) != nrow(meta)) {
  rna_counts <- t(rna_counts)
}
if (is.null(rownames(rna_counts))) {
  rownames(rna_counts) <- as.character(rna_feats)
}
if (is.null(colnames(rna_counts)) && !is.null(rownames(meta))) {
  colnames(rna_counts) <- rownames(meta)
}
cat("  rna_counts: ", paste(dim(rna_counts), collapse = " x "), " (gene x cell)\n")

# parse .gs file (TSV: TRAIT, GENESET; GENESET is "GENE:weight,GENE:weight,...")
gs_lines <- readLines(GS_FILE)
trait_genes <- list()
for (i in 2:length(gs_lines)) {
  parts <- strsplit(gs_lines[i], "\t", fixed = TRUE)[[1]]
  trait <- parts[1]
  pairs <- strsplit(strsplit(parts[2], ",", fixed = TRUE)[[1]], ":", fixed = TRUE)
  trait_genes[[trait]] <- list(
    genes   = vapply(pairs, `[`, character(1), 1),
    weights = as.numeric(vapply(pairs, function(x) if (length(x) >= 2) x[2] else "1",
                                character(1)))
  )
}
cat("  loaded gs traits: ", paste(names(trait_genes), collapse = ", "), "\n\n")

# ---- 3. align cell barcodes (SCDRS dropped some during filter) ----
cat("[3] aligning cell barcodes ...\n")
sample_score <- fread(file.path(SCORE_DIR, paste0(TRAITS["MONO"], ".score.gz")))
score_cells  <- sample_score[[1]]
common <- Reduce(intersect, list(rownames(atac_z_block),
                                  score_cells,
                                  colnames(rna_counts)))
cat("  atac cells:   ", nrow(atac_z_block), "\n")
cat("  scdrs cells:  ", length(score_cells), "\n")
cat("  rna cells:    ", ncol(rna_counts), "\n")
cat("  common cells: ", length(common), "\n")
atac_z_block <- atac_z_block[common, , drop = FALSE]
rna_counts   <- rna_counts[, common, drop = FALSE]
cat("  aligned\n\n")

# ---- 4. precompute size factor (use full-gene total) ----
cell_sums <- colSums(rna_counts)
cat("[4] size factor median: ", median(cell_sums), "\n\n")

# ---- 5. per-trait: RNA PCA + auto-tune ----
fusion_log <- list()
for (tname in names(TRAITS)) {
  trait <- TRAITS[[tname]]
  cat("=== Trait ", tname, " (", trait, ") ===\n", sep = "")

  # 5a. SCDRS norm_score
  cat("  [5a] reading SCDRS score ...\n")
  dt <- fread(file.path(SCORE_DIR, paste0(trait, ".score.gz")))
  setnames(dt, 1, "barcode")
  setkey(dt, barcode)
  norm_score <- dt[common, norm_score]
  names(norm_score) <- common
  stopifnot(!any(is.na(norm_score)))
  cat("    norm_score range: [",
      sprintf("%.2f", min(norm_score)), ", ",
      sprintf("%.2f", max(norm_score)), "]\n", sep = "")

  # 5b. label top 1%
  thr <- as.numeric(quantile(norm_score, 1 - LABEL_TOP_PCT))
  label <- ifelse(norm_score >= thr, "A", "B")
  names(label) <- common
  cat("    label A (top 1%): ", sum(label == "A"),
      "; B: ", sum(label == "B"), "\n")

  # 5c. RNA PCA on SCDRS gene set subset
  trait_g  <- trait_genes[[trait]]$genes
  gene_use <- intersect(trait_g, rownames(rna_counts))
  cat("  [5b] RNA PCA: gs gene ", length(trait_g),
      " -> intersect rna ", length(gene_use), "\n")
  rna_sub      <- rna_counts[gene_use, , drop = FALSE]
  rna_sub_norm <- t(t(rna_sub) / cell_sums) * 1e4
  rna_sub_log  <- log1p(rna_sub_norm)
  cat("    irlba PCA (", length(gene_use), " genes x ", length(common), " cells) ...\n", sep = "")
  t0 <- Sys.time()
  pca <- irlba::irlba(t(rna_sub_log), nv = RNA_PCA_MAX_K, center = TRUE)
  cat("    elapsed: ",
      round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "s\n", sep = "")

  cumvar    <- cumsum(pca$d^2) / sum(pca$d^2)
  k_cumvar  <- which(cumvar >= 0.80)[1]
  K_rna     <- min(max(k_cumvar, RNA_PCA_MIN_K), RNA_PCA_MAX_K)
  cat("    K_rna chosen: ", K_rna,
      " (cumvar80 = ", k_cumvar, ", clipped to [", RNA_PCA_MIN_K, ", ",
      RNA_PCA_MAX_K, "])\n", sep = "")

  rna_emb <- pca$u[, 1:K_rna, drop = FALSE] %*%
             diag(pca$d[1:K_rna], nrow = K_rna)
  rownames(rna_emb) <- common
  rna_z_block <- scale(rna_emb)
  attributes(rna_z_block)$`scaled:center` <- NULL
  attributes(rna_z_block)$`scaled:scale`  <- NULL
  saveRDS(rna_z_block,
          file.path(OUT_DIR, paste0("rna_z_block_v2_", tname, ".rds")))
  cat("    saved rna_z_block_v2_", tname, ".rds (",
      paste(dim(rna_z_block), collapse = " x "), ")\n", sep = "")

  # 5d. auto-tune atac_w
  cat("  [5c] auto-tune atac_w over ",
      paste(ATAC_W_GRID, collapse = "/"), " ...\n", sep = "")
  res <- data.frame(trait = tname, atac_w = ATAC_W_GRID,
                    sil_mean = NA_real_, sil_A_mean = NA_real_)
  for (j in seq_along(ATAC_W_GRID)) {
    w <- ATAC_W_GRID[j]
    fusion <- cbind(sqrt(w)     * atac_z_block,
                    sqrt(1 - w) * rna_z_block)
    sil <- bluster::approxSilhouette(fusion, clusters = label[rownames(fusion)])
    res$sil_mean[j]   <- mean(sil$width)
    res$sil_A_mean[j] <- mean(sil$width[label[rownames(fusion)] == "A"])
    cat("    atac_w=", w,
        " sil_mean=", sprintf("%.4f", res$sil_mean[j]),
        " sil_A_mean=", sprintf("%.4f", res$sil_A_mean[j]), "\n", sep = "")
  }

  best_j <- which.max(res$sil_A_mean)
  best_w <- ATAC_W_GRID[best_j]
  cat("    best atac_w = ", best_w,
      " (sil_A_mean = ", sprintf("%.4f", res$sil_A_mean[best_j]), ")\n", sep = "")

  fusion_best <- cbind(sqrt(best_w)     * atac_z_block,
                       sqrt(1 - best_w) * rna_z_block)
  rownames(fusion_best) <- common
  attr(fusion_best, "atac_w") <- best_w
  attr(fusion_best, "K_rna")  <- K_rna
  saveRDS(fusion_best,
          file.path(OUT_DIR, paste0("fusion_v2_", tname, ".rds")))
  cat("    saved fusion_v2_", tname, ".rds (",
      paste(dim(fusion_best), collapse = " x "), ", atac_w=", best_w, ")\n",
      sep = "")

  fusion_log[[tname]] <- res
  cat("\n")
}

# ---- 6. tune log + checkpoint 2 ----
log_df <- do.call(rbind, fusion_log)
write.csv(log_df, file.path(OUT_DIR, "fusion_tune_log_v2.csv"), row.names = FALSE)
cat("[6] saved fusion_tune_log_v2.csv\n\n")
print(log_df)

cat("\n", strrep("=", 72), "\n", sep = "")
cat("CHECKPOINT 2: at least 2 traits with best atac_w in (0.5, 0.9) interior\n")
cat(strrep("=", 72), "\n")
best_ws <- vapply(fusion_log,
                  function(x) x$atac_w[which.max(x$sil_A_mean)],
                  numeric(1))
cat("  per trait best atac_w:\n")
for (n in names(best_ws)) cat("    ", n, " = ", best_ws[n], "\n", sep = "")
cat("  in (0.5, 0.9) interior: ",
    sum(best_ws > 0.5 & best_ws < 0.9), "/3 traits\n")

cat("\n02_fusion_autotune.R OK.\n")
