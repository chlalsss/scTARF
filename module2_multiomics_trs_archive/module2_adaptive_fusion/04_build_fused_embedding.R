# Build fused ATAC + RNA-reference embedding and reduce by PCA.
#
# Fusion: cbind(atac_block, rna_block) -> centered PCA. MODULE2_FUSION_MODE:
#   - norm_equal_zscore (default): per-dim z-score each block, then norm-equalize
#     Frobenius. Removes per-dim variance asymmetry caused by upstream PCA outputs
#     (e.g. disease latent has variance concentrated in early PCs).
#   - norm_equal: norm-equalize Frobenius without z-score. Backward-compat.
#   - raw: no scaling, fusion_w = 1.

source("D:/RStudio/work/module2_adaptive_fusion/_bootstrap.R")

required_pkgs(c("data.table", "jsonlite"))
suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

paths <- module2_paths()
ensure_module2_dirs(paths)
label <- get_active_label()
cutoff <- if (label == "top1" || startsWith(label, "top1disease")) NA_real_ else get_single_cutoff()
npcs_raw <- Sys.getenv("MODULE2_FUSED_NPCS", "auto")
npcs_auto <- tolower(npcs_raw) == "auto"
npcs_max <- as.integer(Sys.getenv("MODULE2_FUSED_NPCS_MAX", "80"))
npcs <- if (npcs_auto) npcs_max else as.integer(npcs_raw)
if (!is.finite(npcs) || npcs < 1) stop("MODULE2_FUSED_NPCS must be 'auto' or a positive integer")
fusion_mode <- tolower(Sys.getenv("MODULE2_FUSION_MODE", "norm_equal"))

rna_ref_rds <- file.path(paths$out_dir, sprintf("rna_reference_latent_cutoff_%s.rds", label))
if (!file.exists(rna_ref_rds)) {
  stop("Missing RNA reference latent. Run 03_build_rna_reference.R first: ", rna_ref_rds)
}

cat("========== Module 2: build fused embedding ==========\n")
cat("  cutoff       :", cutoff, "\n")
cat("  PCs          :", if (npcs_auto) paste0("auto(max=", npcs, ")") else npcs, "\n")
cat("  fusion mode  :", fusion_mode, "\n")

loaded <- load_shared_embedding(paths$shared_embedding)
emb <- as.data.table(loaded$emb)
latent_cols <- loaded$latent_cols
atac_dt <- emb[modality == "atac"]
rna_ref <- as.data.table(readRDS(rna_ref_rds))
rna_ref_cols <- grep("^rna_ref_", colnames(rna_ref), value = TRUE)

if (!all(atac_dt$cell == rna_ref$atac_cell)) {
  rna_ref <- rna_ref[match(atac_dt$cell, rna_ref$atac_cell)]
}
if (!all(atac_dt$cell == rna_ref$atac_cell)) {
  stop("ATAC cell order mismatch between shared embedding and RNA reference latent")
}

atac_mat <- as.matrix(atac_dt[, ..latent_cols])
rna_ref_mat <- as.matrix(rna_ref[, ..rna_ref_cols])
storage.mode(atac_mat) <- "double"
storage.mode(rna_ref_mat) <- "double"

zscore_block <- function(m) {
  cm <- colMeans(m)
  cs <- apply(m, 2, sd)
  cs[!is.finite(cs) | cs < 1e-8] <- 1
  out <- sweep(sweep(m, 2, cm), 2, cs, "/")
  out[!is.finite(out)] <- 0
  out
}

if (fusion_mode == "norm_equal_zscore") {
  atac_block <- zscore_block(atac_mat)
  rna_block_pre <- zscore_block(rna_ref_mat)
  atac_fro_pre <- sqrt(sum(atac_block^2))
  rna_fro_pre <- sqrt(sum(rna_block_pre^2))
  fusion_w <- if (rna_fro_pre > 1e-12) atac_fro_pre / rna_fro_pre else 0
  rna_block <- rna_block_pre * fusion_w
} else if (fusion_mode == "norm_equal") {
  atac_block <- atac_mat
  atac_fro_pre <- sqrt(sum(atac_block^2))
  rna_fro_pre <- sqrt(sum(rna_ref_mat^2))
  fusion_w <- if (rna_fro_pre > 1e-12) atac_fro_pre / rna_fro_pre else 0
  rna_block <- rna_ref_mat * fusion_w
} else if (fusion_mode == "raw") {
  atac_block <- atac_mat
  rna_block <- rna_ref_mat
  fusion_w <- 1.0
} else {
  stop("Unknown MODULE2_FUSION_MODE: ", fusion_mode)
}

fused_mat <- cbind(atac_block, rna_block)
fused_mat[!is.finite(fused_mat)] <- 0

atac_fro <- sqrt(sum(atac_block^2))
rna_fro <- sqrt(sum(rna_block^2))
cat(sprintf("  pre-PCA dims : %d (atac %d + rna %d)\n",
            ncol(fused_mat), ncol(atac_block), ncol(rna_block)))
cat(sprintf("  block ||F||  : ATAC %.2f, RNA*w %.2f, ratio (RNA/ATAC) %.3f\n",
            atac_fro, rna_fro, rna_fro / max(atac_fro, 1e-12)))

pca_n <- min(npcs, ncol(fused_mat) - 1L, nrow(fused_mat) - 1L)
if (pca_n < 1) stop("Not enough rows/columns for PCA")
if (requireNamespace("irlba", quietly = TRUE)) {
  pca <- irlba::prcomp_irlba(fused_mat, n = pca_n, center = TRUE, scale. = FALSE)
} else {
  warning("Package irlba is not installed; falling back to stats::prcomp.")
  pca <- stats::prcomp(fused_mat, rank. = pca_n, center = TRUE, scale. = FALSE)
}

npcs_use <- if (npcs_auto) auto_select_pca_dim(pca$sdev) else pca_n
cat("  selected PCs :", npcs_use, "\n")

pc_mat <- pca$x[, seq_len(min(npcs_use, ncol(pca$x))), drop = FALSE]
colnames(pc_mat) <- paste0("fusedPC_", seq_len(ncol(pc_mat)))
fused_df <- data.frame(atac_cell = atac_dt$cell, pc_mat, check.names = FALSE)

out_rds <- file.path(paths$out_dir, sprintf("fused_embedding_cutoff_%s_pca%d.rds", label, npcs_use))
saveRDS(fused_df, out_rds, compress = FALSE)
cat("  saved:", out_rds, "\n")

npcs_auto_file <- file.path(paths$out_dir, sprintf("fused_embedding_auto_npcs_cutoff_%s.txt", label))
write_auto_value(npcs_auto_file, npcs_use)

variance <- if (!is.null(pca$sdev)) {
  pca$sdev^2 / sum(pca$sdev^2)
} else {
  rep(NA_real_, ncol(pc_mat))
}
variance_dt <- data.table(
  pc = seq_len(length(variance)),
  variance_fraction = variance
)
variance_out <- file.path(paths$out_dir, sprintf("fused_embedding_variance_cutoff_%s_pca%d.csv", label, npcs_use))
fwrite(variance_dt, variance_out)

jsonlite::write_json(
  list(
    script = "04_build_fused_embedding.R",
    cutoff = cutoff,
    npcs = npcs_use,
    npcs_mode = if (npcs_auto) "auto_elbow" else "manual",
    npcs_requested = npcs_raw,
    npcs_auto_file = npcs_auto_file,
    fusion_mode = fusion_mode,
    fusion_weight = fusion_w,
    fusion = sprintf("mode=%s; per-block z-score=%s; centered PCA",
                     fusion_mode,
                     fusion_mode == "norm_equal_zscore"),
    pca_method = if (npcs_auto) {
      paste0(tolower(Sys.getenv("MODULE2_PCA_AUTO_METHOD", "cumvar")),
             "(threshold=",
             Sys.getenv("MODULE2_PCA_CUMVAR", "0.80"), ")")
    } else "manual",
    block_frobenius_atac = atac_fro,
    block_frobenius_rna_weighted = rna_fro,
    rna_to_atac_frobenius_ratio = rna_fro / max(atac_fro, 1e-12),
    rna_reference = rna_ref_rds,
    output_rds = out_rds,
    variance_csv = variance_out,
    created_at = as.character(Sys.time())
  ),
  file.path(paths$out_dir, sprintf("fused_embedding_manifest_cutoff_%s_pca%d.json", label, npcs_use)),
  pretty = TRUE,
  auto_unbox = TRUE
)

cat("========== Done ==========\n")
