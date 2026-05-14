# module1_method_A_migrate.R
# 模块 1 · 方法 A：把已有 Seurat anchor 产物迁移到统一接口契约格式
#
# 输入（已有，不重算）：
#   multiomics_output_v2/seurat_atac_final.rds       (4.3 GB，含 predicted.id + score)
#   multiomics_output_v2/seurat_rna_normalized.rds   (236 MB)
#   multiomics_output_v2/gene_activity_matrix.rds    (470 MB)
#   multiomics_output_v2/label_transfer_agreement.csv / _highconf.csv
#
# 输出（接口契约）：
#   method_A_seurat/atac_label_transfer.csv         (35,038 ATAC cells, 全量)
#   method_A_seurat/shared_embedding.rds            (70,620 cells × 30 PCs)
#   method_A_seurat/atac_topk_rna_neighbors.rds     (35,038 × top-30)
#   method_A_seurat/label_transfer_agreement.csv    (直接 copy)
#   method_A_seurat/label_transfer_agreement_highconf.csv  (直接 copy)
#   method_A_seurat/timing.json
#
# 为什么不直接用现成的 seurat_coembedding.rds？
#   那个 co-embedding 只采样了 20k 细胞做 Harmony 校正（内存原因），
#   不是全量 70k。为了与方法 C 公平对比，这里重新在全量 70k 上做
#   "RNA PCA + ATAC 投影"，不跑 Harmony（RNA 自身 PCA 就是共同参考系）。
#
# 运行方式：source("D:/RStudio/work/module1_method_A_migrate.R")

# =========================
# 0) 依赖自检（只检查）
# =========================
required_pkgs <- c("Seurat", "Matrix", "data.table", "RANN", "jsonlite")
missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("缺少包: ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(RANN)
  library(jsonlite)
})

set.seed(9527)
script_start <- Sys.time()
timings <- list(
  method = "Seurat anchor (pcaproject + Label Transfer)",
  citation = "Stuart et al. 2019 Cell 177:1888",
  note = "pcaproject shared embedding re-computed on full 70k cells (not the 20k Harmony sample)"
)

# =========================
# 1) 路径
# =========================
data_dir <- "D:/SCAVENGE-reproducibility-main/blood_data"
in_dir   <- file.path(data_dir, "multiomics_output_v2")
out_dir  <- file.path(data_dir, "multiomics_module1_pairing", "method_A_seurat")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("========== 路径配置 ==========\n")
cat("  输入:", in_dir, "\n")
cat("  输出:", out_dir, "\n")

# =========================
# 2) 抽取 label transfer（仅加载 seurat_atac_final 的 metadata）
# =========================
cat("\n========== 抽取 label transfer（全量 35k）==========\n")
t0 <- Sys.time()

# 4.3 GB 文件，只要 meta.data。加载一次，抽完就扔
seurat_atac_final <- readRDS(file.path(in_dir, "seurat_atac_final.rds"))
atac_meta_full <- seurat_atac_final@meta.data

# 与 method C 统一：加 atac_ 前缀
atac_cells <- rownames(atac_meta_full)
atac_cells_prefixed <- paste0("atac_", atac_cells)

label_df <- data.frame(
  atac_cell      = atac_cells_prefixed,
  original.label = as.character(atac_meta_full$BioClassification),
  predicted.id   = as.character(atac_meta_full$predicted.id),
  score          = as.numeric(atac_meta_full$prediction.score.max),
  stringsAsFactors = FALSE
)
fwrite(label_df, file.path(out_dir, "atac_label_transfer.csv"))
cat("  已保存: atac_label_transfer.csv (", nrow(label_df), "cells)\n")

# seurat_atac_final 就这一步要用，立刻释放 4.3 GB
rm(seurat_atac_final); gc()

timings$extract_labels <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("  [计时] 抽 label: %.1f 秒\n", timings$extract_labels))

# =========================
# 3) RNA 侧：ScaleData + RunPCA（全量 35k RNA 细胞）
# =========================
cat("\n========== RNA PCA（共同参考系）==========\n")
t0 <- Sys.time()

seurat_rna <- readRDS(file.path(in_dir, "seurat_rna_normalized.rds"))
hvg <- VariableFeatures(seurat_rna)
cat("  RNA HVG 数:", length(hvg), "\n")

seurat_rna <- ScaleData(seurat_rna, features = hvg, verbose = FALSE)
seurat_rna <- RunPCA(seurat_rna, features = hvg, npcs = 30, verbose = FALSE)

rna_pca_embed <- Embeddings(seurat_rna, "pca")             # 35,582 × 30
pca_loadings  <- Loadings(seurat_rna, "pca")               # hvg × 30
cat("  rna_pca_embed:", nrow(rna_pca_embed), "×", ncol(rna_pca_embed), "\n")

# 记下 RNA 侧每个 HVG 的 mean/sd（用于 ATAC 侧一致缩放）
rna_data <- GetAssayData(seurat_rna, assay = "RNA", layer = "data")[hvg, ]
gene_mean_rna <- Matrix::rowMeans(rna_data)
gene_mean2    <- Matrix::rowMeans(rna_data^2)
gene_sd_rna   <- sqrt(pmax(gene_mean2 - gene_mean_rna^2, 0))
gene_sd_rna[gene_sd_rna == 0] <- 1

timings$rna_pca <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("  [计时] RNA PCA: %.1f 秒\n", timings$rna_pca))

# =========================
# 4) ATAC 侧：Gene Activity → 投影到 RNA PCA 空间
# =========================
cat("\n========== ATAC 投影到 RNA PCA ==========\n")
t0 <- Sys.time()

ga_mat <- readRDS(file.path(in_dir, "gene_activity_matrix.rds"))
shared_hvg <- intersect(hvg, rownames(ga_mat))
cat("  HVG ∩ ATAC genes:", length(shared_hvg), "\n")

ga_sub <- ga_mat[shared_hvg, ]
# Log-normalize per cell (与 RNA 的 LogNormalize scale.factor = 10000 一致)
ga_norm_sp <- log1p(t(t(ga_sub) / Matrix::colSums(ga_sub) * 10000))

# 用 RNA 的 per-gene mean/sd 做 z-score 缩放（保证两模态在同一度量下）
gene_mean_shared <- gene_mean_rna[shared_hvg]
gene_sd_shared   <- gene_sd_rna[shared_hvg]

# 转成稠密矩阵做中心化（会打破 sparsity）
cat("  转稠密矩阵做 z-score（预计内存 ~1 GB）...\n")
ga_dense    <- as.matrix(ga_norm_sp)
rm(ga_norm_sp, ga_sub, ga_mat); gc()
ga_centered <- sweep(ga_dense, 1, gene_mean_shared, "-")
rm(ga_dense); gc()
ga_scaled   <- sweep(ga_centered, 1, gene_sd_shared, "/")
rm(ga_centered); gc()

# 投影到 RNA PCA: atac_pca = ga_scaled^T %*% pca_loadings[shared_hvg, ]
pca_loadings_shared <- pca_loadings[shared_hvg, , drop = FALSE]
atac_pca_embed <- t(ga_scaled) %*% pca_loadings_shared
rm(ga_scaled); gc()
colnames(atac_pca_embed) <- colnames(rna_pca_embed)
cat("  atac_pca_embed:", nrow(atac_pca_embed), "×", ncol(atac_pca_embed), "\n")

timings$atac_project <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("  [计时] ATAC 投影: %.1f 秒\n", timings$atac_project))

# =========================
# 5) 合并得到 shared_embedding（70k × 30）
# =========================
cat("\n========== shared_embedding ==========\n")

rna_cell_ids  <- paste0("rna_",  rownames(rna_pca_embed))
atac_cell_ids <- paste0("atac_", rownames(atac_pca_embed))
rownames(rna_pca_embed)  <- rna_cell_ids
rownames(atac_pca_embed) <- atac_cell_ids

shared_embed <- rbind(rna_pca_embed, atac_pca_embed)
modality <- c(rep("rna",  nrow(rna_pca_embed)),
              rep("atac", nrow(atac_pca_embed)))

shared_embed_df <- data.frame(
  cell     = rownames(shared_embed),
  modality = modality,
  shared_embed,
  check.names = FALSE
)

saveRDS(shared_embed_df, file.path(out_dir, "shared_embedding.rds"),
        compress = FALSE)
cat("  已保存: shared_embedding.rds (", nrow(shared_embed_df), "×",
    ncol(shared_embed_df), ")\n")

# =========================
# 6) ATAC → RNA top-30 kNN（在 30 维 PCA 空间里）
# =========================
cat("\n========== top-30 RNA 邻居查询 ==========\n")
t0 <- Sys.time()

K_TOPK <- 30
nn_res <- RANN::nn2(data = rna_pca_embed, query = atac_pca_embed, k = K_TOPK)

atac_topk_rna <- lapply(seq_len(nrow(atac_pca_embed)), function(i) {
  data.frame(
    rna_cell = rna_cell_ids[nn_res$nn.idx[i, ]],
    distance = nn_res$nn.dists[i, ]
  )
})
names(atac_topk_rna) <- atac_cell_ids

saveRDS(atac_topk_rna, file.path(out_dir, "atac_topk_rna_neighbors.rds"),
        compress = FALSE)
cat("  已保存: atac_topk_rna_neighbors.rds (k =", K_TOPK, ")\n")

timings$knn <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("  [计时] kNN: %.1f 秒\n", timings$knn))

# =========================
# 7) 复制 agreement csv（已有现成的）
# =========================
cat("\n========== 复制 agreement csv ==========\n")

for (fname in c("label_transfer_agreement.csv",
                "label_transfer_agreement_highconf.csv")) {
  src <- file.path(in_dir, fname)
  dst <- file.path(out_dir, fname)
  if (file.exists(src)) {
    file.copy(src, dst, overwrite = TRUE)
    cat("  已复制:", fname, "\n")
  } else {
    cat("  [缺]", fname, "— 源文件不存在\n")
  }
}

# =========================
# 8) 顺手生成 strict 子集的 agreement（与方法 C 的 strict 对齐）
#    方法 A 的 score 是 prediction.score.max，语义不同，
#    但同样可以扫到最高分位看"真正有把握"子集
# =========================
cat("\n========== strict 子集 agreement（score >= 0.95）==========\n")

# 对方法 A，score 分布更两极化，用 0.95 作为严阈值比较有意义
label_df_ok <- label_df[!is.na(label_df$original.label) &
                         !is.na(label_df$predicted.id) &
                         !is.na(label_df$score), ]
strict_thr <- 0.95
strict_mask <- label_df_ok$score >= strict_thr
cat("  score >=", strict_thr, "细胞数:", sum(strict_mask), "/",
    nrow(label_df_ok), "\n")

if (sum(strict_mask) > 0) {
  conf <- table(Original    = label_df_ok$original.label[strict_mask],
                Transferred = label_df_ok$predicted.id[strict_mask])
  orig_t  <- rownames(conf); trans_t <- colnames(conf)
  agr <- sapply(orig_t, function(ct) {
    if (ct %in% trans_t) conf[ct, ct] / sum(conf[ct, ]) else 0
  })
  strict_df <- data.frame(
    celltype         = orig_t,
    n_cells_strict   = as.integer(rowSums(conf)),
    agreement_strict = round(agr * 100, 1)
  )
  strict_df <- strict_df[order(-strict_df$agreement_strict), ]
  fwrite(strict_df, file.path(out_dir, "label_transfer_agreement_strict.csv"))
  cat("  已保存: label_transfer_agreement_strict.csv (thr =", strict_thr, ")\n")
}

# =========================
# 9) timing.json
# =========================
timings$total <- as.numeric(difftime(Sys.time(), script_start, units = "secs"))
write_json(timings, file.path(out_dir, "timing.json"),
           pretty = TRUE, auto_unbox = TRUE)

# =========================
# 10) 总结
# =========================
cat("\n")
cat("======================================================\n")
cat("  方法 A 迁移完成（接口契约格式）\n")
cat("======================================================\n")
cat("  Label transfer  :", nrow(label_df), "cells\n")
cat("  shared_embedding:", nrow(shared_embed_df), "×",
    ncol(shared_embed_df), "\n")
cat("  top-k (k = 30)  :", length(atac_topk_rna), "ATAC cells\n")
cat("  总耗时          :", round(timings$total / 60, 1), "分钟\n")
cat("  输出目录        :", out_dir, "\n")
cat("======================================================\n")
