# module1_method_B_bindSC.R
# 模块 1 · 方法 B：用 bindSC 官方 R 包做 scATAC–scRNA 跨模态整合
#
# 范式：Bi-order canonical correlation analysis
# 学术界成熟方法，不自写算法。
#   - R 包 bindSC (KChen-lab):  https://github.com/KChen-lab/bindSC
#   - 论文: Dou et al. 2022 Genome Biology 23:196
#           https://doi.org/10.1186/s13059-022-02679-x
#
# 运行方式（干净 R session）：
#   1. 先按装包步骤装好 bindSC 及依赖（见项目文档 13）
#   2. Session -> Restart R
#   3. source("D:/RStudio/work/module1_method_B_bindSC.R")
#
# 输入与方法 C 完全一致：
#   gene_activity_matrix.rds + seurat_rna_normalized.rds + atac_metadata.rds
# 输出遵循接口契约，与方法 A / C 可直接对比，可作为模块二输入。

# =========================
# 0) 依赖自检（只检查不安装）
# =========================
required_pkgs <- c("bindSC", "Seurat", "Matrix", "data.table", "dplyr",
                   "ggplot2", "RANN", "irlba", "umap", "jsonlite")
missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("缺少包: ", paste(missing_pkgs, collapse = ", "),
       "\n先在干净 R session 里把装包命令跑完。")
}

suppressPackageStartupMessages({
  library(bindSC)
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(RANN)
  library(irlba)
  library(umap)
  library(jsonlite)
})

cat("bindSC 版本:", as.character(packageVersion("bindSC")), "\n")

set.seed(9527)
script_start <- Sys.time()
timings <- list()

# =========================
# 1) 路径
# =========================
data_dir <- "D:/SCAVENGE-reproducibility-main/blood_data"
in_dir   <- file.path(data_dir, "multiomics_output_v2")
out_dir  <- file.path(data_dir, "multiomics_module1_pairing", "method_B_bindSC")
tmp_dir  <- file.path(out_dir, "bicca_tmp")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)

cat("========== 路径配置 ==========\n")
cat("  输入:", in_dir, "\n")
cat("  输出:", out_dir, "\n")

# =========================
# 2) 读入已有产物
# =========================
cat("\n========== 加载数据 ==========\n")
t0 <- Sys.time()

ga_mat     <- readRDS(file.path(in_dir, "gene_activity_matrix.rds"))
seurat_rna <- readRDS(file.path(in_dir, "seurat_rna_normalized.rds"))
rna_counts <- GetAssayData(seurat_rna, assay = "RNA", layer = "counts")
atac_meta  <- readRDS(file.path(in_dir, "atac_metadata.rds"))
rna_meta   <- as.data.frame(seurat_rna@meta.data)

cat("  Gene Activity:", nrow(ga_mat),     "genes x", ncol(ga_mat),     "cells\n")
cat("  scRNA counts :", nrow(rna_counts), "genes x", ncol(rna_counts), "cells\n")

timings$load <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("  [计时] 加载: %.1f 秒\n", timings$load))

# =========================
# 3) 共享基因子集
# =========================
shared_genes <- intersect(rownames(rna_counts), rownames(ga_mat))
cat("\n  共享基因数:", length(shared_genes), "\n")
rna_counts <- rna_counts[shared_genes, ]
ga_mat     <- ga_mat[shared_genes, ]

# =========================
# 3b) 子采样开关（通路验证用）
#     - SUBSAMPLE_FOR_TEST = TRUE  : 只取 SUBSAMPLE_N 个细胞，5 分钟内出结果
#     - SUBSAMPLE_FOR_TEST = FALSE : 全量 70k 细胞
#     产物会落到 method_B_bindSC/test_subsample_<N>/ 子目录，不污染正式产物目录
# =========================
SUBSAMPLE_FOR_TEST <- TRUE
SUBSAMPLE_N        <- 8000

if (SUBSAMPLE_FOR_TEST) {
  cat("\n[子采样模式] 每个模态随机取", SUBSAMPLE_N, "个细胞\n")
  set.seed(1)
  rna_keep <- sample.int(ncol(rna_counts), min(SUBSAMPLE_N, ncol(rna_counts)))
  set.seed(2)
  ga_keep  <- sample.int(ncol(ga_mat),    min(SUBSAMPLE_N, ncol(ga_mat)))
  rna_counts <- rna_counts[, rna_keep]
  ga_mat     <- ga_mat[,    ga_keep]
  out_dir <- file.path(out_dir, paste0("test_subsample_", SUBSAMPLE_N))
  tmp_dir <- file.path(out_dir, "bicca_tmp")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
  cat("  子采样后 RNA:", ncol(rna_counts), "  ATAC:", ncol(ga_mat), "\n")
  cat("  子采样输出目录:", out_dir, "\n")
}

# =========================
# 4) 归一化 + HVG
#    - RNA: LogNormalize + 2000 HVG
#    - Gene Activity: LogNormalize，复用 RNA 的 HVG（保持基因对齐）
# =========================
cat("\n========== 归一化 + 选 HVG ==========\n")
t0 <- Sys.time()

rna_seurat_tmp <- CreateSeuratObject(counts = rna_counts)
rna_seurat_tmp <- NormalizeData(rna_seurat_tmp, verbose = FALSE)
rna_seurat_tmp <- FindVariableFeatures(rna_seurat_tmp, nfeatures = 2000,
                                       verbose = FALSE)
hvg      <- VariableFeatures(rna_seurat_tmp)
rna_norm <- GetAssayData(rna_seurat_tmp, assay = "RNA", layer = "data")[hvg, ]
rm(rna_seurat_tmp); gc()

ga_seurat_tmp <- CreateSeuratObject(counts = ga_mat)
ga_seurat_tmp <- NormalizeData(ga_seurat_tmp, verbose = FALSE)
ga_norm       <- GetAssayData(ga_seurat_tmp, assay = "RNA", layer = "data")[hvg, ]
rm(ga_seurat_tmp); gc()

cat("  HVG 数      :", length(hvg), "\n")
cat("  rna_norm    :", nrow(rna_norm), "x", ncol(rna_norm), "\n")
cat("  ga_norm     :", nrow(ga_norm),  "x", ncol(ga_norm),  "\n")

# bindSC BiCCA 需要 dense 矩阵（HVG 子集下 ~560MB/侧，可承受）
rna_norm <- as.matrix(rna_norm)
ga_norm  <- as.matrix(ga_norm)

timings$preprocess <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("  [计时] 预处理: %.1f 秒\n", timings$preprocess))

# =========================
# 5) 初始聚类（BiCCA 要求两侧各给 cluster 标签）
#    - 不用 BioClassification（那是评估金标准，不能泄漏给训练）
#    - 各自做快速 PCA + kmeans，k=20 与 latent 维度对齐
# =========================
cat("\n========== 初始聚类 (irlba PCA + kmeans) ==========\n")
t0 <- Sys.time()

K_LATENT <- 20
N_CLST   <- 20

pca_rna <- irlba::irlba(t(rna_norm), nv = K_LATENT)$u
pca_ga  <- irlba::irlba(t(ga_norm),  nv = K_LATENT)$u
rownames(pca_rna) <- colnames(rna_norm)
rownames(pca_ga)  <- colnames(ga_norm)

km_rna <- kmeans(pca_rna, centers = N_CLST, nstart = 10, iter.max = 50)
km_ga  <- kmeans(pca_ga,  centers = N_CLST, nstart = 10, iter.max = 50)

X.clst <- as.character(km_rna$cluster)
Y.clst <- as.character(km_ga$cluster)

cat("  RNA cluster 数:", length(unique(X.clst)), "\n")
cat("  ATAC cluster 数:", length(unique(Y.clst)), "\n")

timings$initial_clust <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("  [计时] 初始聚类: %.1f 秒\n", timings$initial_clust))

# =========================
# 6) BiCCA 整合
#    无 peak 数据，Y 与 Z0 均用 Gene Activity（与方法 C 输入口径一致）
# =========================
cat("\n========== BiCCA 整合 (K=20, alpha=0.5, lambda=0.5) ==========\n")
cat("  预计耗时 10-30 分钟...\n")
t0 <- Sys.time()

res <- bindSC::BiCCA(
  X             = rna_norm,
  Y             = ga_norm,
  Z0            = ga_norm,
  X.clst        = X.clst,
  Y.clst        = Y.clst,
  alpha         = 0.5,
  lambda        = 0.5,
  K             = K_LATENT,
  num.iteration = 20,     # 从 50 降到 20：BiCCA 一般 10–15 轮收敛
  tolerance     = 0.05,   # 从 0.01 放松：允许提前早停
  temp.path     = tmp_dir
)

timings$bicca <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("  [计时] BiCCA: %.1f 秒\n", timings$bicca))

# =========================
# 7) 共享嵌入
#    bindSC 2.x 返回 res$u (RNA cell x K) / res$r (ATAC cell x K)
#    防御性兼容：如果字段名不同，从 res 里挑形状匹配的矩阵
# =========================
cat("\n========== 共享嵌入提取 ==========\n")

pick_latent <- function(res, expected_nrow, expected_ncol) {
  for (nm in names(res)) {
    mt <- res[[nm]]
    if (is.matrix(mt) && nrow(mt) == expected_nrow && ncol(mt) == expected_ncol) {
      return(list(name = nm, mat = mt))
    }
  }
  NULL
}

u_rna  <- if (!is.null(res$u)) res$u else NULL
r_atac <- if (!is.null(res$r)) res$r else NULL

if (is.null(u_rna) || is.null(r_atac)) {
  cat("[诊断] res 字段:\n"); print(names(res))
  if (is.null(u_rna)) {
    hit <- pick_latent(res, ncol(rna_norm), K_LATENT)
    if (!is.null(hit)) { u_rna <- hit$mat; cat("  兜底选中 RNA latent:", hit$name, "\n") }
  }
  if (is.null(r_atac)) {
    hit <- pick_latent(res, ncol(ga_norm), K_LATENT)
    if (!is.null(hit)) { r_atac <- hit$mat; cat("  兜底选中 ATAC latent:", hit$name, "\n") }
  }
}
if (is.null(u_rna) || is.null(r_atac)) {
  stop("BiCCA 返回的 latent 无法识别。把 names(res) 与各字段 dim 贴给我。")
}

rownames(u_rna)  <- colnames(rna_norm)
rownames(r_atac) <- colnames(ga_norm)

cat("  u_rna  (RNA latent):", nrow(u_rna),  "x", ncol(u_rna),  "\n")
cat("  r_atac (ATAC latent):", nrow(r_atac), "x", ncol(r_atac), "\n")

H_shared     <- rbind(u_rna, r_atac)
modality_vec <- c(rep("rna", nrow(u_rna)), rep("atac", nrow(r_atac)))
stopifnot(nrow(H_shared) == (ncol(rna_counts) + ncol(ga_mat)))

# =========================
# 8) ATAC -> top-k RNA 邻居（接口契约：atac_topk_rna_neighbors.rds）
# =========================
cat("\n========== top-k RNA 邻居查询 ==========\n")
t0 <- Sys.time()

K_TOPK <- 30
nn_res <- RANN::nn2(data = u_rna, query = r_atac, k = K_TOPK)

rna_cell_ids  <- rownames(u_rna)
atac_cell_ids <- rownames(r_atac)

atac_topk_rna <- lapply(seq_along(atac_cell_ids), function(i) {
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
# 9) Label Transfer（top-k 多数投票，口径与方法 C 完全一致）
# =========================
cat("\n========== Label Transfer (投票) ==========\n")
t0 <- Sys.time()

rna_labels_by_id  <- setNames(as.character(rna_meta$BioClassification),
                              rownames(rna_meta))
atac_labels_by_id <- setNames(as.character(atac_meta$BioClassification),
                              rownames(atac_meta))

predict_by_vote <- function(rna_cells) {
  lbls <- rna_labels_by_id[rna_cells]
  lbls <- lbls[!is.na(lbls)]
  if (length(lbls) == 0) return(c(predicted = NA_character_, score = NA_real_))
  tab <- sort(table(lbls), decreasing = TRUE)
  c(predicted = names(tab)[1], score = as.numeric(tab[1] / length(lbls)))
}

vote_mat <- t(sapply(atac_topk_rna, function(df) predict_by_vote(df$rna_cell)))

label_df <- data.frame(
  atac_cell      = atac_cell_ids,
  original.label = atac_labels_by_id[atac_cell_ids],
  predicted.id   = vote_mat[, "predicted"],
  score          = as.numeric(vote_mat[, "score"]),
  stringsAsFactors = FALSE
)

fwrite(label_df, file.path(out_dir, "atac_label_transfer.csv"))
cat("  已保存: atac_label_transfer.csv\n")

timings$label_transfer <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

# =========================
# 10) 一致率评估：overall / highconf (>=0.5) / strict (=1.0)
# =========================
cat("\n========== 一致率评估 ==========\n")

label_df_ok <- label_df[!is.na(label_df$original.label) &
                         !is.na(label_df$predicted.id), ]

conf_mat <- table(Original    = label_df_ok$original.label,
                  Transferred = label_df_ok$predicted.id)

orig_types    <- rownames(conf_mat)
trans_types   <- colnames(conf_mat)
shared_labels <- intersect(orig_types, trans_types)

agreement_per_type <- sapply(orig_types, function(ct) {
  if (ct %in% trans_types) conf_mat[ct, ct] / sum(conf_mat[ct, ]) else 0
})

overall_agreement <- sum(sapply(shared_labels,
                                function(ct) conf_mat[ct, ct])) / sum(conf_mat)
cat("  总体一致率:", round(overall_agreement * 100, 1), "%\n")

agreement_df <- data.frame(
  celltype  = orig_types,
  n_cells   = as.integer(rowSums(conf_mat)),
  agreement = round(agreement_per_type * 100, 1)
) %>% arrange(desc(agreement))
fwrite(agreement_df, file.path(out_dir, "label_transfer_agreement.csv"))
cat("  已保存: label_transfer_agreement.csv\n")

compute_subset_agreement <- function(mask, fname_suffix) {
  if (sum(mask) == 0) return(NA_real_)
  conf_mat_s <- table(Original    = label_df_ok$original.label[mask],
                      Transferred = label_df_ok$predicted.id[mask])
  orig_s  <- rownames(conf_mat_s)
  trans_s <- colnames(conf_mat_s)
  agree_s <- sapply(orig_s, function(ct) {
    if (ct %in% trans_s) conf_mat_s[ct, ct] / sum(conf_mat_s[ct, ]) else 0
  })
  agree_df <- data.frame(
    celltype   = orig_s,
    n_cells    = as.integer(rowSums(conf_mat_s)),
    agreement  = round(agree_s * 100, 1)
  ) %>% arrange(desc(agreement))
  fwrite(agree_df, file.path(out_dir,
                             paste0("label_transfer_agreement_", fname_suffix, ".csv")))
  cat("  已保存: label_transfer_agreement_", fname_suffix, ".csv\n", sep = "")
  mean(label_df_ok$original.label[mask] == label_df_ok$predicted.id[mask])
}

hc_mask <- !is.na(label_df_ok$score) & label_df_ok$score >= 0.5
cat("  高置信度 (>= 0.5):", sum(hc_mask), "/", nrow(label_df_ok), "\n")
hc_overall <- compute_subset_agreement(hc_mask, "highconf")
if (!is.na(hc_overall)) cat("  highconf 总体一致率:", round(hc_overall * 100, 1), "%\n")

strict_mask <- !is.na(label_df_ok$score) & label_df_ok$score >= 0.999
cat("  strict (score = 1):", sum(strict_mask), "/", nrow(label_df_ok), "\n")
strict_overall <- compute_subset_agreement(strict_mask, "strict")
if (!is.na(strict_overall)) cat("  strict 总体一致率:", round(strict_overall * 100, 1), "%\n")

# =========================
# 11) 保存共享嵌入（接口契约文件）
# =========================
cat("\n========== 保存共享嵌入 ==========\n")

cells_all <- c(rownames(u_rna), rownames(r_atac))
shared_embed_df <- data.frame(
  cell     = cells_all,
  modality = modality_vec,
  H_shared,
  check.names = FALSE
)
saveRDS(shared_embed_df, file.path(out_dir, "shared_embedding.rds"),
        compress = FALSE)
cat("  已保存: shared_embedding.rds (", nrow(shared_embed_df), "x",
    ncol(shared_embed_df), ")\n")

saveRDS(H_shared, file.path(out_dir, "shared_embedding_Hnorm.rds"),
        compress = FALSE)
cat("  已保存: shared_embedding_Hnorm.rds (与方法 C 命名对齐)\n")

# =========================
# 12) UMAP 可视化（20k sample，全量 70k 走 umap 包太慢）
# =========================
cat("\n========== UMAP 可视化 ==========\n")
t0 <- Sys.time()

n_total  <- nrow(H_shared)
n_sample <- min(20000, n_total)
set.seed(42)
idx_vis <- sort(sample.int(n_total, n_sample))

umap_out <- umap::umap(H_shared[idx_vis, ], n_neighbors = 30, min_dist = 0.3)
umap_df  <- data.frame(
  UMAP1    = umap_out$layout[, 1],
  UMAP2    = umap_out$layout[, 2],
  modality = modality_vec[idx_vis],
  cell     = cells_all[idx_vis],
  stringsAsFactors = FALSE
)

umap_df$celltype <- NA_character_
rna_idx_vis  <- umap_df$modality == "rna"
atac_idx_vis <- umap_df$modality == "atac"
umap_df$celltype[rna_idx_vis]  <- rna_labels_by_id[umap_df$cell[rna_idx_vis]]
umap_df$celltype[atac_idx_vis] <- label_df$predicted.id[
  match(umap_df$cell[atac_idx_vis], label_df$atac_cell)
]

p_mod <- ggplot(umap_df, aes(UMAP1, UMAP2, color = modality)) +
  geom_point(size = 0.3, alpha = 0.5) +
  scale_color_manual(values = c("rna" = "#4DBBD5", "atac" = "#E64B35")) +
  labs(title = "bindSC Co-embedding (BiCCA latent + UMAP)",
       subtitle = paste0("n_sample=", n_sample, " / ", n_total)) +
  theme_classic()
ggsave(file.path(out_dir, "Fig_bindSC_umap_modality.png"),
       plot = p_mod, width = 8, height = 6, dpi = 300)

p_ct <- ggplot(umap_df, aes(UMAP1, UMAP2, color = celltype)) +
  geom_point(size = 0.3, alpha = 0.5) +
  labs(title = "bindSC Co-embedding: Cell Type") +
  theme_classic() + theme(legend.position = "none")
ggsave(file.path(out_dir, "Fig_bindSC_umap_celltype.png"),
       plot = p_ct, width = 8, height = 6, dpi = 300)

p_score <- ggplot(label_df, aes(score)) +
  geom_histogram(bins = 50, fill = "#3C5488", color = "white") +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "red") +
  labs(title = "bindSC: top-k Vote Score Distribution",
       x = "Max-vote proportion", y = "Cells") +
  theme_classic()
ggsave(file.path(out_dir, "Fig_bindSC_prediction_score.png"),
       plot = p_score, width = 6, height = 4, dpi = 300)

cat("  已保存 3 张图 (UMAP on", n_sample, "sample)\n")

timings$umap <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

# =========================
# 13) 保存 BiCCA 对象
# =========================
cat("\n========== 保存 BiCCA 对象 ==========\n")
saveRDS(res, file.path(out_dir, "bindSC_BiCCA_result.rds"), compress = FALSE)
cat("  已保存: bindSC_BiCCA_result.rds\n")

# =========================
# 14) timing.json
# =========================
timings$total          <- as.numeric(difftime(Sys.time(), script_start, units = "secs"))
timings$method         <- "bindSC (BiCCA, K=20, alpha=0.5, lambda=0.5, num.iter=20, tol=0.05)"
timings$citation       <- "Dou et al. 2022 Genome Biology 23:196"
timings$bindSC_version <- as.character(packageVersion("bindSC"))
timings$note           <- paste0(
  "No peak-level ATAC data used; Y = Z0 = gene activity (gene-level alignment only).",
  " HVG=2000 from RNA. Initial clusters via kmeans(k=20) on each modality's irlba PCA.",
  if (SUBSAMPLE_FOR_TEST) paste0(
    " [SUBSAMPLE test run: ", SUBSAMPLE_N, " cells/modality; results NOT for final comparison.]"
  ) else " [FULL run, all cells.]"
)
timings$subsample      <- SUBSAMPLE_FOR_TEST
timings$subsample_n    <- if (SUBSAMPLE_FOR_TEST) SUBSAMPLE_N else NA_integer_
jsonlite::write_json(timings, file.path(out_dir, "timing.json"),
                     pretty = TRUE, auto_unbox = TRUE)

# =========================
# 15) 总结
# =========================
cat("\n")
cat("======================================================\n")
cat("  模块 1 · 方法 B (bindSC) 完成\n")
cat("======================================================\n")
cat("  bindSC 版本    :", as.character(packageVersion("bindSC")), "\n")
cat("  共享嵌入维度   :", dim(H_shared)[1], "x", dim(H_shared)[2], "\n")
cat("  整体一致率     :", round(overall_agreement * 100, 1), "%\n")
if (!is.na(hc_overall))     cat("  highconf 一致率:", round(hc_overall * 100, 1), "%\n")
if (!is.na(strict_overall)) cat("  strict 一致率  :", round(strict_overall * 100, 1), "%\n")
cat("  总耗时         :", round(timings$total / 60, 1), "分钟\n")
cat("  输出目录       :", out_dir, "\n")
cat("======================================================\n")
