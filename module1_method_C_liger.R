# module1_method_C_liger.R
# 模块 1 · 方法 C：用 rliger (LIGER) 官方 API 做 scATAC–scRNA 跨模态整合
#
# 来源（学术界成熟方法，非自写算法）：
#   - R 包 rliger (welch-lab 官方):  https://github.com/welch-lab/liger
#   - 文档:                           https://welch-lab.github.io/liger/
#   - 原始算法论文:
#       * Welch et al. 2019 Cell — iNMF 原文
#       * Gao et al. 2021 Nat Biotechnol — online iNMF
#   - 官方教程: https://welch-lab.github.io/liger/articles/online_iNMF_tutorial.html
#
# 版本适配：
#   本脚本针对 rliger 2.x (>= 1.99) 的 API。旧版参数名（use.raw / dims.use /
#   n_neighbors / min_dist / datasets.use）已全部替换为 2.x 的 camelCase 名。
#
# 运行方式（干净 R session）：
#   1. 先 source("D:/RStudio/work/module1_install_deps.R") 装 rliger + RcppPlanc
#   2. Session → Restart R
#   3. source("D:/RStudio/work/module1_method_C_liger.R")

# =========================
# 0) 依赖自检（只检查不安装）
# =========================
required_pkgs <- c("rliger", "Seurat", "Matrix", "data.table", "dplyr",
                   "ggplot2", "RANN", "jsonlite")
missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("缺少包: ", paste(missing_pkgs, collapse = ", "),
       "\n请先在干净 R session 里 source('D:/RStudio/work/module1_install_deps.R')。")
}
if (!requireNamespace("RcppPlanc", quietly = TRUE)) {
  message("[警告] RcppPlanc 未装，iNMF 会用慢速回退实现。")
}

suppressPackageStartupMessages({
  library(rliger)
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(RANN)
  library(jsonlite)
})

# 打印 rliger 版本，方便排查
cat("rliger 版本:", as.character(packageVersion("rliger")), "\n")

set.seed(9527)
script_start <- Sys.time()
timings <- list()

# =========================
# 1) 路径
# =========================
data_dir <- "D:/SCAVENGE-reproducibility-main/blood_data"
in_dir   <- file.path(data_dir, "multiomics_output_v2")
out_dir  <- file.path(data_dir, "multiomics_module1_pairing", "method_C_liger")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("========== 路径配置 ==========\n")
cat("  输入:", in_dir, "\n")
cat("  输出:", out_dir, "\n")

# =========================
# 2) 读入已有产物（不重算 Gene Activity）
# =========================
cat("\n========== 加载数据 ==========\n")
t0 <- Sys.time()

ga_mat     <- readRDS(file.path(in_dir, "gene_activity_matrix.rds"))
seurat_rna <- readRDS(file.path(in_dir, "seurat_rna_normalized.rds"))
rna_counts <- GetAssayData(seurat_rna, assay = "RNA", layer = "counts")
atac_meta  <- readRDS(file.path(in_dir, "atac_metadata.rds"))
rna_meta   <- as.data.frame(seurat_rna@meta.data)

cat("  Gene Activity:", nrow(ga_mat),     "genes ×", ncol(ga_mat),     "cells\n")
cat("  scRNA counts :", nrow(rna_counts), "genes ×", ncol(rna_counts), "cells\n")

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
# 4) 构建 LIGER 对象
# =========================
cat("\n========== 构建 LIGER 对象 ==========\n")
t0 <- Sys.time()

lig <- createLiger(
  rawData = list(rna = rna_counts, atac = ga_mat),
  modal   = c("rna", "atac"),
  verbose = TRUE
)
cat("  LIGER 对象: 2 datasets,", ncol(rna_counts) + ncol(ga_mat), "cells\n")

timings$create_liger <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

# =========================
# 5) Normalize + selectGenes + scaleNotCenter（2.x 新参数名）
# =========================
cat("\n========== 归一化 + 选变基因 + scale ==========\n")
t0 <- Sys.time()

lig <- normalize(lig)
# 2.x: useDatasets 而不是 datasets.use
lig <- selectGenes(lig, useDatasets = "rna", thresh = 0.1)
cat("  HVG 数 (来自 RNA):", length(varFeatures(lig)), "\n")
lig <- scaleNotCenter(lig)

timings$preprocess <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("  [计时] 预处理: %.1f 秒\n", timings$preprocess))

# =========================
# 6) iNMF 整合
# =========================
cat("\n========== iNMF 整合 (k=20, lambda=5) ==========\n")
cat("  预计耗时 10–30 分钟...\n")
t0 <- Sys.time()

lig <- runIntegration(lig, k = 20, lambda = 5, method = "iNMF")

timings$inmf <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("  [计时] iNMF: %.1f 秒\n", timings$inmf))

# =========================
# 7) Quantile Normalization
# =========================
cat("\n========== Quantile Normalization ==========\n")
t0 <- Sys.time()

lig <- quantileNorm(lig)

H_norm <- getMatrix(lig, slot = "H.norm")
cat("  H.norm 维度:", dim(H_norm)[1], "×", dim(H_norm)[2], "\n")
stopifnot(nrow(H_norm) == (ncol(rna_counts) + ncol(ga_mat)))

timings$quantile_norm <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

# =========================
# 8) UMAP（rliger 2.x 新参数名）
# =========================
cat("\n========== UMAP ==========\n")
t0 <- Sys.time()

# 2.x: useRaw / useDims / nNeighbors / minDist；dimredName 控制存放名字
lig <- runUMAP(lig,
               useRaw      = FALSE,
               useDims     = 1:20,
               nNeighbors  = 30,
               minDist     = 0.3,
               dimredName  = "UMAP")

timings$umap <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("  [计时] UMAP: %.1f 秒\n", timings$umap))

# 断言 UMAP 坐标确实写入
get_umap_coords <- function(obj) {
  if (.hasSlot(obj, "dimReds") && length(obj@dimReds) > 0) {
    dr <- obj@dimReds
    nm <- names(dr)
    pick <- nm[grepl("UMAP", nm, ignore.case = TRUE)]
    if (length(pick) > 0) return(dr[[pick[1]]])
    return(dr[[1]])
  }
  cm <- cellMeta(obj)
  umap_cols <- grep("^UMAP", colnames(cm), value = TRUE)
  if (length(umap_cols) >= 2) {
    mat <- as.matrix(cm[, umap_cols[1:2], drop = FALSE])
    rownames(mat) <- rownames(cm)
    return(mat)
  }
  return(NULL)
}

umap_coords <- get_umap_coords(lig)
if (is.null(umap_coords)) {
  cat("[诊断] slotNames:\n"); print(slotNames(lig))
  cat("[诊断] cellMeta 列:\n"); print(names(cellMeta(lig)))
  stop("UMAP 坐标未找到。请检查上面诊断信息，把输出贴给我。")
}
cat("  UMAP 坐标:", nrow(umap_coords), "×", ncol(umap_coords), "\n")

# =========================
# 9) ATAC → RNA top-k 邻居查询
# =========================
cat("\n========== top-k RNA 邻居查询 ==========\n")
t0 <- Sys.time()

cm <- cellMeta(lig)
is_rna  <- cm$dataset == "rna"
is_atac <- cm$dataset == "atac"

H_rna  <- H_norm[is_rna,  , drop = FALSE]
H_atac <- H_norm[is_atac, , drop = FALSE]
cat("  H_rna :", nrow(H_rna),  "×", ncol(H_rna),  "\n")
cat("  H_atac:", nrow(H_atac), "×", ncol(H_atac), "\n")

K_TOPK <- 30
nn_res <- RANN::nn2(data = H_rna, query = H_atac, k = K_TOPK)
timings$knn <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("  [计时] kNN: %.1f 秒\n", timings$knn))

rna_cell_ids  <- rownames(H_rna)
atac_cell_ids <- rownames(H_atac)

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

# =========================
# 10) Label Transfer：top-k 邻居多数投票
# =========================
cat("\n========== Label Transfer (投票) ==========\n")
t0 <- Sys.time()

# rliger 在 cell 名前加了 "rna_" / "atac_" 前缀
rna_labels_by_id  <- setNames(as.character(rna_meta$BioClassification),
                              paste0("rna_",  rownames(rna_meta)))
atac_labels_by_id <- setNames(as.character(atac_meta$BioClassification),
                              paste0("atac_", rownames(atac_meta)))

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
# 11) 一致率评估（沿用 A 方法的非方阵处理）
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

diag_sum          <- sum(sapply(shared_labels, function(ct) conf_mat[ct, ct]))
overall_agreement <- diag_sum / sum(conf_mat)

cat("  总体一致率:", round(overall_agreement * 100, 1), "%\n")

agreement_df <- data.frame(
  celltype   = orig_types,
  n_cells    = as.integer(rowSums(conf_mat)),
  agreement  = round(agreement_per_type * 100, 1)
) %>% arrange(desc(agreement))

fwrite(agreement_df, file.path(out_dir, "label_transfer_agreement.csv"))
cat("  已保存: label_transfer_agreement.csv\n")

# 高置信度子集
conf_thr <- 0.5
high_conf_mask <- !is.na(label_df_ok$score) & label_df_ok$score >= conf_thr
cat("  高置信度 (>=", conf_thr, "):", sum(high_conf_mask), "/",
    nrow(label_df_ok), "\n")

if (sum(high_conf_mask) > 0) {
  conf_mat_hc <- table(
    Original    = label_df_ok$original.label[high_conf_mask],
    Transferred = label_df_ok$predicted.id[high_conf_mask]
  )
  orig_hc <- rownames(conf_mat_hc); trans_hc <- colnames(conf_mat_hc)
  agreement_hc <- sapply(orig_hc, function(ct) {
    if (ct %in% trans_hc) conf_mat_hc[ct, ct] / sum(conf_mat_hc[ct, ]) else 0
  })
  agreement_hc_df <- data.frame(
    celltype     = orig_hc,
    n_cells_hc   = as.integer(rowSums(conf_mat_hc)),
    agreement_hc = round(agreement_hc * 100, 1)
  ) %>% arrange(desc(agreement_hc))
  fwrite(agreement_hc_df, file.path(out_dir, "label_transfer_agreement_highconf.csv"))
  cat("  已保存: label_transfer_agreement_highconf.csv\n")

  hc_overall <- mean(label_df_ok$original.label[high_conf_mask] ==
                     label_df_ok$predicted.id[high_conf_mask])
  cat("  高置信度总体一致率:", round(hc_overall * 100, 1), "%\n")
}

# =========================
# 12) 保存共享嵌入（接口契约文件）
# =========================
cat("\n========== 保存共享嵌入 ==========\n")

cells_all     <- rownames(H_norm)
modality_all  <- ifelse(is_rna, "rna", "atac")
shared_embed_df <- data.frame(
  cell     = cells_all,
  modality = modality_all,
  H_norm,
  check.names = FALSE
)
saveRDS(shared_embed_df, file.path(out_dir, "shared_embedding.rds"),
        compress = FALSE)
cat("  已保存: shared_embedding.rds (", nrow(shared_embed_df), "×",
    ncol(shared_embed_df), ")\n")

saveRDS(H_norm, file.path(out_dir, "shared_embedding_Hnorm.rds"),
        compress = FALSE)
cat("  已保存: shared_embedding_Hnorm.rds\n")

# =========================
# 13) UMAP 可视化
# =========================
cat("\n========== UMAP 可视化 ==========\n")

umap_df <- data.frame(
  UMAP1   = umap_coords[, 1],
  UMAP2   = umap_coords[, 2],
  dataset = cm$dataset,
  cell    = rownames(cm),
  stringsAsFactors = FALSE
)

# 加细胞类型
umap_df$celltype <- NA_character_
rna_idx  <- umap_df$dataset == "rna"
atac_idx <- umap_df$dataset == "atac"
umap_df$celltype[rna_idx]  <- rna_labels_by_id[umap_df$cell[rna_idx]]
umap_df$celltype[atac_idx] <- label_df$predicted.id[
  match(umap_df$cell[atac_idx], label_df$atac_cell)
]

p_mod <- ggplot(umap_df, aes(UMAP1, UMAP2, color = dataset)) +
  geom_point(size = 0.3, alpha = 0.5) +
  scale_color_manual(values = c("rna" = "#4DBBD5", "atac" = "#E64B35")) +
  labs(title = "LIGER Co-embedding (H.norm + UMAP)",
       subtitle = paste0("n_RNA=", sum(rna_idx), ", n_ATAC=", sum(atac_idx))) +
  theme_classic()
ggsave(file.path(out_dir, "Fig_liger_umap_modality.png"),
       plot = p_mod, width = 8, height = 6, dpi = 300)

p_ct <- ggplot(umap_df, aes(UMAP1, UMAP2, color = celltype)) +
  geom_point(size = 0.3, alpha = 0.5) +
  labs(title = "LIGER Co-embedding: Cell Type") +
  theme_classic() + theme(legend.position = "none")
ggsave(file.path(out_dir, "Fig_liger_umap_celltype.png"),
       plot = p_ct, width = 8, height = 6, dpi = 300)

p_score <- ggplot(label_df, aes(score)) +
  geom_histogram(bins = 50, fill = "#3C5488", color = "white") +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "red") +
  labs(title = "LIGER: top-k Vote Score Distribution",
       x = "Max-vote proportion", y = "Cells") +
  theme_classic()
ggsave(file.path(out_dir, "Fig_liger_prediction_score.png"),
       plot = p_score, width = 6, height = 4, dpi = 300)

cat("  已保存 3 张图\n")

# =========================
# 14) 保存 LIGER 对象
# =========================
cat("\n========== 保存 LIGER 对象 ==========\n")
saveRDS(lig, file.path(out_dir, "liger_object.rds"), compress = FALSE)
cat("  已保存: liger_object.rds\n")

# =========================
# 15) timing.json
# =========================
timings$total    <- as.numeric(difftime(Sys.time(), script_start, units = "secs"))
timings$method   <- "LIGER (rliger::runIntegration iNMF k=20 lambda=5)"
timings$citation <- "Welch et al. 2019 Cell; Gao et al. 2021 Nat Biotechnol"
timings$rliger_version <- as.character(packageVersion("rliger"))
jsonlite::write_json(timings, file.path(out_dir, "timing.json"),
                     pretty = TRUE, auto_unbox = TRUE)

# =========================
# 16) 总结
# =========================
cat("\n")
cat("======================================================\n")
cat("  模块 1 · 方法 C (LIGER) 完成\n")
cat("======================================================\n")
cat("  rliger 版本    :", as.character(packageVersion("rliger")), "\n")
cat("  H.norm 维度    :", dim(H_norm)[1], "×", dim(H_norm)[2], "\n")
cat("  整体一致率     :", round(overall_agreement * 100, 1), "%\n")
if (sum(high_conf_mask) > 0) {
  cat("  高置信度一致率 :", round(hc_overall * 100, 1), "%\n")
}
cat("  总耗时         :", round(timings$total / 60, 1), "分钟\n")
cat("  输出目录       :", out_dir, "\n")
cat("======================================================\n")
