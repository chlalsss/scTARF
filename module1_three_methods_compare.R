# module1_three_methods_compare.R
# 模块 1 · 三方法对比骨架（A: Seurat anchor / B: bindSC / C: LIGER）
#
# 本脚本是**占位 + 接口**，等三种方法产物齐全后再真正运行。
# 对比逻辑不自写机器学习算法，只用标准指标：
#   1. ATAC cell 的 predicted.id 在三方法之间的 Cohen's κ
#   2. 高置信度子集上逐细胞类型一致率差异
#   3. top-k 邻居重合率（Jaccard on top-k RNA cell sets）
#   4. 耗时对比（读 timing.json）
#
# 运行前提（每个目录下都已生成）：
#   multiomics_module1_pairing/method_A_seurat/  - shared_embedding.rds,
#                                                  atac_topk_rna_neighbors.rds,
#                                                  atac_label_transfer.csv, timing.json
#   method_B_bindSC/ - 同上（2026-04-23 从 GLUE 切换为 bindSC）
#   method_C_liger/  - 同上（本次刚产出）

library(data.table)
library(dplyr)
library(ggplot2)
library(jsonlite)

data_dir <- "D:/SCAVENGE-reproducibility-main/blood_data"
base_dir <- file.path(data_dir, "multiomics_module1_pairing")
out_dir  <- file.path(base_dir, "comparison")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

method_dirs <- list(
  A_Seurat = file.path(base_dir, "method_A_seurat"),
  B_bindSC = file.path(base_dir, "method_B_bindSC"),
  C_LIGER  = file.path(base_dir, "method_C_liger")
)

# =========================
# 1) 读入每种方法的 label transfer 结果
# =========================
load_label <- function(d) {
  f <- file.path(d, "atac_label_transfer.csv")
  if (!file.exists(f)) return(NULL)
  fread(f)
}

labels_by_method <- lapply(method_dirs, load_label)
ok_methods <- names(labels_by_method)[!sapply(labels_by_method, is.null)]
cat("可用方法:", paste(ok_methods, collapse = ", "), "\n")

if (length(ok_methods) < 2) {
  stop("至少需要两种方法的产物才能对比。当前仅 ",
       length(ok_methods), " 个。")
}

# =========================
# 2) 两两 Cohen's κ（预测标签一致性）
# =========================
# 占位：kappa 计算 - 等三份数据齐全后实现
# 关键步骤：
#   1) 对齐所有方法共有的 atac_cell
#   2) 对每对方法 (i, j) 计算 Cohen's κ 与 agreement%
#   3) 输出 pairwise 表 + heatmap

# =========================
# 3) top-k 邻居 Jaccard
# =========================
# 占位：读 atac_topk_rna_neighbors.rds，对每个 ATAC 细胞在 A/B/C 下的 top-k RNA 集合取 Jaccard
# 输出每方法对的 Jaccard 分布直方图

# =========================
# 4) 耗时对比
# =========================
timings <- lapply(ok_methods, function(m) {
  f <- file.path(method_dirs[[m]], "timing.json")
  if (file.exists(f)) fromJSON(f) else NULL
})
names(timings) <- ok_methods
saveRDS(timings, file.path(out_dir, "timings_raw.rds"), compress = FALSE)

# =========================
# 5) 共享嵌入在同一 UMAP 上对比（可选）
# =========================
# 占位：对每个方法的 shared_embedding.rds 在 ATAC 子集上做 UMAP，
# 同色点画三联图

cat("\n[占位] 等三方法产物齐全后打开实际计算逻辑。\n")
cat("当前仅建立骨架与接口，未做真实对比。\n")
