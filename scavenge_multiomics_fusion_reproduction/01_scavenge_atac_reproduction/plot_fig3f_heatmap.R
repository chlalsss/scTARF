# plot_fig3f_heatmap.R
# 目的：生成论文 Fig.3f 风格的 heatmap（median TRS，traits 为行，cell types 为列）
# 说明：尽量沿用你现有脚本的结构与数据流，注释为中文，便于直接运行与复查

library(data.table)
library(dplyr)
library(pheatmap)

# =========================
# 0) 运行前检查：你需要先在当前 R session 中准备好 SE_heme_small
#    - SE_heme_small: SummarizedExperiment
#    - colData(SE_heme_small) 必须含 BioClassification
# =========================
if (!exists("SE_heme_small")) {
  stop("当前环境中找不到对象 SE_heme_small。请先 readRDS() 加载小数据集的 SummarizedExperiment。")
}
if (!("BioClassification" %in% colnames(colData(SE_heme_small)))) {
  stop("colData(SE_heme_small) 中不包含 BioClassification。请确认小数据集对象是否正确。")
}

# =========================
# 1) 输出目录（改成你自己的路径）
# =========================
out_dir <- "D:/SCAVENGE-reproducibility-main/blood_data/scavenge_out_small_all_traits"

# =========================
# 2) 读 TRS 全细胞矩阵（fullcells 版本）
#    该矩阵是 Fig.3f 的“每细胞 TRS”来源
# =========================
trs_path <- file.path(out_dir, "TRS_matrix_small_all_traits_fullcells.rds")
if (!file.exists(trs_path)) {
  stop("找不到 TRS 矩阵文件：", trs_path)
}
trs_full <- readRDS(trs_path)

# 兼容性检查：必须有 cell 列
if (!("cell" %in% colnames(trs_full))) {
  stop("TRS 矩阵 trs_full 中缺少 'cell' 列。请确认 fullcells 版本的结构：应包含 cell + traits。")
}

# trait 列名（去掉 cell）
trait_cols <- setdiff(colnames(trs_full), "cell")
if (length(trait_cols) < 1) {
  stop("未在 trs_full 中找到 trait 列（除 cell 外没有其他列）。")
}

# =========================
# 3) 从 SE 里取 cell + BioClassification（这里只需要细胞类型标签）
# =========================
cd <- as.data.frame(colData(SE_heme_small))
# 注意：你的旧代码用 rownames(cd) 当 cell
# 但 SummarizedExperiment 通常 colnames(se) 才是细胞 ID
# 这里做一个稳妥处理：优先用 colnames(SE_heme_small)
cd$cell <- colnames(SE_heme_small)

# 保留 cell 和细胞类型标签
cd <- cd[, c("cell", "BioClassification")]

# =========================
# 4) 合并 cell type + TRS（按 cell 对齐）
#    left_join：以 cd 为主，保证每个细胞都有 BioClassification（TRS 缺失则为 NA）
# =========================
df <- dplyr::left_join(cd, trs_full, by = "cell")

# （可选）如果你想强制把 NA 当作 0（你们 fullcells 理论上已补齐/置0）
# df[, trait_cols] <- lapply(df[, trait_cols, drop=FALSE], function(x) { x[is.na(x)] <- 0; x })

# =========================
# 5) 每个 cell type 取 median TRS（对应论文 Fig.3f 的“median TRS”）
#    输出：每行一个 BioClassification，每列一个 trait 的 median(TRS)
# =========================
med_df <- df %>%
  group_by(BioClassification) %>%
  summarise(across(all_of(trait_cols), ~median(.x, na.rm = TRUE)), .groups = "drop")

# 这里的 mat：行=cell type，列=trait
mat <- as.matrix(med_df[, trait_cols, drop = FALSE])
rownames(mat) <- med_df$BioClassification

# =========================
# 6) 按论文风格调整：
#    - 论文是“traits 为行，cell types 为列”
#    - 论文显示的是 median TRS 的强弱（渐变色），不是 z-score 红蓝图
# =========================
mat_tc <- t(mat)  # 转置后：行=trait，列=cell type

# =========================
# 7) 输出 heatmap（自动创建目录 + 时间戳，避免写入失败/占用）
# =========================
fig_dir <- file.path(out_dir, "fig3_like")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

out_png <- file.path(
  fig_dir,
  paste0("Fig3f_like_medianTRS_heatmap_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
)

# 关键参数：
# - scale="none"：不标准化，更像论文
# - clustering_method="complete"：层次聚类（无监督）
pheatmap::pheatmap(
  mat_tc,
  scale = "none",
  clustering_method = "complete",
  filename = out_png,
  width = 10, height = 6
)

cat("saved to:", out_png, "\n")

# =========================
# 8)（可选）对 traits 做 k-means=4 分组，并按组排序再画一张
#    论文提到 traits 被分成 4 个 cluster（k-means）
#    注意：你当前是 15 traits 版本，与论文 22 traits 不一致属正常
# =========================
set.seed(1)
km <- kmeans(mat_tc, centers = 4)     # 在 trait（行）层面做 kmeans
trait_order <- order(km$cluster)      # 按 cluster 编号排序
mat_tc_km <- mat_tc[trait_order, , drop = FALSE]

out_png_km <- file.path(
  fig_dir,
  paste0("Fig3f_like_medianTRS_heatmap_kmeans4_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
)

pheatmap::pheatmap(
  mat_tc_km,
  scale = "none",
  clustering_method = "complete",
  filename = out_png_km,
  width = 10, height = 6
)

cat("saved to:", out_png_km, "\n")