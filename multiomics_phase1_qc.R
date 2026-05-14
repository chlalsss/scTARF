# multiomics_phase1_qc.R
# 阶段 1：加载 scRNA-seq 数据，与现有 scATAC-seq 做细胞类型体系对比与 QC
# 依赖：SummarizedExperiment, data.table, dplyr, ggplot2

library(SummarizedExperiment)
library(data.table)
library(dplyr)
library(ggplot2)

set.seed(9527)

# =========================
# 0) 路径设置
# =========================
data_dir   <- "D:/SCAVENGE-reproducibility-main/blood_data"
atac_path  <- file.path(data_dir, "scATAC-Healthy-Hematopoiesis-191120.rds")
rna_path   <- file.path(data_dir, "scRNA-Healthy-Hematopoiesis-191120.rds")

# 输出目录
out_dir <- file.path(data_dir, "multiomics_output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# =========================
# 1) 加载两个 SummarizedExperiment
# =========================
cat("========== 加载数据 ==========\n")

cat("[1/2] 加载 scATAC...\n")
SE_atac <- readRDS(atac_path)
cat("  scATAC: ", ncol(SE_atac), " cells x ", nrow(SE_atac), " peaks\n")

cat("[2/2] 加载 scRNA...\n")
SE_rna <- readRDS(rna_path)
cat("  scRNA:  ", ncol(SE_rna), " cells x ", nrow(SE_rna), " genes\n")

# =========================
# 2) 检查 colData 元数据结构
# =========================
cat("\n========== 元数据检查 ==========\n")

cd_atac <- as.data.frame(colData(SE_atac))
cd_rna  <- as.data.frame(colData(SE_rna))

cat("[scATAC colData] 列名:\n")
print(names(cd_atac))
cat("[scRNA  colData] 列名:\n")
print(names(cd_rna))

# 检查两边是否都有 BioClassification
has_bio_atac <- "BioClassification" %in% names(cd_atac)
has_bio_rna  <- "BioClassification" %in% names(cd_rna)
cat("\nBioClassification 存在？  scATAC:", has_bio_atac, " scRNA:", has_bio_rna, "\n")

# =========================
# 3) 细胞类型体系对比
# =========================
cat("\n========== 细胞类型体系对比 ==========\n")

types_atac <- sort(unique(cd_atac$BioClassification))
types_rna  <- sort(unique(cd_rna$BioClassification))

cat("[scATAC] 细胞类型 (n=", length(types_atac), "):\n")
print(types_atac)

cat("\n[scRNA]  细胞类型 (n=", length(types_rna), "):\n")
print(types_rna)

# 交集、ATAC独有、RNA独有
shared_types   <- intersect(types_atac, types_rna)
atac_only      <- setdiff(types_atac, types_rna)
rna_only       <- setdiff(types_rna, types_atac)

cat("\n--- 对比结果 ---\n")
cat("共有类型 (n=", length(shared_types), "):", paste(shared_types, collapse=", "), "\n")
cat("ATAC 独有 (n=", length(atac_only), "):",
    ifelse(length(atac_only) == 0, "无", paste(atac_only, collapse=", ")), "\n")
cat("RNA  独有 (n=", length(rna_only), "):",
    ifelse(length(rna_only) == 0, "无", paste(rna_only, collapse=", ")), "\n")

# =========================
# 4) 细胞类型分布统计
# =========================
cat("\n========== 细胞类型分布统计 ==========\n")

count_atac <- cd_atac %>%
  count(BioClassification) %>%
  arrange(desc(n)) %>%
  mutate(pct = round(n / sum(n) * 100, 1)) %>%
  rename(celltype = BioClassification, n_atac = n, pct_atac = pct)

count_rna <- cd_rna %>%
  count(BioClassification) %>%
  arrange(desc(n)) %>%
  mutate(pct = round(n / sum(n) * 100, 1)) %>%
  rename(celltype = BioClassification, n_rna = n, pct_rna = pct)

# 合并两边分布，方便对比
count_merged <- full_join(count_atac, count_rna, by = "celltype") %>%
  arrange(celltype)

cat("\n细胞类型分布对比表：\n")
print(as.data.frame(count_merged), row.names = FALSE)

# 保存对比表
fwrite(count_merged, file.path(out_dir, "celltype_distribution_comparison.csv"))
cat("\n已保存: celltype_distribution_comparison.csv\n")

# =========================
# 5) 可视化：并排柱状图对比细胞类型分布
# =========================
cat("\n========== 生成可视化 ==========\n")

# 准备长表
plot_df <- bind_rows(
  count_atac %>% select(celltype, n = n_atac) %>% mutate(modality = "scATAC"),
  count_rna  %>% select(celltype, n = n_rna)  %>% mutate(modality = "scRNA")
)

p_bar <- ggplot(plot_df, aes(x = reorder(celltype, -n), y = n, fill = modality)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_manual(values = c("scATAC" = "#E64B35", "scRNA" = "#4DBBD5")) +
  labs(
    title = "Cell Type Distribution: scATAC vs scRNA",
    x = NULL, y = "Number of Cells", fill = "Modality"
  ) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

ggsave(file.path(out_dir, "Fig_phase1_celltype_bar_comparison.png"),
       plot = p_bar, width = 12, height = 6, dpi = 300)
cat("已保存: Fig_phase1_celltype_bar_comparison.png\n")

# =========================
# 6) 可视化：细胞比例散点图（检验比例一致性）
# =========================
pct_merged <- count_merged %>% filter(!is.na(pct_atac) & !is.na(pct_rna))

p_scatter <- ggplot(pct_merged, aes(x = pct_atac, y = pct_rna, label = celltype)) +
  geom_point(size = 3, color = "#3C5488") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_text(vjust = -0.8, size = 3) +
  labs(
    title = "Cell Type Proportion Consistency (ATAC vs RNA)",
    x = "scATAC (%)", y = "scRNA (%)"
  ) +
  theme_classic(base_size = 12) +
  coord_equal()

ggsave(file.path(out_dir, "Fig_phase1_celltype_proportion_scatter.png"),
       plot = p_scatter, width = 7, height = 7, dpi = 300)
cat("已保存: Fig_phase1_celltype_proportion_scatter.png\n")

# =========================
# 7) scRNA 基础 QC 指标
# =========================
cat("\n========== scRNA 基础 QC ==========\n")

rna_mat <- assay(SE_rna)

genes_per_cell <- Matrix::colSums(rna_mat > 0)
umis_per_cell  <- Matrix::colSums(rna_mat)
cells_per_gene <- Matrix::rowSums(rna_mat > 0)

cat("每细胞检测基因数:  median =", median(genes_per_cell),
    " mean =", round(mean(genes_per_cell)),
    " range = [", min(genes_per_cell), ",", max(genes_per_cell), "]\n")
cat("每细胞 UMI 总数:   median =", median(umis_per_cell),
    " mean =", round(mean(umis_per_cell)),
    " range = [", min(umis_per_cell), ",", max(umis_per_cell), "]\n")
cat("每基因被检测细胞数: median =", median(cells_per_gene),
    " mean =", round(mean(cells_per_gene)), "\n")

# QC 直方图
p_genes <- ggplot(data.frame(x = genes_per_cell), aes(x)) +
  geom_histogram(bins = 80, fill = "#4DBBD5", color = "white") +
  labs(title = "scRNA: Genes per Cell", x = "Number of Genes", y = "Cells") +
  theme_classic()

p_umis <- ggplot(data.frame(x = umis_per_cell), aes(x)) +
  geom_histogram(bins = 80, fill = "#4DBBD5", color = "white") +
  labs(title = "scRNA: UMIs per Cell", x = "Total UMIs", y = "Cells") +
  theme_classic()

ggsave(file.path(out_dir, "Fig_phase1_rna_genes_per_cell.png"),
       plot = p_genes, width = 6, height = 4, dpi = 300)
ggsave(file.path(out_dir, "Fig_phase1_rna_umis_per_cell.png"),
       plot = p_umis, width = 6, height = 4, dpi = 300)
cat("已保存: Fig_phase1_rna_genes_per_cell.png\n")
cat("已保存: Fig_phase1_rna_umis_per_cell.png\n")

# =========================
# 8) scATAC 稀疏性回顾（与 scRNA 对比）
# =========================
cat("\n========== 跨模态稀疏性对比 ==========\n")

atac_mat <- assay(SE_atac)
atac_features_per_cell <- Matrix::colSums(atac_mat > 0)

sparsity_compare <- data.frame(
  modality = c("scATAC", "scRNA"),
  n_cells = c(ncol(SE_atac), ncol(SE_rna)),
  n_features = c(nrow(SE_atac), nrow(SE_rna)),
  median_features_per_cell = c(median(atac_features_per_cell), median(genes_per_cell)),
  overall_sparsity = c(
    1 - Matrix::nnzero(atac_mat) / (as.double(nrow(atac_mat)) * ncol(atac_mat)),
    1 - Matrix::nnzero(rna_mat)  / (as.double(nrow(rna_mat))  * ncol(rna_mat))
  )
)
sparsity_compare$overall_sparsity <- round(sparsity_compare$overall_sparsity, 4)

cat("\n跨模态数据规模与稀疏性对比：\n")
print(sparsity_compare, row.names = FALSE)

fwrite(sparsity_compare, file.path(out_dir, "sparsity_comparison_atac_vs_rna.csv"))
cat("\n已保存: sparsity_comparison_atac_vs_rna.csv\n")

# =========================
# 9) 已有 TRS 结果与 scRNA 的可衔接性检查
# =========================
cat("\n========== TRS 结果衔接性检查 ==========\n")

trs_path <- file.path(data_dir, "scavenge_out_small_all_traits",
                       "TRS_matrix_small_all_traits_fullcells.rds")
if (file.exists(trs_path)) {
  trs_full <- readRDS(trs_path)
  cat("TRS 矩阵: ", nrow(trs_full), " cells x ", ncol(trs_full) - 1, " traits\n")

  # TRS 的细胞 ID 和 scATAC 的细胞 ID 一致性
  trs_cells  <- trs_full$cell
  atac_cells <- colnames(SE_atac)
  cat("TRS cell ∩ ATAC cell: ", length(intersect(trs_cells, atac_cells)),
      " / ", length(trs_cells), "\n")

  # scRNA 不直接共享 cell ID（不同实验批次的单细胞），
  # 但通过整合后的 anchor 可以建立对应关系
  rna_cells <- colnames(SE_rna)
  direct_overlap <- length(intersect(trs_cells, rna_cells))
  cat("TRS cell ∩ RNA cell (直接ID): ", direct_overlap, "\n")
  if (direct_overlap == 0) {
    cat("  → 预期结果：scATAC 和 scRNA 来自不同细胞，无直接 barcode 重叠\n")
    cat("  → 后续通过跨模态整合（阶段3）建立 cell-cell 对应关系\n")
  }
} else {
  cat("未找到 TRS 矩阵文件，跳过衔接性检查\n")
}

# =========================
# 10) 阶段 1 总结
# =========================
cat("\n")
cat("══════════════════════════════════════════\n")
cat("  阶段 1 QC 总结\n")
cat("══════════════════════════════════════════\n")
cat("  scATAC:  ", ncol(SE_atac), " cells x ", nrow(SE_atac), " peaks\n")
cat("  scRNA:   ", ncol(SE_rna),  " cells x ", nrow(SE_rna),  " genes\n")
cat("  共有细胞类型: ", length(shared_types), "/",
    max(length(types_atac), length(types_rna)), "\n")
if (length(atac_only) > 0) cat("  ATAC独有: ", paste(atac_only, collapse=", "), "\n")
if (length(rna_only) > 0)  cat("  RNA独有:  ", paste(rna_only, collapse=", "), "\n")
cat("  输出目录: ", out_dir, "\n")
cat("══════════════════════════════════════════\n")
cat("\n阶段 1 完成。下一步：阶段 2（构建共享特征空间）\n")
