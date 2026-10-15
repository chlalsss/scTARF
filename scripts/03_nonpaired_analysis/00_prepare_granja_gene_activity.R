# multiomics_phase2_shared_space.R
# 阶段 2：构建 scATAC 与 scRNA 的共享特征空间（Gene Activity）
#
# 核心思路：
#   scATAC 的原始特征是 peak（染色质开放区），scRNA 的特征是基因表达
#   要跨模态整合，需要将 ATAC peak 转换为"基因活性"（Gene Activity Score）
#   转换逻辑：统计每个基因 body + 上游 promoter 区域内所有 peak 的 accessibility
#   转换后，scATAC 侧也有一个 cell x gene 矩阵，和 scRNA 共享特征空间
#
# 依赖包：Seurat, Signac, GenomeInfoDb, EnsDb.Hsapiens.v75 (hg19 对应)
# 注意：首次运行前需确保已安装上述包

# =========================
# 0) 加载依赖
# =========================
library(Seurat)
library(Signac)
library(SummarizedExperiment)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)  # hg19 对应的基因注释
library(Matrix)
library(data.table)
library(dplyr)
library(ggplot2)

set.seed(9527)

# =========================
# 1) 路径设置
# =========================
source("scripts/_paths.R")
data_dir <- sctarf_blood_root()
out_dir  <- file.path(data_dir, "multiomics_output_v2")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

atac_path <- file.path(data_dir, "scATAC-Healthy-Hematopoiesis-191120.rds")
rna_path  <- file.path(data_dir, "scRNA-Healthy-Hematopoiesis-191120.rds")

# =========================
# 2) 加载 SummarizedExperiment 对象
# =========================
cat("========== 加载数据 ==========\n")
SE_atac <- readRDS(atac_path)
SE_rna  <- readRDS(rna_path)
cat("scATAC:", ncol(SE_atac), "cells x", nrow(SE_atac), "peaks\n")
cat("scRNA: ", ncol(SE_rna),  "cells x", nrow(SE_rna),  "genes\n")

# =========================
# 3) 将 scATAC SE 转换为 Seurat + ChromatinAssay
# =========================
cat("\n========== 构建 ATAC Seurat 对象 ==========\n")

# 3a) 提取 peak x cell 计数矩阵
atac_counts <- assay(SE_atac)

# 3b) 从 rowRanges 获取 peak 坐标，构建 Signac 需要的 GRanges
peak_gr <- rowRanges(SE_atac)

# 3c) 检查 peak 坐标的 seqlevels 风格（需要是 "chr1" 格式）
cat("Peak seqlevels 样例:", head(seqlevels(peak_gr), 5), "\n")

# 3d) 创建 ChromatinAssay
#     ChromatinAssay 是 Signac 专门为 ATAC 数据设计的 assay 类型
#     它知道每行是一个基因组区间（peak），而不仅仅是一个 feature name
chrom_assay <- CreateChromatinAssay(
  counts   = atac_counts,
  ranges   = peak_gr,
  sep      = c(":", "-"),  # peak name 格式分隔符
  min.cells = 0,
  min.features = 0
)

# 3e) 创建 Seurat 对象
seurat_atac <- CreateSeuratObject(
  counts  = chrom_assay,
  assay   = "ATAC",
  meta.data = as.data.frame(colData(SE_atac))
)
cat("ATAC Seurat 对象:", ncol(seurat_atac), "cells x", nrow(seurat_atac), "peaks\n")

# =========================
# 4) 获取 hg19 基因注释（用于 GeneActivity 计算）
# =========================
cat("\n========== 提取 hg19 基因注释 ==========\n")

# EnsDb.Hsapiens.v75 对应 GRCh37 / hg19
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)

# 确保 seqlevel 风格和 ATAC 一致（UCSC 风格: chr1, chr2...）
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "hg19"

# 将注释附加到 ChromatinAssay
Annotation(seurat_atac) <- annotations
cat("基因注释加载完成:", length(annotations), "条记录\n")

# =========================
# 5) 计算 Gene Activity Score（基于 peak overlap，不需要 fragments 文件）
# =========================
cat("\n========== 计算 Gene Activity Score ==========\n")
cat("方法：从 EnsDb 提取基因坐标 → 延伸 promoter 2kb → 与 peak 做 overlap → 求和\n")
cat("预计耗时 1-5 分钟...\n")

step5_start <- Sys.time()

# 5a) 从 EnsDb 提取基因坐标（只取 protein_coding 基因，去重）
gene_coords <- ensembldb::genes(EnsDb.Hsapiens.v75,
                                filter = ~ gene_biotype == "protein_coding",
                                columns = c("gene_name", "gene_biotype"))
seqlevelsStyle(gene_coords) <- "UCSC"
gene_coords <- keepStandardChromosomes(gene_coords, pruning.mode = "coarse")

# 去除重复基因名（保留最长的转录本区间）
gene_coords <- gene_coords[!duplicated(gene_coords$gene_name)]
names(gene_coords) <- gene_coords$gene_name

# 5b) 按转录方向延伸 promoter 2kb
#     正链：start 上游延伸；负链：end 下游延伸
gene_regions <- gene_coords
strand_pos <- as.character(strand(gene_regions)) == "+"
strand_neg <- as.character(strand(gene_regions)) == "-"
strand_star <- as.character(strand(gene_regions)) == "*"

start(gene_regions)[strand_pos | strand_star] <-
  pmax(start(gene_regions)[strand_pos | strand_star] - 2000L, 1L)
end(gene_regions)[strand_neg] <-
  end(gene_regions)[strand_neg] + 2000L

cat("基因区域数:", length(gene_regions), "\n")

# 5c) 找 peak 与基因区域的 overlap
peak_gr <- rowRanges(SE_atac)
overlaps <- GenomicRanges::findOverlaps(
  query   = gene_regions,
  subject = peak_gr
)
cat("Peak-Gene overlaps:", length(overlaps), "\n")

# 5d) 构建 gene x cell 矩阵（对每个基因，sum 所有 overlap 的 peak counts）
atac_counts <- assay(SE_atac)

# 构建 gene x peak 的稀疏映射矩阵
n_genes <- length(gene_regions)
n_peaks <- nrow(atac_counts)

gene_peak_map <- Matrix::sparseMatrix(
  i = queryHits(overlaps),
  j = subjectHits(overlaps),
  x = 1,
  dims = c(n_genes, n_peaks)
)

# 矩阵乘法一步完成：gene_activity = gene_peak_map %*% peak_counts
gene_activities <- gene_peak_map %*% atac_counts
rownames(gene_activities) <- names(gene_regions)

# 移除全零基因
gene_keep <- Matrix::rowSums(gene_activities) > 0
gene_activities <- gene_activities[gene_keep, ]

step5_end <- Sys.time()
cat(sprintf("Gene Activity 计算完成，耗时: %.1f 秒\n",
            as.numeric(difftime(step5_end, step5_start, units = "secs"))))
cat("Gene Activity 矩阵:", nrow(gene_activities), "genes x", ncol(gene_activities), "cells\n")

# 将 Gene Activity 作为新 assay 存入 Seurat 对象
seurat_atac[["RNA_from_ATAC"]] <- CreateAssayObject(counts = gene_activities)

# =========================
# 6) 归一化 Gene Activity（log-normalize，与 scRNA 处理一致）
# =========================
cat("\n========== 归一化 ==========\n")

DefaultAssay(seurat_atac) <- "RNA_from_ATAC"
seurat_atac <- NormalizeData(seurat_atac, normalization.method = "LogNormalize",
                             scale.factor = 10000)
cat("ATAC Gene Activity 归一化完成 (LogNormalize)\n")

# =========================
# 7) 构建 scRNA Seurat 对象并归一化
# =========================
cat("\n========== 构建 RNA Seurat 对象 ==========\n")

rna_counts <- assay(SE_rna)
seurat_rna <- CreateSeuratObject(
  counts    = rna_counts,
  assay     = "RNA",
  meta.data = as.data.frame(colData(SE_rna))
)

seurat_rna <- NormalizeData(seurat_rna, normalization.method = "LogNormalize",
                            scale.factor = 10000)
seurat_rna <- FindVariableFeatures(seurat_rna, selection.method = "vst",
                                   nfeatures = 3000)

cat("RNA Seurat 对象:", ncol(seurat_rna), "cells x", nrow(seurat_rna), "genes\n")
cat("高变基因数:", length(VariableFeatures(seurat_rna)), "\n")

# =========================
# 8) 检查共享基因空间
# =========================
cat("\n========== 共享特征空间检查 ==========\n")

genes_atac_ga <- rownames(gene_activities)
genes_rna     <- rownames(seurat_rna)
shared_genes  <- intersect(genes_atac_ga, genes_rna)

cat("ATAC Gene Activity 基因数:", length(genes_atac_ga), "\n")
cat("scRNA 基因数:             ", length(genes_rna), "\n")
cat("共享基因数:               ", length(shared_genes), "\n")
cat("共享率 (占RNA):           ", round(length(shared_genes)/length(genes_rna)*100, 1), "%\n")

# 高变基因中有多少在 ATAC Gene Activity 中
hvg <- VariableFeatures(seurat_rna)
hvg_in_atac <- sum(hvg %in% genes_atac_ga)
cat("高变基因在 ATAC GA 中:   ", hvg_in_atac, "/", length(hvg),
    "(", round(hvg_in_atac/length(hvg)*100, 1), "%)\n")

# =========================
# 9) Gene Activity 的基本分布检查
# =========================
cat("\n========== Gene Activity 分布检查 ==========\n")

ga_per_cell <- Matrix::colSums(gene_activities)
ga_genes_per_cell <- Matrix::colSums(gene_activities > 0)

cat("每细胞 GA 总和:    median =", round(median(ga_per_cell)),
    " mean =", round(mean(ga_per_cell)), "\n")
cat("每细胞非零基因数:  median =", median(ga_genes_per_cell),
    " mean =", round(mean(ga_genes_per_cell)), "\n")

# GA 分布直方图
p_ga <- ggplot(data.frame(x = ga_genes_per_cell), aes(x)) +
  geom_histogram(bins = 80, fill = "#E64B35", color = "white") +
  labs(title = "scATAC: Genes with Activity per Cell",
       x = "Number of Genes (GA > 0)", y = "Cells") +
  theme_classic()

ggsave(file.path(out_dir, "Fig_phase2_ga_genes_per_cell.png"),
       plot = p_ga, width = 6, height = 4, dpi = 300)
cat("已保存: Fig_phase2_ga_genes_per_cell.png\n")

# =========================
# 10) 保存中间结果
# =========================
cat("\n========== 保存中间结果 ==========\n")

saveRDS(seurat_atac, file.path(out_dir, "seurat_atac_with_GA.rds"))
saveRDS(seurat_rna,  file.path(out_dir, "seurat_rna_normalized.rds"))
saveRDS(gene_activities, file.path(out_dir, "gene_activity_matrix.rds"))
# Metadata extraction retained from multiomics_phase3to6.R.
saveRDS(as.data.frame(colData(SE_atac)), file.path(out_dir, "atac_metadata.rds"), compress = FALSE)

cat("已保存: seurat_atac_with_GA.rds\n")
cat("已保存: seurat_rna_normalized.rds\n")
cat("已保存: gene_activity_matrix.rds\n")

# 共享基因列表
writeLines(shared_genes, file.path(out_dir, "shared_genes_atac_rna.txt"))
cat("已保存: shared_genes_atac_rna.txt\n")

# =========================
# 11) 阶段 2 总结
# =========================
cat("\n")
cat("══════════════════════════════════════════\n")
cat("  阶段 2 总结：共享特征空间已构建\n")
cat("══════════════════════════════════════════\n")
cat("  ATAC Gene Activity: ", nrow(gene_activities), " genes x ",
    ncol(gene_activities), " cells\n")
cat("  scRNA:              ", nrow(seurat_rna), " genes x ",
    ncol(seurat_rna), " cells\n")
cat("  共享基因:           ", length(shared_genes), "\n")
cat("  高变基因覆盖:       ", hvg_in_atac, "/", length(hvg), "\n")
cat("  输出目录:           ", out_dir, "\n")
cat("══════════════════════════════════════════\n")
cat("\n阶段 2 完成。下一步：阶段 3（跨模态整合 + 联合嵌入）\n")
