obj <- readRDS("D:/SCAVENGE-reproducibility-main/blood_data/scATAC-Healthy-Hematopoiesis-191120.rds")
class(obj)



library(SummarizedExperiment)
library(Matrix)

assayNames(obj)          # 看看有哪些矩阵，一般会有 "counts"
M <- assay(obj, "counts")  # 如果这里报错，把 "counts" 换成你看到的名字

dim(M)  # peaks x cells



P <- nrow(M)  # number of peaks
N <- ncol(M)  # number of cells

# 非零计数（>0）的个数
nnz_cell <- Matrix::colSums(M > 0)  # per cell: how many peaks have signal
nnz_peak <- Matrix::rowSums(M > 0)  # per peak: how many cells have signal

# 论文里的 sparsity = absence of signal 的比例 = 0 的比例
cell_sparsity <- 1 - (nnz_cell / P)  # per cell
peak_sparsity <- 1 - (nnz_peak / N)  # per peak

# 全矩阵总体稀疏性（0 的比例）
overall_sparsity <- 1 - (Matrix::nnzero(M) / (as.numeric(P) * as.numeric(N)))

# 输出几个关键数
overall_sparsity
summary(cell_sparsity)
summary(peak_sparsity)

# 常用汇总（更适合写进论文/补充材料）
sparsity_summary <- data.frame(
  overall_sparsity = overall_sparsity,
  cell_sparsity_median = median(cell_sparsity),
  cell_sparsity_IQR = IQR(cell_sparsity),
  cell_sparsity_p90 = as.numeric(quantile(cell_sparsity, 0.90)),
  cell_sparsity_p95 = as.numeric(quantile(cell_sparsity, 0.95)),
  peak_sparsity_median = median(peak_sparsity),
  peak_sparsity_IQR = IQR(peak_sparsity),
  peak_sparsity_p90 = as.numeric(quantile(peak_sparsity, 0.90)),
  peak_sparsity_p95 = as.numeric(quantile(peak_sparsity, 0.95))
)
sparsity_summary



library(ggplot2)

df_cell <- data.frame(cell_sparsity = as.numeric(cell_sparsity))
df_peak <- data.frame(peak_sparsity = as.numeric(peak_sparsity))

p1 <- ggplot(df_cell, aes(cell_sparsity)) +
  geom_histogram(bins = 60) +
  labs(title = "Cell sparsity distribution",
       x = "Cell sparsity (proportion of zeros across peaks)",
       y = "Number of cells")

p2 <- ggplot(df_peak, aes(peak_sparsity)) +
  geom_histogram(bins = 60) +
  labs(title = "Peak sparsity distribution",
       x = "Peak sparsity (proportion of zeros across cells)",
       y = "Number of peaks")

p1; p2




write.csv(df_cell, "cell_sparsity.csv", row.names = FALSE)
write.csv(df_peak, "peak_sparsity.csv", row.names = FALSE)
write.csv(sparsity_summary, "sparsity_summary.csv", row.names = FALSE)

ggsave("D:/SCAVENGE-reproducibility-main/sparsity_output/fig_cell_sparsity_hist.png", p1, width=6, height=4, dpi=300)
ggsave("D:/SCAVENGE-reproducibility-main/sparsity_output/fig_peak_sparsity_hist.png", p2, width=6, height=4, dpi=300)
