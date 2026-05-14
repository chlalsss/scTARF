library(data.table)
library(dplyr)
library(ggplot2)
library(BuenColors)
library(igraph)

# =========================
# [新增] 0) 设置路径 + 选择 trait
# =========================
out_dir <- "D:/SCAVENGE-reproducibility-main/blood_data/scavenge_out_small_all_traits"
trait_to_plot <- "MCV"  # 你可以改成 "LYM" "EOS" "MONO" "NEU" "RBC" "HGB" "HCT" 等

in_file <- file.path(out_dir, paste0("SCAVENGE_small_", trait_to_plot, ".tsv.gz"))
mono_mat_raw <- data.table::fread(in_file)

# =========================
# [新增/改动] 1) 适配 getstarted：把 UMAP1/UMAP2/BioClassification 映射到 x/y/color
# =========================
mono_mat <- mono_mat_raw %>%
  dplyr::rename(
    x = UMAP1,                 # [改动] getstarted 直接用 x，这里用 UMAP1 替代
    y = UMAP2,                 # [改动] getstarted 直接用 y，这里用 UMAP2 替代
    color = BioClassification  # [改动] getstarted 用 color，这里用 BioClassification（cell type label）
  )

# [新增] 如果你要运行 permutation 那段：把 seed 列名补成 seed_idx（getstarted 用 seed_idx）
if ("seed" %in% names(mono_mat) && !"seed_idx" %in% names(mono_mat)) {
  mono_mat$seed_idx <- mono_mat$seed
}

# =========================
# 下面开始：getstarted 原样画图代码（不改）
# =========================

# UMAP plots of cell type annotation and cell-to-cell graph
# Cell type annotation
p <- ggplot(data=mono_mat, aes(x, y, color=color)) + 
  geom_point(size=1, na.rm = TRUE) + 
  pretty_plot() + 
  theme(legend.title = element_blank()) + 
  labs(x="UMAP 1",y="UMAP 2")
p

# Comparsion before and after SCAVENGE analysis
# Z score based visualization
# Scatter plot
p1 <- ggplot(data=mono_mat, aes(x, y, color=z_score)) + 
  geom_point(size=1, na.rm = TRUE, alpha = 0.6) + 
  scale_color_viridis_c() + 
  scale_alpha() +
  pretty_plot() + 
  theme(legend.title = element_blank()) + 
  labs(x="UMAP 1", y="UMAP 2")
p1

# Bar plot
pp1 <- ggplot(data=mono_mat,  aes(x=color, y=z_score))  +
  geom_boxplot(aes(fill=color, color=color), outlier.shape=NA) + 
  guides(fill=FALSE) + 
  pretty_plot(fontsize = 10) +
  stat_summary(geom = "crossbar", width=0.65, fatten=0, color="black", 
               fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) + 
  theme(legend.position = "none")
pp1

# SCAVENGE TRS based visualization
# Scatter plot
p2 <- ggplot(data=mono_mat, aes(x, y, color=TRS)) + 
  geom_point(size=1, na.rm = TRUE, alpha = 0.6) + 
  scale_color_viridis_c() + 
  scale_alpha() + pretty_plot() + 
  theme(legend.title = element_blank()) + 
  labs(x="UMAP 1", y="UMAP 2")
p2

# Bar plot
pp2 <- ggplot(data=mono_mat,  aes(x=color, y=TRS))  +
  geom_boxplot(aes(fill=color, color=color), outlier.shape=NA) + 
  guides(fill=FALSE) + 
  pretty_plot(fontsize = 10) +
  stat_summary(geom = "crossbar", width=0.65, fatten=0, color="black",
               fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) + 
  theme(legend.position = "none")
pp2

# === 保存图片（新增）===
fig_dir <- file.path(out_dir, "fig3_like")
dir.create(fig_dir, showWarnings = FALSE)

ggsave(file.path(fig_dir, paste0("Fig3_", trait_to_plot, "_a_celltype_umap.png")),
       plot = p, width = 7, height = 5, dpi = 300)

ggsave(file.path(fig_dir, paste0("Fig3_", trait_to_plot, "_b_zscore_umap.png")),
       plot = p1, width = 7, height = 5, dpi = 300)

ggsave(file.path(fig_dir, paste0("Fig3_", trait_to_plot, "_c_zscore_box.png")),
       plot = pp1, width = 10, height = 5, dpi = 300)

ggsave(file.path(fig_dir, paste0("Fig3_", trait_to_plot, "_d_trs_umap.png")),
       plot = p2, width = 7, height = 5, dpi = 300)

ggsave(file.path(fig_dir, paste0("Fig3_", trait_to_plot, "_e_trs_box.png")),
       plot = pp2, width = 10, height = 5, dpi = 300)