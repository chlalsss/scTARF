# plot_trs_umap_comparison.R
# 生成 4.4 节所需的 UMAP 对比图：
#   - 细胞类型注释 UMAP
#   - ATAC-only baseline TRS UMAP (每个 trait)
#   - MultiOmicsFusion 融合 TRS UMAP (每个 trait)
# 输出: D:/RStudio/work/scavenge_multiomics_fusion_reproduction/reproduction_outputs/BMMC/04_dual_modality/figures/

# 依赖：uwot, ggplot2, data.table, patchwork
# 如未安装，先在独立 session 运行:
# install.packages(c("uwot", "ggplot2", "data.table", "patchwork"))

library(uwot)
library(ggplot2)
library(data.table)
library(patchwork)

set.seed(42)

BMMC_DIR  <- "D:/RStudio/work/scavenge_multiomics_fusion_reproduction/reproduction_outputs/BMMC"
FIG_DIR   <- file.path(BMMC_DIR, "04_dual_modality", "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# 1. 加载 LSI 嵌入（30 维，来自 baseline）
lsi <- readRDS(file.path(BMMC_DIR, "03_scavenge_baseline", "atac_lsi.rds"))
# lsi 可能是 matrix 或 data.frame；行 = 细胞
if (is.data.frame(lsi)) lsi <- as.matrix(lsi[, -1])  # 去掉可能的 cell 列
cat("LSI dim:", paste(dim(lsi), collapse = " x "), "\n")

# 2. 计算 UMAP（仅做一次，所有 trait 共用）
cat("Computing UMAP from LSI...\n")
umap_coords <- uwot::umap(lsi, n_neighbors = 30, min_dist = 0.3,
                           n_components = 2, n_threads = 4)
umap_dt <- data.table(
  barcode = rownames(lsi),
  UMAP1   = umap_coords[, 1],
  UMAP2   = umap_coords[, 2]
)
saveRDS(umap_dt, file.path(FIG_DIR, "bmmc_atac_umap_coords.rds"))
cat("UMAP done. Rows:", nrow(umap_dt), "\n")

# 3. 加载 metadata（cell_type）
meta <- fread(file.path(BMMC_DIR, "00_precheck", "BMMC_meta_raw.csv"))
setnames(meta, "barcode", "barcode")

umap_dt <- merge(umap_dt, meta[, .(barcode, cell_type, donor)], by = "barcode")

# 4. 加载 baseline TRS
baseline_mat <- readRDS(file.path(BMMC_DIR, "03_scavenge_baseline", "TRS_matrix_baseline.rds"))
if (inherits(baseline_mat, "data.table")) baseline_mat <- as.data.frame(baseline_mat)
if ("cell" %in% colnames(baseline_mat)) {
  rownames(baseline_mat) <- baseline_mat$cell
  baseline_mat$cell <- NULL
}

# 5. 辅助函数：拼接 TRS 并生成 UMAP 散点图
make_trs_umap <- function(dt, trs_vec, title, low = "navy", high = "gold") {
  dt2 <- merge(dt, data.table(barcode = names(trs_vec), trs = trs_vec),
               by = "barcode", all.x = TRUE)
  dt2[is.na(trs), trs := 0]
  q99 <- quantile(dt2$trs, 0.99, na.rm = TRUE)
  dt2[trs > q99, trs := q99]

  ggplot(dt2[sample(.N)], aes(UMAP1, UMAP2, color = trs)) +
    geom_point(size = 0.6, alpha = 0.6) +
    scale_color_gradient(low = low, high = high, name = "TRS") +
    labs(title = title) +
    theme_void(base_size = 14) +
    theme(plot.title    = element_text(size = 14, hjust = 0.5, face = "bold"),
          legend.title  = element_text(size = 12),
          legend.text   = element_text(size = 11),
          legend.key.size = unit(0.5, "cm"))
}

make_celltype_umap <- function(dt, title = "Cell type") {
  dt2 <- dt[sample(.N)]
  ggplot(dt2, aes(UMAP1, UMAP2, color = cell_type)) +
    geom_point(size = 0.6, alpha = 0.6) +
    labs(title = title) +
    theme_void(base_size = 14) +
    theme(plot.title    = element_text(size = 14, hjust = 0.5, face = "bold"),
          legend.text   = element_text(size = 10),
          legend.title  = element_blank(),
          legend.key.size = unit(0.4, "cm")) +
    guides(color = guide_legend(override.aes = list(size = 3), ncol = 2))
}

# 6. 每个 trait 生成三联图
traits <- c("MONO", "PLT", "LYM")

for (tr in traits) {
  cat("Plotting", tr, "...\n")

  # baseline TRS
  if (!tr %in% colnames(baseline_mat)) {
    cat("  baseline TRS for", tr, "not found; skip\n"); next
  }
  bline <- setNames(baseline_mat[, tr], rownames(baseline_mat))

  # MultiOmicsFusion TRS (atac_w=0.7)
  MultiOmicsFusion_file <- file.path(BMMC_DIR, "04_dual_modality",
                       paste0("TRS_MultiOmicsFusion_", tr, "_w0.7.tsv.gz"))
  if (!file.exists(MultiOmicsFusion_file)) {
    cat("  MultiOmicsFusion TRS file not found:", MultiOmicsFusion_file, "; skip\n"); next
  }
  MultiOmicsFusion_dt  <- fread(MultiOmicsFusion_file)
  MultiOmicsFusion_trs <- setNames(MultiOmicsFusion_dt$TRS, MultiOmicsFusion_dt$barcode)

  p_ct   <- make_celltype_umap(umap_dt, "Cell type")
  p_base <- make_trs_umap(umap_dt, bline,   paste0(tr, " – ATAC-only TRS"))
  p_MultiOmicsFusion   <- make_trs_umap(umap_dt, MultiOmicsFusion_trs,  paste0(tr, " – MultiOmicsFusion fusion TRS"))

  combined <- p_ct | p_base | p_MultiOmicsFusion
  out_png <- file.path(FIG_DIR, paste0("trs_umap_comparison_", tr, ".png"))
  ggsave(out_png, combined, width = 18, height = 6, dpi = 300)
  cat("  saved:", out_png, "\n")
}

# 7. TRS 箱线图（baseline vs MultiOmicsFusion）
# baseline 和 MultiOmicsFusion 的 TRS 使用不同 scale factor，绝对值不可比。
# 在每种方法内部分别做 min-max 归一化到 [0,1]，再合并比较相对富集程度。
minmax <- function(x) {
  r <- range(x, na.rm = TRUE)
  if (diff(r) == 0) return(rep(0, length(x)))
  (x - r[1]) / diff(r)
}

for (tr in traits) {
  cat("Boxplot", tr, "...\n")
  if (!tr %in% colnames(baseline_mat)) next
  MultiOmicsFusion_file <- file.path(BMMC_DIR, "04_dual_modality",
                       paste0("TRS_MultiOmicsFusion_", tr, "_w0.7.tsv.gz"))
  if (!file.exists(MultiOmicsFusion_file)) next

  bline <- data.table(barcode = rownames(baseline_mat),
                      trs = minmax(baseline_mat[, tr]),
                      method = "ATAC-only")
  MultiOmicsFusion_raw <- fread(MultiOmicsFusion_file)
  MultiOmicsFusion_dt  <- data.table(barcode = MultiOmicsFusion_raw$barcode,
                       trs = minmax(MultiOmicsFusion_raw$TRS),
                       method = "MultiOmicsFusion fusion")

  comb <- rbind(bline, MultiOmicsFusion_dt)
  comb <- merge(comb, meta[, .(barcode, cell_type)], by = "barcode")

  # 按 MultiOmicsFusion fusion 各细胞类型中位数降序排列 x 轴
  ord <- comb[method == "MultiOmicsFusion fusion",
              .(med = median(trs, na.rm = TRUE)), by = cell_type][order(-med)]$cell_type
  comb[, cell_type := factor(cell_type, levels = ord)]

  p <- ggplot(comb, aes(cell_type, trs, fill = method)) +
    geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.4, linewidth = 0.4) +
    scale_fill_manual(values = c("ATAC-only" = "#aac4de", "MultiOmicsFusion fusion" = "#2c7fb8")) +
    labs(title = paste0(tr, ": normalized TRS by cell type (ATAC-only vs MultiOmicsFusion fusion)"),
         x = NULL, y = "Normalized TRS (min-max within method)") +
    theme_bw(base_size = 13) +
    theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 11),
          axis.text.y  = element_text(size = 11),
          axis.title.y = element_text(size = 12),
          plot.title   = element_text(size = 14, face = "bold"),
          legend.text  = element_text(size = 12),
          legend.title = element_blank(),
          legend.position = "top")

  out_png <- file.path(FIG_DIR, paste0("trs_boxplot_comparison_", tr, ".png"))
  ggsave(out_png, p, width = 16, height = 6, dpi = 300)
  cat("  saved:", out_png, "\n")
}

cat("\nAll plots done.\n")
