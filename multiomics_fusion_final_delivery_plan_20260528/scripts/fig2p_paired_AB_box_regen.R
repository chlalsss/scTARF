#!/usr/bin/env Rscript
# 重新生成 fig2P 配对 BMMC 主要性能图：
#   A  三方法分组柱状图 (ATAC-only baseline / scDRS-seed w=1.0 / paired fusion w=0.7)
#   B  拆分 Delta 堆叠柱 (scDRS 基因集效应 = w1.0-baseline ; 配对增量 = w0.7-w1.0)
#   3  横向 cell-type TRS 箱线图 (ATAC-only vs v2 fusion), MONO / PLT / LYM
# 只展示三个主要性状: MONO / PLT / LYM_wide(主报靶, 图中标 LYM)
# 数据源: D:/SCAVENGE-reproducibility-main/BMMC (配对实验真实输出) + 已验证 roc_table45

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

BMMC      <- "D:/SCAVENGE-reproducibility-main/BMMC"
ROC_TAB   <- "D:/RStudio/work/ppt_dataset_figures/tables/roc_table45_auc_check.csv"
OUT       <- "D:/RStudio/work/figure_generation_work_20260604/fig2p_assets"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

## ---- 配色与方法标签 ----
m_base <- "ATAC-only baseline"
m_scdrs <- "scDRS-seed (no RNA, w=1.0)"
m_fuse <- "Paired fusion (w=0.7)"
method_levels <- c(m_base, m_scdrs, m_fuse)
method_colors <- c(
  "#9AA0A6",  # baseline 灰
  "#4C78A8",  # scDRS-seed 蓝
  "#8E63CE"   # paired fusion 紫
)
names(method_colors) <- method_levels

trait_levels <- c("MONO", "PLT", "LYM")

## ================= A : 三方法分组柱状图 =================
roc <- fread(ROC_TAB)
# roc method 列: "ATAC-only baseline" / "v2 no RNA (w=1.0)" / "v2 fused (w=0.7)"
map_method <- c(
  "ATAC-only baseline" = m_base,
  "v2 no RNA (w=1.0)"  = m_scdrs,
  "v2 fused (w=0.7)"   = m_fuse
)
map_trait <- c("MONO" = "MONO", "PLT" = "PLT", "LYM_wide" = "LYM")
roc <- roc[target %in% names(map_trait)]
roc[, trait := factor(map_trait[target], levels = trait_levels)]
roc[, method := factor(map_method[method], levels = method_levels)]

pA <- ggplot(roc, aes(trait, auroc, fill = method)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.72, color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.3f", auroc)),
            position = position_dodge(width = 0.78),
            vjust = -0.4, size = 3.4, color = "#333333") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "#B0B0B0", linewidth = 0.4) +
  scale_fill_manual(values = method_colors, name = NULL) +
  scale_y_continuous(limits = c(0, 1.04), breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
  labs(x = NULL, y = "AUROC (cell-type recognition)",
       title = "Paired BMMC: per-trait AUROC") +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(face = "bold", size = 15),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "top",
        legend.text = element_text(size = 11),
        axis.text = element_text(size = 12),
        axis.title.y = element_text(size = 12))
ggsave(file.path(OUT, "fig2p_A_auroc.png"), pA, width = 8.4, height = 5.0, dpi = 300, bg = "white")

## ================= B : 拆分 Delta 堆叠柱 =================
wide <- dcast(roc, trait ~ method, value.var = "auroc")
setnames(wide, c(m_base, m_scdrs, m_fuse), c("base", "scdrs", "fuse"))
delta <- rbindlist(list(
  data.table(trait = wide$trait, component = "scDRS gene-set effect (w1.0 - baseline)",
             value = wide$scdrs - wide$base),
  data.table(trait = wide$trait, component = "Pairing increment (w0.7 - w1.0)",
             value = wide$fuse - wide$scdrs)
))
comp_levels <- c("scDRS gene-set effect (w1.0 - baseline)", "Pairing increment (w0.7 - w1.0)")
delta[, component := factor(component, levels = comp_levels)]
delta[, trait := factor(trait, levels = trait_levels)]
comp_colors <- setNames(c("#4C78A8", "#8E63CE"), comp_levels)

# 顶部两行标签: 总增益 + 拆分(基因集 / 配对); 段太小不放段内白字避免重叠
tot <- data.table(trait = wide$trait,
                  set = wide$scdrs - wide$base,
                  pair = wide$fuse - wide$scdrs)
tot[, total := set + pair]
tot[, trait := factor(trait, levels = trait_levels)]

pB <- ggplot(delta, aes(trait, value, fill = component)) +
  geom_col(width = 0.62, color = "white", linewidth = 0.3) +
  geom_text(data = tot,
            aes(trait, total, label = sprintf("Δ %.3f\n(set +%.3f · pair +%.3f)", total, set, pair)),
            inherit.aes = FALSE, vjust = -0.35, size = 3.5, fontface = "bold",
            color = "#222222", lineheight = 0.9) +
  scale_fill_manual(values = comp_colors, name = NULL) +
  scale_y_continuous(limits = c(0, max(tot$total) * 1.24), expand = c(0, 0)) +
  labs(x = NULL, y = "Δ AUROC vs baseline (decomposed)",
       title = "Paired BMMC: improvement attribution") +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(face = "bold", size = 15),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "top",
        legend.direction = "vertical",
        legend.text = element_text(size = 10),
        axis.text = element_text(size = 12),
        axis.title.y = element_text(size = 12))
ggsave(file.path(OUT, "fig2p_B_delta.png"), pB, width = 8.4, height = 5.0, dpi = 300, bg = "white")

## ================= 3 : 横向 cell-type TRS 箱线图 =================
meta <- fread(file.path(BMMC, "00_precheck", "BMMC_meta_raw.csv"))
baseline_mat <- readRDS(file.path(BMMC, "03_scavenge_baseline", "TRS_matrix_baseline.rds"))
if (inherits(baseline_mat, "data.table")) baseline_mat <- as.data.frame(baseline_mat)
if ("cell" %in% colnames(baseline_mat)) {
  rownames(baseline_mat) <- baseline_mat$cell; baseline_mat$cell <- NULL
}
minmax <- function(x) { r <- range(x, na.rm = TRUE); if (diff(r) == 0) return(rep(0, length(x))); (x - r[1]) / diff(r) }

box_traits <- c("MONO" = "MONO", "PLT" = "PLT", "LYM" = "LYM")  # baseline_mat 列名
for (tr in names(box_traits)) {
  col <- box_traits[[tr]]
  v2_file <- file.path(BMMC, "04_dual_modality", paste0("TRS_v2_", col, "_w0.7.tsv.gz"))
  if (!col %in% colnames(baseline_mat) || !file.exists(v2_file)) {
    cat("skip", tr, "\n"); next
  }
  bline <- data.table(barcode = rownames(baseline_mat), trs = minmax(baseline_mat[, col]), method = "ATAC-only")
  v2raw <- fread(v2_file)
  v2dt  <- data.table(barcode = v2raw$barcode, trs = minmax(v2raw$TRS), method = "v2 fusion")
  comb  <- merge(rbind(bline, v2dt), meta[, .(barcode, cell_type)], by = "barcode")

  # 按 v2 fusion 中位数排序; coord_flip 后高值在顶
  ord <- comb[method == "v2 fusion", .(med = median(trs, na.rm = TRUE)), by = cell_type][order(med)]$cell_type
  comb[, cell_type := factor(cell_type, levels = ord)]
  comb[, method := factor(method, levels = c("ATAC-only", "v2 fusion"))]

  p <- ggplot(comb, aes(cell_type, trs, fill = method)) +
    geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.25, linewidth = 0.35, width = 0.7) +
    scale_fill_manual(values = c("ATAC-only" = "#AAC4DE", "v2 fusion" = "#8E63CE"), name = NULL) +
    coord_flip() +
    labs(title = paste0(tr, " — cell-type TRS"), x = NULL,
         y = "Normalized TRS (min-max within method)") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 14),
          axis.text.y = element_text(size = 9.5),
          axis.text.x = element_text(size = 10),
          axis.title.x = element_text(size = 10.5),
          legend.position = "top",
          legend.text = element_text(size = 10),
          panel.grid.minor = element_blank())
  ggsave(file.path(OUT, paste0("fig2p_box_", tr, ".png")), p, width = 6.2, height = 7.6, dpi = 220, bg = "white")
  cat("saved box", tr, "\n")
}

cat("ALL DONE ->", OUT, "\n")
