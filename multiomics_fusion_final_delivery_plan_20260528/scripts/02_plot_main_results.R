# Plot main AUROC comparison figures from the unified table.
# Run after scripts/01_build_master_result_tables.R.
#
# 现行标准 (2026-06):
#   - 目标性状: paired = MONO/PLT/LYM(LYM_wide 标为 LYM); non-paired = MONO/RDW/LYM
#   - 版本命名: scTARF x.x-关键词 (Baseline / scDRS 作为参考基线单独命名)
#   - non-paired 的 scDRS-only(RNA 基线)取自 RNA reference scDRS 同性状结果
#   - 方法"做法"留到流程图, 这里只展示名称与结果

source("D:/RStudio/work/multiomics_fusion_final_delivery_plan_20260528/scripts/_final_delivery_common.R")
ensure_final_dirs()
need_packages(c("ggplot2"))
library(ggplot2)

in_file <- file.path(FINAL_ROOT, "01_master_tables/main_result_comparison_table.csv")
if (!file.exists(in_file)) {
  stop("Missing main table. Run scripts/01_build_master_result_tables.R first.", call. = FALSE)
}
out_dir  <- file.path(FINAL_ROOT, "02_main_figures")
fig4_dir <- "D:/RStudio/work/figure_generation_work_20260604/fig4_main_assets"
dir.create(fig4_dir, recursive = TRUE, showWarnings = FALSE)

x <- read_csv_if_exists(in_file)
x <- x[x$metric == "AUROC" & is.finite(x$value), , drop = FALSE]
if (!"final_version" %in% names(x)) {
  x$final_version <- classify_final_version(
    x$method_version, x$method, x$feature_mode, x$method_group, x$pairing_setting
  )
}

## ---------- 版本命名与配色 (scTARF x.x-关键词) ----------
keep_versions <- c("Baseline", "scDRS-only", "Paired-Fusion",
                   "v2.1", "v2.2", "v3.1", "v3.2", "v3.3")
new_label <- c(
  "Baseline"      = "Baseline (ATAC-only)",
  "scDRS-only"    = "scDRS-only (RNA)",
  "Paired-Fusion" = "scTARF1.0-paired",
  "v2.1"          = "scTARF2.1-pseudopair",
  "v2.2"          = "scTARF2.2-nbFDR",
  "v3.1"          = "scTARF3.1-top1",
  "v3.2"          = "scTARF3.2-top5sig",
  "v3.3"          = "scTARF3.3-sparse"
)
label_levels <- unname(new_label[keep_versions])
version_colors <- c(
  "Baseline (ATAC-only)"   = "#9AA0A6",
  "scDRS-only (RNA)"       = "#E8A33D",
  "scTARF1.0-paired"       = "#8E63CE",
  "scTARF2.1-pseudopair"   = "#4E9C8F",
  "scTARF2.2-nbFDR"        = "#4C78A8",
  "scTARF3.1-top1"         = "#6FB1E0",
  "scTARF3.2-top5sig"      = "#D17DB0",
  "scTARF3.3-sparse"       = "#C44E4E"
)

PAIR_LAB <- "Paired BMMC multiome"
NONP_LAB <- "Non-paired Granja ATAC + RNA reference"
cols <- c("pairing_setting", "dataset", "trait", "final_version", "value",
          "method_group", "source_file", "note")
sub  <- function(d) d[, cols, drop = FALSE]

## ---------- 目标性状筛选 ----------
# paired: MONO / PLT / LYM_wide(显示 LYM)
pr <- x[x$pairing_setting == "paired" &
        x$trait %in% c("MONO", "PLT", "LYM_wide") &
        x$final_version %in% keep_versions, , drop = FALSE]
pr <- sub(pr)
pr$trait_disp  <- ifelse(pr$trait == "LYM_wide", "LYM", pr$trait)
pr$setting_lab <- PAIR_LAB

# non-paired fusion: MONO / RDW / LYM
nf <- x[x$pairing_setting == "non_paired" &
        x$trait %in% c("MONO", "RDW", "LYM") &
        x$final_version %in% keep_versions, , drop = FALSE]
nf <- sub(nf)
# non-paired RNA baseline: 取 RNA reference scDRS 同性状结果, 归到 non-paired
ns <- x[x$pairing_setting == "rna_reference_only" &
        x$final_version == "scDRS-only" &
        x$trait %in% c("MONO", "RDW", "LYM"), , drop = FALSE]
ns <- sub(ns)
np <- rbind(nf, ns)
np$trait_disp  <- np$trait
np$setting_lab <- NONP_LAB

plot_df <- rbind(pr, np)
plot_df <- best_by_key(plot_df,
                       keys = c("setting_lab", "trait_disp", "final_version"))
plot_df$display_method <- factor(unname(new_label[plot_df$final_version]),
                                 levels = label_levels)
plot_df$trait_disp  <- factor(plot_df$trait_disp,
                              levels = c("MONO", "PLT", "LYM", "RDW"))
plot_df$setting_lab <- factor(plot_df$setting_lab, levels = c(PAIR_LAB, NONP_LAB))

## ================= 图1: 各性状 AUROC 分组柱 =================
p1 <- ggplot(plot_df, aes(trait_disp, value, fill = display_method)) +
  geom_col(position = position_dodge2(width = 0.85, preserve = "single"),
           width = 0.78, color = "white", linewidth = 0.25) +
  geom_hline(yintercept = 0.5, linetype = "dashed",
             color = "#B0B0B0", linewidth = 0.4) +
  facet_wrap(~ setting_lab, scales = "free_x") +
  scale_fill_manual(values = version_colors, name = NULL, drop = FALSE) +
  scale_y_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.2),
                     expand = c(0, 0)) +
  labs(x = NULL, y = "AUROC", title = "Main AUROC comparison by trait") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 2))
ggsave(file.path(out_dir, "fig_main_auroc_by_trait.png"), p1,
       width = 11, height = 5.4, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "fig_main_auroc_by_trait.pdf"), p1,
       width = 11, height = 5.4)

## ================= 图2: 相对 Baseline 的 Delta =================
base_df <- plot_df[plot_df$final_version == "Baseline",
                   c("setting_lab", "trait_disp", "value")]
names(base_df)[3] <- "baseline_auroc"
delta_df <- merge(plot_df[plot_df$final_version != "Baseline", ],
                  base_df, by = c("setting_lab", "trait_disp"), all.x = TRUE)
delta_df$delta <- delta_df$value - delta_df$baseline_auroc
delta_df <- delta_df[is.finite(delta_df$delta), , drop = FALSE]

p2 <- ggplot(delta_df, aes(trait_disp, delta, fill = display_method)) +
  geom_hline(yintercept = 0, color = "grey40", linewidth = 0.4) +
  geom_col(position = position_dodge2(width = 0.85, preserve = "single"),
           width = 0.78, color = "white", linewidth = 0.25) +
  facet_wrap(~ setting_lab, scales = "free_x") +
  scale_fill_manual(values = version_colors, name = NULL, drop = FALSE) +
  labs(x = NULL, y = "AUROC delta vs Baseline (ATAC-only)",
       title = "Change relative to ATAC-only baseline") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 2))
ggsave(file.path(out_dir, "fig_main_delta_vs_baseline.png"), p2,
       width = 11, height = 5.4, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "fig_main_delta_vs_baseline.pdf"), p2,
       width = 11, height = 5.4)

## ================= 图3: 版本 x 性状 AUROC 热图 =================
heat <- plot_df
heat$panel_trait <- paste(ifelse(heat$setting_lab == PAIR_LAB, "paired", "non-paired"),
                          heat$trait_disp, sep = " / ")
trait_row_levels <- c("paired / MONO", "paired / PLT", "paired / LYM",
                      "non-paired / MONO", "non-paired / RDW", "non-paired / LYM")
grid <- expand.grid(display_method = label_levels,
                    panel_trait = trait_row_levels,
                    stringsAsFactors = FALSE)
heat <- merge(grid, heat[c("display_method", "panel_trait", "value")],
              by = c("display_method", "panel_trait"), all.x = TRUE, sort = FALSE)
heat$has_auroc <- is.finite(heat$value)
heat$display_method <- factor(heat$display_method, levels = label_levels)
heat$panel_trait <- factor(heat$panel_trait, levels = rev(trait_row_levels))

p3 <- ggplot(heat, aes(display_method, panel_trait, fill = value)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = ifelse(has_auroc, sprintf("%.3f", value), "")),
            size = 3.1) +
  scale_fill_gradient(low = "#F7FBFF", high = "#08519C", limits = c(0, 1),
                      na.value = "grey94", name = "AUROC") +
  labs(x = NULL, y = NULL, title = "Method-by-trait AUROC heatmap") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 35, hjust = 1),
        panel.grid = element_blank())
ggsave(file.path(out_dir, "fig_main_method_heatmap.png"), p3,
       width = 9, height = 5.4, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "fig_main_method_heatmap.pdf"), p3,
       width = 9, height = 5.4)

## ---------- 同步到 fig4 资产目录 + 落数据 ----------
for (f in c("fig_main_auroc_by_trait.png", "fig_main_delta_vs_baseline.png",
            "fig_main_method_heatmap.png")) {
  file.copy(file.path(out_dir, f), file.path(fig4_dir, f), overwrite = TRUE)
}
write_csv(plot_df,  file.path(out_dir, "main_auroc_plot_data.csv"))
write_csv(delta_df, file.path(out_dir, "main_delta_vs_baseline_plot_data.csv"))
write_csv(heat,     file.path(out_dir, "main_method_heatmap_plot_data.csv"))

cat("Wrote standardized main figures to:\n  ", out_dir, "\n  ", fig4_dir, "\n", sep = "")
