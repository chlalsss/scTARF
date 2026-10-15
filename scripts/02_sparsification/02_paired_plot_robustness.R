# 04_plot_robustness.R
# §3 阶段 3 主图: ATAC sparsification AUROC sensitivity curves.
# 输入: atac_sparsify_eval.csv (由 01_atac_sparsify.R 写出)
# 输出: figures/atac_sparsify.png
#       figures/atac_fusion_advantage.png
#       atac_sparsify_summary.csv (聚合表)
#
# 脚本功能：
# 读取 ATAC sparsification 压力测试的逐次评估结果。
# 按 dropout 比例、trait、target 和 atac_w 聚合 AUROC。
# 输出 AUROC 曲线图，以及融合方法相对 ATAC-only 的 AUROC 差值图。

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

source(file.path("scripts", "_paths.R"))
OUT_DIR  <- file.path(sctarf_bmmc_root(), "05_robustness")
FIG_DIR  <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

EVAL_FILE <- file.path(OUT_DIR, "atac_sparsify_eval.csv")

# ============================================================================
# 1. load + filter + dedup
# ============================================================================
dt <- fread(EVAL_FILE)
cat("loaded:", nrow(dt), "rows\n")

# 主报靶 + 敏感性分析靶
# - LYM_wide (NK+T+B 全淋巴谱系): 教科书定义, 独立于 SCDRS, 无自证 -> 主报
# - LYM_TNK (NK+T): SCDRS group_analysis 衍生定义, 含自证 -> 敏感性分析
# - MONO/PLT 评估靶基于 SCDRS group_analysis 命中, 但跟教科书一致 (单核、巨核-红祖)
keep_pair <- paste(dt$trait, dt$target, sep = "/")
dt <- dt[keep_pair %in% c("MONO/MONO", "PLT/PLT", "LYM/LYM_wide", "LYM/LYM_TNK")]
cat("after filter (4 targets: 3 main + 1 sensitivity):", nrow(dt), "rows\n")

# dedup (csv unique 没去掉的浮点重复)
# 同一个 p_drop/seed/trait/target/atac_w 组合只保留一条记录，
# 避免重复行影响均值和置信区间。
dt <- unique(dt, by = c("p_drop", "seed", "trait", "target", "atac_w"))
cat("after dedup:", nrow(dt), "rows\n")

# ============================================================================
# 2. aggregate: mean / sd / se / 95%CI per (p_drop, target, atac_w)
# ============================================================================
# 对多个 seed 重复实验取均值，并计算标准误和 95% 置信区间。
agg <- dt[, .(
  mean_auroc = mean(auroc),
  sd_auroc   = sd(auroc),
  n          = .N
), by = .(p_drop, trait, target, atac_w)]
agg[, se_auroc := sd_auroc / sqrt(n)]
agg[, lo95 := mean_auroc - 1.96 * se_auroc]
agg[, hi95 := mean_auroc + 1.96 * se_auroc]

# 加好看的 facet 标签
# panel 用于图中的分面标题；atac_w_label 用于图例。
target_label_map <- c(
  "MONO/MONO"    = "MONO  (CD14+/CD16+ Mono)",
  "PLT/PLT"      = "PLT  (MK/E prog)",
  "LYM/LYM_wide" = "LYM  (NK+T+B all lymphoid)",
  "LYM/LYM_TNK"  = "LYM_TNK  (NK+T, SCDRS-implicated)"
)
agg[, panel := target_label_map[paste(trait, target, sep = "/")]]
agg[, panel := factor(panel, levels = unname(target_label_map))]
agg[, atac_w_label := factor(
  ifelse(atac_w == 0.7, "Fused (atac_w=0.7)", "ATAC-only (atac_w=1.0)"),
  levels = c("Fused (atac_w=0.7)", "ATAC-only (atac_w=1.0)")
)]

write.csv(agg, file.path(OUT_DIR, "atac_sparsify_summary.csv"), row.names = FALSE)
cat("saved atac_sparsify_summary.csv (", nrow(agg), " rows)\n", sep = "")

# 拆分主图 (3 panel) vs 敏感性分析 (1 panel)
# 主图使用 MONO、PLT、LYM_wide；LYM_TNK 单独作为敏感性分析。
agg_main <- agg[paste(trait, target, sep = "/") %in%
                c("MONO/MONO", "PLT/PLT", "LYM/LYM_wide")]
agg_sens <- agg[paste(trait, target, sep = "/") == "LYM/LYM_TNK"]

# ============================================================================
# 3. main plot — AUROC vs p_drop, 3 panels (评估靶独立于 SCDRS, 无自证)
# ============================================================================
# 输入 d 是聚合后的 AUROC 表。
# 曲线表示不同 p_drop 下的平均 AUROC，阴影表示 95% CI。
make_main_plot <- function(d, title_suffix, subtitle_suffix) {
  ggplot(d, aes(x = p_drop, y = mean_auroc,
                color = atac_w_label, group = atac_w_label)) +
    geom_ribbon(aes(ymin = lo95, ymax = hi95, fill = atac_w_label),
                alpha = 0.18, color = NA) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.0) +
    scale_x_continuous(breaks = seq(0, 0.9, by = 0.1),
                       labels = scales::percent_format(accuracy = 1)) +
    scale_color_manual(values = c("Fused (atac_w=0.7)"     = "#8E5BA6",
                                  "ATAC-only (atac_w=1.0)" = "#3B6FB6")) +
    scale_fill_manual(values  = c("Fused (atac_w=0.7)"     = "#8E5BA6",
                                  "ATAC-only (atac_w=1.0)" = "#3B6FB6")) +
    facet_wrap(~ panel, nrow = 1, scales = "free_y") +
    labs(
      title    = paste0("BMMC ATAC sparsification — ", title_suffix),
      subtitle = paste0("3-seed mean (ribbon = 95% CI). Fused (purple) vs ATAC-only (blue). ",
                        subtitle_suffix),
      x = "ATAC read drop-out proportion",
      y = "AUROC",
      color = NULL, fill = NULL
    ) +
    theme_bw(base_size = 14) +
    theme(
      legend.position  = "bottom",
      legend.text      = element_text(size = 13),
      axis.text.x      = element_text(size = 11),
      axis.text.y      = element_text(size = 12),
      axis.title       = element_text(size = 13),
      plot.title       = element_text(size = 14, face = "bold"),
      plot.subtitle    = element_text(size = 11),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey92", color = NA),
      strip.text       = element_text(face = "bold", size = 12)
    )
}

p_main <- make_main_plot(agg_main,
  title_suffix    = "robustness of paired fusion (main)",
  subtitle_suffix = "Targets defined independently of SCDRS (textbook cell-type definitions).")
ggsave(file.path(FIG_DIR, "atac_sparsify.png"),
       p_main, width = 13, height = 5.0, dpi = 300)
cat("saved figures/atac_sparsify.png (3 panel main, no self-confirmation)\n")

# 敏感性分析图: LYM_TNK (含自证, 仅作补充)
p_sens <- make_main_plot(agg_sens,
  title_suffix    = "sensitivity analysis on SCDRS-implicated target",
  subtitle_suffix = "LYM_TNK target derived from SCDRS group_analysis (contains self-confirmation, supplementary).")
ggsave(file.path(FIG_DIR, "atac_sparsify_lym_tnk_sensitivity.png"),
       p_sens, width = 6, height = 5.0, dpi = 300)
cat("saved figures/atac_sparsify_lym_tnk_sensitivity.png (sensitivity, contains self-confirmation)\n")

# ============================================================================
# 4. fusion-advantage plot — Δ = (w=0.7) - (w=1.0) vs p_drop
# ============================================================================
# 将长表转换为宽表，使 w=0.7 和 w=1.0 的 AUROC 出现在同一行。
# delta > 0 表示融合方法优于 ATAC-only。
adv <- dcast(agg, p_drop + trait + target + panel ~ atac_w,
             value.var = c("mean_auroc", "sd_auroc", "n"))
adv[, delta := mean_auroc_0.7 - mean_auroc_1]
# pooled SE for the difference (assumes independence between w=0.7 and w=1.0,
# which is rough but conservative — they share the same KNN)
# se_delta 是差值的标准误，用于画 delta 的置信区间。
adv[, se_delta := sqrt((sd_auroc_0.7 / sqrt(n_0.7))^2 +
                       (sd_auroc_1   / sqrt(n_1))^2)]
adv[, lo95 := delta - 1.96 * se_delta]
adv[, hi95 := delta + 1.96 * se_delta]

# adv_d 是包含 delta 和置信区间的表。
# y=0 的虚线表示融合和 ATAC-only 持平。
make_adv_plot <- function(adv_d, title_suffix) {
  ggplot(adv_d, aes(x = p_drop, y = delta)) +
    geom_hline(yintercept = 0, color = "grey60", linetype = 2) +
    geom_ribbon(aes(ymin = lo95, ymax = hi95), fill = "#8E5BA6", alpha = 0.18) +
    geom_line(color = "#8E5BA6", linewidth = 0.9) +
    geom_point(color = "#8E5BA6", size = 2.0) +
    scale_x_continuous(breaks = seq(0, 0.9, by = 0.1),
                       labels = scales::percent_format(accuracy = 1)) +
    facet_wrap(~ panel, nrow = 1, scales = "free_y") +
    labs(
      title    = paste0("Fusion advantage under ATAC stress — ", title_suffix),
      subtitle = "Δ AUROC = Fused (atac_w=0.7) - ATAC-only (atac_w=1.0). Positive = fusion helps.",
      x        = "ATAC read drop-out proportion",
      y        = "Δ AUROC (Fused − ATAC-only)"
    ) +
    theme_bw(base_size = 14) +
    theme(
      axis.text.x      = element_text(size = 11),
      axis.text.y      = element_text(size = 12),
      axis.title       = element_text(size = 13),
      plot.title       = element_text(size = 14, face = "bold"),
      plot.subtitle    = element_text(size = 11),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey92", color = NA),
      strip.text       = element_text(face = "bold", size = 12)
    )
}

# 主图 Δ (3 panel, 主报)
adv_main <- adv[paste(trait, target, sep = "/") %in%
                c("MONO/MONO", "PLT/PLT", "LYM/LYM_wide")]
p_adv <- make_adv_plot(adv_main, "main targets (no self-confirmation)")
ggsave(file.path(FIG_DIR, "atac_fusion_advantage.png"),
       p_adv, width = 13, height = 4.5, dpi = 300)
cat("saved figures/atac_fusion_advantage.png (3 panel main)\n")

# 敏感性 Δ (LYM_TNK)
adv_sens <- adv[paste(trait, target, sep = "/") == "LYM/LYM_TNK"]
p_adv_sens <- make_adv_plot(adv_sens, "LYM_TNK sensitivity (contains self-confirmation)")
ggsave(file.path(FIG_DIR, "atac_fusion_advantage_lym_tnk_sensitivity.png"),
       p_adv_sens, width = 6, height = 4.5, dpi = 300)
cat("saved figures/atac_fusion_advantage_lym_tnk_sensitivity.png (sensitivity)\n")

# ============================================================================
# 5. headline numbers (console)
# ============================================================================
# 控制台打印关键数值，便于不打开图片时直接查看每个 p_drop 下的差异。
cat("\n", strrep("=", 72), "\n", sep = "")
cat("HEADLINE — fusion advantage at key p_drop levels\n")
cat(strrep("=", 72), "\n", sep = "")
print(adv[order(panel, p_drop),
          .(panel, p_drop,
            Fused    = round(mean_auroc_0.7, 4),
            ATAConly = round(mean_auroc_1,   4),
            Delta    = round(delta,          4),
            CI_lo    = round(lo95,           4),
            CI_hi    = round(hi95,           4))])

cat("\n04_plot_robustness.R OK.\n")
