#!/usr/bin/env Rscript
# 按 module4 (non-paired Granja) 真实结果生成 fig2N 的 A-F 面板, 与 fig2P 对应。
# 性状 MONO / RDW / LYM。三方法: ATAC-only baseline / trait-aware seed(w=1.0) / +pseudo-pair RNA(w=0.7)
# E/F 鲁棒性只画已完成的 RDW / LYM; MONO fusion sweep 计算中, 标注占位。

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

M4  <- "D:/RStudio/work/module4_traitaware_pseudopair"
OUT <- "D:/RStudio/work/figure_generation_work_20260604/fig2n_assets"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

traits <- c("MONO", "RDW", "LYM")
m_base <- "ATAC-only baseline"
m_seed <- "Trait-aware seed (w=1.0)"
m_pair <- "+ pseudo-pair RNA (w=0.7)"
mlev <- c(m_base, m_seed, m_pair)
mcol <- setNames(c("#9AA0A6", "#4C78A8", "#8E63CE"), mlev)

## ---------- 读 AUROC (A/B) ----------
fus_files <- list(
  MONO = file.path(M4, "04_evaluation/v3_sparse_mapping_top5significant/mono_neighborhood_fdr_ge0p80/v3_sparse_mapping_top5significant_mono_neighborhood_fdr_ge0p80_auroc.csv"),
  RDW  = file.path(M4, "04_evaluation/v3_sparse_other_traits/neighborhood_fdr_ge0p80/v3_sparse_other_traits_neighborhood_fdr_ge0p80_auroc.csv"),
  LYM  = file.path(M4, "04_evaluation/v3_sparse_lym_granja_rna/neighborhood_fdr_ge0p80/v3_sparse_lym_granja_rna_neighborhood_fdr_ge0p80_auroc.csv")
)
base_au <- fread(file.path(M4, "04_evaluation/granja_baseline_auroc.csv"))

get_au <- function(tr) {
  f <- fread(fus_files[[tr]])
  # baseline 优先取融合文件内的 baseline_atac 行(同一评估口径); 缺则回退总表
  b <- f[method == "baseline_atac" & target == tr, auroc][1]
  if (is.na(b)) b <- base_au[method == "baseline_atac" & target == tr, auroc][1]
  w07 <- f[target == tr & abs(atac_w - 0.7) < 1e-6, auroc][1]
  w10 <- f[target == tr & abs(atac_w - 1.0) < 1e-6, auroc][1]
  c(base = b, seed = w10, pair = w07)
}
A <- rbindlist(lapply(traits, function(tr) {
  v <- get_au(tr)
  data.table(trait = tr,
             method = factor(c(m_base, m_seed, m_pair), levels = mlev),
             auroc = c(v["base"], v["seed"], v["pair"]))
}))
A[, trait := factor(trait, levels = traits)]

pA <- ggplot(A, aes(trait, auroc, fill = method)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.72, color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.3f", auroc)), position = position_dodge(width = 0.78),
            vjust = -0.4, size = 3.3, color = "#333333") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "#B0B0B0", linewidth = 0.4) +
  scale_fill_manual(values = mcol, name = NULL) +
  scale_y_continuous(limits = c(0, 1.04), breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
  labs(x = NULL, y = "AUROC (cell-type recognition)", title = "Non-paired Granja: per-trait AUROC") +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(face = "bold", size = 15), panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(), legend.position = "top",
        legend.text = element_text(size = 10.5), axis.text = element_text(size = 12))
ggsave(file.path(OUT, "fig2n_A_auroc.png"), pA, width = 8.4, height = 5.0, dpi = 300, bg = "white")

## ---------- B 拆分Δ (seed效应 + 伪配对增量) ----------
wide <- dcast(A, trait ~ method, value.var = "auroc")
setnames(wide, c(m_base, m_seed, m_pair), c("base", "seed", "pair"))
comp_levels <- c("Trait-aware seed effect (w1.0 - baseline)", "Pseudo-pairing increment (w0.7 - w1.0)")
B <- rbindlist(list(
  data.table(trait = wide$trait, component = comp_levels[1], value = wide$seed - wide$base),
  data.table(trait = wide$trait, component = comp_levels[2], value = wide$pair - wide$seed)
))
B[, component := factor(component, levels = comp_levels)]
B[, trait := factor(trait, levels = traits)]
ccol <- setNames(c("#4C78A8", "#8E63CE"), comp_levels)
net <- data.table(trait = wide$trait, net = wide$pair - wide$base)
net[, trait := factor(trait, levels = traits)]
B[, vj := ifelse(value >= 0, -0.4, 1.3)]

pB <- ggplot(B, aes(trait, value, fill = component)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.62, color = "white", linewidth = 0.3) +
  geom_hline(yintercept = 0, color = "#444444", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%+.3f", value), vjust = vj), position = position_dodge(width = 0.7),
            size = 3.1, color = "#333333", show.legend = FALSE) +
  scale_fill_manual(values = ccol, name = NULL) +
  labs(x = NULL, y = "Δ AUROC component", title = "Non-paired: improvement attribution") +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(face = "bold", size = 15), panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(), legend.position = "top", legend.direction = "vertical",
        legend.text = element_text(size = 10), axis.text = element_text(size = 12))
ggsave(file.path(OUT, "fig2n_B_delta.png"), pB, width = 8.4, height = 5.0, dpi = 300, bg = "white")

## ---------- C/D 需要 per-cell TRS ----------
TRS_DIR <- file.path(M4, "04_evaluation/baseline_original_scavenge_robustness",
                     "module4_original_scavenge_robust_remaining_traits_seed9527_traitfilter_fixed/trs")
fus_trs <- list(
  MONO = file.path(M4, "03_pseudopair_pipeline/output/v3_sparse_mapping_top5significant/module3_mono_neighborhood_fdr_ge0p80/TRS_v2_MONO_w1.0.tsv.gz"),
  RDW  = file.path(M4, "03_pseudopair_pipeline/output/v3_sparse_other_traits/module3_neighborhood_fdr_ge0p80/TRS_v2_RDW_w1.0.tsv.gz"),
  LYM  = file.path(M4, "03_pseudopair_pipeline/output/v3_sparse_lym_granja_rna/module3_neighborhood_fdr_ge0p80/TRS_v2_LYM_w1.0.tsv.gz")
)
minmax <- function(x) { r <- range(x, na.rm = TRUE); if (diff(r) == 0) return(rep(0, length(x))); (x - r[1]) / diff(r) }

trs_umap <- function(dt, col, title, low = "navy", high = "gold") {
  q <- quantile(dt[[col]], 0.99, na.rm = TRUE)
  d2 <- copy(dt); d2[, capval := pmin(get(col), q)]
  ggplot(d2[sample(.N)], aes(UMAP1, UMAP2, color = capval)) +
    geom_point(size = 0.45, alpha = 0.55) +
    scale_color_gradient(low = low, high = high, name = "TRS") +
    labs(title = title) + theme_void(base_size = 13) +
    theme(plot.title = element_text(size = 13, hjust = 0.5, face = "bold"),
          legend.key.size = unit(0.45, "cm"))
}
ct_umap <- function(dt) {
  ggplot(dt[sample(.N)], aes(UMAP1, UMAP2, color = cell_type)) +
    geom_point(size = 0.45, alpha = 0.55) + labs(title = "Cell type") +
    theme_void(base_size = 13) +
    theme(plot.title = element_text(size = 13, hjust = 0.5, face = "bold"),
          legend.title = element_blank(), legend.text = element_text(size = 7.5),
          legend.key.size = unit(0.32, "cm")) +
    guides(color = guide_legend(override.aes = list(size = 2.4), ncol = 2))
}

for (tr in traits) {
  base <- fread(file.path(TRS_DIR, sprintf("SCAVENGE_ATAC_p0.00_seed1_%s.tsv.gz", tr)))
  fus  <- fread(fus_trs[[tr]])
  setnames(fus, c("barcode", "TRS"), c("cell", "TRS_fus"))
  dt <- merge(base[, .(cell, cell_type, UMAP1, UMAP2, TRS_base = TRS)], fus[, .(cell, TRS_fus)],
              by = "cell", all.x = TRUE)
  dt[is.na(TRS_fus), TRS_fus := 0]
  # C: 三联
  pc <- ct_umap(dt) | trs_umap(dt, "TRS_base", paste0(tr, " - ATAC-only TRS")) |
        trs_umap(dt, "TRS_fus", paste0(tr, " - trait-aware fusion TRS"))
  ggsave(file.path(OUT, sprintf("fig2n_C_umap_%s.png", tr)), pc, width = 18, height = 6, dpi = 200, bg = "white")
  # D: 横向箱线 (ATAC-only vs fusion, minmax 各方法内)
  comb <- rbind(
    data.table(cell_type = dt$cell_type, trs = minmax(dt$TRS_base), method = "ATAC-only"),
    data.table(cell_type = dt$cell_type, trs = minmax(dt$TRS_fus),  method = "fusion (w=1.0)")
  )
  ord <- comb[method == "fusion (w=1.0)", .(med = median(trs, na.rm = TRUE)), by = cell_type][order(med)]$cell_type
  comb[, cell_type := factor(cell_type, levels = ord)]
  comb[, method := factor(method, levels = c("ATAC-only", "fusion (w=1.0)"))]
  pd <- ggplot(comb, aes(cell_type, trs, fill = method)) +
    geom_boxplot(outlier.size = 0.25, outlier.alpha = 0.2, linewidth = 0.32, width = 0.7) +
    scale_fill_manual(values = c("ATAC-only" = "#AAC4DE", "fusion (w=1.0)" = "#8E63CE"), name = NULL) +
    coord_flip() + labs(title = paste0(tr, " - cell-type TRS"), x = NULL,
                        y = "Normalized TRS (min-max within method)") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 13), axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 9.5), legend.position = "top",
          legend.text = element_text(size = 9.5), panel.grid.minor = element_blank())
  ggsave(file.path(OUT, sprintf("fig2n_D_box_%s.png", tr)), pd, width = 6.2, height = 7.6, dpi = 220, bg = "white")
  cat("C/D done", tr, "\n")
}

## ---------- E/F 鲁棒性 (RDW/LYM 已完成; MONO 计算中) ----------
base_drop <- fread(file.path(M4, "04_evaluation/baseline_original_scavenge_robustness",
                  "module4_original_scavenge_robust_remaining_traits_seed9527_traitfilter_fixed/scavenge_atac_dropout_summary.csv"))
base_drop <- base_drop[trait %in% traits, .(trait, p_drop, method = m_base, auroc = mean_auroc, lo95, hi95)]

rdw_fd <- fread(file.path(M4, "05_robustness/output/v3_sparse_cite_erythroid_raw_atac_full_pipeline_dropout/v3_sparse_cite_erythroid_raw_atac_full_pipeline_dropout_plot_data.csv"))
lym_fd <- fread(file.path(M4, "05_robustness/output/v3_sparse_lym_granja_rna_raw_atac_full_pipeline_dropout/v3_sparse_lym_granja_rna_raw_atac_full_pipeline_dropout_plot_data.csv"))
fus_fd <- rbind(rdw_fd[target == "RDW"], lym_fd[target == "LYM"], fill = TRUE)
fus_fd[, trait := target]
fus_fd[, method := ifelse(abs(atac_w - 1.0) < 1e-6, "Fused (w=1.0)", "Fused (w=0.7)")]
fus_fd <- fus_fd[, .(trait, p_drop, method, auroc, lo95, hi95)]

robust <- rbind(base_drop, fus_fd)
robust[, trait := factor(trait, levels = traits)]
rlev <- c(m_base, "Fused (w=1.0)", "Fused (w=0.7)")
rcol <- setNames(c("#4C78A8", "#C81D25", "#E8943A"), rlev)
robust[, method := factor(method, levels = rlev)]

# MONO 计算中标注
mono_note <- data.table(trait = factor("MONO", levels = traits), p_drop = 0.45, auroc = NA)

pE <- ggplot(robust, aes(p_drop, auroc, color = method, fill = method)) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.0) + geom_point(size = 1.6) +
  facet_wrap(~ trait, nrow = 1, scales = "free_y") +
  scale_color_manual(values = rcol, name = NULL, drop = FALSE) +
  scale_fill_manual(values = rcol, guide = "none", drop = FALSE) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "ATAC read drop-out proportion", y = "AUROC",
       title = "Non-paired ATAC sparsification robustness (MONO fusion sweep in progress)") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 14), legend.position = "top",
        strip.text = element_text(face = "bold", size = 12), panel.grid.minor = element_blank())
# 在 MONO 分面标注 fusion 计算中
pE <- pE + geom_text(data = data.table(trait = factor("MONO", levels = traits)),
                     aes(x = 0.45, y = Inf, label = "fusion sweep 计算中"),
                     inherit.aes = FALSE, vjust = 1.6, size = 3.4, color = "#999999")
ggsave(file.path(OUT, "fig2n_E_robustness.png"), pE, width = 13, height = 4.8, dpi = 300, bg = "white")

## F: 融合优势 ΔAUROC = Fused(w=1.0) - baseline, vs dropout (RDW/LYM)
adv <- merge(
  fus_fd[method == "Fused (w=1.0)", .(trait, p_drop, f_auroc = auroc, f_lo = lo95, f_hi = hi95)],
  base_drop[, .(trait, p_drop, b_auroc = auroc, b_lo = lo95, b_hi = hi95)],
  by = c("trait", "p_drop"))
adv[, delta := f_auroc - b_auroc]
adv[, dse := (abs(f_hi - f_lo) / 2 / 1.96)]
adv[, bse := (abs(b_hi - b_lo) / 2 / 1.96)]
adv[, dlo := delta - 1.96 * sqrt(dse^2 + bse^2)]
adv[, dhi := delta + 1.96 * sqrt(dse^2 + bse^2)]
adv[, trait := factor(trait, levels = traits)]

pF <- ggplot(adv, aes(p_drop, delta)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#888888") +
  geom_ribbon(aes(ymin = dlo, ymax = dhi), fill = "#C81D25", alpha = 0.15) +
  geom_line(color = "#C81D25", linewidth = 1.0) + geom_point(color = "#C81D25", size = 1.6) +
  facet_wrap(~ trait, nrow = 1, drop = FALSE) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "ATAC read drop-out proportion", y = "Δ AUROC (Fused w=1.0 - baseline)",
       title = "Fusion advantage under ATAC stress (MONO in progress)") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 14),
        strip.text = element_text(face = "bold", size = 12), panel.grid.minor = element_blank())
pF <- pF + geom_text(data = data.table(trait = factor("MONO", levels = traits), p_drop = 0.45, delta = 0),
                     aes(label = "计算中"), color = "#999999", size = 3.6)
ggsave(file.path(OUT, "fig2n_F_advantage.png"), pF, width = 13, height = 4.5, dpi = 300, bg = "white")

cat("ALL DONE ->", OUT, "\n")
print(A)
print(net)
