#!/usr/bin/env Rscript
# fig5: 用 MONO 演示 scTARF 全流程的中间变量面板。
# Track P 配对 BMMC multiome (module3): P1 G_trait / P2 z_rna / P3 seed / P4 TRS / P5 AUROC
# Track N 非配对 Granja (module4):    N1 RNA选择 / N2 mapping / N3 伪配对 / N4 seed / N5 TRS / N6 AUROC
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

B3 <- "D:/SCAVENGE-reproducibility-main/BMMC/04_dual_modality"
M4 <- "D:/RStudio/work/module4_traitaware_pseudopair"
OUT <- "D:/RStudio/work/figure_generation_work_20260604/fig5_assets"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
W <- 5; H <- 4.4; DPI <- 200
COL_base <- "#9AA0A6"; COL_fus <- "#8E63CE"; COL_atac <- "#4C78A8"; COL_grn <- "#59A14F"
th  <- theme_bw(base_size = 12) + theme(plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 9.5, color = "#555555"), panel.grid.minor = element_blank())
thv <- theme_void(base_size = 12) + theme(plot.title = element_text(face = "bold", size = 14, hjust = 0),
        plot.subtitle = element_text(size = 9.3, color = "#555555", hjust = 0),
        legend.key.size = unit(0.4, "cm"), legend.title = element_text(size = 10), legend.text = element_text(size = 8))
sv <- function(p, n) ggsave(file.path(OUT, n), p, width = W, height = H, dpi = DPI, bg = "white")

## ===================== Track P (paired) =====================
# P1 G_trait
g <- as.data.table(readRDS(file.path(B3, "Gtrait_MONO.rds")))
top5 <- g[rank <= 5]
pP1 <- ggplot(g, aes(rank, gene_score)) +
  geom_line(color = COL_atac, linewidth = 0.8) + geom_point(color = COL_atac, size = 0.7) +
  geom_point(data = top5, color = COL_fus, size = 2.3) +
  geom_text(data = g[rank == 1], aes(label = gene_symbol), nudge_x = 7, hjust = 0, color = COL_fus, fontface = "bold", size = 3.4) +
  labs(title = "P1 · G_trait 构建", subtitle = "100 genes (kneedle); top5: CNN2,NCAPD2,MRPL51,TMEM259,ENDOD1",
       x = "gene rank", y = "gene score (norm)") + th
sv(pP1, "P1_gtrait.png")

# P2 scDRS z_rna
z <- as.numeric(readRDS(file.path(B3, "zrna_MONO.rds"))); dz <- data.table(z = z)
qz <- quantile(dz$z, c(.95, .99)); mz <- mean(dz$z)
pP2 <- ggplot(dz, aes(z)) + geom_histogram(bins = 80, fill = COL_atac, alpha = 0.85) +
  geom_vline(xintercept = mz, color = "#333333", linewidth = 0.5) +
  geom_vline(xintercept = qz, linetype = "dashed", color = "#C81D25", linewidth = 0.5) +
  annotate("text", x = qz[2], y = Inf, label = "q99", vjust = 1.6, hjust = -0.1, size = 3, color = "#C81D25") +
  coord_cartesian(xlim = c(-4, 8)) +
  labs(title = "P2 · scDRS z_rna", subtitle = "69,249 cells; mean≈0.14, heavy right tail (max 13.6)",
       x = "scDRS z (MONO gene set)", y = "cells") + th
sv(pP2, "P2_zrna.png")

# P3 seed UMAP
um <- as.data.table(readRDS(file.path(B3, "figures/bmmc_atac_umap_coords.rds")))
sd <- readRDS(file.path(B3, "seed_v2_MONO.rds")); sdt <- data.table(barcode = names(sd), seed = as.integer(sd))
d3 <- merge(um, sdt, by = "barcode")
pP3 <- ggplot() +
  geom_point(data = d3[seed == 0], aes(UMAP1, UMAP2), color = "#DDDDDD", size = 0.3, alpha = 0.5) +
  geom_point(data = d3[seed == 1], aes(UMAP1, UMAP2), color = COL_fus, size = 0.8) +
  labs(title = "P3 · seed 选择", subtitle = "SCDRS norm_score top 1% → 684 seeds (thr=3.044)") + thv
sv(pP3, "P3_seed.png")

# P4 RWR -> TRS UMAP
trs <- fread(file.path(B3, "TRS_v2_MONO_w0.7.tsv.gz")); setnames(trs, c("barcode", "TRS"))
d4 <- merge(um, trs, by = "barcode"); q4 <- quantile(d4$TRS, .99); d4[, TRSc := pmin(TRS, q4)]
pP4 <- ggplot(d4[sample(.N)], aes(UMAP1, UMAP2, color = TRSc)) + geom_point(size = 0.4, alpha = 0.6) +
  scale_color_gradient(low = "navy", high = "gold", name = "TRS") +
  labs(title = "P4 · RWR → TRS", subtitle = "paired fusion TRS (atac_w=0.7)") + thv
sv(pP4, "P4_trs.png")

# P5 AUROC
ev <- fread(file.path(B3, "eval_v2_auroc.csv"))
d5 <- data.table(method = c("ATAC-only\nbaseline", "paired fusion\n(v2, w=0.7)"),
                 auroc = c(ev[target == "MONO" & method == "baseline_atac", auroc][1],
                           ev[target == "MONO" & method == "v2" & abs(atac_w - 0.7) < 1e-6, auroc][1]))
d5[, method := factor(method, levels = method)]
pP5 <- ggplot(d5, aes(method, auroc, fill = method)) + geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%.3f", auroc)), vjust = -0.4, size = 4.2) +
  scale_fill_manual(values = c(COL_base, COL_fus), guide = "none") +
  scale_y_continuous(limits = c(0, 1.06), breaks = seq(0, 1, .25), expand = c(0, 0)) +
  labs(title = "P5 · AUROC", subtitle = "MONO 0.532 → 0.976 (Δ+0.444)", x = NULL, y = "AUROC") +
  th + theme(axis.text.x = element_text(size = 10))
sv(pP5, "P5_auroc.png")

## ===================== Track N (non-paired) =====================
# N1 RNA reference selection
rs <- fread(file.path(M4, "01_data_prep_rna_reference_screening/output/rna_suitability_summary.csv"))
rsm <- rs[trait == "MONO", .(reference, gene_coverage_score, rna_signal_score, rna_suitability_score)]
rl <- melt(rsm, id.vars = "reference", variable.name = "metric", value.name = "score")
rl[, metric := factor(metric, levels = c("gene_coverage_score", "rna_signal_score", "rna_suitability_score"),
                      labels = c("gene\ncoverage", "RNA\nsignal", "suitability\n(final)"))]
rl[, reference := factor(ifelse(reference == "granja_rna", "Granja", "CITE (用于本次 MONO)"),
                         levels = c("Granja", "CITE (用于本次 MONO)"))]
pN1 <- ggplot(rl, aes(metric, score, fill = reference)) +
  geom_col(position = position_dodge(0.75), width = 0.7) +
  geom_text(aes(label = sprintf("%.1f", score)), position = position_dodge(0.75), vjust = -0.3, size = 3) +
  geom_hline(yintercept = 70, linetype = "dotted", color = "#999999") +
  scale_fill_manual(values = c("Granja" = "#C7C7C7", "CITE (用于本次 MONO)" = COL_grn), name = NULL) +
  scale_y_continuous(limits = c(0, 108), expand = c(0, 0)) +
  labs(title = "N1 · RNA 参考选择", subtitle = "Granja/CITE 均过适宜性门槛；MONO 映射采用 CITE RNA (90k cells)", x = NULL, y = "score") +
  th + theme(axis.text.x = element_text(size = 9), legend.position = "top", legend.text = element_text(size = 9))
sv(pN1, "N1_rna_select.png")

# N2 / N3 pseudo transfer
pz <- fread(file.path(M4, "03_pseudopair_pipeline/output/v3_sparse_mapping_top5significant/transfer_mono/granja_pseudo_zrna_MONO.csv"))
pN2 <- ggplot(pz, aes(top1_cos_dist)) + geom_histogram(bins = 60, fill = COL_atac, alpha = 0.85) +
  labs(title = "N2 · 稀疏 trait mapping", subtitle = "top-1 RNA 邻居 cosine 距离; k=5, 保留基因中位20 (密度3.2%)",
       x = "top-1 neighbor cosine distance", y = "ATAC cells") + th
sv(pN2, "N2_mapping.png")

pN3 <- ggplot(pz, aes(pseudo_z_rna)) + geom_histogram(bins = 70, fill = COL_fus, alpha = 0.85) +
  geom_vline(xintercept = median(pz$pseudo_z_rna), linetype = "dashed", color = "#333333", linewidth = 0.5) +
  labs(title = "N3 · 伪配对 (pseudo-pair)", subtitle = sprintf("伪 scDRS z = 5-NN 加权; 中位 %.2f, 右长尾", median(pz$pseudo_z_rna)),
       x = "pseudo z_rna", y = "ATAC cells") + th
sv(pN3, "N3_pseudopair.png")

# N4 neighborhood-FDR seed
d4n <- data.table(method = c("旧: top-1%", "邻域-FDR (≤0.1)"), n = c(350, 3861))
d4n[, method := factor(method, levels = method)]
pN4 <- ggplot(d4n, aes(method, n, fill = method)) + geom_col(width = 0.6) +
  geom_text(aes(label = n), vjust = -0.3, size = 4.2) +
  scale_fill_manual(values = c("#C7C7C7", COL_fus), guide = "none") +
  scale_y_continuous(limits = c(0, 4300), expand = c(0, 0)) +
  labs(title = "N4 · 邻域-FDR 种子", subtitle = "seed 350 → 3861 (覆盖更充分)", x = NULL, y = "n seed") +
  th + theme(axis.text.x = element_text(size = 10))
sv(pN4, "N4_seed.png")

# N5 TRS UMAP (Granja)
TRS_DIR <- file.path(M4, "04_evaluation/baseline_original_scavenge_robustness/module4_original_scavenge_robust_remaining_traits_seed9527_traitfilter_fixed/trs")
base <- fread(file.path(TRS_DIR, "SCAVENGE_ATAC_p0.00_seed1_MONO.tsv.gz"))
fus <- fread(file.path(M4, "03_pseudopair_pipeline/output/v3_sparse_mapping_top5significant/module3_mono_neighborhood_fdr_ge0p80/TRS_v2_MONO_w1.0.tsv.gz"))
setnames(fus, c("barcode", "TRS_fus"))
d5n <- merge(base[, .(cell, UMAP1, UMAP2)], fus, by.x = "cell", by.y = "barcode")
q5 <- quantile(d5n$TRS_fus, .99); d5n[, TRSc := pmin(TRS_fus, q5)]
pN5 <- ggplot(d5n[sample(.N)], aes(UMAP1, UMAP2, color = TRSc)) + geom_point(size = 0.4, alpha = 0.6) +
  scale_color_gradient(low = "navy", high = "gold", name = "TRS") +
  labs(title = "N5 · RWR → TRS", subtitle = "非配对 fusion TRS (Granja UMAP, w=1.0)") + thv
sv(pN5, "N5_trs.png")

# N6 AUROC
au <- fread(file.path(M4, "04_evaluation/v3_sparse_mapping_top5significant/mono_neighborhood_fdr_ge0p80/v3_sparse_mapping_top5significant_mono_neighborhood_fdr_ge0p80_auroc.csv"))
d6 <- data.table(method = c("ATAC-only\nbaseline", "trait-aware\nfusion (w=1.0)"),
                 auroc = c(au[target == "MONO" & method == "baseline_atac", auroc][1],
                           au[target == "MONO" & abs(atac_w - 1) < 1e-6, auroc][1]))
d6[, method := factor(method, levels = method)]
pN6 <- ggplot(d6, aes(method, auroc, fill = method)) + geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%.3f", auroc)), vjust = -0.4, size = 4.2) +
  scale_fill_manual(values = c(COL_base, COL_fus), guide = "none") +
  scale_y_continuous(limits = c(0, 1.06), breaks = seq(0, 1, .25), expand = c(0, 0)) +
  labs(title = "N6 · AUROC", subtitle = "MONO 0.863 → 0.918 (Δ+0.055)", x = NULL, y = "AUROC") +
  th + theme(axis.text.x = element_text(size = 10))
sv(pN6, "N6_auroc.png")

cat("ALL DONE ->", OUT, "\n")
