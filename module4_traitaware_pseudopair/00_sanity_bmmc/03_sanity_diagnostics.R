# module4 / Phase 0.3: sanity diagnostics + gate alpha
# ------------------------------------------------------------------
# 目的: 检查 Phase 0.2 mapping 算法没出严重 bug. 三件事:
#   (1) 输出格式 sanity (无 NaN, 权重和=1, isolated cell < 5%)
#   (2) cell-type 一致率 (核心): 三 trait 都 ≥ 50% → 关卡 α 通过
#   (3) 距离分布合理性 (Euclidean → cosine_dist 换算后中位在 [0.2, 0.8])
#
# 输入:
#   00_sanity_bmmc/output/bmmc_mapping_{MONO,PLT,LYM}.rds    (Phase 0.2)
#   00_sanity_bmmc/output/bmmc_meta_shuffled.rds              (Phase 0.1)
#   D:/SCAVENGE-reproducibility-main/BMMC/01_data_qc/BMMC_cell_meta.rds  (原始, ATAC barcode -> 真 cell_type)
#
# 输出:
#   00_sanity_bmmc/output/sanity_report.csv         主诊断 (3 trait × 检查项)
#   00_sanity_bmmc/output/mapping_quality_by_celltype.csv  按 cell_type 一致率
#   00_sanity_bmmc/output/figures/sanity_distance_distribution.png
#   00_sanity_bmmc/output/figures/sanity_celltype_consistency_heatmap.png
# ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(ggplot2)
})

OUT_DIR  <- "D:/RStudio/work/module4_traitaware_pseudopair/00_sanity_bmmc/output"
FIG_DIR  <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- run variant ----
# 跟 02_trait_mapping_bmmc.R 顶部 USE_WG 联动:
#   USE_WG = TRUE  -> SUFFIX = ""           (原 weighted 方案)
#   USE_WG = FALSE -> SUFFIX = "_noweight"  (A0 对照)
SUFFIX <- "_noweight"
cat(">>> Phase 0.3 reading variant: SUFFIX = '", SUFFIX, "'\n", sep = "")

META_ORIG_PATH <- "D:/SCAVENGE-reproducibility-main/BMMC/01_data_qc/BMMC_cell_meta.rds"
META_SHUF_PATH <- file.path(OUT_DIR, "bmmc_meta_shuffled.rds")

TRAITS         <- c("MONO", "PLT", "LYM")
ISOLATED_QUANT <- 0.95   # cell is "isolated" if top1 distance > Q95 of all top1
COS_LO <- 0.2; COS_HI <- 0.8

# ------------------------------------------------------------------
# load metas
# ------------------------------------------------------------------
cat("[1/5] Loading metas ...\n")

# Original meta: gives true cell_type per ATAC barcode (BMMC paired, ATAC bc = RNA bc)
meta_orig <- readRDS(META_ORIG_PATH)
meta_orig <- as.data.frame(meta_orig)
stopifnot("barcode" %in% colnames(meta_orig))
stopifnot("cell_type" %in% colnames(meta_orig))
rownames(meta_orig) <- meta_orig$barcode
cat("       meta_orig: ", paste(dim(meta_orig), collapse = " x "), "\n")

# Shuffled meta: rows aligned to rna_shuffled columns; gives neighbor cell_type
meta_shuf <- readRDS(META_SHUF_PATH)
meta_shuf <- as.data.frame(meta_shuf)
cat("       meta_shuf: ", paste(dim(meta_shuf), collapse = " x "), "\n")
shuf_ct  <- meta_shuf$cell_type      # length = N_rna_cells_shuffled

# Sanity: all cell_type values are character (drop factor if any)
shuf_ct <- as.character(shuf_ct)

all_celltypes <- sort(unique(c(as.character(meta_orig$cell_type), shuf_ct)))
cat("       union cell_types: ", length(all_celltypes), " types\n")

# ------------------------------------------------------------------
# helper: weighted-majority vote per ATAC cell
# ------------------------------------------------------------------
# nbr_types : N x K character matrix of neighbor cell_type
# nbr_w     : N x K numeric weight matrix
# returns   : character N (predicted cell_type)
weighted_majority <- function(nbr_types, nbr_w) {
  N <- nrow(nbr_types); K <- ncol(nbr_types)
  vapply(seq_len(N), function(i) {
    t <- nbr_types[i, ]; w <- nbr_w[i, ]
    # tapply sums weights by type; pick argmax (ties: first appearance)
    s <- tapply(w, t, sum)
    names(s)[which.max(s)]
  }, character(1))
}

# ------------------------------------------------------------------
# per-trait diagnostics
# ------------------------------------------------------------------
report_rows  <- list()
ct_long      <- list()       # for by-celltype CSV
dist_long    <- list()       # for distance density plot
confusion_dt <- list()       # for confusion heatmap

for (tname in TRAITS) {
  cat("\n[2/5] === Trait ", tname, " ===\n", sep = "")
  map_file <- file.path(OUT_DIR,
                        paste0("bmmc_mapping_", tname, SUFFIX, ".rds"))
  m <- readRDS(map_file)

  atac_bc       <- m$atac_barcodes
  nbr_idx       <- m$rna_neighbor_idx          # N × K integer
  nbr_dist      <- m$neighbor_distance         # N × K Euclidean
  nbr_w         <- m$neighbor_weight           # N × K softmax
  K             <- m$k
  N             <- length(atac_bc)

  # ---- (1) format sanity ----
  n_nan_dist  <- sum(!is.finite(nbr_dist))
  n_nan_w     <- sum(!is.finite(nbr_w))
  w_rowsum    <- rowSums(nbr_w)
  rowsum_ok   <- all(abs(w_rowsum - 1) < 1e-6)
  d_top1      <- nbr_dist[, 1]
  iso_thr     <- as.numeric(quantile(d_top1, ISOLATED_QUANT))
  n_iso       <- sum(d_top1 > iso_thr * 1.05)   # heuristic
  pct_iso     <- 100 * n_iso / N

  cat(sprintf("       format: NaN_dist=%d  NaN_w=%d  rowsum_ok=%s  isolated %.2f%%\n",
              n_nan_dist, n_nan_w, rowsum_ok, pct_iso))

  # ---- (2) distance distribution (Euclidean & cosine_dist) ----
  d_med      <- median(nbr_dist)
  d_q05      <- quantile(nbr_dist, 0.05); names(d_q05) <- NULL
  d_q95      <- quantile(nbr_dist, 0.95); names(d_q95) <- NULL
  cos_med    <- d_med^2 / 2
  cos_q05    <- d_q05^2 / 2
  cos_q95    <- d_q95^2 / 2
  dist_in_range <- (cos_med >= COS_LO) && (cos_med <= COS_HI)
  cat(sprintf("       cos_dist (top-5 all): median=%.3f  q05=%.3f  q95=%.3f  in [0.2, 0.8]: %s\n",
              cos_med, cos_q05, cos_q95, dist_in_range))

  # gather for density plot
  dist_long[[tname]] <- data.frame(
    trait      = tname,
    cos_dist   = as.vector(nbr_dist)^2 / 2,
    stringsAsFactors = FALSE
  )

  # ---- (3) cell-type consistency ----
  true_ct <- meta_orig[atac_bc, "cell_type"]      # length N
  stopifnot(!any(is.na(true_ct)))                  # ATAC bc all present in meta_orig

  # neighbor cell_types via shuffled meta
  nbr_types <- matrix(shuf_ct[as.integer(nbr_idx)],
                      nrow = N, ncol = K)

  cat("       running weighted majority vote (N = ", N, ") ...\n", sep = "")
  t0 <- Sys.time()
  pred_ct <- weighted_majority(nbr_types, nbr_w)
  cat("       vote elapsed: ",
      round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "s\n", sep = "")

  agreement_total <- mean(pred_ct == true_ct)
  # high-confidence subset: top1 distance < median(top1)
  hi_mask <- d_top1 < median(d_top1)
  agreement_hi    <- mean(pred_ct[hi_mask] == true_ct[hi_mask])

  cat(sprintf("       cell-type agreement: overall = %.3f  high-conf subset = %.3f\n",
              agreement_total, agreement_hi))

  # by-cell_type breakdown
  ct_dt <- data.table(
    trait     = tname,
    cell_type = true_ct,
    correct   = pred_ct == true_ct
  )
  bt <- ct_dt[, .(n_cells = .N,
                  agreement_rate = mean(correct)),
              by = .(cell_type)]
  bt[, trait := tname]
  setcolorder(bt, c("trait", "cell_type", "n_cells", "agreement_rate"))
  ct_long[[tname]] <- bt

  # confusion: agreement_rate per cell_type for heatmap (long format)
  conf_dt <- data.table(
    trait     = tname,
    true_ct   = true_ct,
    pred_ct   = pred_ct
  )[, .N, by = .(trait, true_ct, pred_ct)]
  confusion_dt[[tname]] <- conf_dt

  # ---- (4) report row ----
  pass_format     <- (n_nan_dist == 0) && (n_nan_w == 0) && rowsum_ok && (pct_iso < 5)
  pass_distance   <- dist_in_range
  pass_consistency <- agreement_total >= 0.5
  pass_all <- pass_format && pass_distance && pass_consistency

  report_rows[[tname]] <- data.frame(
    trait              = tname,
    n_atac             = N,
    n_g_trait_used     = m$n_g_trait_used,
    n_nan_dist         = n_nan_dist,
    n_nan_weight       = n_nan_w,
    weight_rowsum_ok   = rowsum_ok,
    pct_isolated       = round(pct_iso, 3),
    cos_dist_median    = round(cos_med, 4),
    cos_dist_q05       = round(cos_q05, 4),
    cos_dist_q95       = round(cos_q95, 4),
    distance_in_range  = pass_distance,
    agreement_overall  = round(agreement_total, 4),
    agreement_hi_conf  = round(agreement_hi,    4),
    pass_format        = pass_format,
    pass_distance      = pass_distance,
    pass_consistency   = pass_consistency,
    pass_alpha_gate    = pass_all,
    stringsAsFactors   = FALSE
  )
}

# ------------------------------------------------------------------
# write CSVs
# ------------------------------------------------------------------
cat("\n[3/5] Writing CSVs ...\n")
report_df <- do.call(rbind, report_rows)
report_df$variant <- if (SUFFIX == "") "weighted" else sub("^_", "", SUFFIX)
write.csv(report_df,
          file.path(OUT_DIR, paste0("sanity_report", SUFFIX, ".csv")),
          row.names = FALSE)

ct_df <- data.table::rbindlist(ct_long)
ct_df[, variant := if (SUFFIX == "") "weighted" else sub("^_", "", SUFFIX)]
data.table::fwrite(ct_df,
                   file.path(OUT_DIR,
                             paste0("mapping_quality_by_celltype",
                                     SUFFIX, ".csv")))

confusion_all <- data.table::rbindlist(confusion_dt)

# ------------------------------------------------------------------
# figures
# ------------------------------------------------------------------
cat("[4/5] Drawing figures ...\n")

# (a) distance distribution (cos_dist density per trait)
dist_df <- do.call(rbind, dist_long)
p1 <- ggplot(dist_df, aes(x = cos_dist, fill = trait, color = trait)) +
  geom_density(alpha = 0.3) +
  geom_vline(xintercept = c(COS_LO, COS_HI), linetype = 2, color = "gray40") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  labs(title = "Phase 0.2 top-5 neighbor cosine_dist (= d_euc^2 / 2)",
       subtitle = "dashed lines: gate alpha range [0.2, 0.8]",
       x = "cosine distance", y = "density") +
  theme_bw()
ggsave(file.path(FIG_DIR,
                  paste0("sanity_distance_distribution", SUFFIX, ".png")),
       p1, width = 7, height = 4.5, dpi = 150)

# (b) by-cell_type agreement bar chart
bt_all <- ct_df[order(trait, -agreement_rate)]
# keep cell_types with >= 50 cells
bt_keep <- bt_all[n_cells >= 50]
p2 <- ggplot(bt_keep, aes(x = reorder(cell_type, agreement_rate),
                          y = agreement_rate, fill = trait)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 0.5, linetype = 2, color = "red") +
  coord_flip() +
  facet_wrap(~ trait, ncol = 3) +
  labs(title = "Cell-type-level mapping agreement (red dashed: 0.5)",
       x = NULL, y = "agreement rate") +
  theme_bw() +
  theme(legend.position = "none")
ggsave(file.path(FIG_DIR,
                  paste0("sanity_celltype_consistency_by_type", SUFFIX, ".png")),
       p2, width = 11, height = 6, dpi = 150)

# (c) confusion heatmap: true cell_type x predicted cell_type, per trait
# Normalize per row (true) to get fraction
conf_norm <- confusion_all[, .(N = sum(N)),
                            by = .(trait, true_ct, pred_ct)]
conf_norm[, frac := N / sum(N), by = .(trait, true_ct)]

p3 <- ggplot(conf_norm, aes(x = pred_ct, y = true_ct, fill = frac)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "steelblue", limits = c(0, 1)) +
  facet_wrap(~ trait, ncol = 1) +
  labs(title = "Confusion: true vs predicted cell_type (row-normalized)",
       x = "predicted", y = "true", fill = "fraction") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        axis.text.y = element_text(size = 7))
ggsave(file.path(FIG_DIR,
                  paste0("sanity_celltype_consistency_heatmap", SUFFIX, ".png")),
       p3, width = 11, height = 14, dpi = 150)

# ------------------------------------------------------------------
# final report
# ------------------------------------------------------------------
cat("\n[5/5] === Phase 0.3 sanity report ===\n")
print(report_df)

cat("\n--- Gate alpha summary ---\n")
all_pass <- all(report_df$pass_alpha_gate)
for (i in seq_len(nrow(report_df))) {
  r <- report_df[i, ]
  cat(sprintf("  %s : format=%s  distance=%s  consistency=%s (%.3f)  ALPHA = %s\n",
              r$trait, r$pass_format, r$pass_distance,
              r$pass_consistency, r$agreement_overall, r$pass_alpha_gate))
}

cat("\nGATE ALPHA ", if (all_pass) "PASS" else "FAIL", "\n", sep = "")
cat("Outputs in: ", OUT_DIR, "\n", sep = "")
