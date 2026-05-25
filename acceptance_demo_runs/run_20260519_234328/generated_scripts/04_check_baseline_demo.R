# 03_scavenge_baseline/04_check_baseline.R
#
# Sanity-check: does the ATAC-only baseline TRS actually concentrate on the
# expected target cell types per trait? If MONO TRS top5 isn't enriched for
# CD14+/CD16+ Mono, something is wrong upstream and there's no point comparing
# new methods against this baseline.
#
# Inputs:
#   <out_dir>/03_scavenge_baseline/TRS_matrix_baseline.rds
#   <out_dir>/01_data_qc/BMMC_cell_meta.rds
#   bmmc_target_types() from _config.R
#
# Outputs:
#   <out_dir>/03_scavenge_baseline/baseline_enrichment.csv
#   <out_dir>/03_scavenge_baseline/figures/baseline_top5_celltype.png

source("D:/RStudio/work/module3_multiome/_bootstrap.R")
required_pkgs(c("data.table", "ggplot2"))

paths <- m3_paths()
out_dir <- file.path(paths$out_dir, "03_scavenge_baseline")
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

trs_mat <- readRDS(file.path(out_dir, "TRS_matrix_baseline.rds"))
meta    <- readRDS(file.path(paths$out_dir, "01_data_qc", "BMMC_cell_meta.rds"))
targets <- bmmc_target_types()

trs_dt <- data.table::as.data.table(trs_mat)
meta_dt <- data.table::as.data.table(meta)
merged <- merge(trs_dt, meta_dt[, .(cell = barcode, cell_type)], by = "cell")
message("Merged cells with TRS + cell_type: ", nrow(merged))

# AUROC: rank-based, immune to target-frequency imbalance.
# label = is_target_celltype, predictor = TRS. Computed via the
# Mann-Whitney U identity: AUROC = (sum(rank(TRS in pos)) - n_pos*(n_pos+1)/2)
#                                  / (n_pos * n_neg)
auroc_from_ranks <- function(scores, labels) {
  ok <- !is.na(scores) & !is.na(labels)
  scores <- scores[ok]; labels <- as.logical(labels[ok])
  n_pos <- sum(labels); n_neg <- sum(!labels)
  if (n_pos == 0 || n_neg == 0) return(NA_real_)
  r <- rank(scores)
  (sum(r[labels]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

rows <- list()
for (tr in names(targets)) {
  if (!tr %in% colnames(merged)) {
    message("Trait ", tr, " missing from TRS matrix; skip"); next
  }
  ct_target <- targets[[tr]]

  vals    <- merged[[tr]]
  is_tgt  <- merged$cell_type %in% ct_target
  ord     <- order(vals, decreasing = TRUE, na.last = NA)
  n_total <- length(ord)
  n_top5  <- ceiling(0.05 * n_total)
  top5_idx <- ord[seq_len(n_top5)]

  top5_target <- sum(merged$cell_type[top5_idx] %in% ct_target)
  bg_target   <- sum(is_tgt)
  bg_total    <- nrow(merged)

  observed_frac <- top5_target / n_top5
  expected_frac <- bg_target / bg_total
  enrichment    <- observed_frac / expected_frac
  auroc         <- auroc_from_ranks(vals, is_tgt)

  rows[[tr]] <- data.table::data.table(
    trait             = tr,
    n_top5            = n_top5,
    n_top5_in_target  = top5_target,
    target_frac_top5  = observed_frac,
    target_frac_bg    = expected_frac,
    enrichment_ratio  = enrichment,
    AUROC             = auroc
  )
  message(sprintf(
    "  %s: top5 in target = %d/%d (%.1f%%) | bg target = %.1f%% | enr = %.2fx | AUROC = %.3f",
    tr, top5_target, n_top5, observed_frac * 100,
    expected_frac * 100, enrichment, auroc
  ))
}

summary_dt <- data.table::rbindlist(rows)
data.table::fwrite(summary_dt, file.path(out_dir, "baseline_enrichment.csv"))
print(summary_dt)

# Plot: stacked celltype distribution within top5 per trait
plot_rows <- list()
for (tr in names(targets)) {
  if (!tr %in% colnames(merged)) next
  vals <- merged[[tr]]
  ord  <- order(vals, decreasing = TRUE, na.last = NA)
  top5_idx <- ord[seq_len(ceiling(0.05 * length(ord)))]
  ct_counts <- merged[top5_idx, .N, by = cell_type]
  ct_counts[, trait := tr]
  ct_counts[, is_target := cell_type %in% targets[[tr]]]
  plot_rows[[tr]] <- ct_counts
}
plot_dt <- data.table::rbindlist(plot_rows)
p <- ggplot2::ggplot(
  plot_dt,
  ggplot2::aes(x = trait, y = N, fill = is_target)
) +
  ggplot2::geom_col(position = "stack") +
  ggplot2::scale_fill_manual(
    values = c(`TRUE` = "#2c7fb8", `FALSE` = "#cccccc"),
    name = "in target set"
  ) +
  ggplot2::labs(
    title = "Baseline TRS top 5% — cell type composition",
    x = "Trait", y = "Cells in TRS top 5%"
  ) +
  ggplot2::theme_bw()
ggplot2::ggsave(
  file.path(fig_dir, "baseline_top5_celltype.png"),
  p, width = 5, height = 4, dpi = 150
)
message("\nWrote: ",
        file.path(out_dir, "baseline_enrichment.csv"), "\n",
        "       ",
        file.path(fig_dir, "baseline_top5_celltype.png"))

# --- Ablation: alternative LYM target (B-lineage only) -------------------
# Diagnostic showed LYM TRS pulls B-lineage to mid-top while T/NK stays at
# bottom — i.e. SCAVENGE on this peak set treats LYM as a B-cell-related
# signal rather than pan-lymphocyte. Test that hypothesis directly.

if ("LYM" %in% colnames(merged)) {
  b_lineage <- c("Naive CD20+ B", "Transitional B", "B1 B",
                 "Plasma cell", "Lymph prog")
  is_b   <- merged$cell_type %in% b_lineage
  auroc_b <- auroc_from_ranks(merged$LYM, is_b)
  message(sprintf(
    "  LYM (B-lineage target only, %d cells = %.1f%% of bg): AUROC = %.3f",
    sum(is_b), 100 * sum(is_b) / nrow(merged), auroc_b
  ))
  # also: T/NK only
  t_nk <- c("CD4+ T activated", "CD4+ T naive", "CD8+ T",
            "CD8+ T naive", "NK", "ILC")
  is_tnk   <- merged$cell_type %in% t_nk
  auroc_tnk <- auroc_from_ranks(merged$LYM, is_tnk)
  message(sprintf(
    "  LYM (T/NK target only, %d cells = %.1f%% of bg): AUROC = %.3f",
    sum(is_tnk), 100 * sum(is_tnk) / nrow(merged), auroc_tnk
  ))
}

# Verdict — use AUROC as primary, enrichment_ratio as secondary
# AUROC interpretation: 0.5 random, >0.7 OK, >0.85 strong. Robust to imbalanced targets.
ok_traits   <- summary_dt[AUROC >= 0.70, trait]
weak_traits <- summary_dt[AUROC <  0.70, trait]
if (length(weak_traits) == 0) {
  message("\nBASELINE OK — all traits AUROC >= 0.70.")
} else {
  message("\nBASELINE WEAK on (AUROC < 0.70): ", paste(weak_traits, collapse = ", "),
          ".")
  message("This is exactly the gap dual-modality methods aim to close. ",
          "Weak ATAC-only baseline does not block stage 2 — it widens the headroom.")
}
