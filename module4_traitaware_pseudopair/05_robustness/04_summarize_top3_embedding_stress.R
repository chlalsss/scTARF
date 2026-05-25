# module4 robustness: summarize and plot top3 embedding stress results
# -----------------------------------------------------------------------------
# Input:
#   05_robustness/output/top3_embedding_stress/top3_embedding_stress_auroc_all.csv
#
# Outputs:
#   top3_embedding_stress_summary.csv
#   top3_embedding_stress_fusion_advantage.csv
#   figures/top3_embedding_stress_auroc.png
#   figures/top3_embedding_stress_delta_vs_baseline.png
#   figures/top3_embedding_stress_fusion_advantage.png
#
# Run in RStudio:
#   setwd("D:/RStudio/work/module4_traitaware_pseudopair")
#   source("05_robustness/04_summarize_top3_embedding_stress.R")
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

ROOT <- "D:/RStudio/work/module4_traitaware_pseudopair"
ROBUST_DIR <- file.path(ROOT, "05_robustness/output/top3_embedding_stress")
FIG_DIR <- file.path(ROBUST_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

IN_FILE <- file.path(ROBUST_DIR, "top3_embedding_stress_auroc_all.csv")
SUMMARY_FILE <- file.path(ROBUST_DIR, "top3_embedding_stress_summary.csv")
ADV_FILE <- file.path(ROBUST_DIR, "top3_embedding_stress_fusion_advantage.csv")

cat(">>> Summarize top3 ATAC embedding stress robustness\n\n")
if (!file.exists(IN_FILE)) {
  stop("Missing input file: ", IN_FILE, call. = FALSE)
}

dt <- fread(IN_FILE)
dt <- unique(dt, by = c("run_id", "noise_sd", "stress_seed", "method", "target", "atac_w"))

baseline <- dt[method == "baseline_atac", .(
  baseline_auroc = mean(auroc, na.rm = TRUE),
  baseline_n_positive = n_positive[1],
  baseline_n_total = n_total[1]
), by = target]

v2 <- dt[method == "v2"]
summary_dt <- v2[, .(
  mean_auroc = mean(auroc, na.rm = TRUE),
  sd_auroc = stats::sd(auroc, na.rm = TRUE),
  n_runs = .N
), by = .(noise_sd, target, atac_w)]
summary_dt[, se_auroc := sd_auroc / sqrt(n_runs)]
summary_dt[, lo95 := mean_auroc - 1.96 * se_auroc]
summary_dt[, hi95 := mean_auroc + 1.96 * se_auroc]
summary_dt <- merge(summary_dt, baseline, by = "target", all.x = TRUE)
summary_dt[, delta_vs_baseline := mean_auroc - baseline_auroc]
setorder(summary_dt, target, noise_sd, atac_w)

wide <- dcast(
  summary_dt,
  noise_sd + target + baseline_auroc ~ atac_w,
  value.var = c("mean_auroc", "sd_auroc", "n_runs", "delta_vs_baseline")
)
if (all(c("mean_auroc_0.7", "mean_auroc_1") %in% names(wide))) {
  wide[, fusion_advantage_w0p7_minus_w1p0 := mean_auroc_0.7 - mean_auroc_1]
}

fwrite(summary_dt, SUMMARY_FILE)
fwrite(wide, ADV_FILE)

cat("Summary:\n")
print(summary_dt[, .(
  target,
  noise_sd,
  atac_w,
  mean_auroc = round(mean_auroc, 4),
  delta_vs_baseline = round(delta_vs_baseline, 4),
  n_runs
)])

cat("\nFusion advantage table:\n")
cols_show <- intersect(
  c("target", "noise_sd", "baseline_auroc", "mean_auroc_0.7", "mean_auroc_1",
    "fusion_advantage_w0p7_minus_w1p0"),
  names(wide)
)
print(wide[, ..cols_show])

if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)

  summary_dt[, atac_w_label := factor(
    ifelse(atac_w == 0.7, "w=0.7", "w=1.0"),
    levels = c("w=0.7", "w=1.0")
  )]

  p1 <- ggplot(summary_dt, aes(
    x = noise_sd,
    y = mean_auroc,
    color = atac_w_label,
    fill = atac_w_label,
    group = atac_w_label
  )) +
    geom_ribbon(aes(ymin = lo95, ymax = hi95), alpha = 0.15, color = NA) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    geom_hline(aes(yintercept = baseline_auroc), linetype = 2, color = "grey45") +
    facet_wrap(~ target, scales = "free_y") +
    theme_bw(base_size = 13) +
    labs(
      title = "Top3 trait robustness under ATAC embedding stress",
      subtitle = "Dashed line = baseline ATAC AUROC. Ribbon = 95% CI across stress seeds.",
      x = "ATAC embedding noise SD",
      y = "AUROC",
      color = NULL,
      fill = NULL
    )
  ggplot2::ggsave(file.path(FIG_DIR, "top3_embedding_stress_auroc.png"),
                  p1, width = 10, height = 5, dpi = 300)

  p2 <- ggplot(summary_dt, aes(
    x = noise_sd,
    y = delta_vs_baseline,
    color = atac_w_label,
    fill = atac_w_label,
    group = atac_w_label
  )) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey45") +
    geom_ribbon(aes(
      ymin = delta_vs_baseline - 1.96 * se_auroc,
      ymax = delta_vs_baseline + 1.96 * se_auroc
    ), alpha = 0.15, color = NA) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    facet_wrap(~ target, scales = "free_y") +
    theme_bw(base_size = 13) +
    labs(
      title = "Delta versus baseline under ATAC embedding stress",
      x = "ATAC embedding noise SD",
      y = "Delta AUROC versus baseline",
      color = NULL,
      fill = NULL
    )
  ggplot2::ggsave(file.path(FIG_DIR, "top3_embedding_stress_delta_vs_baseline.png"),
                  p2, width = 10, height = 5, dpi = 300)

  if ("fusion_advantage_w0p7_minus_w1p0" %in% names(wide)) {
    p3 <- ggplot(wide, aes(
      x = noise_sd,
      y = fusion_advantage_w0p7_minus_w1p0
    )) +
      geom_hline(yintercept = 0, linetype = 2, color = "grey45") +
      geom_line(color = "#D7263D", linewidth = 0.8) +
      geom_point(color = "#D7263D", size = 2) +
      facet_wrap(~ target, scales = "free_y") +
      theme_bw(base_size = 13) +
      labs(
        title = "Fusion advantage under ATAC embedding stress",
        subtitle = "Delta = AUROC(w=0.7) - AUROC(w=1.0). Positive means RNA embedding helps the graph.",
        x = "ATAC embedding noise SD",
        y = "Delta AUROC: w=0.7 minus w=1.0"
      )
    ggplot2::ggsave(file.path(FIG_DIR, "top3_embedding_stress_fusion_advantage.png"),
                    p3, width = 10, height = 5, dpi = 300)
  }

  cat("\nFigures are in: ", FIG_DIR, "\n", sep = "")
} else {
  cat("\nPackage ggplot2 is not installed; skipped figures.\n")
}

cat("\nWrote:\n")
cat("  ", SUMMARY_FILE, "\n", sep = "")
cat("  ", ADV_FILE, "\n", sep = "")
