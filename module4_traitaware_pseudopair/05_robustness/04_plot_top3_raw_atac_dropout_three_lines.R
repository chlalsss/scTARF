# module4 robustness: clear three-line AUROC plots
# -----------------------------------------------------------------------------
# Plot baseline, w=0.7, and w=1.0 in the same figure.
#
# Input:
#   05_robustness/output/top3_raw_atac_dropout/top3_raw_atac_dropout_summary.csv
#
# Output:
#   05_robustness/output/top3_raw_atac_dropout/figures/top3_raw_atac_dropout_three_lines.png
#   05_robustness/output/top3_raw_atac_dropout/figures/top3_raw_atac_dropout_three_lines.pdf
#
# Run in RStudio:
#   setwd("D:/RStudio/work/module4_traitaware_pseudopair")
#   source("05_robustness/04_plot_top3_raw_atac_dropout_three_lines.R")
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

ROOT <- "D:/RStudio/work/module4_traitaware_pseudopair"
ROBUST_DIR <- file.path(ROOT, "05_robustness/output/top3_raw_atac_dropout")
FIG_DIR <- file.path(ROBUST_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

SUMMARY_FILE <- file.path(ROBUST_DIR, "top3_raw_atac_dropout_summary.csv")
OUT_PNG <- file.path(FIG_DIR, "top3_raw_atac_dropout_three_lines.png")
OUT_PDF <- file.path(FIG_DIR, "top3_raw_atac_dropout_three_lines.pdf")

if (!file.exists(SUMMARY_FILE)) {
  stop("Missing summary file: ", SUMMARY_FILE, call. = FALSE)
}

summary_dt <- fread(SUMMARY_FILE)
summary_dt <- summary_dt[trait %in% c("MONO", "RDW", "RBC")]

v2_dt <- summary_dt[, .(
  trait,
  target,
  p_drop,
  auroc = mean_auroc,
  lo95,
  hi95,
  method = fifelse(abs(atac_w - 0.7) < 1e-8, "w=0.7", "w=1.0")
)]

baseline_dt <- unique(summary_dt[, .(trait, target, p_drop, baseline_auroc)])
baseline_dt <- baseline_dt[, .(
  trait,
  target,
  p_drop,
  auroc = baseline_auroc,
  lo95 = baseline_auroc,
  hi95 = baseline_auroc,
  method = "baseline"
)]

plot_dt <- rbindlist(list(baseline_dt, v2_dt), use.names = TRUE)
plot_dt[, method := factor(method, levels = c("baseline", "w=0.7", "w=1.0"))]
plot_dt[, panel := factor(target, levels = c("MONO", "RDW", "RBC"))]

p <- ggplot(plot_dt, aes(
  x = p_drop,
  y = auroc,
  color = method,
  linetype = method,
  group = method
)) +
  geom_ribbon(
    data = plot_dt[method != "baseline"],
    aes(ymin = lo95, ymax = hi95, fill = method),
    alpha = 0.14,
    color = NA,
    show.legend = FALSE
  ) +
  geom_line(linewidth = 1.05) +
  geom_point(size = 2.2) +
  facet_wrap(~ panel, nrow = 1, scales = "free_y") +
  scale_x_continuous(
    breaks = seq(0, 0.9, by = 0.1),
    labels = scales::percent_format(accuracy = 1)
  ) +
  scale_color_manual(
    values = c(
      baseline = "#4D4D4D",
      `w=0.7` = "#2C7FB8",
      `w=1.0` = "#D95F02"
    )
  ) +
  scale_fill_manual(
    values = c(
      `w=0.7` = "#2C7FB8",
      `w=1.0` = "#D95F02"
    )
  ) +
  scale_linetype_manual(
    values = c(
      baseline = "dashed",
      `w=0.7` = "solid",
      `w=1.0` = "solid"
    )
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "top",
    legend.title = element_blank(),
    strip.background = element_rect(fill = "grey92", color = "grey60"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    title = "Raw ATAC read drop-out robustness",
    subtitle = "Baseline, w=0.7, and w=1.0 are shown together. Ribbon = 95% CI across 3 stress seeds.",
    x = "ATAC read drop-out proportion",
    y = "AUROC"
  )

ggsave(OUT_PNG, p, width = 11.5, height = 4.8, dpi = 300)
ggsave(OUT_PDF, p, width = 11.5, height = 4.8)

cat("Saved:\n")
cat("  ", OUT_PNG, "\n", sep = "")
cat("  ", OUT_PDF, "\n", sep = "")
