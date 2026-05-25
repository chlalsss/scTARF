# module4 robustness: plot v3 sparse raw ATAC drop-out
# -----------------------------------------------------------------------------
# One plot per trait. Each plot contains three curves:
#   baseline ATAC, w=0.7, w=1.0
#
# Inputs:
#   05_robustness/output/v3_sparse_mono_raw_atac_dropout/atac_sparsify_eval.csv
#   05_robustness/output/v3_sparse_mono_raw_atac_dropout/dynamic_baseline/
#     dynamic_baseline_summary.csv
#
# Usually this is called by:
#   source("05_robustness/07_run_v3_sparse_mono_raw_atac_dropout.R")
#
# Direct debug run in RStudio:
#   setwd("D:/RStudio/work/module4_traitaware_pseudopair")
#   source("05_robustness/08_plot_v3_sparse_mono_raw_atac_dropout.R")
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

ROOT <- "D:/RStudio/work/module4_traitaware_pseudopair"
ROBUST_DIR <- file.path(ROOT, "05_robustness/output/v3_sparse_mono_raw_atac_dropout")
FIG_DIR <- file.path(ROBUST_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

EVAL_FILE <- file.path(ROBUST_DIR, "atac_sparsify_eval.csv")
SUMMARY_FILE <- file.path(ROBUST_DIR, "v3_sparse_mono_raw_atac_dropout_summary.csv")
PLOT_DATA_FILE <- file.path(ROBUST_DIR, "v3_sparse_mono_raw_atac_dropout_plot_data.csv")
BASELINE_SUMMARY_FILE <- file.path(
  ROBUST_DIR,
  "dynamic_baseline",
  "dynamic_baseline_summary.csv"
)

EXPECTED_P_DROP <- round(seq(0.0, 0.9, by = 0.1), 2)

stop_missing <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0L) {
    stop("Missing required file(s):\n  ", paste(missing, collapse = "\n  "),
         call. = FALSE)
  }
}

cat(">>> Plot v3 sparse MONO raw ATAC drop-out robustness\n\n")
stop_missing(c(EVAL_FILE, BASELINE_SUMMARY_FILE))

eval_dt <- fread(EVAL_FILE)
eval_dt[, p_drop := round(as.numeric(p_drop), 2)]
eval_dt[, atac_w := round(as.numeric(atac_w), 1)]
eval_dt <- eval_dt[trait == target]
eval_dt <- eval_dt[p_drop %in% EXPECTED_P_DROP]
eval_dt <- unique(eval_dt, by = c("p_drop", "seed", "trait", "target", "atac_w"))

missing_grid <- eval_dt[, .N, by = .(trait, p_drop, atac_w)]
missing_grid <- missing_grid[!p_drop %in% EXPECTED_P_DROP | !(atac_w %in% c(0.7, 1.0))]
if (nrow(missing_grid) > 0L) {
  stop("Unexpected p_drop / atac_w values in eval file.", call. = FALSE)
}

summary_dt <- eval_dt[, .(
  auroc = mean(auroc, na.rm = TRUE),
  sd_auroc = stats::sd(auroc, na.rm = TRUE),
  n_runs = .N
), by = .(trait, target, p_drop, atac_w)]
summary_dt[, se_auroc := sd_auroc / sqrt(n_runs)]
summary_dt[, lo95 := auroc - 1.96 * se_auroc]
summary_dt[, hi95 := auroc + 1.96 * se_auroc]
summary_dt[, method := fifelse(abs(atac_w - 0.7) < 1e-8, "w=0.7", "w=1.0")]

observed_grid <- unique(summary_dt[, .(trait, method, p_drop)])
expected_grid <- CJ(
  trait = unique(summary_dt$trait),
  method = c("w=0.7", "w=1.0"),
  p_drop = EXPECTED_P_DROP
)
missing_grid <- fsetdiff(expected_grid, observed_grid)
if (nrow(missing_grid) > 0L) {
  stop("Not all p_drop levels 0.0-0.9 are present for every trait/method:\n",
       paste(capture.output(print(missing_grid)), collapse = "\n"),
       call. = FALSE)
}

baseline_summary <- fread(BASELINE_SUMMARY_FILE)
baseline_summary[, p_drop := round(as.numeric(p_drop), 2)]
baseline_summary <- baseline_summary[
  trait %in% unique(summary_dt$trait) &
    target %in% unique(summary_dt$target) &
    p_drop %in% EXPECTED_P_DROP
]

observed_baseline <- unique(baseline_summary[, .(trait, method = "baseline", p_drop)])
expected_baseline <- CJ(
  trait = unique(summary_dt$trait),
  method = "baseline",
  p_drop = EXPECTED_P_DROP
)
missing_baseline <- fsetdiff(expected_baseline, observed_baseline)
if (nrow(missing_baseline) > 0L) {
  stop("Dynamic baseline is missing p_drop levels 0.0-0.9:\n",
       paste(capture.output(print(missing_baseline)), collapse = "\n"),
       call. = FALSE)
}

baseline_plot <- baseline_summary[, .(
  trait,
  target,
  p_drop,
  auroc = baseline_mean_auroc,
  lo95 = baseline_lo95,
  hi95 = baseline_hi95,
  method = "baseline"
)]

method_plot <- summary_dt[, .(
  trait,
  target,
  p_drop,
  auroc,
  lo95,
  hi95,
  method
)]

plot_dt <- rbindlist(list(baseline_plot, method_plot), use.names = TRUE)
plot_dt[, method := factor(method, levels = c("baseline", "w=0.7", "w=1.0"))]
setorder(plot_dt, trait, method, p_drop)

fwrite(summary_dt, SUMMARY_FILE)
fwrite(plot_dt, PLOT_DATA_FILE)

plot_one_trait <- function(trait_name) {
  d <- plot_dt[trait == trait_name]
  out_png <- file.path(FIG_DIR, paste0("v3_sparse_mono_raw_atac_dropout_", trait_name, "_three_lines.png"))
  out_pdf <- file.path(FIG_DIR, paste0("v3_sparse_mono_raw_atac_dropout_", trait_name, "_three_lines.pdf"))

  p <- ggplot(d, aes(
    x = p_drop,
    y = auroc,
    color = method,
    fill = method,
    linetype = method,
    group = method
  )) +
    geom_ribbon(aes(ymin = lo95, ymax = hi95), alpha = 0.12, color = NA) +
    geom_line(linewidth = 1.05) +
    geom_point(size = 2.2) +
    scale_x_continuous(
      breaks = EXPECTED_P_DROP,
      labels = scales::percent_format(accuracy = 1)
    ) +
    scale_color_manual(values = c(
      baseline = "#4D4D4D",
      `w=0.7` = "#2C7FB8",
      `w=1.0` = "#D95F02"
    )) +
    scale_fill_manual(values = c(
      baseline = "#4D4D4D",
      `w=0.7` = "#2C7FB8",
      `w=1.0` = "#D95F02"
    )) +
    scale_linetype_manual(values = c(
      baseline = "dashed",
      `w=0.7` = "solid",
      `w=1.0` = "solid"
    )) +
    theme_bw(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(
      title = paste0(trait_name, " raw ATAC drop-out robustness"),
      subtitle = "Dynamic ATAC-only baseline, w=0.7, and w=1.0. Ribbon = 95% interval across stress seeds.",
      x = "ATAC read drop-out proportion",
      y = "AUROC"
    )

  ggsave(out_png, p, width = 7.5, height = 5.2, dpi = 300)
  ggsave(out_pdf, p, width = 7.5, height = 5.2)
  c(out_png, out_pdf)
}

saved <- unlist(lapply(unique(plot_dt$trait), plot_one_trait), use.names = FALSE)

cat("Saved summary:\n")
cat("  ", SUMMARY_FILE, "\n", sep = "")
cat("  ", PLOT_DATA_FILE, "\n", sep = "")
cat("Saved figures:\n")
cat(paste0("  ", saved, collapse = "\n"), "\n", sep = "")
