# Plot pure-ATAC original SCAVENGE drop-out robustness curves.
#
# Run after 01_run_scavenge_atac_dropout.R:
#   setwd("D:/RStudio/work/scavenge_atac_baseline_robustness")
#   source("02_plot_robustness.R")

this_file <- tryCatch(normalizePath(sys.frame(1)$ofile, winslash = "/"),
                      error = function(e) "")
candidate_dirs <- character()
if (nzchar(this_file)) candidate_dirs <- c(candidate_dirs, dirname(this_file))
if (exists("ROBUST_SCRIPT_DIR")) candidate_dirs <- c(candidate_dirs, ROBUST_SCRIPT_DIR)
candidate_dirs <- c(
  candidate_dirs,
  "D:/RStudio/work/scavenge_atac_baseline_robustness",
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
)
candidate_dirs <- unique(normalizePath(candidate_dirs, winslash = "/", mustWork = FALSE))
config_hits <- candidate_dirs[file.exists(file.path(candidate_dirs, "00_config.R"))]
if (length(config_hits) == 0L) {
  stop("Cannot locate 00_config.R. Candidate dirs:\n  ",
       paste(candidate_dirs, collapse = "\n  "), call. = FALSE)
}
ROBUST_SCRIPT_DIR <- config_hits[1L]
source(file.path(ROBUST_SCRIPT_DIR, "00_config.R"), local = environment())

required_pkgs <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) {
    stop("Missing R package(s): ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
}
required_pkgs(c("data.table", "ggplot2"))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

if (!file.exists(EVAL_FILE)) {
  stop("Missing eval file: ", EVAL_FILE,
       "\nRun 01_run_scavenge_atac_dropout.R first.", call. = FALSE)
}

cat(">>> Plot pure-ATAC original SCAVENGE drop-out robustness\n\n")
cat("Input : ", EVAL_FILE, "\n", sep = "")
cat("Output: ", FIG_DIR, "\n\n", sep = "")

eval_dt <- fread(EVAL_FILE)
eval_dt <- eval_dt[
  method == "original_scavenge_atac" &
    trait %in% TRAITS &
    p_drop %in% P_DROP_LIST
]
eval_dt <- unique(eval_dt, by = c("method", "p_drop", "seed", "trait", "target"))

if (!"auroc" %in% names(eval_dt) || all(!is.finite(eval_dt$auroc))) {
  stop("No finite AUROC values found. Check META_FILE, META_CELLTYPE_COL, and TARGETS.",
       call. = FALSE)
}

expected_grid <- CJ(
  trait = TRAITS,
  p_drop = P_DROP_LIST,
  seed = SEED_LIST
)
observed_grid <- unique(eval_dt[, .(trait, p_drop, seed)])
missing_grid <- fsetdiff(expected_grid, observed_grid)
if (nrow(missing_grid) > 0L) {
  cat("Warning: missing p_drop x seed rows; plotting available data only.\n")
  print(missing_grid)
  cat("\n")
}

summary_dt <- eval_dt[is.finite(auroc), .(
  mean_auroc = mean(auroc, na.rm = TRUE),
  sd_auroc = stats::sd(auroc, na.rm = TRUE),
  n_runs = .N,
  n_positive = n_positive[1],
  n_eval = n_eval[1]
), by = .(trait, target, p_drop)]
summary_dt[, se_auroc := sd_auroc / sqrt(n_runs)]
summary_dt[, lo95 := mean_auroc - 1.96 * se_auroc]
summary_dt[, hi95 := mean_auroc + 1.96 * se_auroc]
summary_dt[, trait := factor(trait, levels = TRAITS)]
setorder(summary_dt, trait, p_drop)

fwrite(summary_dt, SUMMARY_FILE)
plot_data_file <- file.path(OUTPUT_DIR, "scavenge_atac_dropout_plot_data.csv")
fwrite(summary_dt, plot_data_file)

p <- ggplot(summary_dt, aes(
  x = p_drop,
  y = mean_auroc,
  group = trait
)) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95), fill = "#6B7280", alpha = 0.16) +
  geom_line(color = "#1F77B4", linewidth = 0.9) +
  geom_point(color = "#1F77B4", size = 2.1) +
  facet_wrap(~ trait, nrow = 1, scales = "free_y") +
  scale_x_continuous(
    breaks = P_DROP_LIST,
    labels = function(x) paste0(round(100 * x), "%")
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "grey92", color = "grey65"),
    strip.text = element_text(face = "bold")
  ) +
  labs(
    title = "Original SCAVENGE pure-ATAC robustness",
    subtitle = "Raw ATAC read drop-out; line is mean AUROC and ribbon is 95% interval across configured seeds.",
    x = "ATAC read drop-out proportion",
    y = "AUROC"
  )

out_png <- file.path(FIG_DIR, "scavenge_atac_dropout_auroc.png")
out_pdf <- file.path(FIG_DIR, "scavenge_atac_dropout_auroc.pdf")
ggsave(out_png, p, width = max(6.5, 3.5 * length(unique(summary_dt$trait))),
       height = 4.8, dpi = 300)
ggsave(out_pdf, p, width = max(6.5, 3.5 * length(unique(summary_dt$trait))),
       height = 4.8)

cat("Summary:\n")
print(summary_dt[, .(
  trait,
  p_drop,
  mean_auroc = round(mean_auroc, 4),
  lo95 = round(lo95, 4),
  hi95 = round(hi95, 4),
  n_runs
)])

cat("\nWrote:\n")
cat("  ", SUMMARY_FILE, "\n", sep = "")
cat("  ", plot_data_file, "\n", sep = "")
cat("  ", out_png, "\n", sep = "")
cat("  ", out_pdf, "\n", sep = "")
