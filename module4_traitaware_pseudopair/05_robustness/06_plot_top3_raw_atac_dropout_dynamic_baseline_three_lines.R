# module4 robustness: plot dynamic baseline + w=0.7 + w=1.0
# -----------------------------------------------------------------------------
# This is the corrected robustness plot.
#
# Why:
#   If x-axis is ATAC read drop-out, baseline must also be recomputed at every
#   p_drop level. A fixed horizontal baseline line is only a reference line and
#   should not be used as the main robustness comparison.
#
# Input:
#   dynamic_baseline/baseline_TRS_p<drop>_seed<seed>_<trait>.tsv.gz
#   top3_raw_atac_dropout_summary.csv
#
# Output:
#   dynamic_baseline/dynamic_baseline_eval_full.csv
#   dynamic_baseline/dynamic_baseline_summary_full.csv
#   figures/top3_raw_atac_dropout_dynamic_baseline_three_lines.png
#   figures/top3_raw_atac_dropout_dynamic_baseline_three_lines.pdf
#
# Run in RStudio:
#   setwd("D:/RStudio/work/module4_traitaware_pseudopair")
#   source("05_robustness/06_plot_top3_raw_atac_dropout_dynamic_baseline_three_lines.R")
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(pROC)
  library(ggplot2)
})

ROOT <- "D:/RStudio/work/module4_traitaware_pseudopair"
ROBUST_DIR <- file.path(ROOT, "05_robustness/output/top3_raw_atac_dropout")
BASELINE_DIR <- file.path(ROBUST_DIR, "dynamic_baseline")
FIG_DIR <- file.path(ROBUST_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

META_FILE <- file.path(ROOT, "04_evaluation/inputs/granja_cell_meta.rds")
V2_SUMMARY_FILE <- file.path(ROBUST_DIR, "top3_raw_atac_dropout_summary.csv")

BASELINE_EVAL_FULL <- file.path(BASELINE_DIR, "dynamic_baseline_eval_full.csv")
BASELINE_SUMMARY_FULL <- file.path(BASELINE_DIR, "dynamic_baseline_summary_full.csv")
PLOT_DATA_FILE <- file.path(ROBUST_DIR, "top3_raw_atac_dropout_dynamic_three_line_plot_data.csv")

OUT_PNG <- file.path(FIG_DIR, "top3_raw_atac_dropout_dynamic_baseline_three_lines.png")
OUT_PDF <- file.path(FIG_DIR, "top3_raw_atac_dropout_dynamic_baseline_three_lines.pdf")

TRAITS <- c("MONO", "RDW", "RBC")
TARGETS <- list(
  MONO = c("11_CD14.Mono.1", "12_CD14.Mono.2"),
  RDW = c("02_Early.Eryth", "03_Late.Eryth", "04_Early.Baso"),
  RBC = c("02_Early.Eryth", "03_Late.Eryth", "04_Early.Baso")
)

stop_missing <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0L) {
    stop("Missing required file(s):\n  ", paste(missing, collapse = "\n  "),
         call. = FALSE)
  }
}

load_meta <- function(path) {
  meta <- readRDS(path)
  if (inherits(meta, "data.table")) meta <- as.data.frame(meta)
  if (is.null(rownames(meta)) ||
      identical(rownames(meta), as.character(seq_len(nrow(meta))))) {
    if ("barcode" %in% colnames(meta)) rownames(meta) <- meta$barcode
  }
  meta
}

compute_auroc <- function(predictor, labels) {
  keep <- is.finite(predictor) & !is.na(labels)
  predictor <- predictor[keep]
  labels <- labels[keep]
  if (length(predictor) < 20L) return(NA_real_)
  if (sum(labels == 1L) < 10L || sum(labels == 0L) < 10L) return(NA_real_)
  v <- stats::var(predictor)
  if (!is.finite(v) || v < 1e-12) return(NA_real_)
  roc_obj <- suppressMessages(pROC::roc(
    labels, predictor,
    quiet = TRUE,
    direction = "<",
    levels = c(0L, 1L)
  ))
  as.numeric(pROC::auc(roc_obj))
}

rebuild_dynamic_baseline <- function() {
  trs_files <- list.files(
    BASELINE_DIR,
    pattern = "^baseline_TRS_p[0-9.]+_seed[0-9]+_(MONO|RDW|RBC)\\.tsv\\.gz$",
    full.names = TRUE
  )
  if (length(trs_files) == 0L) {
    stop("No dynamic baseline TRS files found in: ", BASELINE_DIR, call. = FALSE)
  }

  meta <- load_meta(META_FILE)
  rows <- vector("list", length(trs_files))

  for (i in seq_along(trs_files)) {
    fn <- basename(trs_files[i])
    m <- regexec(
      "^baseline_TRS_p([0-9.]+)_seed([0-9]+)_([^_]+)\\.tsv\\.gz$",
      fn
    )
    parts <- regmatches(fn, m)[[1]]
    if (length(parts) != 4L) {
      stop("Could not parse baseline TRS file name: ", fn, call. = FALSE)
    }

    p_drop <- as.numeric(parts[2])
    seed <- as.integer(parts[3])
    trait <- parts[4]

    dt <- fread(trs_files[i])
    if (!all(c("barcode", "TRS") %in% names(dt))) {
      stop("TRS file lacks barcode/TRS columns: ", trs_files[i], call. = FALSE)
    }

    common_eval <- intersect(dt$barcode, rownames(meta))
    lbl <- as.integer(as.character(meta[common_eval, "cell_type"]) %in% TARGETS[[trait]])
    score <- dt$TRS[match(common_eval, dt$barcode)]

    rows[[i]] <- data.table(
      p_drop = p_drop,
      seed = seed,
      trait = trait,
      target = trait,
      method = "dynamic_baseline",
      auroc = compute_auroc(score, lbl),
      n_positive = sum(lbl),
      n_total = length(lbl)
    )
  }

  eval_dt <- rbindlist(rows)
  eval_dt <- unique(eval_dt, by = c("p_drop", "seed", "trait", "target", "method"))
  setorder(eval_dt, trait, p_drop, seed)

  summary_dt <- eval_dt[, .(
    auroc = mean(auroc, na.rm = TRUE),
    sd_auroc = stats::sd(auroc, na.rm = TRUE),
    n_runs = .N,
    n_positive = n_positive[1],
    n_total = n_total[1]
  ), by = .(p_drop, trait, target, method)]
  summary_dt[, se_auroc := sd_auroc / sqrt(n_runs)]
  summary_dt[, lo95 := auroc - 1.96 * se_auroc]
  summary_dt[, hi95 := auroc + 1.96 * se_auroc]
  setorder(summary_dt, trait, p_drop)

  fwrite(eval_dt, BASELINE_EVAL_FULL)
  fwrite(summary_dt, BASELINE_SUMMARY_FULL)

  list(eval = eval_dt, summary = summary_dt)
}

cat(">>> Build corrected dynamic-baseline three-line robustness plot\n\n")
stop_missing(c(META_FILE, V2_SUMMARY_FILE))

cat("[1/3] Rebuilding full dynamic baseline AUROC from TRS files ...\n")
baseline <- rebuild_dynamic_baseline()
baseline_summary <- baseline$summary
cat("  dynamic baseline rows: ", nrow(baseline$eval), " eval; ",
    nrow(baseline_summary), " summary\n", sep = "")

cat("[2/3] Loading w=0.7 and w=1.0 summary ...\n")
v2 <- fread(V2_SUMMARY_FILE)
v2 <- v2[trait %in% TRAITS]
v2_plot <- v2[, .(
  p_drop,
  trait,
  target,
  method = fifelse(abs(atac_w - 0.7) < 1e-8, "w=0.7", "w=1.0"),
  auroc = mean_auroc,
  sd_auroc,
  n_runs,
  lo95,
  hi95
)]

baseline_plot <- baseline_summary[, .(
  p_drop,
  trait,
  target,
  method,
  auroc,
  sd_auroc,
  n_runs,
  lo95,
  hi95
)]

plot_dt <- rbindlist(list(baseline_plot, v2_plot), use.names = TRUE, fill = TRUE)
plot_dt[, method := factor(method, levels = c("dynamic_baseline", "w=0.7", "w=1.0"))]
plot_dt[, method_label := fcase(
  method == "dynamic_baseline", "dynamic baseline",
  method == "w=0.7", "w=0.7",
  method == "w=1.0", "w=1.0"
)]
plot_dt[, method_label := factor(method_label, levels = c("dynamic baseline", "w=0.7", "w=1.0"))]
plot_dt[, panel := factor(target, levels = TRAITS)]
setorder(plot_dt, trait, p_drop, method)
fwrite(plot_dt, PLOT_DATA_FILE)

cat("[3/3] Plotting ...\n")
p <- ggplot(plot_dt, aes(
  x = p_drop,
  y = auroc,
  color = method_label,
  fill = method_label,
  linetype = method_label,
  group = method_label
)) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95), alpha = 0.12, color = NA) +
  geom_line(linewidth = 1.05) +
  geom_point(size = 2.2) +
  facet_wrap(~ panel, nrow = 1, scales = "free_y") +
  scale_x_continuous(
    breaks = seq(0, 0.9, by = 0.1),
    labels = scales::percent_format(accuracy = 1)
  ) +
  scale_color_manual(
    values = c(
      `dynamic baseline` = "#4D4D4D",
      `w=0.7` = "#2C7FB8",
      `w=1.0` = "#D95F02"
    )
  ) +
  scale_fill_manual(
    values = c(
      `dynamic baseline` = "#4D4D4D",
      `w=0.7` = "#2C7FB8",
      `w=1.0` = "#D95F02"
    )
  ) +
  scale_linetype_manual(
    values = c(
      `dynamic baseline` = "dashed",
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
    subtitle = "Dynamic ATAC-only baseline is recomputed at each drop-out level. Ribbon = 95% CI across 3 seeds.",
    x = "ATAC read drop-out proportion",
    y = "AUROC"
  )

ggsave(OUT_PNG, p, width = 11.5, height = 4.8, dpi = 300)
ggsave(OUT_PDF, p, width = 11.5, height = 4.8)

cat("\nDynamic baseline summary:\n")
print(baseline_summary[, .(
  trait,
  p_drop,
  auroc = round(auroc, 4),
  n_runs
)])

cat("\nSaved:\n")
cat("  ", BASELINE_EVAL_FULL, "\n", sep = "")
cat("  ", BASELINE_SUMMARY_FULL, "\n", sep = "")
cat("  ", PLOT_DATA_FILE, "\n", sep = "")
cat("  ", OUT_PNG, "\n", sep = "")
cat("  ", OUT_PDF, "\n", sep = "")
