# module4 robustness: summarize raw ATAC read drop-out results
# -----------------------------------------------------------------------------
# Input:
#   05_robustness/output/top3_raw_atac_dropout/atac_sparsify_eval.csv
#
# This file is written by the patched module3 robustness script. It contains v2
# AUROC rows for each:
#   p_drop x random seed x trait x target x atac_w
#
# This summary adds baseline AUROC from the completed module4 evaluations, then
# writes aggregate tables and figures.
#
# Run in RStudio after 02_run_top3_raw_atac_dropout.R finishes:
#   setwd("D:/RStudio/work/module4_traitaware_pseudopair")
#   source("05_robustness/03_summarize_top3_raw_atac_dropout.R")
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(pROC)
})

if (!exists("ROOT")) ROOT <- "D:/RStudio/work/module4_traitaware_pseudopair"
if (!exists("ROBUST_DIR")) {
  ROBUST_DIR <- file.path(ROOT, "05_robustness/output/top3_raw_atac_dropout")
}
if (!exists("ROBUST_TRAITS")) ROBUST_TRAITS <- c("MONO", "RDW", "RBC")
if (!exists("FIG_DIR")) FIG_DIR <- file.path(ROBUST_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

EVAL_FILE <- file.path(ROBUST_DIR, "atac_sparsify_eval.csv")
SUMMARY_FILE <- file.path(ROBUST_DIR, "top3_raw_atac_dropout_summary.csv")
ADV_FILE <- file.path(ROBUST_DIR, "top3_raw_atac_dropout_fusion_advantage.csv")

if (!exists("BASELINE_FILES")) {
  BASELINE_FILES <- c(
    file.path(
      ROOT,
      "04_evaluation/module4_v2_neighborhood_fdr_ge0p80/module4_v2_neighborhood_fdr_ge0p80_auroc.csv"
    ),
    file.path(
      ROOT,
      "04_evaluation/trait_screening_ge0p80/trait_screening_ge0p80_auroc.csv"
    )
  )
}

stop_missing <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0L) {
    stop("Missing required file(s):\n  ", paste(missing, collapse = "\n  "),
         call. = FALSE)
  }
}

cat(">>> Summarize module4 top3 raw ATAC read drop-out robustness\n\n")
stop_missing(BASELINE_FILES)

META_FILE <- file.path(ROOT, "04_evaluation/inputs/granja_cell_meta.rds")
stop_missing(META_FILE)

TARGETS_ALL <- list(
  MONO = c("11_CD14.Mono.1", "12_CD14.Mono.2"),
  RDW = c("02_Early.Eryth", "03_Late.Eryth", "04_Early.Baso"),
  RBC = c("02_Early.Eryth", "03_Late.Eryth", "04_Early.Baso")
)
missing_targets <- setdiff(ROBUST_TRAITS, names(TARGETS_ALL))
if (length(missing_targets) > 0L) {
  stop("No target definition for trait(s): ", paste(missing_targets, collapse = ", "),
       call. = FALSE)
}
TARGETS <- TARGETS_ALL[ROBUST_TRAITS]

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

rebuild_eval_from_trs <- function() {
  trait_pattern <- paste(ROBUST_TRAITS, collapse = "|")
  trs_files <- list.files(
    ROBUST_DIR,
    pattern = paste0("^TRS_p[0-9.]+_seed[0-9]+_(", trait_pattern, ")_w[0-9.]+\\.tsv\\.gz$"),
    full.names = TRUE
  )
  if (length(trs_files) == 0L) {
    stop("No TRS files found in: ", ROBUST_DIR, call. = FALSE)
  }

  meta <- readRDS(META_FILE)
  if (inherits(meta, "data.table")) meta <- as.data.frame(meta)
  if (is.null(rownames(meta)) ||
      identical(rownames(meta), as.character(seq_len(nrow(meta))))) {
    if ("barcode" %in% colnames(meta)) rownames(meta) <- meta$barcode
  }
  if (!"cell_type" %in% colnames(meta)) {
    stop("META_FILE does not contain a cell_type column: ", META_FILE, call. = FALSE)
  }

  rows <- vector("list", length(trs_files))
  for (i in seq_along(trs_files)) {
    fn <- basename(trs_files[i])
    m <- regexec(
      "^TRS_p([0-9.]+)_seed([0-9]+)_([^_]+)_w([0-9.]+)\\.tsv\\.gz$",
      fn
    )
    parts <- regmatches(fn, m)[[1]]
    if (length(parts) != 5L) {
      stop("Could not parse TRS file name: ", fn, call. = FALSE)
    }

    p_drop <- as.numeric(parts[2])
    seed <- as.integer(parts[3])
    trait <- parts[4]
    atac_w <- as.numeric(parts[5])

    v2 <- fread(trs_files[i])
    if (!all(c("barcode", "TRS") %in% names(v2))) {
      stop("TRS file lacks barcode/TRS columns: ", trs_files[i], call. = FALSE)
    }
    common_eval <- intersect(v2$barcode, rownames(meta))
    lbl <- as.integer(as.character(meta[common_eval, "cell_type"]) %in% TARGETS[[trait]])
    score <- v2$TRS[match(common_eval, v2$barcode)]
    a <- compute_auroc(score, lbl)

    rows[[i]] <- data.table(
      p_drop = p_drop,
      seed = seed,
      trait = trait,
      target = trait,
      atac_w = atac_w,
      auroc = a
    )
  }

  out <- rbindlist(rows)
  out <- unique(out, by = c("p_drop", "seed", "trait", "target", "atac_w"))
  setorder(out, p_drop, seed, trait, atac_w)
  fwrite(out, EVAL_FILE)
  out
}

cat("[1/3] Rebuilding AUROC table from TRS files ...\n")
dt <- rebuild_eval_from_trs()
cat("  rebuilt rows: ", nrow(dt), "\n", sep = "")

dt <- unique(dt, by = c("p_drop", "seed", "trait", "target", "atac_w"))
dt <- dt[trait %in% ROBUST_TRAITS]

base_dt <- rbindlist(lapply(BASELINE_FILES, fread), use.names = TRUE, fill = TRUE)
baseline <- base_dt[
  method == "baseline_atac" & target %in% ROBUST_TRAITS,
  .(
    baseline_auroc = auroc[1],
    baseline_ci_lo = auroc_ci_lo[1],
    baseline_ci_hi = auroc_ci_hi[1],
    n_positive = n_positive[1],
    n_total = n_total[1]
  ),
  by = target
]

summary_dt <- dt[, .(
  mean_auroc = mean(auroc, na.rm = TRUE),
  sd_auroc = stats::sd(auroc, na.rm = TRUE),
  n_runs = .N
), by = .(p_drop, trait, target, atac_w)]
summary_dt[, se_auroc := sd_auroc / sqrt(n_runs)]
summary_dt[, lo95 := mean_auroc - 1.96 * se_auroc]
summary_dt[, hi95 := mean_auroc + 1.96 * se_auroc]
summary_dt <- merge(summary_dt, baseline, by = "target", all.x = TRUE)
summary_dt[, delta_vs_baseline := mean_auroc - baseline_auroc]
setorder(summary_dt, trait, target, p_drop, atac_w)

adv <- dcast(
  summary_dt,
  p_drop + trait + target + baseline_auroc ~ atac_w,
  value.var = c("mean_auroc", "sd_auroc", "n_runs", "delta_vs_baseline")
)
if (all(c("mean_auroc_0.7", "mean_auroc_1") %in% names(adv))) {
  adv[, fusion_advantage_w0p7_minus_w1p0 := mean_auroc_0.7 - mean_auroc_1]
}

fwrite(summary_dt, SUMMARY_FILE)
fwrite(adv, ADV_FILE)

cat("Summary:\n")
print(summary_dt[, .(
  trait,
  target,
  p_drop,
  atac_w,
  mean_auroc = round(mean_auroc, 4),
  delta_vs_baseline = round(delta_vs_baseline, 4),
  n_runs
)])

cat("\nFusion advantage:\n")
cols_show <- intersect(
  c("trait", "target", "p_drop", "baseline_auroc", "mean_auroc_0.7",
    "mean_auroc_1", "fusion_advantage_w0p7_minus_w1p0"),
  names(adv)
)
print(adv[, ..cols_show])

if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)

  summary_dt[, atac_w_label := factor(
    ifelse(atac_w == 0.7, "w=0.7", "w=1.0"),
    levels = c("w=0.7", "w=1.0")
  )]
  summary_dt[, panel := factor(target, levels = ROBUST_TRAITS)]

  p1 <- ggplot(summary_dt, aes(
    x = p_drop,
    y = mean_auroc,
    color = atac_w_label,
    fill = atac_w_label,
    group = atac_w_label
  )) +
    geom_hline(aes(yintercept = baseline_auroc), linetype = 2, color = "grey45") +
    geom_ribbon(aes(ymin = lo95, ymax = hi95), alpha = 0.15, color = NA) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    facet_wrap(~ panel, scales = "free_y") +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    theme_bw(base_size = 13) +
    labs(
      title = "Module4 robustness under raw ATAC read drop-out",
      subtitle = "Patched from module3 robustness. Dashed line = baseline ATAC AUROC.",
      x = "ATAC read drop-out proportion",
      y = "AUROC",
      color = NULL,
      fill = NULL
    )
  ggplot2::ggsave(
    file.path(FIG_DIR, "top3_raw_atac_dropout_auroc.png"),
    p1, width = 10, height = 5, dpi = 300
  )

  p2 <- ggplot(summary_dt, aes(
    x = p_drop,
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
    facet_wrap(~ panel, scales = "free_y") +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    theme_bw(base_size = 13) +
    labs(
      title = "Delta versus baseline under raw ATAC read drop-out",
      x = "ATAC read drop-out proportion",
      y = "Delta AUROC versus baseline",
      color = NULL,
      fill = NULL
    )
  ggplot2::ggsave(
    file.path(FIG_DIR, "top3_raw_atac_dropout_delta_vs_baseline.png"),
    p2, width = 10, height = 5, dpi = 300
  )

  if ("fusion_advantage_w0p7_minus_w1p0" %in% names(adv)) {
      adv[, panel := factor(target, levels = ROBUST_TRAITS)]
    p3 <- ggplot(adv, aes(
      x = p_drop,
      y = fusion_advantage_w0p7_minus_w1p0
    )) +
      geom_hline(yintercept = 0, linetype = 2, color = "grey45") +
      geom_line(color = "#D7263D", linewidth = 0.8) +
      geom_point(color = "#D7263D", size = 2) +
      facet_wrap(~ panel, scales = "free_y") +
      scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
      theme_bw(base_size = 13) +
      labs(
        title = "Fusion advantage under raw ATAC read drop-out",
        subtitle = "Delta = AUROC(w=0.7) - AUROC(w=1.0). Positive means RNA embedding helps the graph.",
        x = "ATAC read drop-out proportion",
        y = "Delta AUROC: w=0.7 minus w=1.0"
      )
    ggplot2::ggsave(
      file.path(FIG_DIR, "top3_raw_atac_dropout_fusion_advantage.png"),
      p3, width = 10, height = 5, dpi = 300
    )
  }

  cat("\nFigures are in: ", FIG_DIR, "\n", sep = "")
} else {
  cat("\nPackage ggplot2 is not installed; skipped figures.\n")
}

cat("\nWrote:\n")
cat("  ", SUMMARY_FILE, "\n", sep = "")
cat("  ", ADV_FILE, "\n", sep = "")
