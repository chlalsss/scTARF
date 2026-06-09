suppressPackageStartupMessages({
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting.", call. = FALSE)
  }
})

library(ggplot2)

ROOT <- normalizePath("D:/RStudio/work/module4_traitaware_pseudopair", mustWork = TRUE)

ORIG_DIR <- file.path(
  ROOT,
  "04_evaluation/baseline_original_scavenge_robustness",
  "module4_original_scavenge_robust_remaining_traits_seed9527_traitfilter_fixed"
)
ORIG_SUMMARY <- file.path(ORIG_DIR, "scavenge_atac_dropout_summary.csv")
ORIG_IMPORT_MANIFEST <- file.path(ORIG_DIR, "trait_import_manifest.csv")

IMPROVED_FILES <- c(
  file.path(
    ROOT,
    "05_robustness/output/v3_sparse_mono_raw_atac_full_pipeline_dropout",
    "v3_sparse_mono_raw_atac_full_pipeline_dropout_plot_data.csv"
  ),
  file.path(
    ROOT,
    "05_robustness/output/v3_sparse_cite_erythroid_raw_atac_full_pipeline_dropout",
    "v3_sparse_cite_erythroid_raw_atac_full_pipeline_dropout_plot_data.csv"
  ),
  file.path(
    ROOT,
    "05_robustness/output/v3_sparse_lym_granja_rna_raw_atac_full_pipeline_dropout",
    "v3_sparse_lym_granja_rna_raw_atac_full_pipeline_dropout_plot_data.csv"
  )
)

OUT_DIR <- file.path(
  ROOT,
  "05_robustness/output/original_scavenge_vs_improved_w_dropout_curves"
)
FIG_DIR <- file.path(OUT_DIR, "figures")
WARNING_FILE <- file.path(OUT_DIR, "WARNING_INVALID_ORIGINAL_SCAVENGE_SOURCE.txt")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(ORIG_SUMMARY)) {
  stop(
    "Missing fixed original SCAVENGE summary:\n  ", ORIG_SUMMARY,
    "\nRun RUN_MODULE4_ORIGINAL_SCAVENGE_BASELINE.R first after the trait-filter fix.",
    call. = FALSE
  )
}
if (file.exists(ORIG_IMPORT_MANIFEST)) {
  import_manifest <- read.csv(ORIG_IMPORT_MANIFEST, check.names = FALSE)
  if (nrow(import_manifest) > 1L &&
      "n_gwas_rows" %in% names(import_manifest) &&
      length(unique(import_manifest[["n_gwas_rows"]])) == 1L) {
    stop(
      "Original SCAVENGE trait import manifest looks invalid: every trait has ",
      "the same n_gwas_rows. Re-run the fixed original baseline in a clean output tag.",
      call. = FALSE
    )
  }
}
# Plot whatever improved result files already exist; skip the ones still running
# (e.g. MONO before its full-pipeline sweep finishes) instead of failing.
missing_improved <- IMPROVED_FILES[!file.exists(IMPROVED_FILES)]
if (length(missing_improved) > 0L) {
  cat("Skipping improved result file(s) not yet produced:\n  ",
      paste(missing_improved, collapse = "\n  "), "\n", sep = "")
}
IMPROVED_FILES <- IMPROVED_FILES[file.exists(IMPROVED_FILES)]
if (length(IMPROVED_FILES) == 0L) {
  stop("No improved result files exist yet; run the improved sweeps first.", call. = FALSE)
}

orig <- read.csv(ORIG_SUMMARY, check.names = FALSE)
improved <- do.call(
  rbind,
  lapply(IMPROVED_FILES, read.csv, check.names = FALSE)
)

orig_plot <- data.frame(
  trait = orig[["trait"]],
  target = orig[["target"]],
  p_drop = as.numeric(orig[["p_drop"]]),
  method = "Original SCAVENGE",
  auroc = as.numeric(orig[["mean_auroc"]]),
  lo95 = as.numeric(orig[["lo95"]]),
  hi95 = as.numeric(orig[["hi95"]]),
  sd_auroc = as.numeric(orig[["sd_auroc"]]),
  n_runs = as.integer(orig[["n_runs"]]),
  n_positive = as.integer(orig[["n_positive"]]),
  n_total = as.integer(orig[["n_eval"]]),
  source_file = ORIG_SUMMARY,
  stringsAsFactors = FALSE
)

improved <- improved[improved[["atac_w"]] %in% c(0.7, 1.0), , drop = FALSE]
improved_plot <- data.frame(
  trait = improved[["trait"]],
  target = improved[["target"]],
  p_drop = as.numeric(improved[["p_drop"]]),
  method = ifelse(
    abs(as.numeric(improved[["atac_w"]]) - 1.0) < 1e-8,
    "Improved w=1",
    "Improved w=0.7"
  ),
  auroc = as.numeric(improved[["auroc"]]),
  lo95 = as.numeric(improved[["lo95"]]),
  hi95 = as.numeric(improved[["hi95"]]),
  sd_auroc = as.numeric(improved[["sd_auroc"]]),
  n_runs = as.integer(improved[["n_runs"]]),
  n_positive = as.integer(improved[["n_positive"]]),
  n_total = as.integer(improved[["n_total"]]),
  source_file = NA_character_,
  stringsAsFactors = FALSE
)

complete_improved <- split(improved_plot[["method"]], improved_plot[["trait"]])
complete_improved <- names(complete_improved)[vapply(
  complete_improved,
  function(x) all(c("Improved w=1", "Improved w=0.7") %in% unique(x)),
  logical(1)
)]

traits <- intersect(unique(orig_plot[["trait"]]), complete_improved)
trait_order <- c("MONO", "MCH", "RDW", "RBC", "LYM")
traits <- c(trait_order[trait_order %in% traits], sort(setdiff(traits, trait_order)))

plot_data <- rbind(
  orig_plot[orig_plot[["trait"]] %in% traits, ],
  improved_plot[improved_plot[["trait"]] %in% traits, ]
)
plot_data[["method"]] <- factor(
  plot_data[["method"]],
  levels = c("Original SCAVENGE", "Improved w=1", "Improved w=0.7")
)
plot_data <- plot_data[order(plot_data[["trait"]], plot_data[["method"]], plot_data[["p_drop"]]), ]

write.csv(
  plot_data,
  file.path(OUT_DIR, "original_scavenge_vs_improved_w_plot_data.csv"),
  row.names = FALSE
)

manifest <- data.frame(
  key = c(
    "original_summary",
    "original_trait_import_manifest",
    paste0("improved_source_", seq_along(IMPROVED_FILES)),
    "output_dir",
    "traits_plotted"
  ),
  value = c(ORIG_SUMMARY, ORIG_IMPORT_MANIFEST, IMPROVED_FILES, OUT_DIR, paste(traits, collapse = ",")),
  stringsAsFactors = FALSE
)
write.csv(manifest, file.path(OUT_DIR, "plot_manifest.csv"), row.names = FALSE)

palette <- c(
  "Original SCAVENGE" = "#333333",
  "Improved w=1" = "#0072B2",
  "Improved w=0.7" = "#D55E00"
)
line_types <- c(
  "Original SCAVENGE" = "solid",
  "Improved w=1" = "longdash",
  "Improved w=0.7" = "dotdash"
)

for (trait_id in traits) {
  dt <- plot_data[plot_data[["trait"]] == trait_id, , drop = FALSE]
  out_base <- file.path(
    FIG_DIR,
    paste0("trait_", trait_id, "_original_scavenge_vs_improved_w")
  )
  y_min <- max(0, floor((min(dt[["auroc"]], na.rm = TRUE) - 0.03) * 20) / 20)
  y_max <- min(1, ceiling((max(dt[["auroc"]], na.rm = TRUE) + 0.03) * 20) / 20)
  if (!is.finite(y_min) || !is.finite(y_max) || y_min >= y_max) {
    y_min <- 0
    y_max <- 1
  }

  p <- ggplot(dt, aes(x = p_drop, y = auroc, color = method, linetype = method)) +
    geom_line(linewidth = 1.05) +
    geom_point(size = 2.2) +
    scale_color_manual(values = palette, drop = FALSE) +
    scale_linetype_manual(values = line_types, drop = FALSE) +
    scale_x_continuous(
      breaks = seq(0, 0.9, by = 0.1),
      limits = c(0, 0.9),
      labels = function(x) paste0(round(x * 100), "%")
    ) +
    scale_y_continuous(breaks = pretty(c(y_min, y_max), n = 6)) +
    coord_cartesian(ylim = c(y_min, y_max)) +
    labs(
      title = paste0(trait_id, " dropout robustness"),
      x = "Raw ATAC dropout",
      y = "AUROC",
      color = NULL,
      linetype = NULL
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      legend.position = "top",
      legend.justification = "left",
      legend.key.width = unit(1.5, "cm"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major.y = element_line(color = "grey88", linewidth = 0.3),
      panel.grid.minor = element_blank()
    )

  ggsave(paste0(out_base, ".png"), p, width = 6.8, height = 4.4, dpi = 320)
  ggsave(paste0(out_base, ".pdf"), p, width = 6.8, height = 4.4, device = cairo_pdf)
}

if (file.exists(WARNING_FILE)) {
  unlink(WARNING_FILE)
}
cat("Wrote plot data: ", file.path(OUT_DIR, "original_scavenge_vs_improved_w_plot_data.csv"), "\n", sep = "")
cat("Wrote figures to: ", FIG_DIR, "\n", sep = "")
cat("Traits plotted: ", paste(traits, collapse = ", "), "\n", sep = "")
