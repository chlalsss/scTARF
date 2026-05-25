#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(pROC)
})

base_dir <- Sys.getenv(
  "BMMC_RESULT_DIR",
  "D:/SCAVENGE-reproducibility-main/BMMC"
)

out_root <- "D:/RStudio/work/ppt_dataset_figures"
fig_dir <- file.path(out_root, "figures")
tab_dir <- file.path(out_root, "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

meta_file <- file.path(base_dir, "01_data_qc", "BMMC_cell_meta.rds")
baseline_file <- file.path(base_dir, "03_scavenge_baseline", "TRS_matrix_baseline.rds")
dual_dir <- file.path(base_dir, "04_dual_modality")

stopifnot(file.exists(meta_file), file.exists(baseline_file), dir.exists(dual_dir))

meta <- readRDS(meta_file)
if (inherits(meta, "data.table")) meta <- as.data.frame(meta)
if (is.null(rownames(meta)) ||
    identical(rownames(meta), as.character(seq_len(nrow(meta))))) {
  if (!"barcode" %in% colnames(meta)) {
    stop("BMMC metadata must contain a barcode column.")
  }
  rownames(meta) <- meta$barcode
}
if (!"cell_type" %in% colnames(meta)) {
  stop("BMMC metadata must contain a cell_type column.")
}

baseline <- readRDS(baseline_file)
if (inherits(baseline, "data.table")) baseline <- as.data.frame(baseline)
if ("cell" %in% colnames(baseline)) {
  rownames(baseline) <- baseline$cell
  baseline$cell <- NULL
}

targets <- list(
  MONO = list(
    trait = "MONO",
    panel = "MONO\nCD14+/CD16+ monocytes",
    positives = c("CD14+ Mono", "CD16+ Mono")
  ),
  PLT = list(
    trait = "PLT",
    panel = "PLT\nMK/E progenitors",
    positives = c("MK/E prog")
  ),
  LYM_wide = list(
    trait = "LYM",
    panel = "LYM_wide\npan-lymphoid cells",
    positives = c(
      "NK", "ILC", "CD4+ T activated", "CD4+ T naive",
      "CD8+ T", "CD8+ T naive", "Lymph prog",
      "B1 B", "Naive CD20+ B", "Transitional B", "Plasma cell"
    )
  )
)

method_levels <- c(
  "ATAC-only baseline",
  "v2 no RNA (w=1.0)",
  "v2 fused (w=0.7)"
)

method_colors <- c(
  "ATAC-only baseline" = "#2B6CB0",
  "v2 no RNA (w=1.0)" = "#6B7280",
  "v2 fused (w=0.7)" = "#C81D25"
)

load_v2_scores <- function(trait, weight) {
  f <- file.path(dual_dir, sprintf("TRS_v2_%s_w%.1f.tsv.gz", trait, weight))
  if (!file.exists(f)) stop("Missing TRS file: ", f)
  x <- data.table::fread(f)
  if (!all(c("barcode", "TRS") %in% colnames(x))) {
    stop("Unexpected columns in ", f, ": ", paste(colnames(x), collapse = ", "))
  }
  scores <- x$TRS
  names(scores) <- x$barcode
  scores
}

make_roc <- function(scores, positive_celltypes, target_name, panel, method_label) {
  common <- intersect(names(scores), rownames(meta))
  if (length(common) < 100) {
    stop("Too few cells after barcode alignment for ", target_name, " / ", method_label)
  }
  labels <- as.integer(as.character(meta[common, "cell_type"]) %in% positive_celltypes)
  if (sum(labels == 1) < 10 || sum(labels == 0) < 10) {
    stop("Insufficient positive or negative cells for ", target_name)
  }

  roc_obj <- suppressMessages(pROC::roc(labels, as.numeric(scores[common]), quiet = TRUE))
  auc_value <- as.numeric(pROC::auc(roc_obj))
  curve <- data.frame(
    target = target_name,
    panel = panel,
    method = method_label,
    fpr = 1 - roc_obj$specificities,
    tpr = roc_obj$sensitivities,
    auroc = auc_value,
    n_positive = sum(labels == 1),
    n_total = length(labels),
    stringsAsFactors = FALSE
  )
  curve[order(curve$fpr, curve$tpr), ]
}

roc_list <- list()
for (target_name in names(targets)) {
  spec <- targets[[target_name]]
  trait <- spec$trait
  if (!trait %in% colnames(baseline)) {
    stop("Baseline TRS matrix lacks trait column: ", trait)
  }

  base_scores <- baseline[[trait]]
  names(base_scores) <- rownames(baseline)

  roc_list[[paste(target_name, "baseline", sep = "_")]] <- make_roc(
    base_scores, spec$positives, target_name, spec$panel, "ATAC-only baseline"
  )
  roc_list[[paste(target_name, "v2_w1", sep = "_")]] <- make_roc(
    load_v2_scores(trait, 1.0), spec$positives, target_name, spec$panel,
    "v2 no RNA (w=1.0)"
  )
  roc_list[[paste(target_name, "v2_w07", sep = "_")]] <- make_roc(
    load_v2_scores(trait, 0.7), spec$positives, target_name, spec$panel,
    "v2 fused (w=0.7)"
  )
}

roc_df <- do.call(rbind, roc_list)
roc_df$method <- factor(roc_df$method, levels = method_levels)
roc_df$panel <- factor(
  roc_df$panel,
  levels = vapply(targets, `[[`, character(1), "panel")
)

auc_df <- unique(roc_df[c("target", "panel", "method", "auroc", "n_positive", "n_total")])
auc_df <- auc_df[order(auc_df$panel, auc_df$method), ]
short_method <- c(
  "ATAC-only baseline" = "Baseline",
  "v2 no RNA (w=1.0)" = "No RNA",
  "v2 fused (w=0.7)" = "Fused"
)
auc_df$label <- sprintf(
  "%s AUC %.3f",
  unname(short_method[as.character(auc_df$method)]),
  auc_df$auroc
)

label_y <- c(
  "ATAC-only baseline" = 0.27,
  "v2 no RNA (w=1.0)" = 0.19,
  "v2 fused (w=0.7)" = 0.11
)
auc_df$x <- 0.06
auc_df$y <- unname(label_y[as.character(auc_df$method)])

write.csv(
  auc_df[c("target", "method", "auroc", "n_positive", "n_total")],
  file.path(tab_dir, "roc_table45_auc_check.csv"),
  row.names = FALSE
)

write.csv(
  roc_df[c("target", "method", "fpr", "tpr", "auroc")],
  file.path(tab_dir, "roc_table45_curve_points.csv"),
  row.names = FALSE
)

p <- ggplot(roc_df, aes(x = fpr, y = tpr, color = method)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "#B8B8B8", linewidth = 0.45) +
  geom_path(linewidth = 1.1, lineend = "round") +
  geom_label(
    data = auc_df,
    aes(x = x, y = y, label = label, color = method),
    inherit.aes = FALSE,
    hjust = 0,
    size = 3.05,
    fill = "white",
    linewidth = 0,
    label.padding = unit(0.10, "lines"),
    show.legend = FALSE
  ) +
  facet_wrap(~ panel, nrow = 1) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  scale_color_manual(values = method_colors, drop = FALSE) +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  scale_y_continuous(breaks = seq(0, 1, 0.25)) +
  labs(
    title = "ROC curves for cell-type recognition",
    subtitle = "Existing per-cell TRS outputs; dashed diagonal indicates random ranking",
    x = "False positive rate",
    y = "True positive rate",
    color = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 16, margin = margin(b = 4)),
    plot.subtitle = element_text(size = 10, color = "#555555", margin = margin(b = 10)),
    strip.text = element_text(face = "bold", size = 11, lineheight = 0.95),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#E3E3E3", linewidth = 0.35),
    axis.title = element_text(face = "bold"),
    legend.position = "bottom",
    legend.key.width = unit(1.5, "cm"),
    panel.spacing.x = unit(2.0, "lines"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(10, 14, 8, 10)
  )

png_file <- file.path(fig_dir, "fig_roc_table45.png")
pdf_file <- file.path(fig_dir, "fig_roc_table45.pdf")

ggsave(png_file, p, width = 13, height = 4.8, dpi = 320, bg = "white")
ggsave(pdf_file, p, width = 13, height = 4.8, bg = "white")

cat("Saved ROC figure:\n")
cat("  ", png_file, "\n", sep = "")
cat("  ", pdf_file, "\n", sep = "")
cat("Saved AUROC check table:\n")
cat("  ", file.path(tab_dir, "roc_table45_auc_check.csv"), "\n", sep = "")
