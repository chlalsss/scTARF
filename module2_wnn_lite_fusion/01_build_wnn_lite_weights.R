# Build per-cell RNA/ATAC weights for WNN-lite graph fusion.
#
# Only unsupervised RNA-reference quality metrics are used for graph weights:
#   - mean_similarity
#   - n_selected
#
# Cell-type agreement columns, when present, are copied only for diagnostics.

source("D:/RStudio/work/module2_wnn_lite_fusion/_bootstrap.R")

required_pkgs(c("data.table", "ggplot2", "jsonlite"))
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(jsonlite)
})

paths <- wnn_lite_paths()
ensure_wnn_lite_dirs(paths)

cutoff <- get_single_cutoff()
label <- cutoff_label(cutoff)
min_weight <- env_numeric("WNN_LITE_RNA_WEIGHT_MIN", 0.20)
max_weight <- env_numeric("WNN_LITE_RNA_WEIGHT_MAX", 0.80)
lower_q <- env_numeric("WNN_LITE_SCALE_LOWER_Q", 0)
upper_q <- env_numeric("WNN_LITE_SCALE_UPPER_Q", 1)

if (min_weight < 0 || max_weight > 1 || min_weight >= max_weight) {
  stop("WNN_LITE_RNA_WEIGHT_MIN/MAX must satisfy 0 <= min < max <= 1", call. = FALSE)
}
if (lower_q < 0 || upper_q > 1 || lower_q >= upper_q) {
  stop("WNN_LITE_SCALE_LOWER_Q/UPPER_Q must satisfy 0 <= lower < upper <= 1", call. = FALSE)
}

by_cell_file <- adaptive_by_cell_path(paths, cutoff)
if (!file.exists(by_cell_file)) {
  stop("Missing adaptive by-cell diagnostics. Run adaptive fusion 02 first: ", by_cell_file)
}

cat("========== Module 2 WNN-lite: build weights ==========\n")
cat("  cutoff       :", cutoff, "\n")
cat("  input        :", by_cell_file, "\n")
cat("  weight range :", sprintf("[%.2f, %.2f]", min_weight, max_weight), "\n")
cat("  scale q      :", sprintf("[%.3f, %.3f]", lower_q, upper_q), "\n")

by_cell <- fread(by_cell_file)
needed <- c("atac_cell", "n_selected", "mean_similarity")
if (!all(needed %in% colnames(by_cell))) {
  stop("By-cell table must contain: ", paste(needed, collapse = ", "), call. = FALSE)
}
if (anyDuplicated(by_cell$atac_cell)) {
  stop("Duplicated atac_cell values in ", by_cell_file, call. = FALSE)
}

weights <- copy(by_cell)
weights[, n_selected_log1p := log1p(n_selected)]
weights[, mean_similarity_score := scale01(mean_similarity, lower_q = lower_q, upper_q = upper_q)]
weights[, n_selected_score := scale01(n_selected_log1p, lower_q = lower_q, upper_q = upper_q)]
weights[, rna_reliability := mean_similarity_score * n_selected_score]
weights[, rna_weight := min_weight + (max_weight - min_weight) * rna_reliability]
weights[, atac_weight := 1 - rna_weight]

front_cols <- c(
  "atac_cell",
  "rna_weight",
  "atac_weight",
  "rna_reliability",
  "mean_similarity_score",
  "n_selected_score",
  "mean_similarity",
  "n_selected",
  "n_selected_log1p"
)
weights <- weights[, c(front_cols, setdiff(colnames(weights), front_cols)), with = FALSE]

out_csv <- wnn_weights_path(paths, cutoff)
fwrite(weights, out_csv)

quantile_row <- function(x, metric) {
  qs <- as.numeric(stats::quantile(x, c(0, 0.10, 0.25, 0.50, 0.75, 0.90, 1), na.rm = TRUE))
  data.table(
    metric = metric,
    min = qs[1],
    q10 = qs[2],
    q25 = qs[3],
    median = qs[4],
    q75 = qs[5],
    q90 = qs[6],
    max = qs[7],
    mean = mean(x, na.rm = TRUE)
  )
}
summary_dt <- rbind(
  quantile_row(weights$rna_weight, "rna_weight"),
  quantile_row(weights$atac_weight, "atac_weight"),
  quantile_row(weights$rna_reliability, "rna_reliability"),
  quantile_row(weights$mean_similarity, "mean_similarity"),
  quantile_row(weights$n_selected, "n_selected")
)
summary_out <- file.path(paths$out_dir, sprintf("wnn_lite_weight_summary_cutoff_%s.csv", label))
fwrite(summary_dt, summary_out)

p_weight <- ggplot(weights, aes(x = rna_weight)) +
  geom_histogram(bins = 80, fill = "#4C78A8", color = "white", linewidth = 0.15) +
  geom_vline(
    xintercept = stats::median(weights$rna_weight, na.rm = TRUE),
    color = "#D62728",
    linewidth = 0.7
  ) +
  labs(
    title = sprintf("WNN-lite RNA weight distribution (cutoff=%s)", label),
    x = "RNA reference weight",
    y = "ATAC cell count"
  ) +
  theme_bw(base_size = 12)
fig_weight <- file.path(paths$fig_dir, sprintf("Fig_wnn_lite_weights_cutoff_%s.png", label))
ggsave(fig_weight, p_weight, width = 6.5, height = 4.6, dpi = 300)

p_components <- ggplot(
  weights,
  aes(x = mean_similarity_score, y = n_selected_score, color = rna_weight)
) +
  geom_point(alpha = 0.35, size = 0.35) +
  scale_color_viridis_c(option = "C", end = 0.95) +
  labs(
    title = "WNN-lite reliability components",
    x = "Scaled mean_similarity",
    y = "Scaled log1p(n_selected)",
    color = "RNA weight"
  ) +
  theme_bw(base_size = 12)
fig_components <- file.path(paths$fig_dir, sprintf("Fig_wnn_lite_weight_components_cutoff_%s.png", label))
ggsave(fig_components, p_components, width = 6.2, height = 5.2, dpi = 300)

manifest <- list(
  script = "01_build_wnn_lite_weights.R",
  cutoff = cutoff,
  input_by_cell_csv = by_cell_file,
  formula = "rna_weight = min + (max - min) * scale01(mean_similarity) * scale01(log1p(n_selected)); atac_weight = 1 - rna_weight",
  cell_type_columns_used_for_weighting = FALSE,
  rna_weight_min = min_weight,
  rna_weight_max = max_weight,
  scale_lower_q = lower_q,
  scale_upper_q = upper_q,
  n_atac = nrow(weights),
  output_weights_csv = out_csv,
  output_summary_csv = summary_out,
  weight_figure = fig_weight,
  component_figure = fig_components,
  created_at = as.character(Sys.time())
)
jsonlite::write_json(
  manifest,
  file.path(paths$out_dir, sprintf("wnn_lite_weight_manifest_cutoff_%s.json", label)),
  pretty = TRUE,
  auto_unbox = TRUE
)

cat("  saved weights :", out_csv, "\n")
cat("  saved summary :", summary_out, "\n")
cat("  saved figure  :", fig_weight, "\n")
cat("\n  weight summary\n")
print(summary_dt)
cat("========== Done ==========\n")

