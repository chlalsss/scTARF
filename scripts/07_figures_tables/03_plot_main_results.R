# Plot main AUROC comparison figures from the unified table.
# Run after scripts/07_figures_tables/02_build_master_result_tables.R.

current_script_dir <- function() {
  frames <- sys.frames()
  for (i in rev(seq_along(frames))) {
    if (!is.null(frames[[i]]$ofile)) {
      return(dirname(normalizePath(frames[[i]]$ofile, mustWork = FALSE)))
    }
  }
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE)))
  }
  normalizePath("scripts/07_figures_tables", mustWork = FALSE)
}

source(file.path(current_script_dir(), "00_final_delivery_common.R"))
ensure_final_dirs()
need_packages(c("ggplot2"))
library(ggplot2)

in_file <- file.path(FINAL_ROOT, "01_master_tables/main_result_comparison_table.csv")
if (!file.exists(in_file)) {
  stop("Missing main table. Run scripts/07_figures_tables/02_build_master_result_tables.R first.", call. = FALSE)
}

out_dir <- file.path(FINAL_ROOT, "02_main_figures")
fig4_dir <- file.path(out_dir, "fig4_main_assets")
dir.create(fig4_dir, recursive = TRUE, showWarnings = FALSE)

x <- read_csv_if_exists(in_file)
x <- x[x$metric == "AUROC" & is.finite(x$value), , drop = FALSE]
if (!"final_version" %in% names(x)) {
  x$final_version <- classify_final_version(
    x$method_version, x$method, x$feature_mode, x$method_group, x$pairing_setting
  )
}

keep_versions <- c("Baseline", "scDRS-only", "Paired-Fusion", "v3.3")
new_label <- c(
  "Baseline" = "Baseline (ATAC-only)",
  "scDRS-only" = "scDRS-only (RNA)",
  "Paired-Fusion" = "scTARF1.0-paired",
  "v3.3" = "scTARF3.3-sparse"
)
label_levels <- unname(new_label[keep_versions])
version_colors <- c(
  "Baseline (ATAC-only)" = "#9AA0A6",
  "scDRS-only (RNA)" = "#E8A33D",
  "scTARF1.0-paired" = "#8E63CE",
  "scTARF3.3-sparse" = "#C44E4E"
)

paired_label <- "Paired BMMC multiome"
nonpaired_label <- "Non-paired Granja ATAC + RNA reference"
cols <- c("pairing_setting", "dataset", "trait", "final_version", "value",
          "method_group", "source_file", "note")
select_cols <- function(d) d[, cols, drop = FALSE]

paired <- x[x$pairing_setting == "paired" &
              x$trait %in% c("MONO", "PLT", "LYM_wide") &
              x$final_version %in% keep_versions, , drop = FALSE]
paired <- select_cols(paired)
paired$trait_disp <- ifelse(paired$trait == "LYM_wide", "LYM", paired$trait)
paired$setting_lab <- paired_label

nonpaired_fusion <- x[x$pairing_setting == "non_paired" &
                       x$trait %in% c("MONO", "RDW", "LYM") &
                       x$final_version %in% keep_versions, , drop = FALSE]
nonpaired_fusion <- select_cols(nonpaired_fusion)

nonpaired_scdrs <- x[x$pairing_setting == "rna_reference_only" &
                       x$final_version == "scDRS-only" &
                       x$trait %in% c("MONO", "RDW", "LYM"), , drop = FALSE]
nonpaired_scdrs <- select_cols(nonpaired_scdrs)

nonpaired <- rbind(nonpaired_fusion, nonpaired_scdrs)
nonpaired$trait_disp <- nonpaired$trait
nonpaired$setting_lab <- nonpaired_label

plot_df <- rbind(paired, nonpaired)
plot_df <- best_by_key(plot_df, keys = c("setting_lab", "trait_disp", "final_version"))
plot_df$display_method <- factor(unname(new_label[plot_df$final_version]), levels = label_levels)
plot_df$trait_disp <- factor(plot_df$trait_disp, levels = c("MONO", "PLT", "LYM", "RDW"))
plot_df$setting_lab <- factor(plot_df$setting_lab, levels = c(paired_label, nonpaired_label))

p1 <- ggplot(plot_df, aes(trait_disp, value, fill = display_method)) +
  geom_col(position = position_dodge2(width = 0.85, preserve = "single"),
           width = 0.78, color = "white", linewidth = 0.25) +
  geom_hline(yintercept = 0.5, linetype = "dashed",
             color = "#B0B0B0", linewidth = 0.4) +
  facet_wrap(~ setting_lab, scales = "free_x") +
  scale_fill_manual(values = version_colors, name = NULL, drop = FALSE) +
  scale_y_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.2),
                     expand = c(0, 0)) +
  labs(x = NULL, y = "AUROC", title = "Main AUROC comparison by trait") +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "bottom"
  ) +
  guides(fill = guide_legend(nrow = 1))
ggsave(file.path(out_dir, "fig_main_auroc_by_trait.png"), p1,
       width = 11, height = 5.4, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "fig_main_auroc_by_trait.pdf"), p1,
       width = 11, height = 5.4)

base_df <- plot_df[plot_df$final_version == "Baseline",
                   c("setting_lab", "trait_disp", "value")]
names(base_df)[3] <- "baseline_auroc"
delta_df <- merge(plot_df[plot_df$final_version != "Baseline", ],
                  base_df, by = c("setting_lab", "trait_disp"), all.x = TRUE)
delta_df$delta <- delta_df$value - delta_df$baseline_auroc
delta_df <- delta_df[is.finite(delta_df$delta), , drop = FALSE]

p2 <- ggplot(delta_df, aes(trait_disp, delta, fill = display_method)) +
  geom_hline(yintercept = 0, color = "grey40", linewidth = 0.4) +
  geom_col(position = position_dodge2(width = 0.85, preserve = "single"),
           width = 0.78, color = "white", linewidth = 0.25) +
  facet_wrap(~ setting_lab, scales = "free_x") +
  scale_fill_manual(values = version_colors, name = NULL, drop = FALSE) +
  labs(x = NULL, y = "AUROC delta vs Baseline (ATAC-only)",
       title = "Change relative to ATAC-only baseline") +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "bottom"
  ) +
  guides(fill = guide_legend(nrow = 1))
ggsave(file.path(out_dir, "fig_main_delta_vs_baseline.png"), p2,
       width = 11, height = 5.4, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "fig_main_delta_vs_baseline.pdf"), p2,
       width = 11, height = 5.4)

heat <- plot_df
heat$panel_trait <- paste(ifelse(heat$setting_lab == paired_label, "paired", "non-paired"),
                          heat$trait_disp, sep = " / ")
trait_row_levels <- c("paired / MONO", "paired / PLT", "paired / LYM",
                      "non-paired / MONO", "non-paired / RDW", "non-paired / LYM")
grid <- expand.grid(
  display_method = label_levels,
  panel_trait = trait_row_levels,
  stringsAsFactors = FALSE
)
heat <- merge(grid, heat[c("display_method", "panel_trait", "value")],
              by = c("display_method", "panel_trait"), all.x = TRUE, sort = FALSE)
heat$has_auroc <- is.finite(heat$value)
heat$display_method <- factor(heat$display_method, levels = label_levels)
heat$panel_trait <- factor(heat$panel_trait, levels = rev(trait_row_levels))

p3 <- ggplot(heat, aes(display_method, panel_trait, fill = value)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = ifelse(has_auroc, sprintf("%.3f", value), "")), size = 3.1) +
  scale_fill_gradient(low = "#F7FBFF", high = "#08519C", limits = c(0, 1),
                      na.value = "grey94", name = "AUROC") +
  labs(x = NULL, y = NULL, title = "Method-by-trait AUROC heatmap") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 35, hjust = 1),
    panel.grid = element_blank()
  )
ggsave(file.path(out_dir, "fig_main_method_heatmap.png"), p3,
       width = 9, height = 5.4, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "fig_main_method_heatmap.pdf"), p3,
       width = 9, height = 5.4)

for (f in c("fig_main_auroc_by_trait.png", "fig_main_delta_vs_baseline.png",
            "fig_main_method_heatmap.png")) {
  file.copy(file.path(out_dir, f), file.path(fig4_dir, f), overwrite = TRUE)
}
write_csv(plot_df, file.path(out_dir, "main_auroc_plot_data.csv"))
write_csv(delta_df, file.path(out_dir, "main_delta_vs_baseline_plot_data.csv"))
write_csv(heat, file.path(out_dir, "main_method_heatmap_plot_data.csv"))

cat("Wrote standardized main figures to:\n  ", out_dir, "\n", sep = "")
