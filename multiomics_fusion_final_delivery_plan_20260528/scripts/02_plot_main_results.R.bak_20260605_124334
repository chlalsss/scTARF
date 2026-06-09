# Plot main AUROC comparison figures from the unified table.
# Run after scripts/01_build_master_result_tables.R.

source("D:/RStudio/work/multiomics_fusion_final_delivery_plan_20260528/scripts/_final_delivery_common.R")
ensure_final_dirs()
need_packages(c("ggplot2"))

library(ggplot2)

in_file <- file.path(FINAL_ROOT, "01_master_tables/main_result_comparison_table.csv")
if (!file.exists(in_file)) {
  stop("Missing main table. Run scripts/01_build_master_result_tables.R first.", call. = FALSE)
}

out_dir <- file.path(FINAL_ROOT, "02_main_figures")
x <- read_csv_if_exists(in_file)
x <- x[x$metric == "AUROC" & is.finite(x$value), , drop = FALSE]
if (!"final_version" %in% names(x)) {
  x$final_version <- classify_final_version(
    x$method_version, x$method, x$feature_mode, x$method_group, x$pairing_setting
  )
}
x$display_method <- final_version_label(x$final_version)

keep_versions <- c("Baseline", "scDRS-only", "Paired-Fusion", "v2.1", "v2.2", "v3.1", "v3.2", "v3.3")
plot_df <- x[x$final_version %in% keep_versions, , drop = FALSE]
plot_df <- best_by_key(plot_df, keys = c("dataset", "pairing_setting", "trait", "final_version"))

plot_df$trait <- factor(plot_df$trait, levels = sort(unique(plot_df$trait)))
plot_df$final_version <- factor(plot_df$final_version, levels = keep_versions)
plot_df$display_method <- factor(plot_df$display_method, levels = final_version_label(keep_versions))

p1 <- ggplot(plot_df, aes(x = trait, y = value, fill = display_method)) +
  geom_col(position = position_dodge(width = 0.82), width = 0.75) +
  facet_wrap(~ pairing_setting + dataset, scales = "free_x") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(x = NULL, y = "AUROC", fill = NULL, title = "Main AUROC comparison by trait") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")

ggsave(file.path(out_dir, "fig_main_auroc_by_trait.png"), p1, width = 12, height = 7, dpi = 300)
ggsave(file.path(out_dir, "fig_main_auroc_by_trait.pdf"), p1, width = 12, height = 7)

baseline <- x[x$final_version == "Baseline", c("dataset", "pairing_setting", "trait", "value")]
names(baseline)[names(baseline) == "value"] <- "baseline_auroc"
delta_df <- merge(plot_df[plot_df$method_group != "original_scavenge", ],
                  baseline, by = c("dataset", "pairing_setting", "trait"), all.x = TRUE)
delta_df$delta_vs_baseline <- delta_df$value - delta_df$baseline_auroc
delta_df <- delta_df[is.finite(delta_df$delta_vs_baseline), , drop = FALSE]

p2 <- ggplot(delta_df, aes(x = trait, y = delta_vs_baseline, fill = display_method)) +
  geom_hline(yintercept = 0, color = "grey40", linewidth = 0.4) +
  geom_col(position = position_dodge(width = 0.82), width = 0.75) +
  facet_wrap(~ pairing_setting + dataset, scales = "free_x") +
  labs(x = NULL, y = "AUROC delta vs original SCAVENGE", fill = NULL,
       title = "Change relative to original SCAVENGE") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")

ggsave(file.path(out_dir, "fig_main_delta_vs_baseline.png"), p2, width = 12, height = 7, dpi = 300)
ggsave(file.path(out_dir, "fig_main_delta_vs_baseline.pdf"), p2, width = 12, height = 7)

heat <- best_by_key(plot_df, keys = c("dataset", "pairing_setting", "trait", "display_method"))
heat$panel_trait <- paste(heat$pairing_setting, heat$trait, sep = " / ")
panel_levels <- unique(heat$panel_trait[order(heat$pairing_setting, heat$trait)])
method_levels <- final_version_label(keep_versions)
heat_grid <- expand.grid(
  display_method = method_levels,
  panel_trait = panel_levels,
  stringsAsFactors = FALSE
)
heat <- merge(
  heat_grid,
  heat[c("display_method", "panel_trait", "value", "dataset", "pairing_setting", "trait", "source_file", "note")],
  by = c("display_method", "panel_trait"),
  all.x = TRUE,
  sort = FALSE
)
heat$has_auroc <- is.finite(heat$value)
panel_parts <- strsplit(as.character(heat$panel_trait), " / ", fixed = TRUE)
panel_pairing <- vapply(panel_parts, `[`, character(1), 1)
panel_trait <- vapply(panel_parts, `[`, character(1), 2)
missing_pairing <- is.na(heat$pairing_setting) | !nzchar(as.character(heat$pairing_setting))
missing_trait <- is.na(heat$trait) | !nzchar(as.character(heat$trait))
heat$pairing_setting[missing_pairing] <- panel_pairing[missing_pairing]
heat$trait[missing_trait] <- panel_trait[missing_trait]
heat$dataset[is.na(heat$dataset) & heat$pairing_setting == "paired"] <- "BMMC paired multiome"
heat$dataset[is.na(heat$dataset) & heat$pairing_setting == "non_paired"] <- "Granja non-paired ATAC + RNA reference"
heat$dataset[is.na(heat$dataset) & heat$pairing_setting == "rna_reference_only"] <- "RNA reference scDRS"
version_by_label <- setNames(keep_versions, method_levels)
heat$final_version <- unname(version_by_label[as.character(heat$display_method)])
heat$display_method <- factor(heat$display_method, levels = method_levels)
heat$panel_trait <- factor(heat$panel_trait, levels = rev(panel_levels))

p3 <- ggplot(heat, aes(x = display_method, y = panel_trait, fill = value)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(has_auroc, sprintf("%.3f", value), "")), size = 3) +
  scale_fill_gradient(low = "#F7FBFF", high = "#08519C", limits = c(0, 1), na.value = "grey95", name = "AUROC") +
  labs(x = NULL, y = NULL, title = "Method-by-trait AUROC heatmap") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())

ggsave(file.path(out_dir, "fig_main_method_heatmap.png"), p3, width = 11, height = 8, dpi = 300)
ggsave(file.path(out_dir, "fig_main_method_heatmap.pdf"), p3, width = 11, height = 8)

write_csv(delta_df, file.path(out_dir, "main_delta_vs_baseline_plot_data.csv"))
write_csv(heat, file.path(out_dir, "main_method_heatmap_plot_data.csv"))

cat("Wrote main figures under: ", out_dir, "\n", sep = "")
