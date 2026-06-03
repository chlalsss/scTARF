# Build ablation and parameter-analysis tables/figures.
# Run after scripts/01_build_master_result_tables.R.

source("D:/RStudio/work/multiomics_fusion_final_delivery_plan_20260528/scripts/_final_delivery_common.R")
ensure_final_dirs()
need_packages(c("ggplot2"))

library(ggplot2)

main_file <- file.path(FINAL_ROOT, "01_master_tables/main_result_comparison_table.csv")
if (!file.exists(main_file)) {
  stop("Missing main result table. Run scripts/01_build_master_result_tables.R first.", call. = FALSE)
}

out_dir <- file.path(FINAL_ROOT, "03_ablation_parameter")
main <- read_csv_if_exists(main_file)
main <- main[main$metric == "AUROC" & is.finite(main$value), , drop = FALSE]
if (!"final_version" %in% names(main)) {
  main$final_version <- classify_final_version(
    main$method_version, main$method, main$feature_mode, main$method_group, main$pairing_setting
  )
}
main$final_version_label <- final_version_label(main$final_version)
main$ablation <- main$final_version_label

baseline <- main[main$final_version == "Baseline", c("dataset", "pairing_setting", "trait", "value")]
names(baseline)[names(baseline) == "value"] <- "baseline_auroc"
abl <- merge(main, baseline, by = c("dataset", "pairing_setting", "trait"), all.x = TRUE)
abl$delta_vs_baseline <- abl$value - abl$baseline_auroc

abl$final_version <- factor(abl$final_version, levels = final_version_order())
abl <- abl[order(abl$dataset, abl$trait, abl$final_version, abl$atac_w), ]
abl$final_version <- as.character(abl$final_version)
write_csv(abl, file.path(out_dir, "ablation_result_table.csv"))

best_abl <- best_by_key(abl, keys = c("dataset", "pairing_setting", "trait", "final_version", "ablation"))
write_csv(best_abl, file.path(out_dir, "ablation_best_by_trait.csv"))

weight_sweep <- abl[abl$method_group == "fusion" & is.finite(abl$atac_w), , drop = FALSE]
write_csv(weight_sweep, file.path(out_dir, "parameter_weight_sweep_table.csv"))

p1 <- ggplot(best_abl, aes(x = trait, y = value, fill = ablation)) +
  geom_col(position = position_dodge(width = 0.82), width = 0.75) +
  facet_wrap(~ pairing_setting + dataset, scales = "free_x") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(x = NULL, y = "AUROC", fill = NULL, title = "Ablation comparison") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")

ggsave(file.path(out_dir, "fig_ablation_by_trait.png"), p1, width = 12, height = 7, dpi = 300)
ggsave(file.path(out_dir, "fig_ablation_by_trait.pdf"), p1, width = 12, height = 7)

if (nrow(weight_sweep) > 0L) {
  p2 <- ggplot(weight_sweep, aes(x = atac_w, y = value, color = method_version, group = method_version)) +
    geom_point(size = 2) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~ trait + pairing_setting, scales = "free_y") +
    labs(x = "ATAC weight", y = "AUROC", color = NULL, title = "ATAC/RNA weight sweep") +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())

  ggsave(file.path(out_dir, "fig_weight_sweep.png"), p2, width = 12, height = 7, dpi = 300)
  ggsave(file.path(out_dir, "fig_weight_sweep.pdf"), p2, width = 12, height = 7)
}

feature_sources <- c(
  file.path(MODULE4_ROOT, "04_evaluation/v3_top5_significant_diagnostics/v3_top5_significant_feature_summary.csv"),
  file.path(MODULE4_ROOT, "03_pseudopair_pipeline/output/v3_top5_significant/transfer_mono/top5_significant_transfer_summary.csv"),
  file.path(MODULE4_ROOT, "03_pseudopair_pipeline/output/v3_sparse_other_traits/transfer_trait_screening/top5_significant_transfer_summary.csv")
)
feature_rows <- lapply(feature_sources, function(f) {
  x <- read_csv_if_exists(f)
  if (is.null(x)) return(NULL)
  x$source_file <- f
  x
})
feature_summary <- bind_rows_fill(feature_rows)
if (nrow(feature_summary) > 0L) {
  write_csv(feature_summary, file.path(out_dir, "top5_significant_feature_filter_summary.csv"))
}

notes <- c(
  "# Ablation Interpretation Notes",
  "",
  "- Use `ablation_result_table.csv` for the complete long-format evidence table.",
  "- Use `ablation_best_by_trait.csv` to discuss the best setting per trait.",
  "- Use `parameter_weight_sweep_table.csv` to decide whether `atac_w = 0.7` or `atac_w = 1.0` is better for each version and trait.",
  "- If `atac_w = 1.0` is strongest, interpret RNA as helping mainly through seed/completion diagnostics rather than direct graph weighting.",
  "- If top5 significant filtering improves stability but not mean AUROC, report it as a denoising/stability module rather than a universal AUROC booster."
)
writeLines(notes, file.path(out_dir, "ablation_interpretation_template.md"), useBytes = TRUE)

cat("Wrote ablation outputs under: ", out_dir, "\n", sep = "")
