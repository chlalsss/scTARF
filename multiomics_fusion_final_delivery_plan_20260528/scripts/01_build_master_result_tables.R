# Build unified main result tables from existing paired and non-paired outputs.
# This script reads existing result CSVs only; it does not rerun original analyses.

source("D:/RStudio/work/multiomics_fusion_final_delivery_plan_20260528/scripts/_final_delivery_common.R")
ensure_final_dirs()

out_dir <- file.path(FINAL_ROOT, "01_master_tables")

result_specs <- data.frame(
  source_file = c(
    file.path(PPT_FIG_ROOT, "tables/roc_table45_auc_check.csv"),
    file.path(MODULE4_ROOT, "04_evaluation/granja_baseline_auroc.csv"),
    file.path(MODULE4_ROOT, "04_evaluation/module4_v2/module4_v2_auroc.csv"),
    file.path(MODULE4_ROOT, "04_evaluation/module4_v2_neighborhood_fdr_ge0p80/module4_v2_neighborhood_fdr_ge0p80_auroc.csv"),
    file.path(MODULE4_ROOT, "04_evaluation/trait_screening_ge0p80/trait_screening_ge0p80_auroc.csv"),
    file.path(MODULE4_ROOT, "04_evaluation/trait_screening_remaining_ge0p80/trait_screening_remaining_ge0p80_auroc.csv"),
    file.path(MODULE4_ROOT, "04_evaluation/v3_dense_top1/mono_neighborhood_fdr_ge0p80/v3_dense_top1_mono_neighborhood_fdr_ge0p80_auroc.csv"),
    file.path(MODULE4_ROOT, "04_evaluation/v3_dense_top1/trait_screening_ge0p80/v3_dense_top1_trait_screening_ge0p80_auroc.csv"),
    file.path(MODULE4_ROOT, "04_evaluation/v3_top5_significant/mono_neighborhood_fdr_ge0p80/v3_top5_significant_mono_neighborhood_fdr_ge0p80_auroc.csv"),
    file.path(MODULE4_ROOT, "04_evaluation/v3_top5_significant_p005/mono_neighborhood_fdr_ge0p80/v3_top5_significant_p005_mono_neighborhood_fdr_ge0p80_auroc.csv"),
    file.path(MODULE4_ROOT, "04_evaluation/v3_sparse_mapping_top5significant/mono_neighborhood_fdr_ge0p80/v3_sparse_mapping_top5significant_mono_neighborhood_fdr_ge0p80_auroc.csv"),
    file.path(MODULE4_ROOT, "04_evaluation/v3_sparse_other_traits/neighborhood_fdr_ge0p80/v3_sparse_other_traits_neighborhood_fdr_ge0p80_auroc.csv"),
    file.path(MODULE4_ROOT, "04_evaluation/v3_sparse_lym_granja_rna/neighborhood_fdr_ge0p80/v3_sparse_lym_granja_rna_neighborhood_fdr_ge0p80_auroc.csv"),
    file.path(MODULE4_ROOT, "04_evaluation/v3_top5_significant_diagnostics/v3_top5_significant_auroc_summary.csv")
  ),
  dataset = c(
    "BMMC paired multiome",
    rep("Granja non-paired ATAC + RNA reference", 13)
  ),
  pairing_setting = c("paired", rep("non_paired", 13)),
  method_version = c(
    "module3_paired_v2",
    "baseline_atac",
    "module4_v2",
    "module4_v2_neighborhood_fdr",
    "module4_v2_trait_screening",
    "module4_v2_trait_screening_remaining",
    "v3_dense_top1",
    "v3_dense_top1",
    "v3_top5_significant",
    "v3_top5_significant_p005",
    "v3_sparse_mapping_top5significant",
    "v3_sparse_other_traits",
    "v3_sparse_lym_granja_rna",
    "v3_top5_significant_summary"
  ),
  rna_reference = c(
    "paired BMMC RNA",
    "none",
    rep("CITE / Granja RNA reference", 10),
    "Granja RNA reference",
    "CITE / Granja RNA reference"
  ),
  top_k = c(NA, NA, NA, NA, NA, NA, 1, 1, 5, 5, 5, 5, 5, 5),
  feature_mode = c(
    "paired_fused_trs",
    "atac_only",
    "dense_completion",
    "dense_completion_neighborhood_fdr",
    "trait_screening",
    "trait_screening_remaining",
    "top1_dense",
    "top1_dense",
    "top5_significant",
    "top5_significant_p005",
    "sparse_mapping_top5_significant",
    "sparse_other_traits_top5_significant",
    "sparse_lym_granja_rna",
    "top5_significant_summary"
  ),
  note = c(
    "Prepared by ppt_dataset_figures/03_plot_paper_roc_curves.R",
    "Original ATAC-only SCAVENGE baseline",
    "Earlier non-paired v2 output",
    "v2 with neighborhood-FDR seed selection",
    "trait screening output",
    "remaining trait screening output",
    "top1 mono output",
    "top1 trait-screening output",
    "top5 significant mono output",
    "top5 significant mono p<=0.05 output",
    "sparse mapping + top5 significant mono output",
    "sparse workflow for non-MONO traits",
    "LYM v3.3 using Granja RNA reference (RNA-reference selection improvement)",
    "aggregated diagnostic summary"
  ),
  stringsAsFactors = FALSE
)

rows <- list()
for (i in seq_len(nrow(result_specs))) {
  spec <- result_specs[i, ]
  x <- read_csv_if_exists(spec$source_file)
  if (is.null(x)) next

  if (identical(spec$method_version, "module3_paired_v2")) {
    if ("target" %in% names(x) && !"trait" %in% names(x)) x$trait <- x$target
    x$method <- as.character(x$method)
    x$atac_w <- ifelse(grepl("0.7", x$method), 0.7,
                       ifelse(grepl("1.0", x$method), 1.0, NA_real_))
  }

  rows[[length(rows) + 1L]] <- normalize_auroc_table(
    x = x,
    source_file = spec$source_file,
    dataset = spec$dataset,
    pairing_setting = spec$pairing_setting,
    method_version = spec$method_version,
    rna_reference = spec$rna_reference,
    top_k = spec$top_k,
    feature_mode = spec$feature_mode,
    note = spec$note
  )
}

scdrs_specs <- data.frame(
  source_file = c(
    file.path(MODULE4_ROOT, "01_data_prep/output/granja_scdrs_target_enrichment.csv"),
    file.path(MODULE4_ROOT, "01_data_prep_external_rna/output/cite_scdrs_target_enrichment.csv"),
    file.path(MODULE4_ROOT, "01_data_prep_external_rna/output/trait_screening/cite_scdrs_target_enrichment_trait_screening.csv"),
    file.path(MODULE4_ROOT, "04_evaluation/module4_v2_diagnostics/rna_sources/rna_source_diagnostics_summary.csv")
  ),
  dataset = "RNA reference scDRS",
  pairing_setting = "rna_reference_only",
  rna_reference = c("Granja RNA", "NeurIPS CITE RNA", "NeurIPS CITE RNA trait-screening", "NeurIPS CITE RNA diagnostics"),
  note = c(
    "Generated by 01_data_prep/03_rna_scdrs_python/summarize_scdrs_granja.R",
    "Generated by 01_data_prep_external_rna/05_summarize_cite_scdrs.R",
    "Generated by 01_data_prep_external_rna/08_summarize_cite_scdrs_trait_screening.R",
    "CITE score target AUROC from RNA source diagnostics"
  ),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(scdrs_specs))) {
  spec <- scdrs_specs[i, ]
  x <- read_csv_if_exists(spec$source_file)
  if (is.null(x)) next
  rows[[length(rows) + 1L]] <- normalize_scdrs_table(
    x = x,
    source_file = spec$source_file,
    dataset = spec$dataset,
    pairing_setting = spec$pairing_setting,
    rna_reference = spec$rna_reference,
    note = spec$note
  )
}

supplemental_scdrs_files <- c(
  file.path(FINAL_ROOT, "09_missing_master_table/bmmc_paired_scdrs_only_auroc.csv"),
  file.path(FINAL_ROOT, "09_missing_master_table/nonpaired_trait_screening_scdrs_auroc.csv")
)
for (path in supplemental_scdrs_files) {
  x <- read_csv_if_exists(path)
  if (is.null(x)) next
  rows[[length(rows) + 1L]] <- x
}

main <- bind_rows_fill(rows)
if (nrow(main) == 0L) stop("No result files were found.", call. = FALSE)

main <- main[is.finite(main$value), , drop = FALSE]
main$final_version <- classify_final_version(
  method_version = main$method_version,
  method = main$method,
  feature_mode = main$feature_mode,
  method_group = main$method_group,
  pairing_setting = main$pairing_setting
)
main$final_version_label <- final_version_label(main$final_version)
main$final_version <- factor(main$final_version, levels = final_version_order())
main <- main[order(main$dataset, main$trait, main$final_version, main$method_group, main$method_version, main$atac_w), ]
main$final_version <- as.character(main$final_version)

main_file <- file.path(out_dir, "main_result_comparison_table.csv")
write_csv(main, main_file)

best <- best_by_key(
  main,
  keys = c("dataset", "pairing_setting", "trait", "method_group", "method_version", "feature_mode")
)
best_file <- file.path(out_dir, "main_result_best_by_method_version.csv")
write_csv(best, best_file)

best_final <- best_by_key(
  main,
  keys = c("dataset", "pairing_setting", "trait", "final_version")
)
best_final_file <- file.path(out_dir, "main_result_best_by_final_version.csv")
write_csv(best_final, best_final_file)

source_manifest <- unique(main[c(
  "source_file", "dataset", "pairing_setting", "trait", "method_group",
  "method_version", "method", "feature_mode", "final_version", "note"
)])
source_manifest$final_file <- main_file
source_manifest$script <- "scripts/01_build_master_result_tables.R"
source_exists <- function(path) {
  parts <- strsplit(as.character(path), ";", fixed = TRUE)
  vapply(parts, function(x) {
    x <- trimws(x)
    length(x) > 0L && all(file.exists(x))
  }, logical(1))
}
source_manifest$status <- ifelse(source_exists(source_manifest$source_file), "found", "missing")
source_manifest <- source_manifest[c(
  "final_file", "source_file", "script", "dataset", "pairing_setting",
  "trait", "method_group", "method_version", "method", "feature_mode",
  "final_version", "status", "note"
)]
manifest_file <- file.path(FINAL_ROOT, "99_source_manifest/source_manifest_master_results.csv")
write_csv(source_manifest, manifest_file)
append_manifest(source_manifest)

version_file <- file.path(out_dir, "method_version_definitions.csv")
write_csv(version_definition_table(), version_file)

cat("Wrote:\n")
cat("  ", main_file, "\n", sep = "")
cat("  ", best_file, "\n", sep = "")
cat("  ", best_final_file, "\n", sep = "")
cat("  ", manifest_file, "\n", sep = "")
cat("  ", version_file, "\n", sep = "")
