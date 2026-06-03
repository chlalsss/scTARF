# Compile manuscript/report text templates from final-output tables.

source("D:/RStudio/work/multiomics_fusion_final_delivery_plan_20260528/scripts/_final_delivery_common.R")
ensure_final_dirs()

out_dir <- file.path(FINAL_ROOT, "07_manuscript_text")

main_file <- file.path(FINAL_ROOT, "01_master_tables/main_result_comparison_table.csv")
ablation_file <- file.path(FINAL_ROOT, "03_ablation_parameter/ablation_best_by_trait.csv")
rna_file <- file.path(FINAL_ROOT, "05_rna_reference_dependency/rna_reference_dependency_table.csv")

write_section <- function(path, lines) {
  writeLines(lines, path, useBytes = TRUE)
  invisible(path)
}

method_lines <- c(
  "# Methods: Trait-aware RNA-ATAC Fusion SCAVENGE",
  "",
  "We evaluated paired and non-paired multi-omics settings using a unified downstream SCAVENGE/RWR/TRS framework.",
  "",
  "In paired multiome data, ATAC and RNA feature blocks were constructed for the same cells and used to build fused KNN graphs before trait seed selection and RWR propagation.",
  "",
  "In non-paired data, scATAC cells were matched to scRNA reference cells. RNA reference cells were scored by scDRS, top-k RNA neighbors were used to complete trait-aware RNA signals for each ATAC cell, and top5-completed signals were filtered against random top5 backgrounds before graph construction.",
  "",
  "The final TRS was evaluated against post-hoc trait-relevant cell-type definitions using AUROC and related metrics."
)
write_section(file.path(out_dir, "methods_final_workflow.md"), method_lines)

version_lines <- c(
  "# Method Version Definitions",
  "",
  "Final reporting uses method-level version names. If the same algorithm was run on different traits in separate batches, those runs are merged into the same version.",
  "",
  capture.output(print(version_definition_table(), row.names = FALSE))
)
write_section(file.path(out_dir, "method_version_definitions.md"), version_lines)

main_lines <- c(
  "# Results: Main Method Comparison",
  "",
  "Input table expected:",
  paste0("`", main_file, "`"),
  "",
  "Report SCAVENGE, scDRS-only, Paired-Fusion, v2.1, v2.2, v3.1, v3.2, and v3.3 performance across paired and non-paired settings. Include both improved and non-improved traits."
)
if (file.exists(main_file)) {
  x <- read_csv_if_exists(main_file)
  x <- x[x$metric == "AUROC" & is.finite(x$value), , drop = FALSE]
  if (!"final_version" %in% names(x)) {
    x$final_version <- classify_final_version(
      x$method_version, x$method, x$feature_mode, x$method_group, x$pairing_setting
    )
  }
  x$final_version_label <- final_version_label(x$final_version)
  top <- head(x[order(-x$value), c("dataset", "trait", "final_version", "final_version_label", "method_version", "atac_w", "value")], 12)
  main_lines <- c(main_lines, "", "Current top AUROC rows:", "", capture.output(print(top, row.names = FALSE)))
}
write_section(file.path(out_dir, "results_main_comparison.md"), main_lines)

ablation_lines <- c(
  "# Results: Ablation and Parameter Analysis",
  "",
  paste0("Input table expected: `", ablation_file, "`"),
  "",
  "Discuss the contribution of RNA reference scoring, v2.1 basic pseudo-pair fusion, v2.2 neighborhood-FDR seed selection, v3.1 top1 matching, v3.2 top5 significant filtering, v3.3 sparse mapping, and ATAC/RNA weight choices."
)
if (file.exists(ablation_file)) {
  x <- read_csv_if_exists(ablation_file)
  keep <- intersect(c("dataset", "trait", "ablation", "atac_w", "value", "delta_vs_baseline"), names(x))
  ablation_lines <- c(ablation_lines, "", "Current ablation rows:", "", capture.output(print(head(x[, keep, drop = FALSE], 20), row.names = FALSE)))
}
write_section(file.path(out_dir, "results_ablation.md"), ablation_lines)

eff_lines <- c(
  "# Results: Efficiency",
  "",
  "Use `final_outputs/04_robustness_efficiency/efficiency/` for runtime and memory tables.",
  "",
  "This text intentionally excludes robustness interpretation, which should be prepared separately."
)
write_section(file.path(out_dir, "results_efficiency.md"), eff_lines)

rna_lines <- c(
  "# Results: RNA Reference Dependency",
  "",
  paste0("Input table expected: `", rna_file, "`"),
  "",
  "Compare whether scDRS target enrichment and RNA-source diagnostics remain stable across RNA references."
)
if (file.exists(rna_file)) {
  x <- read_csv_if_exists(rna_file)
  rna_lines <- c(rna_lines, "", "Current RNA reference rows:", "", capture.output(print(x, row.names = FALSE)))
}
write_section(file.path(out_dir, "results_rna_reference_dependency.md"), rna_lines)

case_lines <- c(
  "# Case Study Template",
  "",
  "Choose one trait with either clear AUROC gain or stronger biological interpretability. Suggested candidates: MONO, RDW/RBC, EOS, or WBC.",
  "",
  "Required figures: high-risk cell-type bar plot, original-vs-fusion TRS scatter/rank-shift plot, and UMAP if coordinates are available.",
  "",
  "Required interpretation: data source, disease/trait, high-risk cell types, differences from original SCAVENGE, and biological or precision-medicine meaning."
)
write_section(file.path(out_dir, "case_study.md"), case_lines)

cat("Wrote manuscript text templates under: ", out_dir, "\n", sep = "")
