# Build RNA-reference dependency tables and plots from existing scDRS/source diagnostics.

source("D:/RStudio/work/multiomics_fusion_final_delivery_plan_20260528/scripts/_final_delivery_common.R")
ensure_final_dirs()
need_packages(c("ggplot2"))

library(ggplot2)

out_dir <- file.path(FINAL_ROOT, "05_rna_reference_dependency")

ref_specs <- data.frame(
  source_file = c(
    file.path(MODULE4_ROOT, "01_data_prep/output/granja_scdrs_target_enrichment.csv"),
    file.path(MODULE4_ROOT, "01_data_prep_external_rna/output/cite_scdrs_target_enrichment.csv"),
    file.path(MODULE4_ROOT, "04_evaluation/module4_v2_diagnostics/rna_sources/rna_source_diagnostics_summary.csv")
  ),
  rna_reference = c("Granja RNA", "NeurIPS CITE RNA", "NeurIPS CITE RNA source diagnostics"),
  stringsAsFactors = FALSE
)

rows <- list()
for (i in seq_len(nrow(ref_specs))) {
  spec <- ref_specs[i, ]
  x <- read_csv_if_exists(spec$source_file)
  if (is.null(x)) next
  value_col <- intersect(c("target_auroc", "cite_score_target_auroc", "pseudo_z_target_auroc"), names(x))[1]
  if (is.na(value_col) || !"trait" %in% names(x)) next
  rows[[length(rows) + 1L]] <- data.frame(
    atac_dataset = "Granja ATAC",
    rna_reference = spec$rna_reference,
    trait = as.character(x$trait),
    score_type = value_col,
    auroc = suppressWarnings(as.numeric(x[[value_col]])),
    source_file = spec$source_file,
    stringsAsFactors = FALSE
  )
}

ref_table <- bind_rows_fill(rows)
if (nrow(ref_table) > 0L) {
  write_csv(ref_table, file.path(out_dir, "rna_reference_dependency_table.csv"))

  p1 <- ggplot(ref_table, aes(x = trait, y = auroc, fill = rna_reference)) +
    geom_col(position = position_dodge(width = 0.82), width = 0.75) +
    facet_wrap(~ score_type) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(x = NULL, y = "AUROC", fill = NULL, title = "RNA reference dependency") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.minor = element_blank(),
          legend.position = "bottom")
  ggsave(file.path(out_dir, "fig_rna_reference_dependency.png"), p1, width = 10, height = 6, dpi = 300)
}

source_diag_file <- file.path(MODULE4_ROOT, "04_evaluation/module4_v2_diagnostics/rna_sources/rna_source_diagnostics_summary.csv")
source_diag <- read_csv_if_exists(source_diag_file)
if (!is.null(source_diag) && nrow(source_diag) > 0L) {
  keep <- intersect(c(
    "trait",
    "cite_score_target_auroc",
    "pseudo_z_target_auroc",
    "source_target_weight_auroc",
    "mean_source_target_weight_seed",
    "mean_source_myeloid_weight_seed",
    "top_seed_source_celltypes_by_weight"
  ), names(source_diag))
  write_csv(source_diag[, keep, drop = FALSE],
            file.path(out_dir, "rna_source_diagnostics_compact.csv"))

  long <- data.frame()
  for (metric in intersect(c("cite_score_target_auroc", "pseudo_z_target_auroc", "source_target_weight_auroc"), names(source_diag))) {
    long <- rbind(long, data.frame(
      trait = source_diag$trait,
      metric = metric,
      value = suppressWarnings(as.numeric(source_diag[[metric]])),
      stringsAsFactors = FALSE
    ))
  }
  if (nrow(long) > 0L) {
    p2 <- ggplot(long, aes(x = trait, y = value, fill = metric)) +
      geom_col(position = position_dodge(width = 0.82), width = 0.75) +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
      labs(x = NULL, y = "AUROC", fill = NULL, title = "RNA source diagnostic AUROC") +
      theme_bw(base_size = 11) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "bottom",
            panel.grid.minor = element_blank())
    ggsave(file.path(out_dir, "fig_rna_source_diagnostic_auroc.png"), p2, width = 9, height = 5, dpi = 300)
  }
}

notes <- c(
  "# RNA Reference Dependency Interpretation Template",
  "",
  "- Compare Granja RNA and CITE RNA scDRS-only AUROC by trait.",
  "- Use source diagnostics to check whether top RNA neighbors come from expected target cell types.",
  "- A stable method should preserve both AUROC and high-risk cell-type localization after replacing the RNA reference.",
  "- If one trait depends strongly on a specific reference, report it as a boundary condition rather than hiding it."
)
writeLines(notes, file.path(out_dir, "rna_reference_dependency_interpretation.md"), useBytes = TRUE)

cat("Wrote RNA-reference dependency outputs under: ", out_dir, "\n", sep = "")
