# Summarize v3.2 top5-significant diagnostic runs across output tags.
# Run after one or more RUN_V3_TOP5_SIGNIFICANT.R executions.

ROOT <- "D:/RStudio/work/module4_traitaware_pseudopair"

if (!exists("V3_COMPARE_TAGS")) {
  output_root <- file.path(ROOT, "02_trait_mapping/output")
  tags <- basename(list.dirs(output_root, recursive = FALSE, full.names = TRUE))
  V3_COMPARE_TAGS <- sort(tags[grepl("^v3_.*top5.*significant", tags)])
}

if (!exists("V3_SUMMARY_OUT_DIR")) {
  V3_SUMMARY_OUT_DIR <- file.path(ROOT, "04_evaluation/v3_top5_significant_diagnostics")
}

dir.create(V3_SUMMARY_OUT_DIR, recursive = TRUE, showWarnings = FALSE)

read_csv_if_exists <- function(path) {
  if (!file.exists(path)) return(NULL)
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

bind_rows_or_empty <- function(rows) {
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0) return(data.frame())
  all_names <- unique(unlist(lapply(rows, names), use.names = FALSE))
  rows <- lapply(rows, function(x) {
    missing_names <- setdiff(all_names, names(x))
    for (nm in missing_names) x[[nm]] <- NA
    x[, all_names, drop = FALSE]
  })
  do.call(rbind, rows)
}

feature_rows <- list()
auroc_rows <- list()

for (tag in V3_COMPARE_TAGS) {
  feature_specs <- list(
    mono = file.path(
      ROOT,
      "02_trait_mapping/output",
      tag,
      "mono_completion_diagnostics/top5_significant_transfer_summary.csv"
    ),
    trait_screening = file.path(
      ROOT,
      "02_trait_mapping/output",
      tag,
      "trait_screening_completion_diagnostics/top5_significant_transfer_summary.csv"
    )
  )

  for (analysis_name in names(feature_specs)) {
    x <- read_csv_if_exists(feature_specs[[analysis_name]])
    if (is.null(x)) next
    x$output_tag <- tag
    x$analysis <- analysis_name
    x$source_file <- feature_specs[[analysis_name]]
    feature_rows[[length(feature_rows) + 1L]] <- x
  }

  auroc_specs <- list(
    mono = file.path(
      ROOT,
      "04_evaluation",
      tag,
      "mono_neighborhood_fdr_ge0p80",
      paste0(tag, "_mono_neighborhood_fdr_ge0p80_auroc.csv")
    ),
    trait_screening = file.path(
      ROOT,
      "04_evaluation",
      tag,
      "trait_screening_ge0p80",
      paste0(tag, "_trait_screening_ge0p80_auroc.csv")
    )
  )

  for (analysis_name in names(auroc_specs)) {
    x <- read_csv_if_exists(auroc_specs[[analysis_name]])
    if (is.null(x)) next
    x$output_tag <- tag
    x$analysis <- analysis_name
    x$source_file <- auroc_specs[[analysis_name]]
    auroc_rows[[length(auroc_rows) + 1L]] <- x
  }
}

feature_summary <- bind_rows_or_empty(feature_rows)
auroc_summary <- bind_rows_or_empty(auroc_rows)

feature_out <- file.path(V3_SUMMARY_OUT_DIR, "v3_top5_significant_feature_summary.csv")
auroc_out <- file.path(V3_SUMMARY_OUT_DIR, "v3_top5_significant_auroc_summary.csv")
best_out <- file.path(V3_SUMMARY_OUT_DIR, "v3_top5_significant_best_by_target.csv")

if (nrow(feature_summary) > 0) {
  feature_keep <- intersect(
    c(
      "output_tag", "analysis", "trait", "significance_method", "b_null_top5",
      "empirical_p_threshold", "filter_q_threshold", "retained_gene_median",
      "retained_gene_mean", "retained_gene_q95", "retained_gene_max",
      "empty_cell_rate", "pseudo_rna_density", "source_file"
    ),
    names(feature_summary)
  )
  write.csv(feature_summary[, feature_keep, drop = FALSE], feature_out, row.names = FALSE)
}

if (nrow(auroc_summary) > 0) {
  write.csv(auroc_summary, auroc_out, row.names = FALSE)

  best_rows <- list()
  keys <- unique(paste(auroc_summary$output_tag, auroc_summary$analysis, auroc_summary$target, sep = "\r"))
  for (key in keys) {
    parts <- strsplit(key, "\r", fixed = TRUE)[[1]]
    y <- auroc_summary[
      auroc_summary$output_tag == parts[1] &
        auroc_summary$analysis == parts[2] &
        auroc_summary$target == parts[3],
      ,
      drop = FALSE
    ]
    baseline <- y[y$method == "baseline_atac", , drop = FALSE]
    candidate <- y[y$method != "baseline_atac", , drop = FALSE]
    if (nrow(candidate) == 0) next
    candidate <- candidate[order(-candidate$auroc), , drop = FALSE]
    best <- candidate[1, , drop = FALSE]
    baseline_auroc <- if (nrow(baseline) > 0) baseline$auroc[1] else NA_real_
    best_rows[[length(best_rows) + 1L]] <- data.frame(
      output_tag = parts[1],
      analysis = parts[2],
      target = parts[3],
      baseline_auroc = baseline_auroc,
      best_method = best$method[1],
      best_atac_w = best$atac_w[1],
      best_auroc = best$auroc[1],
      delta_vs_baseline = best$auroc[1] - baseline_auroc,
      stringsAsFactors = FALSE
    )
  }

  best_summary <- bind_rows_or_empty(best_rows)
  if (nrow(best_summary) > 0) {
    write.csv(best_summary, best_out, row.names = FALSE)
  }
}

cat("\n=== v3.2 top5-significant run summary ===\n")
cat("Tags scanned: ", paste(V3_COMPARE_TAGS, collapse = ", "), "\n", sep = "")
if (file.exists(feature_out)) cat("Feature summary: ", feature_out, "\n", sep = "")
if (file.exists(auroc_out)) cat("AUROC summary: ", auroc_out, "\n", sep = "")
if (file.exists(best_out)) cat("Best-by-target summary: ", best_out, "\n", sep = "")

if (exists("best_summary") && nrow(best_summary) > 0) {
  print(best_summary)
}
