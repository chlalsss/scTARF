# Shared helpers for final-delivery scripts.
# These scripts only write under multiomics_fusion_final_delivery_plan_20260528/final_outputs.

WORK_ROOT <- Sys.getenv("SCTARF_WORK_ROOT", normalizePath(".", winslash = "/", mustWork = FALSE))
DELIVERY_ROOT <- Sys.getenv("SCTARF_DELIVERY_ROOT", file.path(WORK_ROOT, "outputs", "final_delivery"))
SCRIPT_ROOT <- file.path(DELIVERY_ROOT, "scripts")
FINAL_ROOT <- file.path(DELIVERY_ROOT, "final_outputs")

MODULE3_ROOT <- Sys.getenv("M3_OUT_DIR", file.path(WORK_ROOT, "local_data", "BMMC"))
MODULE4_ROOT <- Sys.getenv("SCTARF_NONPAIRED_ROOT", file.path(WORK_ROOT, "outputs", "nonpaired"))
PPT_FIG_ROOT <- Sys.getenv("SCTARF_PPT_FIG_ROOT", file.path(WORK_ROOT, "outputs", "paired_figures"))
PAPER_METRIC_ROOT <- Sys.getenv("SCTARF_PAPER_METRIC_ROOT", file.path(WORK_ROOT, "outputs", "paper_metrics"))
ACCEPTANCE_ROOT <- Sys.getenv("SCTARF_ACCEPTANCE_ROOT", file.path(WORK_ROOT, "outputs", "acceptance_demo_runs"))

dir_ensure <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

ensure_final_dirs <- function() {
  dirs <- file.path(FINAL_ROOT, c(
    "00_method_flow",
    "01_master_tables",
    "02_main_figures",
    "03_ablation_parameter",
    "04_robustness_efficiency",
    "04_robustness_efficiency/efficiency",
    "05_rna_reference_dependency",
    "06_case_study",
    "07_manuscript_text",
    "09_missing_master_table",
    "99_source_manifest"
  ))
  invisible(lapply(dirs, dir_ensure))
}

need_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) {
    stop("Missing R package(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

read_csv_if_exists <- function(path) {
  if (!file.exists(path)) return(NULL)
  if (requireNamespace("data.table", quietly = TRUE)) {
    return(data.table::fread(path, data.table = FALSE, check.names = FALSE))
  }
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

write_csv <- function(x, path) {
  dir_ensure(dirname(path))
  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::fwrite(x, path)
  } else {
    write.csv(x, path, row.names = FALSE, quote = TRUE)
  }
  invisible(path)
}

bind_rows_fill <- function(rows) {
  rows <- Filter(Negate(is.null), rows)
  rows <- rows[vapply(rows, nrow, integer(1)) > 0L]
  if (length(rows) == 0L) return(data.frame())
  all_names <- unique(unlist(lapply(rows, names), use.names = FALSE))
  rows <- lapply(rows, function(x) {
    miss <- setdiff(all_names, names(x))
    for (nm in miss) x[[nm]] <- NA
    x[, all_names, drop = FALSE]
  })
  do.call(rbind, rows)
}

first_existing <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0L) return(NA_character_)
  hit[1]
}

metric_group_from_method <- function(method) {
  m <- tolower(as.character(method))
  out <- rep("fusion", length(m))
  out[grepl("baseline|original|atac-only", m)] <- "original_scavenge"
  out[grepl("scdrs", m)] <- "scdrs"
  out[grepl("top1|top5|v2|v3|fusion|fused", m)] <- "fusion"
  out
}

classify_final_version <- function(method_version, method = NA_character_,
                                   feature_mode = NA_character_,
                                   method_group = NA_character_,
                                   pairing_setting = NA_character_) {
  mv <- tolower(as.character(method_version))
  mt <- tolower(as.character(method))
  fm <- tolower(as.character(feature_mode))
  mg <- tolower(as.character(method_group))
  ps <- tolower(as.character(pairing_setting))

  out <- rep(NA_character_, length(mv))
  is_baseline <- mg == "original_scavenge" | mt == "baseline_atac" | grepl("baseline|original|atac_only", mv)
  is_scdrs <- mg == "scdrs" | grepl("scdrs", mv)
  is_candidate <- !is_baseline & !is_scdrs

  out[is_baseline] <- "Baseline"
  out[is_scdrs] <- "scDRS-only"
  out[is_candidate & grepl("module3_paired_v2", mv)] <- "Paired-Fusion"
  out[is_candidate & grepl("v3_sparse_mapping_top5significant|v3_sparse_other_traits|v3_sparse_lym_granja_rna", mv)] <- "v3.3"
  out[is.na(out) & is_candidate & mg == "fusion" & ps == "paired"] <- "Paired-Fusion"
  out[is.na(out) & is_candidate & mg == "fusion"] <- "Unclassified fusion"
  out
}

final_version_label <- function(version) {
  labels <- c(
    "Baseline" = "Baseline: Original SCAVENGE",
    "scDRS-only" = "scDRS-only: RNA baseline",
    "Paired-Fusion" = "Paired-Fusion: paired RNA-ATAC fusion",
    "v3.3" = "v3.3: sparse mapping + top5 significant fusion",
    "Final-best" = "Final-best: best fusion summary"
  )
  out <- labels[as.character(version)]
  out[is.na(out)] <- as.character(version)[is.na(out)]
  unname(out)
}

final_version_order <- function() {
  c("Baseline", "scDRS-only", "Paired-Fusion", "v3.3", "Final-best")
}

version_definition_table <- function() {
  data.frame(
    version = final_version_order(),
    experiment = c(
      "Original SCAVENGE",
      "RNA-only scDRS baseline",
      "Paired multiome fusion",
      "Sparse mapping + top5 significant fusion",
      "Best fusion summary"
    ),
    raw_results = c(
      "baseline_atac; granja_baseline_auroc.csv; BMMC ATAC-only rows",
      "granja_scdrs_target_enrichment.csv; cite_scdrs_target_enrichment.csv; rna_source_diagnostics_summary.csv; bmmc_paired_scdrs_only_auroc.csv; nonpaired_trait_screening_scdrs_auroc.csv",
      "04_dual_modality/eval_v2_auroc.csv",
      "v3_sparse_mapping_top5significant; v3_sparse_other_traits; v3_sparse_lym_granja_rna",
      "Best AUROC selected from fusion versions"
    ),
    data_setting = c(
      "paired + non-paired",
      "RNA reference only",
      "paired multiome",
      "non-paired",
      "paired + non-paired"
    ),
    traits = c(
      "BMMC: MONO, PLT, LYM_wide; Granja: MONO, PLT, LYM, MCH, RDW, RBC, EOS, WBC",
      "BMMC: MONO, PLT, LYM_wide; RNA reference: MONO, PLT, LYM, EOS, MCH, RBC, RDW, WBC",
      "MONO, PLT, LYM_wide",
      "MONO, RDW, LYM",
      "all evaluated traits"
    ),
    main_method = c(
      "ATAC-only SCAVENGE/RWR/TRS.",
      "Evaluate RNA-reference scDRS scores directly, without ATAC or RWR.",
      "Same-cell RNA+ATAC fusion followed by KNN/RWR/TRS.",
      "Sparse mapping plus top5-significant RNA completion, then fused KNN/RWR/TRS. Other-trait runs are the same method applied to more traits.",
      "Reporting summary only; not a new algorithm."
    ),
    stringsAsFactors = FALSE
  )
}

normalize_auroc_table <- function(x, source_file, dataset, pairing_setting,
                                  method_version = NA_character_,
                                  rna_reference = NA_character_,
                                  top_k = NA,
                                  feature_mode = NA_character_,
                                  note = NA_character_) {
  if (is.null(x) || nrow(x) == 0L) return(NULL)

  if (!"method" %in% names(x)) x$method <- method_version
  if (!"target" %in% names(x) && "trait" %in% names(x)) x$target <- x$trait
  if (!"trait" %in% names(x)) x$trait <- x$target
  if (!"atac_w" %in% names(x)) x$atac_w <- NA_real_
  if (!"auroc" %in% names(x) && "AUROC" %in% names(x)) x$auroc <- x$AUROC
  if (!"auroc_ci_lo" %in% names(x)) x$auroc_ci_lo <- NA_real_
  if (!"auroc_ci_hi" %in% names(x)) x$auroc_ci_hi <- NA_real_
  if (!"n_positive" %in% names(x)) x$n_positive <- NA_integer_
  if (!"n_total" %in% names(x)) x$n_total <- NA_integer_

  out <- data.frame(
    dataset = dataset,
    pairing_setting = pairing_setting,
    trait = as.character(x$trait),
    target_definition = as.character(x$target),
    method_group = metric_group_from_method(x$method),
    method_version = ifelse(is.na(method_version), as.character(x$method), method_version),
    method = as.character(x$method),
    rna_reference = rna_reference,
    top_k = top_k,
    feature_mode = feature_mode,
    atac_w = suppressWarnings(as.numeric(x$atac_w)),
    metric = "AUROC",
    value = suppressWarnings(as.numeric(x$auroc)),
    ci_low = suppressWarnings(as.numeric(x$auroc_ci_lo)),
    ci_high = suppressWarnings(as.numeric(x$auroc_ci_hi)),
    n_positive = suppressWarnings(as.integer(x$n_positive)),
    n_total = suppressWarnings(as.integer(x$n_total)),
    source_file = source_file,
    note = note,
    stringsAsFactors = FALSE
  )
  out$final_version <- classify_final_version(
    method_version = out$method_version,
    method = out$method,
    feature_mode = out$feature_mode,
    method_group = out$method_group,
    pairing_setting = out$pairing_setting
  )
  out$final_version_label <- final_version_label(out$final_version)
  out
}

normalize_scdrs_table <- function(x, source_file, dataset, pairing_setting,
                                  rna_reference, note = NA_character_) {
  if (is.null(x) || nrow(x) == 0L) return(NULL)
  if (!"trait" %in% names(x)) return(NULL)

  value_col <- intersect(c("target_auroc", "cite_score_target_auroc", "auroc"), names(x))[1]
  if (is.na(value_col)) return(NULL)
  n_col <- intersect(c("n_cells", "cite_n_cells"), names(x))[1]
  pos_col <- intersect(c("n_positive_target", "n_positive"), names(x))[1]
  target_col <- intersect(c("target_def", "target", "trait"), names(x))[1]

  out <- data.frame(
    dataset = dataset,
    pairing_setting = pairing_setting,
    trait = as.character(x$trait),
    target_definition = as.character(x[[target_col]]),
    method_group = "scdrs",
    method_version = "scdrs_only",
    method = "scDRS-only",
    rna_reference = rna_reference,
    top_k = NA,
    feature_mode = "rna_score_only",
    atac_w = NA_real_,
    metric = "AUROC",
    value = suppressWarnings(as.numeric(x[[value_col]])),
    ci_low = NA_real_,
    ci_high = NA_real_,
    n_positive = if (!is.na(pos_col)) suppressWarnings(as.integer(x[[pos_col]])) else NA_integer_,
    n_total = if (!is.na(n_col)) suppressWarnings(as.integer(x[[n_col]])) else NA_integer_,
    source_file = source_file,
    note = note,
    stringsAsFactors = FALSE
  )
  out$final_version <- "scDRS-only"
  out$final_version_label <- final_version_label(out$final_version)
  out
}

best_by_key <- function(x, keys, value_col = "value") {
  if (nrow(x) == 0L) return(x)
  split_key <- do.call(paste, c(x[keys], sep = "\r"))
  idx <- unlist(lapply(split(seq_len(nrow(x)), split_key), function(ii) {
    yy <- x[ii, , drop = FALSE]
    ii[which.max(yy[[value_col]])]
  }), use.names = FALSE)
  x[idx, , drop = FALSE]
}

append_manifest <- function(rows, path = file.path(FINAL_ROOT, "99_source_manifest", "result_source_manifest.csv")) {
  rows <- as.data.frame(rows, stringsAsFactors = FALSE)
  if (!"final_version" %in% names(rows)) rows$final_version <- "not_applicable"
  rows$final_version[is.na(rows$final_version) | !nzchar(rows$final_version)] <- "not_applicable"
  if (file.exists(path)) {
    old <- read_csv_if_exists(path)
    if (!"final_version" %in% names(old)) old$final_version <- "not_applicable"
    old$final_version[is.na(old$final_version) | !nzchar(old$final_version)] <- "not_applicable"
    if ("script" %in% names(old) && "script" %in% names(rows)) {
      old <- old[!(old$script %in% unique(rows$script)), , drop = FALSE]
    }
    rows <- bind_rows_fill(list(old, rows))
    rows <- unique(rows)
  }
  write_csv(rows, path)
  invisible(path)
}
