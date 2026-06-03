# Update the master long tables and final readable wide table from new
# missing-value result files.
#
# This script does not run analysis. It only reads already-produced CSVs.

WORK_ROOT <- "D:/RStudio/work"
DELIVERY_ROOT <- file.path(WORK_ROOT, "multiomics_fusion_final_delivery_plan_20260528")
MODULE4_ROOT <- file.path(WORK_ROOT, "module4_traitaware_pseudopair")

FINAL_ROOT <- file.path(DELIVERY_ROOT, "final_outputs")
MASTER_DIR <- file.path(FINAL_ROOT, "01_master_tables")
READABLE_TABLE_DIR <- file.path(FINAL_ROOT, "最终输出/01_主结果表")
RESULT_DIR <- file.path(FINAL_ROOT, "09_missing_master_table")

MAIN_FILE <- file.path(MASTER_DIR, "main_result_comparison_table.csv")
BEST_FILE <- file.path(MASTER_DIR, "main_result_best_by_method_version.csv")
WIDE_FILE <- file.path(READABLE_TABLE_DIR, "最终主结果宽表.csv")

if (!exists("BMMC_SCDRS_RESULT")) {
  BMMC_SCDRS_RESULT <- file.path(RESULT_DIR, "bmmc_paired_scdrs_only_auroc.csv")
}
if (!exists("NONPAIRED_SCDRS_RESULT")) {
  NONPAIRED_SCDRS_RESULT <- file.path(RESULT_DIR, "nonpaired_trait_screening_scdrs_auroc.csv")
}
if (!exists("TOP1_AUROC_RESULT")) {
  TOP1_AUROC_RESULT <- file.path(
    MODULE4_ROOT,
    "04_evaluation/v3_dense_top1_missing_traits/neighborhood_fdr_ge0p80",
    "v3_dense_top1_missing_traits_neighborhood_fdr_ge0p80_auroc.csv"
  )
}
if (!exists("OVERWRITE_EXISTING")) OVERWRITE_EXISTING <- FALSE

read_csv_if_exists <- function(path) {
  if (!file.exists(path)) return(NULL)
  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE,
                  na.strings = c("", "NA"))
}

write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, path, row.names = FALSE, quote = TRUE, na = "")
  invisible(path)
}

write_final_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  txt <- utils::capture.output(utils::write.csv(x, row.names = FALSE, na = ""))
  con <- file(path, open = "wb")
  on.exit(close(con), add = TRUE)
  writeBin(as.raw(c(0xEF, 0xBB, 0xBF)), con)
  writeBin(charToRaw(paste0(paste(enc2utf8(txt), collapse = "\n"), "\n")), con)
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

best_by_key <- function(x, keys, value_col = "value") {
  if (nrow(x) == 0L) return(x)
  split_key <- do.call(paste, c(x[keys], sep = "\r"))
  idx <- unlist(lapply(split(seq_len(nrow(x)), split_key), function(ii) {
    yy <- x[ii, , drop = FALSE]
    ii[which.max(suppressWarnings(as.numeric(yy[[value_col]])))]
  }), use.names = FALSE)
  x[idx, , drop = FALSE]
}

standardize_top1 <- function(path) {
  x <- read_csv_if_exists(path)
  if (is.null(x) || nrow(x) == 0L) return(NULL)
  if (!all(c("method", "target", "atac_w", "auroc") %in% names(x))) {
    stop("Top1 AUROC file has unexpected columns: ", path, call. = FALSE)
  }
  method_version <- as.character(x$method)
  method_version[method_version == "v2"] <- "v3_dense_top1_missing_traits"
  data.frame(
    dataset = "Granja non-paired ATAC + RNA reference",
    pairing_setting = "non_paired",
    trait = as.character(x$target),
    target_definition = as.character(x$target),
    method_group = ifelse(grepl("baseline|atac", x$method, ignore.case = TRUE),
                          "original_scavenge", "fusion"),
    method_version = method_version,
    method = as.character(x$method),
    rna_reference = "CITE / Granja RNA reference",
    top_k = ifelse(grepl("baseline|atac", x$method, ignore.case = TRUE), NA, 1),
    feature_mode = "top1_dense_missing_traits",
    atac_w = suppressWarnings(as.numeric(x$atac_w)),
    metric = "AUROC",
    value = suppressWarnings(as.numeric(x$auroc)),
    ci_low = if ("auroc_ci_lo" %in% names(x)) suppressWarnings(as.numeric(x$auroc_ci_lo)) else NA_real_,
    ci_high = if ("auroc_ci_hi" %in% names(x)) suppressWarnings(as.numeric(x$auroc_ci_hi)) else NA_real_,
    n_positive = if ("n_positive" %in% names(x)) suppressWarnings(as.integer(x$n_positive)) else NA_integer_,
    n_total = if ("n_total" %in% names(x)) suppressWarnings(as.integer(x$n_total)) else NA_integer_,
    source_file = path,
    note = "Dense top1 missing-traits non-paired run",
    stringsAsFactors = FALSE
  )
}

is_blank <- function(x) is.na(x) || !nzchar(trimws(as.character(x)))

set_if_needed <- function(df, i, col, value) {
  if (length(value) == 0L || is.na(value[1])) return(df)
  if (OVERWRITE_EXISTING || is_blank(df[[col]][i])) df[[col]][i] <- value[1]
  df
}

pick_best <- function(x) {
  if (is.null(x) || nrow(x) == 0L) return(NULL)
  x <- x[is.finite(suppressWarnings(as.numeric(x$value))), , drop = FALSE]
  if (nrow(x) == 0L) return(NULL)
  x[which.max(as.numeric(x$value)), , drop = FALSE]
}

cat(">>> Update master/wide tables from missing-result CSVs\n\n")
main <- read_csv_if_exists(MAIN_FILE)
if (is.null(main)) stop("Missing main long table: ", MAIN_FILE, call. = FALSE)

new_rows <- bind_rows_fill(list(
  read_csv_if_exists(BMMC_SCDRS_RESULT),
  read_csv_if_exists(NONPAIRED_SCDRS_RESULT),
  standardize_top1(TOP1_AUROC_RESULT)
))

if (nrow(new_rows) > 0L) {
  main <- bind_rows_fill(list(main, new_rows))
  main <- unique(main)
  write_csv(main, MAIN_FILE)
  best <- best_by_key(
    main,
    keys = c("dataset", "pairing_setting", "trait", "method_group",
             "method_version", "feature_mode")
  )
  write_csv(best, BEST_FILE)
  cat("Appended/merged new rows: ", nrow(new_rows), "\n", sep = "")
} else {
  cat("No new result CSVs found; updating wide table from existing long table only.\n")
}

wide <- read_csv_if_exists(WIDE_FILE)
if (is.null(wide)) stop("Missing final wide table: ", WIDE_FILE, call. = FALSE)

for (i in seq_len(nrow(wide))) {
  dataset <- wide[["数据集"]][i]
  pairing <- wide[["数据类型"]][i]
  trait <- wide[["性状"]][i]

  if (identical(pairing, "paired")) {
    scdrs <- pick_best(main[main$dataset == dataset &
                              main$pairing_setting == pairing &
                              main$trait == trait &
                              main$method_group == "scdrs", , drop = FALSE])
    fused <- pick_best(main[main$dataset == dataset &
                              main$pairing_setting == pairing &
                              main$trait == trait &
                              main$method_group == "fusion" &
                              main$method_version == "module3_paired_v2" &
                              suppressWarnings(as.numeric(main$atac_w)) == 0.7, , drop = FALSE])
    if (!is.null(scdrs)) {
      wide <- set_if_needed(wide, i, "纯_scDRS_AUROC", scdrs$value)
      wide <- set_if_needed(wide, i, "纯scDRS来源", scdrs$source_file)
    }
    if (!is.null(fused)) {
      wide <- set_if_needed(wide, i, "v2早期版本_AUROC", fused$value)
      wide <- set_if_needed(wide, i, "top5或最终融合_AUROC", fused$value)
    }
  }

  if (identical(pairing, "non_paired")) {
    scdrs <- pick_best(main[main$method_group == "scdrs" &
                              main$trait == trait, , drop = FALSE])
    top1 <- pick_best(main[main$dataset == dataset &
                             main$pairing_setting == pairing &
                             main$trait == trait &
                             main$method_group == "fusion" &
                             grepl("top1", main$method_version, ignore.case = TRUE), , drop = FALSE])
    if (!is.null(scdrs)) {
      wide <- set_if_needed(wide, i, "纯_scDRS_AUROC", scdrs$value)
      wide <- set_if_needed(wide, i, "纯scDRS来源", scdrs$source_file)
    }
    if (!is.null(top1)) {
      wide <- set_if_needed(wide, i, "top1_AUROC", top1$value)
    }
  }
}

write_final_csv(wide, WIDE_FILE)
cat("Wrote wide table: ", WIDE_FILE, "\n", sep = "")
