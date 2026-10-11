# Diagnostic evaluation of existing BMMC scores; no model stages are rerun.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop('Usage: Rscript recompute_common_cells.R <BMMC_ROOT> <OUTPUT_DIR>')
}
if (!requireNamespace('pROC', quietly = TRUE)) stop('Install the pROC package.')
bmmc_root <- normalizePath(args[[1]], mustWork = TRUE)
out_dir <- args[[2]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

meta_path <- file.path(bmmc_root, '01_data_qc', 'BMMC_cell_meta.rds')
baseline_path <- file.path(bmmc_root, '03_scavenge_baseline', 'TRS_matrix_baseline.rds')
meta <- as.data.frame(readRDS(meta_path))
if (is.null(rownames(meta)) || identical(rownames(meta), as.character(seq_len(nrow(meta))))) {
  stopifnot('barcode' %in% names(meta))
  rownames(meta) <- meta$barcode
}
stopifnot('cell_type' %in% names(meta), !anyDuplicated(rownames(meta)))
baseline <- as.data.frame(readRDS(baseline_path))
if ('cell' %in% names(baseline)) {
  rownames(baseline) <- baseline$cell
  baseline$cell <- NULL
}
stopifnot(!anyDuplicated(rownames(baseline)))

targets <- list(
  MONO = c('CD14+ Mono', 'CD16+ Mono'),
  PLT = 'MK/E prog',
  LYM = c('NK', 'ILC', 'CD4+ T activated', 'CD4+ T naive', 'CD8+ T',
          'CD8+ T naive', 'Lymph prog', 'B1 B', 'Naive CD20+ B',
          'Transitional B', 'Plasma cell')
)
inputs <- c(meta_path, baseline_path)
rows <- list()
load_scores <- function(trait, weight) {
  path <- file.path(bmmc_root, '04_dual_modality',
                    sprintf('TRS_v2_%s_w%.1f.tsv.gz', trait, weight))
  x <- read.delim(gzfile(path), check.names = FALSE)
  stopifnot(all(c('barcode', 'TRS') %in% names(x)), !anyDuplicated(x$barcode))
  setNames(as.numeric(x$TRS), x$barcode)
}
evaluate <- function(scores, cells, positives) {
  scores <- scores[cells]
  labels <- as.integer(as.character(meta[cells, 'cell_type']) %in% positives)
  stopifnot(all(is.finite(scores)), length(unique(labels)) == 2L)
  fixed <- pROC::roc(labels, as.numeric(scores), levels = c(0, 1), direction = '<', quiet = TRUE)
  auto <- pROC::roc(labels, as.numeric(scores), levels = c(0, 1), direction = 'auto', quiet = TRUE)
  list(n = length(cells), n_positive = sum(labels), fixed = as.numeric(pROC::auc(fixed)),
       auto = as.numeric(pROC::auc(auto)), direction = auto$direction)
}
for (trait in names(targets)) {
  stopifnot(trait %in% names(baseline))
  scores <- list(baseline = setNames(as.numeric(baseline[[trait]]), rownames(baseline)),
                 'w0.7' = load_scores(trait, 0.7), 'w1.0' = load_scores(trait, 1.0))
  paths <- file.path(bmmc_root, '04_dual_modality', sprintf('TRS_v2_%s_w%.1f.tsv.gz', trait, c(0.7, 1.0)))
  inputs <- c(inputs, paths)
  valid <- lapply(scores, function(x) intersect(names(x)[is.finite(x)], rownames(meta)))
  common <- Reduce(intersect, valid)
  for (method in names(scores)) {
    original <- evaluate(scores[[method]], valid[[method]], targets[[trait]])
    shared <- evaluate(scores[[method]], common, targets[[trait]])
    rows[[length(rows) + 1L]] <- data.frame(
      trait = trait, method = method, n_evaluated = original$n, n_positive = original$n_positive,
      auc_original_higher_is_positive = original$fixed, auc_original_proc_auto = original$auto,
      original_proc_direction = original$direction, n_shared = shared$n, n_shared_positive = shared$n_positive,
      auc_shared_higher_is_positive_diagnostic = shared$fixed,
      auc_shared_proc_auto_diagnostic = shared$auto, shared_proc_direction = shared$direction,
      n_excluded_for_shared = original$n - shared$n,
      stringsAsFactors = FALSE
    )
  }
}
result <- do.call(rbind, rows)
write.csv(result, file.path(out_dir, 'common_cell_auroc.csv'), row.names = FALSE)
# Hash the saved inputs to identify the exact score/annotation versions.
hashes <- data.frame(input = normalizePath(inputs), md5 = unname(tools::md5sum(inputs)))
write.csv(hashes, file.path(out_dir, 'input_checksums.csv'), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(out_dir, 'session_info.txt'))
print(result[, c('trait', 'method', 'n_evaluated', 'n_shared', 'auc_original_proc_auto',
                 'auc_shared_proc_auto_diagnostic')], row.names = FALSE)
cat('Maximum absolute AUROC change (historical auto-direction):',
    max(abs(result$auc_shared_proc_auto_diagnostic - result$auc_original_proc_auto)), '\n')
