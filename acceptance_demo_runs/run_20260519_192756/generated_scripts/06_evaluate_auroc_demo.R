# 06_evaluate_auroc.R
# §2.F: AUROC eval (baseline / v1 / v2 sweep), 3-way + ablation targets.
# 直接 source(). 自动扫描 OUT_DIR 里所有 TRS_v2_*_w*.tsv.gz, 评估全部.
# 需要 pROC: install.packages("pROC")  (在单独 session 装)

if (!requireNamespace("pROC", quietly = TRUE)) {
  stop("Missing pROC. In a separate R session run: install.packages('pROC')")
}

suppressPackageStartupMessages({
  library(data.table)
  library(pROC)
})

OUT_DIR      <- "D:/RStudio/work/acceptance_demo_runs/run_20260519_192756/outputs/algorithm_core_mono"
META_FILE    <- "D:/SCAVENGE-reproducibility-main/BMMC/01_data_qc/BMMC_cell_meta.rds"
BASELINE_TRS <- "D:/SCAVENGE-reproducibility-main/BMMC/03_scavenge_baseline/TRS_matrix_baseline.rds"

# Target definitions (apple-to-apple v1 风格 + ablation)
TARGETS <- list(
  MONO     = c("CD14+ Mono", "CD16+ Mono"),
  PLT      = c("MK/E prog"),
  LYM      = c("B1 B", "Naive CD20+ B", "Transitional B", "Plasma cell"),
  LYM_TNK  = c("NK", "ILC", "CD4+ T activated", "CD4+ T naive",
               "CD8+ T", "CD8+ T naive"),
  LYM_wide = c("NK", "ILC", "CD4+ T activated", "CD4+ T naive",
               "CD8+ T", "CD8+ T naive", "Lymph prog",
               "B1 B", "Naive CD20+ B", "Transitional B", "Plasma cell")
)
PRED_TRAIT <- list(
  MONO = "MONO", PLT = "PLT",
  LYM = "LYM", LYM_TNK = "LYM", LYM_wide = "LYM"
)

# meta — coerce data.table -> data.frame and set barcode as rownames
meta <- readRDS(META_FILE)
if (inherits(meta, "data.table")) meta <- as.data.frame(meta)
if (is.null(rownames(meta)) ||
    identical(rownames(meta), as.character(seq_len(nrow(meta))))) {
  if ("barcode" %in% colnames(meta)) rownames(meta) <- meta$barcode
  else stop("meta has no rownames and no barcode column")
}
cat("meta dim: ", paste(dim(meta), collapse = " x "), "\n\n")

compute_auroc <- function(predictor, labels) {
  if (sum(labels == 1) < 10 || sum(labels == 0) < 10) return(NA_real_)
  roc_obj <- suppressMessages(pROC::roc(labels, predictor, quiet = TRUE))
  as.numeric(pROC::auc(roc_obj))
}

# load baseline TRS matrix; promote 'cell' column to rownames if present
baseline_mat <- if (file.exists(BASELINE_TRS)) readRDS(BASELINE_TRS) else NULL
if (!is.null(baseline_mat)) {
  if (inherits(baseline_mat, "data.table"))
    baseline_mat <- as.data.frame(baseline_mat)
  if ("cell" %in% colnames(baseline_mat)) {
    rownames(baseline_mat) <- baseline_mat$cell
    baseline_mat$cell <- NULL
  }
  cat("baseline_mat: ",
      paste(dim(baseline_mat), collapse = " x "),
      ", cols: ", paste(colnames(baseline_mat), collapse = ", "), "\n")
}

# scan v2 TRS files
v2_files <- list.files(OUT_DIR,
  pattern = "^TRS_v2_(MONO|PLT|LYM)_w[0-9.]+\\.tsv\\.gz$",
  full.names = TRUE)
cat("v2 files: ", length(v2_files), "\n")
parse_v2 <- function(f) {
  m <- regmatches(basename(f),
                  regexec("^TRS_v2_([A-Z]+)_w([0-9.]+)\\.tsv\\.gz$",
                          basename(f)))[[1]]
  list(trait = m[2], atac_w = as.numeric(m[3]), file = f)
}
v2_meta <- lapply(v2_files, parse_v2)

results <- data.frame()

for (tn in names(TARGETS)) {
  ct_pos <- TARGETS[[tn]]
  pt <- PRED_TRAIT[[tn]]
  cat("\n=== target ", tn, " (positive: ",
      paste(ct_pos, collapse = "/"), ") ===\n", sep = "")

  # baseline
  if (!is.null(baseline_mat) && pt %in% colnames(baseline_mat)) {
    bc <- rownames(baseline_mat)
    common <- intersect(bc, rownames(meta))
    pred <- baseline_mat[common, pt]
    lbl <- as.integer(as.character(meta[common, "cell_type"]) %in% ct_pos)
    a <- compute_auroc(pred, lbl)
    cat(sprintf("  baseline_atac %-10s AUROC = %.3f  (pos=%d / total=%d)\n",
                pt, a, sum(lbl), length(lbl)))
    results <- rbind(results, data.frame(
      method = "baseline_atac", target = tn, atac_w = NA_real_, auroc = a
    ))
  }

  # v1
  v1_file <- file.path(OUT_DIR, paste0("TRS_dual_", pt, ".tsv.gz"))
  if (file.exists(v1_file)) {
    v1 <- fread(v1_file)
    bc_col  <- intersect(c("barcode", "V1", "rn"),  colnames(v1))[1]
    trs_col <- intersect(c("TRS", "trs", "score"), colnames(v1))[1]
    if (is.na(bc_col) || is.na(trs_col)) {
      cat("  v1 file unrecognized cols: ",
          paste(colnames(v1), collapse = ","),
          "; assume col1=barcode col2=TRS\n")
      bc_col  <- colnames(v1)[1]
      trs_col <- colnames(v1)[2]
    }
    pred_v1 <- v1[[trs_col]]
    names(pred_v1) <- v1[[bc_col]]
    common <- intersect(names(pred_v1), rownames(meta))
    lbl <- as.integer(as.character(meta[common, "cell_type"]) %in% ct_pos)
    a <- compute_auroc(pred_v1[common], lbl)
    cat(sprintf("  v1_paired     %-10s AUROC = %.3f\n", pt, a))
    results <- rbind(results, data.frame(
      method = "v1_paired", target = tn, atac_w = NA_real_, auroc = a
    ))
  }

  # v2 sweep
  v2_for <- Filter(function(x) x$trait == pt, v2_meta)
  v2_for <- v2_for[order(sapply(v2_for, `[[`, "atac_w"))]
  for (vm in v2_for) {
    v2 <- fread(vm$file)
    pred_v2 <- v2$TRS
    names(pred_v2) <- v2$barcode
    common <- intersect(names(pred_v2), rownames(meta))
    lbl <- as.integer(as.character(meta[common, "cell_type"]) %in% ct_pos)
    a <- compute_auroc(pred_v2[common], lbl)
    cat(sprintf("  v2 atac_w=%.1f  %-10s AUROC = %.3f\n", vm$atac_w, pt, a))
    results <- rbind(results, data.frame(
      method = "v2", target = tn, atac_w = vm$atac_w, auroc = a
    ))
  }
}

# save full table + headline
out_csv <- file.path(OUT_DIR, "eval_v2_auroc.csv")
write.csv(results, out_csv, row.names = FALSE)
cat("\nsaved: ", basename(out_csv), "\n")

cat("\n", strrep("=", 78), "\n", sep = "")
cat("HEADLINE: best v2 per target vs baseline / v1\n")
cat(strrep("=", 78), "\n")
for (tn in names(TARGETS)) {
  sub <- results[results$target == tn, ]
  base <- sub$auroc[sub$method == "baseline_atac"]
  v1   <- sub$auroc[sub$method == "v1_paired"]
  v2s  <- sub[sub$method == "v2", ]
  if (nrow(v2s) == 0) next
  best <- v2s[which.max(v2s$auroc), ]
  cat(sprintf("  %-9s  baseline=%s  v1=%s  v2_best=%.3f (atac_w=%.1f)\n",
              tn,
              ifelse(length(base) > 0, sprintf("%.3f", base), "    NA"),
              ifelse(length(v1)   > 0, sprintf("%.3f", v1),   "    NA"),
              best$auroc, best$atac_w))
}

cat("\n06_evaluate_auroc.R OK.\n")
