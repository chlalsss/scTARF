# 06b_evaluate_auroc_dualseed.R
# §2.F' AUROC eval — 4-way 对比 (baseline / v1 / v2 SCDRS-seed / v2 dual-seed).
# fork 自 06_evaluate_auroc.R, 加扫描 TRS_dualseed_* 文件, 输出 eval_dualseed_auroc.csv.
# 不动 06_evaluate_auroc.R / eval_v2_auroc.csv, 保留阶段 2 主结果.

if (!requireNamespace("pROC", quietly = TRUE)) {
  stop("Missing pROC. In a separate R session run: install.packages('pROC')")
}

suppressPackageStartupMessages({
  library(data.table)
  library(pROC)
})

OUT_DIR      <- "D:/SCAVENGE-reproducibility-main/BMMC/04_dual_modality"
META_FILE    <- "D:/SCAVENGE-reproducibility-main/BMMC/01_data_qc/BMMC_cell_meta.rds"
BASELINE_TRS <- "D:/SCAVENGE-reproducibility-main/BMMC/03_scavenge_baseline/TRS_matrix_baseline.rds"

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

# baseline TRS
baseline_mat <- if (file.exists(BASELINE_TRS)) readRDS(BASELINE_TRS) else NULL
if (!is.null(baseline_mat)) {
  if (inherits(baseline_mat, "data.table"))
    baseline_mat <- as.data.frame(baseline_mat)
  if ("cell" %in% colnames(baseline_mat)) {
    rownames(baseline_mat) <- baseline_mat$cell
    baseline_mat$cell <- NULL
  }
}

# 扫所有 v2 TRS 文件 (SCDRS-seed)
v2_files <- list.files(OUT_DIR,
  pattern = "^TRS_v2_(MONO|PLT|LYM)_w[0-9.]+\\.tsv\\.gz$",
  full.names = TRUE)
parse_file <- function(f, prefix) {
  m <- regmatches(basename(f),
                  regexec(paste0("^", prefix, "_([A-Z]+)_w([0-9.]+)\\.tsv\\.gz$"),
                          basename(f)))[[1]]
  list(trait = m[2], atac_w = as.numeric(m[3]), file = f)
}
v2_meta  <- lapply(v2_files, parse_file, prefix = "TRS_v2")

# 扫所有 dual-seed TRS 文件
ds_files <- list.files(OUT_DIR,
  pattern = "^TRS_dualseed_(MONO|PLT|LYM)_w[0-9.]+\\.tsv\\.gz$",
  full.names = TRUE)
ds_meta  <- lapply(ds_files, parse_file, prefix = "TRS_dualseed")
cat("v2 SCDRS-seed files:  ", length(v2_files), "\n")
cat("v2 dual-seed files:   ", length(ds_files), "\n")

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
    cat(sprintf("  baseline      %-10s AUROC = %.3f\n", pt, a))
    results <- rbind(results, data.frame(
      method = "baseline_atac", target = tn,
      atac_w = NA_real_, auroc = a
    ))
  }

  # v1
  v1_file <- file.path(OUT_DIR, paste0("TRS_dual_", pt, ".tsv.gz"))
  if (file.exists(v1_file)) {
    v1 <- fread(v1_file)
    bc_col  <- intersect(c("barcode", "V1", "rn"),  colnames(v1))[1]
    trs_col <- intersect(c("TRS", "trs", "score"), colnames(v1))[1]
    if (is.na(bc_col) || is.na(trs_col)) {
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
      method = "v1_paired", target = tn,
      atac_w = NA_real_, auroc = a
    ))
  }

  # v2 (SCDRS-seed only)
  v2_for <- Filter(function(x) x$trait == pt, v2_meta)
  v2_for <- v2_for[order(sapply(v2_for, `[[`, "atac_w"))]
  for (vm in v2_for) {
    v2 <- fread(vm$file)
    pred_v2 <- v2$TRS
    names(pred_v2) <- v2$barcode
    common <- intersect(names(pred_v2), rownames(meta))
    lbl <- as.integer(as.character(meta[common, "cell_type"]) %in% ct_pos)
    a <- compute_auroc(pred_v2[common], lbl)
    cat(sprintf("  v2_scdrs_w%.1f %-10s AUROC = %.3f\n",
                vm$atac_w, pt, a))
    results <- rbind(results, data.frame(
      method = "v2_scdrs_seed", target = tn,
      atac_w = vm$atac_w, auroc = a
    ))
  }

  # v2 dual-seed (本脚本核心新增)
  ds_for <- Filter(function(x) x$trait == pt, ds_meta)
  ds_for <- ds_for[order(sapply(ds_for, `[[`, "atac_w"))]
  for (vm in ds_for) {
    ds <- fread(vm$file)
    pred_ds <- ds$TRS
    names(pred_ds) <- ds$barcode
    common <- intersect(names(pred_ds), rownames(meta))
    lbl <- as.integer(as.character(meta[common, "cell_type"]) %in% ct_pos)
    a <- compute_auroc(pred_ds[common], lbl)
    cat(sprintf("  v2_dual_w%.1f  %-10s AUROC = %.3f\n",
                vm$atac_w, pt, a))
    results <- rbind(results, data.frame(
      method = "v2_dual_seed", target = tn,
      atac_w = vm$atac_w, auroc = a
    ))
  }
}

# save
out_csv <- file.path(OUT_DIR, "eval_dualseed_auroc.csv")
write.csv(results, out_csv, row.names = FALSE)
cat("\nsaved: ", basename(out_csv), "\n")

# 4-way headline
cat("\n", strrep("=", 84), "\n", sep = "")
cat("4-WAY HEADLINE — baseline / v1 / v2 SCDRS-seed / v2 dual-seed (per target)\n")
cat(strrep("=", 84), "\n", sep = "")

main_targets <- c("MONO", "PLT", "LYM_TNK")
for (tn in main_targets) {
  sub <- results[results$target == tn, ]
  base   <- sub$auroc[sub$method == "baseline_atac"]
  v1     <- sub$auroc[sub$method == "v1_paired"]
  v2_07  <- sub$auroc[sub$method == "v2_scdrs_seed" &
                      abs(sub$atac_w - 0.7) < 1e-6]
  v2_10  <- sub$auroc[sub$method == "v2_scdrs_seed" &
                      abs(sub$atac_w - 1.0) < 1e-6]
  ds_07  <- sub$auroc[sub$method == "v2_dual_seed" &
                      abs(sub$atac_w - 0.7) < 1e-6]
  ds_10  <- sub$auroc[sub$method == "v2_dual_seed" &
                      abs(sub$atac_w - 1.0) < 1e-6]
  cat(sprintf("\n  [%s]\n", tn))
  cat(sprintf("    baseline           = %s\n",
              ifelse(length(base) > 0, sprintf("%.4f", base), "NA")))
  cat(sprintf("    v1_paired          = %s\n",
              ifelse(length(v1)   > 0, sprintf("%.4f", v1),   "NA")))
  cat(sprintf("    v2_scdrs  w=0.7    = %s\n",
              ifelse(length(v2_07) > 0, sprintf("%.4f", v2_07), "NA")))
  cat(sprintf("    v2_scdrs  w=1.0    = %s\n",
              ifelse(length(v2_10) > 0, sprintf("%.4f", v2_10), "NA")))
  cat(sprintf("    v2_dual   w=0.7    = %s",
              ifelse(length(ds_07) > 0, sprintf("%.4f", ds_07), "NA")))
  if (length(ds_07) > 0 && length(v2_07) > 0)
    cat(sprintf("  (Δ vs scdrs = %+.4f)", ds_07 - v2_07))
  cat("\n")
  cat(sprintf("    v2_dual   w=1.0    = %s",
              ifelse(length(ds_10) > 0, sprintf("%.4f", ds_10), "NA")))
  if (length(ds_10) > 0 && length(v2_10) > 0)
    cat(sprintf("  (Δ vs scdrs = %+.4f)", ds_10 - v2_10))
  cat("\n")
}

cat("\n06b_evaluate_auroc_dualseed.R OK.\n")
