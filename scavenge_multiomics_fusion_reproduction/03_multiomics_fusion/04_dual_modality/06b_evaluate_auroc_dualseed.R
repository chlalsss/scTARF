# 06b_evaluate_auroc_dualseed.R
# §2.F' AUROC eval — 4-way 对比 (baseline / v1 / MultiOmicsFusion SCDRS-seed / MultiOmicsFusion dual-seed).
# fork 自 06_evaluate_auroc.R, 加扫描 TRS_dualseed_* 文件, 输出 eval_dualseed_auroc.csv.
# 不动 06_evaluate_auroc.R / eval_MultiOmicsFusion_auroc.csv, 保留阶段 2 主结果.

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

# 扫所有 MultiOmicsFusion TRS 文件 (SCDRS-seed)
MultiOmicsFusion_files <- list.files(OUT_DIR,
  pattern = "^TRS_MultiOmicsFusion_(MONO|PLT|LYM)_w[0-9.]+\\.tsv\\.gz$",
  full.names = TRUE)
parse_file <- function(f, prefix) {
  m <- regmatches(basename(f),
                  regexec(paste0("^", prefix, "_([A-Z]+)_w([0-9.]+)\\.tsv\\.gz$"),
                          basename(f)))[[1]]
  list(trait = m[2], atac_w = as.numeric(m[3]), file = f)
}
MultiOmicsFusion_meta  <- lapply(MultiOmicsFusion_files, parse_file, prefix = "TRS_MultiOmicsFusion")

# 扫所有 dual-seed TRS 文件
ds_files <- list.files(OUT_DIR,
  pattern = "^TRS_dualseed_(MONO|PLT|LYM)_w[0-9.]+\\.tsv\\.gz$",
  full.names = TRUE)
ds_meta  <- lapply(ds_files, parse_file, prefix = "TRS_dualseed")
cat("MultiOmicsFusion SCDRS-seed files:  ", length(MultiOmicsFusion_files), "\n")
cat("MultiOmicsFusion dual-seed files:   ", length(ds_files), "\n")

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

  # MultiOmicsFusion (SCDRS-seed only)
  MultiOmicsFusion_for <- Filter(function(x) x$trait == pt, MultiOmicsFusion_meta)
  MultiOmicsFusion_for <- MultiOmicsFusion_for[order(sapply(MultiOmicsFusion_for, `[[`, "atac_w"))]
  for (vm in MultiOmicsFusion_for) {
    MultiOmicsFusion <- fread(vm$file)
    pred_MultiOmicsFusion <- MultiOmicsFusion$TRS
    names(pred_MultiOmicsFusion) <- MultiOmicsFusion$barcode
    common <- intersect(names(pred_MultiOmicsFusion), rownames(meta))
    lbl <- as.integer(as.character(meta[common, "cell_type"]) %in% ct_pos)
    a <- compute_auroc(pred_MultiOmicsFusion[common], lbl)
    cat(sprintf("  MultiOmicsFusion_scdrs_w%.1f %-10s AUROC = %.3f\n",
                vm$atac_w, pt, a))
    results <- rbind(results, data.frame(
      method = "MultiOmicsFusion_scdrs_seed", target = tn,
      atac_w = vm$atac_w, auroc = a
    ))
  }

  # MultiOmicsFusion dual-seed (本脚本核心新增)
  ds_for <- Filter(function(x) x$trait == pt, ds_meta)
  ds_for <- ds_for[order(sapply(ds_for, `[[`, "atac_w"))]
  for (vm in ds_for) {
    ds <- fread(vm$file)
    pred_ds <- ds$TRS
    names(pred_ds) <- ds$barcode
    common <- intersect(names(pred_ds), rownames(meta))
    lbl <- as.integer(as.character(meta[common, "cell_type"]) %in% ct_pos)
    a <- compute_auroc(pred_ds[common], lbl)
    cat(sprintf("  MultiOmicsFusion_dual_w%.1f  %-10s AUROC = %.3f\n",
                vm$atac_w, pt, a))
    results <- rbind(results, data.frame(
      method = "MultiOmicsFusion_dual_seed", target = tn,
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
cat("4-WAY HEADLINE — baseline / v1 / MultiOmicsFusion SCDRS-seed / MultiOmicsFusion dual-seed (per target)\n")
cat(strrep("=", 84), "\n", sep = "")

main_targets <- c("MONO", "PLT", "LYM_TNK")
for (tn in main_targets) {
  sub <- results[results$target == tn, ]
  base   <- sub$auroc[sub$method == "baseline_atac"]
  v1     <- sub$auroc[sub$method == "v1_paired"]
  MultiOmicsFusion_07  <- sub$auroc[sub$method == "MultiOmicsFusion_scdrs_seed" &
                      abs(sub$atac_w - 0.7) < 1e-6]
  MultiOmicsFusion_10  <- sub$auroc[sub$method == "MultiOmicsFusion_scdrs_seed" &
                      abs(sub$atac_w - 1.0) < 1e-6]
  ds_07  <- sub$auroc[sub$method == "MultiOmicsFusion_dual_seed" &
                      abs(sub$atac_w - 0.7) < 1e-6]
  ds_10  <- sub$auroc[sub$method == "MultiOmicsFusion_dual_seed" &
                      abs(sub$atac_w - 1.0) < 1e-6]
  cat(sprintf("\n  [%s]\n", tn))
  cat(sprintf("    baseline           = %s\n",
              ifelse(length(base) > 0, sprintf("%.4f", base), "NA")))
  cat(sprintf("    v1_paired          = %s\n",
              ifelse(length(v1)   > 0, sprintf("%.4f", v1),   "NA")))
  cat(sprintf("    MultiOmicsFusion_scdrs  w=0.7    = %s\n",
              ifelse(length(MultiOmicsFusion_07) > 0, sprintf("%.4f", MultiOmicsFusion_07), "NA")))
  cat(sprintf("    MultiOmicsFusion_scdrs  w=1.0    = %s\n",
              ifelse(length(MultiOmicsFusion_10) > 0, sprintf("%.4f", MultiOmicsFusion_10), "NA")))
  cat(sprintf("    MultiOmicsFusion_dual   w=0.7    = %s",
              ifelse(length(ds_07) > 0, sprintf("%.4f", ds_07), "NA")))
  if (length(ds_07) > 0 && length(MultiOmicsFusion_07) > 0)
    cat(sprintf("  (Δ vs scdrs = %+.4f)", ds_07 - MultiOmicsFusion_07))
  cat("\n")
  cat(sprintf("    MultiOmicsFusion_dual   w=1.0    = %s",
              ifelse(length(ds_10) > 0, sprintf("%.4f", ds_10), "NA")))
  if (length(ds_10) > 0 && length(MultiOmicsFusion_10) > 0)
    cat(sprintf("  (Δ vs scdrs = %+.4f)", ds_10 - MultiOmicsFusion_10))
  cat("\n")
}

cat("\n06b_evaluate_auroc_dualseed.R OK.\n")
