# module4 / 04_evaluation / eval_auroc.R
# ------------------------------------------------------------------
# 通用 AUROC 评估脚本. 算法核心继承自 module3 06_evaluate_auroc.R,
# 接口参数化让同一份脚本能跑 BMMC / Granja / 任意数据集.
#
# 输入预期格式 (用 00_granja_prepare_inputs.R 这种预处理脚本产出):
#   META_FILE     data.frame, rowname = barcode, 含 cell_type 列
#   BASELINE_TRS  data.frame, rowname = barcode, 列 = trait names (MONO/PLT/LYM/...)
#   V2_FILES_DIR  目录, 含 TRS_v2_{TRAIT}_w{atac_w}.tsv.gz 形式的 module4 流水线输出
#                 (空目录或不存在, 仅 baseline 评估时正常使用)
#
# 输出:
#   <OUT_DIR>/eval_auroc.csv  长表 (method, target, atac_w, auroc[, auroc_ci_lo, ci_hi])
#
# 用法 1 (直接 source 跑 Granja baseline):
#   配置顶部常量后 source(__FILE__).
#
# 用法 2 (做 library 用):
#   source 之后 evaluate_auroc(...) 函数可独立调用.
# ------------------------------------------------------------------

if (!requireNamespace("pROC", quietly = TRUE)) {
  stop("Missing pROC. In a separate R session run: install.packages('pROC')")
}

suppressPackageStartupMessages({
  library(data.table)
  library(pROC)
})

# ============================================================
# CONFIG (强制赋值, 想覆盖请直接编辑本节, 不要在 session 里 assign)
# 之前用 if (!exists()) 模式会被 session 里同名变量污染 (踩过 module3
# 那边的 META_FILE / BASELINE_TRS 遗留), 改为强制覆盖更稳.
# ============================================================

META_FILE    <- "D:/RStudio/work/module4_traitaware_pseudopair/04_evaluation/inputs/granja_cell_meta.rds"
BASELINE_TRS <- "D:/RStudio/work/module4_traitaware_pseudopair/04_evaluation/inputs/granja_baseline_TRS_matrix.rds"
V2_FILES_DIR <- "D:/RStudio/work/module4_traitaware_pseudopair/03_pseudopair_pipeline/output/v3_top5_significant_p005/module3_trait_screening_ge0p80"
OUT_DIR <- "D:/RStudio/work/module4_traitaware_pseudopair/04_evaluation/v3_top5_significant_p005/trait_screening_ge0p80"
OUT_CSV_NAME <- "v3_top5_significant_p005_trait_screening_ge0p80_auroc.csv"

# Target cell_type definitions for Granja 体系 (26 类, 命名 "<NN>_<lineage>")
# 由 00_granja_prepare_inputs.R [3/4] 实测确认. 跟 WORKFLOW §B.5 比有两处修正:
#   - 不存在 "13_CD16.Mono" (13_/14_/26_ 是 Unk), MONO 只含 CD14.Mono.1/2
#   - 不存在显式 MK 类, PLT 用红系 + 嗜碱性 (跟 module3 BMMC PLT lineage 一致)
TARGETS <- list(
  RDW = c("02_Early.Eryth", "03_Late.Eryth", "04_Early.Baso"),
  RBC = c("02_Early.Eryth", "03_Late.Eryth", "04_Early.Baso")
)
PRED_TRAIT <- list(RDW = "RDW", RBC = "RBC")

BOOTSTRAP_CI <- TRUE
N_BOOT       <- 200L
SEED         <- 2026L

# ============================================================
# CORE FUNCTIONS
# ============================================================

# load + normalize meta: data.frame with barcode rownames + cell_type col
load_meta <- function(meta_file) {
  m <- readRDS(meta_file)
  if (inherits(m, "data.table")) m <- as.data.frame(m)
  if (is.null(rownames(m)) ||
      identical(rownames(m), as.character(seq_len(nrow(m))))) {
    if ("barcode" %in% colnames(m)) {
      rownames(m) <- m$barcode
    } else {
      stop("meta has no rownames and no barcode column")
    }
  }
  if (!"cell_type" %in% colnames(m)) {
    stop("meta missing cell_type column. Available: ",
         paste(colnames(m), collapse = ", "))
  }
  m
}

# load + normalize baseline TRS: data.frame with barcode rownames
load_baseline_trs <- function(trs_file) {
  if (!file.exists(trs_file)) return(NULL)
  x <- readRDS(trs_file)
  if (inherits(x, "data.table")) x <- as.data.frame(x)
  bc_col_cands <- c("cell", "cell_id", "barcode")
  bc_col       <- intersect(bc_col_cands, colnames(x))
  if (length(bc_col) > 0) {
    rownames(x) <- as.character(x[[bc_col[1]]])
    x[[bc_col[1]]] <- NULL
  }
  x
}

# core AUROC + optional bootstrap CI
compute_auroc <- function(predictor, labels,
                          bootstrap = FALSE, n_boot = 200L, seed = 2026L) {
  # drop NA in either predictor or labels
  keep <- is.finite(predictor) & !is.na(labels)
  predictor <- predictor[keep]
  labels    <- labels[keep]
  if (length(predictor) < 20L) {
    return(list(auroc = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_,
                n_used = length(predictor)))
  }
  if (sum(labels == 1L) < 10L || sum(labels == 0L) < 10L) {
    return(list(auroc = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_,
                n_used = length(predictor)))
  }
  v <- stats::var(predictor)
  if (!is.finite(v) || v < 1e-12) {
    return(list(auroc = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_,
                n_used = length(predictor)))
  }
  roc_obj <- suppressMessages(pROC::roc(labels, predictor,
                                         quiet = TRUE,
                                         direction = "<",
                                         levels = c(0L, 1L)))
  a <- as.numeric(pROC::auc(roc_obj))
  if (!bootstrap) {
    return(list(auroc = a, ci_lo = NA_real_, ci_hi = NA_real_,
                n_used = length(predictor)))
  }
  set.seed(seed)
  ci_obj <- pROC::ci.auc(roc_obj, method = "bootstrap",
                         boot.n = n_boot)
  list(auroc = a,
       ci_lo  = as.numeric(ci_obj[1]),
       ci_hi  = as.numeric(ci_obj[3]),
       n_used = length(predictor))
}

# parse module4 流水线 TRS 文件名:  TRS_v2_<TRAIT>_w<atac_w>.tsv.gz
parse_v2_filename <- function(f) {
  bn <- basename(f)
  m  <- regmatches(bn, regexec(
    "^TRS_v2_([A-Z_]+)_w([0-9.]+)\\.tsv\\.gz$", bn))[[1]]
  if (length(m) < 3) return(NULL)
  list(trait = m[2], atac_w = as.numeric(m[3]), file = f)
}

# main evaluator
evaluate_auroc <- function(meta, baseline_mat = NULL, v2_files = NULL,
                            targets = TARGETS, pred_trait = PRED_TRAIT,
                            bootstrap = BOOTSTRAP_CI,
                            n_boot = N_BOOT, seed = SEED) {
  results <- data.frame()
  v2_meta <- if (is.null(v2_files)) list() else Filter(Negate(is.null),
                                                        lapply(v2_files, parse_v2_filename))

  for (tn in names(targets)) {
    ct_pos <- targets[[tn]]
    pt     <- pred_trait[[tn]]
    cat("\n=== target ", tn, " (positive: ",
        paste(ct_pos, collapse = "/"), ") ===\n", sep = "")

    # baseline
    if (!is.null(baseline_mat) && pt %in% colnames(baseline_mat)) {
      bc     <- rownames(baseline_mat)
      common <- intersect(bc, rownames(meta))
      pred   <- baseline_mat[common, pt]
      lbl    <- as.integer(as.character(meta[common, "cell_type"]) %in% ct_pos)
      r      <- compute_auroc(pred, lbl, bootstrap, n_boot, seed)
      cat(sprintf("  baseline_atac  %-10s AUROC = %.4f%s  (pos=%d / total=%d)\n",
                  pt, r$auroc,
                  if (bootstrap) sprintf(" [%.4f, %.4f]", r$ci_lo, r$ci_hi) else "",
                  sum(lbl), length(lbl)))
      results <- rbind(results, data.frame(
        method = "baseline_atac", target = tn, atac_w = NA_real_,
        auroc = r$auroc, auroc_ci_lo = r$ci_lo, auroc_ci_hi = r$ci_hi,
        n_positive = sum(lbl), n_total = length(lbl),
        stringsAsFactors = FALSE
      ))
    } else if (!is.null(baseline_mat)) {
      cat("  baseline_atac  -- trait '", pt, "' not in baseline cols\n", sep = "")
    }

    # v2 sweep (module4 主测 / 改进路线)
    v2_for <- Filter(function(x) x$trait == pt, v2_meta)
    v2_for <- v2_for[order(sapply(v2_for, `[[`, "atac_w"))]
    for (vm in v2_for) {
      v2 <- fread(vm$file)
      bc_col  <- intersect(c("barcode", "cell", "cell_id", "V1", "rn"),
                            colnames(v2))[1]
      trs_col <- intersect(c("TRS", "trs", "score"), colnames(v2))[1]
      if (is.na(bc_col) || is.na(trs_col)) {
        cat("  v2 file ", basename(vm$file),
            ": unrecognized cols, assume col1=barcode col2=TRS\n", sep = "")
        bc_col  <- colnames(v2)[1]; trs_col <- colnames(v2)[2]
      }
      pred_v2 <- v2[[trs_col]]
      names(pred_v2) <- v2[[bc_col]]
      common <- intersect(names(pred_v2), rownames(meta))
      lbl    <- as.integer(as.character(meta[common, "cell_type"]) %in% ct_pos)
      r      <- compute_auroc(pred_v2[common], lbl, bootstrap, n_boot, seed)
      cat(sprintf("  v2 atac_w=%.2f %-10s AUROC = %.4f%s\n",
                  vm$atac_w, pt, r$auroc,
                  if (bootstrap) sprintf(" [%.4f, %.4f]", r$ci_lo, r$ci_hi) else ""))
      results <- rbind(results, data.frame(
        method = "v2", target = tn, atac_w = vm$atac_w,
        auroc = r$auroc, auroc_ci_lo = r$ci_lo, auroc_ci_hi = r$ci_hi,
        n_positive = sum(lbl), n_total = length(lbl),
        stringsAsFactors = FALSE
      ))
    }
  }
  results
}

# ============================================================
# MAIN  (直接 source 跑, 用 CONFIG 段顶部常量)
# ============================================================

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat(">>> eval_auroc.R running\n")
cat("    META_FILE   : ", META_FILE,    "\n", sep = "")
cat("    BASELINE_TRS: ", BASELINE_TRS, "\n", sep = "")
cat("    V2_FILES_DIR: ", V2_FILES_DIR, "\n", sep = "")
cat("    OUT_DIR     : ", OUT_DIR,      "\n", sep = "")
cat("    BOOTSTRAP_CI: ", BOOTSTRAP_CI, " (n_boot = ", N_BOOT,
    ", seed = ", SEED, ")\n\n", sep = "")

# load
meta <- load_meta(META_FILE)
cat("loaded meta: ", paste(dim(meta), collapse = " x "),
    "; cell_type col present\n", sep = "")

baseline_mat <- load_baseline_trs(BASELINE_TRS)
if (!is.null(baseline_mat)) {
  cat("loaded baseline TRS: ",
      paste(dim(baseline_mat), collapse = " x "),
      "; cols (first 8): ",
      paste(head(colnames(baseline_mat), 8), collapse = ", "), "\n", sep = "")
} else {
  cat("baseline TRS file not found, skipping baseline.\n")
}

# scan v2 files (optional)
v2_files <- character(0)
if (dir.exists(V2_FILES_DIR)) {
  v2_files <- list.files(V2_FILES_DIR,
    pattern = "^TRS_v2_[A-Z_]+_w[0-9.]+\\.tsv\\.gz$",
    full.names = TRUE)
}
cat("scanned v2 TRS files: ", length(v2_files), "\n", sep = "")

# run
results <- evaluate_auroc(meta, baseline_mat, v2_files,
                           targets = TARGETS, pred_trait = PRED_TRAIT,
                           bootstrap = BOOTSTRAP_CI,
                           n_boot = N_BOOT, seed = SEED)

# save
out_csv <- file.path(OUT_DIR, OUT_CSV_NAME)
write.csv(results, out_csv, row.names = FALSE)
cat("\nsaved: ", out_csv, "\n", sep = "")

# headline
cat("\n", strrep("=", 78), "\n", sep = "")
cat("HEADLINE: AUROC per target\n")
cat(strrep("=", 78), "\n", sep = "")
for (tn in names(TARGETS)) {
  sub <- results[results$target == tn, ]
  if (nrow(sub) == 0) next
  base <- sub$auroc[sub$method == "baseline_atac"]
  base_ci_lo <- sub$auroc_ci_lo[sub$method == "baseline_atac"]
  base_ci_hi <- sub$auroc_ci_hi[sub$method == "baseline_atac"]
  v2s  <- sub[sub$method == "v2", ]
  base_str <- if (length(base) > 0 && !is.na(base))
    sprintf("%.4f [%.4f, %.4f]", base, base_ci_lo, base_ci_hi) else "  NA"
  v2_str <- if (nrow(v2s) > 0) {
    best <- v2s[which.max(v2s$auroc), ]
    sprintf("v2_best=%.4f (atac_w=%.2f)", best$auroc, best$atac_w)
  } else "v2_best=NA"
  cat(sprintf("  %-9s  baseline=%s  %s\n", tn, base_str, v2_str))
}

cat("\neval_auroc.R OK.\n")
