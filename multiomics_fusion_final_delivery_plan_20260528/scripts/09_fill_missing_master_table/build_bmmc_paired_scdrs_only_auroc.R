# Compute paired BMMC scDRS-only AUROC from existing scDRS score files.
#
# This script only uses paired BMMC inputs. It does not run or call any
# non-paired module4 workflow.

WORK_ROOT <- "D:/RStudio/work"
DELIVERY_ROOT <- file.path(WORK_ROOT, "multiomics_fusion_final_delivery_plan_20260528")
OUT_DIR <- file.path(DELIVERY_ROOT, "final_outputs/09_missing_master_table")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

if (!exists("BMMC_META_FILE")) {
  BMMC_META_FILE <- "D:/SCAVENGE-reproducibility-main/BMMC/01_data_qc/BMMC_cell_meta.rds"
}
if (!exists("BMMC_SCDRS_SCORE_DIRS")) {
  BMMC_SCDRS_SCORE_DIRS <- c(
    file.path(WORK_ROOT, "module3_multiome/04_dual_modality/01_scdrs_python/output"),
    "D:/SCAVENGE-reproducibility-main/BMMC/04_dual_modality/01_scdrs_python/output"
  )
}
if (!exists("OUT_CSV")) {
  OUT_CSV <- file.path(OUT_DIR, "bmmc_paired_scdrs_only_auroc.csv")
}

TRAIT_SCORE_FILES <- c(
  MONO = "UKB_460K.blood_MONOCYTE_COUNT.score.gz",
  PLT = "UKB_460K.blood_PLATELET_COUNT.score.gz",
  LYM_wide = "UKB_460K.blood_LYMPHOCYTE_COUNT.score.gz"
)

TARGETS <- list(
  MONO = c("CD14+ Mono", "CD16+ Mono"),
  PLT = c("MK/E prog"),
  LYM_wide = c(
    "NK", "ILC", "CD4+ T activated", "CD4+ T naive",
    "CD8+ T", "CD8+ T naive", "Lymph prog",
    "B1 B", "Naive CD20+ B", "Transitional B", "Plasma cell"
  )
)

stop_missing <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0L) {
    stop("Missing required file(s):\n  ", paste(missing, collapse = "\n  "),
         call. = FALSE)
  }
}

first_existing_dir <- function(paths) {
  hit <- paths[dir.exists(paths)]
  if (length(hit) == 0L) {
    stop("No existing scDRS score directory found among:\n  ",
         paste(paths, collapse = "\n  "), call. = FALSE)
  }
  hit[1]
}

load_meta <- function(path) {
  x <- readRDS(path)
  if (inherits(x, "data.table")) x <- as.data.frame(x)
  if (is.null(rownames(x)) ||
      identical(rownames(x), as.character(seq_len(nrow(x))))) {
    if ("barcode" %in% names(x)) rownames(x) <- as.character(x$barcode)
  }
  if (!"cell_type" %in% names(x)) {
    stop("BMMC metadata must contain cell_type: ", path, call. = FALSE)
  }
  x
}

score_col <- function(x) {
  hit <- intersect(c("norm_score", "zscore", "raw_score", "score"), names(x))[1]
  if (is.na(hit)) stop("No score column found in: ", paste(names(x), collapse = ", "))
  hit
}

cell_col <- function(x) {
  hit <- intersect(c("barcode", "cell", "cell_id", "V1"), names(x))[1]
  if (is.na(hit)) names(x)[1] else hit
}

compute_auc <- function(score, label) {
  keep <- is.finite(score) & !is.na(label)
  score <- score[keep]
  label <- as.integer(label[keep])
  if (length(score) < 20L || sum(label == 1L) < 10L || sum(label == 0L) < 10L) {
    return(NA_real_)
  }
  if (requireNamespace("pROC", quietly = TRUE)) {
    return(as.numeric(pROC::auc(pROC::roc(label, score, quiet = TRUE,
                                           direction = "<",
                                           levels = c(0L, 1L)))))
  }
  r <- rank(score, ties.method = "average")
  n_pos <- sum(label == 1L)
  n_neg <- sum(label == 0L)
  (sum(r[label == 1L]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

cat(">>> BMMC paired scDRS-only AUROC\n\n")
score_dir <- first_existing_dir(BMMC_SCDRS_SCORE_DIRS)
stop_missing(c(BMMC_META_FILE, file.path(score_dir, TRAIT_SCORE_FILES)))

meta <- load_meta(BMMC_META_FILE)
rows <- list()
for (trait in names(TRAIT_SCORE_FILES)) {
  path <- file.path(score_dir, TRAIT_SCORE_FILES[[trait]])
  cat("[trait] ", trait, " <- ", path, "\n", sep = "")
  sc <- utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  cc <- cell_col(sc)
  vc <- score_col(sc)
  score <- as.numeric(sc[[vc]])
  names(score) <- as.character(sc[[cc]])

  common <- intersect(names(score), rownames(meta))
  label <- as.integer(as.character(meta[common, "cell_type"]) %in% TARGETS[[trait]])
  auc <- compute_auc(score[common], label)

  rows[[trait]] <- data.frame(
    dataset = "BMMC paired multiome",
    pairing_setting = "paired",
    trait = trait,
    target_definition = trait,
    method_group = "scdrs",
    method_version = "scdrs_only",
    method = "scDRS-only",
    rna_reference = "paired BMMC RNA",
    top_k = NA_real_,
    feature_mode = "rna_score_only",
    atac_w = NA_real_,
    metric = "AUROC",
    value = auc,
    ci_low = NA_real_,
    ci_high = NA_real_,
    n_positive = sum(label == 1L),
    n_total = length(label),
    source_file = path,
    note = "Computed from paired BMMC scDRS score.gz and BMMC cell_type metadata",
    stringsAsFactors = FALSE
  )
}

out <- do.call(rbind, rows)
utils::write.csv(out, OUT_CSV, row.names = FALSE, quote = TRUE)
cat("\nWrote: ", OUT_CSV, "\n", sep = "")
