# module4 / Phase 1.5.2: prepare module3 v2-compatible inputs
# -----------------------------------------------------------------------------
# 目的:
#   mapping + transfer 之后, 下游尽量直接复用 module3 的 v2 流水线。
#   module3 v2 的 KNN/RWR/TRS 脚本期望看到这些标准文件:
#
#     atac_z_block_v2.rds
#     rna_z_block_v2_MONO.rds
#     rna_z_block_v2_PLT.rds
#     rna_z_block_v2_LYM.rds
#     seed_v2_MONO.rds
#     seed_v2_PLT.rds
#     seed_v2_LYM.rds
#
#   本脚本只做“适配输入格式”, 不重写 module3 后续算法。
#
# 输入:
#   1. Granja ATAC LSI: blood_lsi.npy + atac_cells.csv
#   2. module4 transfer 输出:
#      - granja_pseudo_rna_profile_<trait>.rds
#      - granja_pseudo_zrna_<trait>.rds
#
# 输出:
#   03_pseudopair_pipeline/output/module3_v2_compat/
#     atac_z_block_v2.rds
#     rna_z_block_v2_<trait>.rds
#     seed_v2_<trait>.rds
#     fusion_summary.csv
#     seed_summary_v2.csv
#     module3_v2_compat_manifest.csv
#
# 关键原则:
#   - 不使用 cell_type。
#   - ATAC 是主空间: 使用 ATAC LSI 前 30 维, 并按列 z-score。
#   - RNA 是辅助空间: 使用 pseudo RNA profile 做 PCA, 再按列 z-score。
#   - seed 沿用 module3 v2 规则: pseudo_z_rna top 1%。
#
# 运行:
#   setwd("D:/RStudio/work/module4_traitaware_pseudopair")
#   source("03_pseudopair_pipeline/02_prepare_module3_v2_inputs.R")
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(irlba)
  library(matrixStats)
})

ROOT <- "D:/RStudio/work/module4_traitaware_pseudopair"
TRANSFER_DIR <- file.path(ROOT, "03_pseudopair_pipeline/output/transfer")
OUT_DIR <- file.path(ROOT, "03_pseudopair_pipeline/output/module3_v2_compat")
TMP_DIR <- file.path(OUT_DIR, "_tmp")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TMP_DIR, recursive = TRUE, showWarnings = FALSE)

ATAC_LSI_NPY <- "D:/SCAVENGE-reproducibility-main/blood_data/multiomics_module1_pairing/method_B_scSHEFT/scsheft_input/blood_lsi.npy"
ATAC_CELLS_CSV <- "D:/SCAVENGE-reproducibility-main/blood_data/multiomics_module1_pairing/method_B_scSHEFT/scsheft_input/atac_cells.csv"
ATAC_LSI_CSV <- file.path(TMP_DIR, "granja_atac_lsi_from_npy.csv")
PY_HELPER <- file.path(ROOT, "03_pseudopair_pipeline/export_lsi_npy.py")
if (!exists("PYTHON_EXE")) PYTHON_EXE <- "D:/python/python.exe"

if (!exists("TRAITS_TO_RUN")) TRAITS_TO_RUN <- c("MONO", "PLT", "LYM")
if (!exists("ATAC_N_DIMS")) ATAC_N_DIMS <- 30L
if (!exists("RNA_PCA_MIN_K")) RNA_PCA_MIN_K <- 10L
if (!exists("RNA_PCA_MAX_K")) RNA_PCA_MAX_K <- 30L
if (!exists("RNA_CUMVAR_TGT")) RNA_CUMVAR_TGT <- 0.80
if (!exists("SEED_TOP_PCT")) SEED_TOP_PCT <- 0.01
if (!exists("FORCE_EXPORT_LSI")) FORCE_EXPORT_LSI <- FALSE

stop_missing <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0L) {
    stop("Missing required file(s):\n  ", paste(missing, collapse = "\n  "),
         call. = FALSE)
  }
}

zscore_columns <- function(x) {
  x <- as.matrix(x)
  mu <- matrixStats::colMeans2(x)
  sd <- matrixStats::colSds(x)
  sd[sd == 0 | !is.finite(sd)] <- 1
  x <- sweep(x, 2L, mu, "-")
  x <- sweep(x, 2L, sd, "/")
  x
}

export_lsi_if_needed <- function() {
  stop_missing(c(PYTHON_EXE, PY_HELPER, ATAC_LSI_NPY, ATAC_CELLS_CSV))
  if (!FORCE_EXPORT_LSI && file.exists(ATAC_LSI_CSV)) {
    cat("[ATAC] Reusing exported LSI CSV: ", ATAC_LSI_CSV, "\n", sep = "")
    return(invisible(NULL))
  }
  cat("[ATAC] Exporting blood_lsi.npy to barcode-attached CSV via Python ...\n")
  args <- c(
    PY_HELPER,
    "--npy", ATAC_LSI_NPY,
    "--cells", ATAC_CELLS_CSV,
    "--out-csv", ATAC_LSI_CSV
  )
  status <- system2(PYTHON_EXE, args = args)
  if (!identical(status, 0L)) {
    stop("Python LSI export failed with status ", status, call. = FALSE)
  }
}

load_atac_z_block <- function(reference_cells) {
  export_lsi_if_needed()
  lsi_dt <- data.table::fread(ATAC_LSI_CSV)
  if (!"barcode" %in% colnames(lsi_dt)) {
    stop("Exported ATAC LSI CSV must contain barcode column.", call. = FALSE)
  }
  lsi_cols <- grep("^LSI_", colnames(lsi_dt), value = TRUE)
  if (length(lsi_cols) < ATAC_N_DIMS) {
    stop("ATAC LSI has only ", length(lsi_cols), " dimensions; need ",
         ATAC_N_DIMS, call. = FALSE)
  }
  lsi_dt <- as.data.frame(lsi_dt)
  rownames(lsi_dt) <- lsi_dt$barcode
  missing_cells <- setdiff(reference_cells, rownames(lsi_dt))
  if (length(missing_cells) > 0L) {
    stop("ATAC LSI is missing ", length(missing_cells),
         " reference cells. Example: ", paste(head(missing_cells), collapse = ", "),
         call. = FALSE)
  }
  atac_lsi <- as.matrix(lsi_dt[reference_cells, lsi_cols[seq_len(ATAC_N_DIMS)], drop = FALSE])
  atac_z <- zscore_columns(atac_lsi)
  rownames(atac_z) <- reference_cells
  colnames(atac_z) <- paste0("ATAC_LSI_", seq_len(ncol(atac_z)))
  atac_z
}

choose_k <- function(cumvar, min_k, max_k, target) {
  hit <- which(cumvar >= target)[1]
  if (is.na(hit)) hit <- max_k
  max(min_k, min(max_k, hit))
}

build_rna_z_block <- function(pseudo_rna_obj, reference_cells) {
  mat <- pseudo_rna_obj$matrix
  if (!inherits(mat, "Matrix")) mat <- as(mat, "CsparseMatrix")
  missing_cells <- setdiff(reference_cells, colnames(mat))
  if (length(missing_cells) > 0L) {
    stop("pseudo RNA matrix is missing ", length(missing_cells),
         " reference cells.", call. = FALSE)
  }
  mat <- mat[, reference_cells, drop = FALSE]  # gene x ATAC cell

  # module3 v2 的思路: RNA 端单独 PCA, 只保留前 10-30 个 PC。
  # 这里输入已经是 pseudo RNA profile, 所以不再重新 normalize。
  x <- Matrix::t(mat)  # ATAC cell x gene
  k_max <- min(RNA_PCA_MAX_K, ncol(x) - 1L, nrow(x) - 1L)
  if (k_max < RNA_PCA_MIN_K) {
    stop("Too few dimensions for RNA PCA: k_max = ", k_max, call. = FALSE)
  }
  cat("    irlba PCA on pseudo RNA: ", nrow(x), " cells x ", ncol(x),
      " genes; nv=", k_max, "\n", sep = "")
  pca <- irlba::irlba(x, nv = k_max, center = TRUE)
  scores <- pca$u %*% diag(pca$d, nrow = length(pca$d))

  # 与 module3 v2 保持一致: 在已计算出的前 max_k 个 PC 内算累计解释度。
  cumvar <- cumsum(pca$d^2) / sum(pca$d^2)
  k_rna <- choose_k(cumvar, RNA_PCA_MIN_K, RNA_PCA_MAX_K, RNA_CUMVAR_TGT)
  scores <- scores[, seq_len(k_rna), drop = FALSE]
  rna_z <- zscore_columns(scores)
  rownames(rna_z) <- reference_cells
  colnames(rna_z) <- paste0("RNA_PC_", seq_len(ncol(rna_z)))

  list(
    rna_z = rna_z,
    k_rna = k_rna,
    cumvar_at_k = cumvar[k_rna],
    cumvar = cumvar,
    singular_values = pca$d
  )
}

make_seed <- function(pseudo_z_obj, reference_cells, top_pct = 0.01) {
  dt <- as.data.frame(pseudo_z_obj$scores)
  if (!all(c("barcode", "pseudo_z_rna") %in% colnames(dt))) {
    stop("pseudo zRNA object must contain scores$barcode and scores$pseudo_z_rna.",
         call. = FALSE)
  }
  score <- dt$pseudo_z_rna[match(reference_cells, dt$barcode)]
  if (anyNA(score) || any(!is.finite(score))) {
    stop("pseudo_z_rna has missing values after barcode alignment.", call. = FALSE)
  }
  names(score) <- reference_cells
  thr <- as.numeric(quantile(score, 1 - top_pct, na.rm = TRUE))
  seed <- numeric(length(score))
  names(seed) <- reference_cells
  seed[score >= thr] <- 1
  list(seed = seed, threshold = thr, n_seed = sum(seed > 0), score = score)
}

cat(">>> Prepare module3 v2-compatible inputs for module4\n\n")
cat("Output directory: ", OUT_DIR, "\n\n", sep = "")

transfer_needed <- c(
  file.path(TRANSFER_DIR, paste0("granja_pseudo_rna_profile_", TRAITS_TO_RUN, ".rds")),
  file.path(TRANSFER_DIR, paste0("granja_pseudo_zrna_", TRAITS_TO_RUN, ".rds"))
)
stop_missing(transfer_needed)

cat("[1/4] Defining reference ATAC cells from transfer outputs ...\n")
first_obj <- readRDS(file.path(TRANSFER_DIR, paste0("granja_pseudo_rna_profile_", TRAITS_TO_RUN[1], ".rds")))
reference_cells <- first_obj$atac_barcodes
if (is.null(reference_cells)) reference_cells <- colnames(first_obj$matrix)
if (is.null(reference_cells)) {
  stop("Cannot find ATAC barcodes in pseudo RNA object.", call. = FALSE)
}
for (trait in TRAITS_TO_RUN) {
  obj <- readRDS(file.path(TRANSFER_DIR, paste0("granja_pseudo_rna_profile_", trait, ".rds")))
  cells <- obj$atac_barcodes
  if (is.null(cells)) cells <- colnames(obj$matrix)
  if (!identical(reference_cells, cells)) {
    stop("ATAC cell order differs across pseudo RNA profiles. Trait: ", trait,
         call. = FALSE)
  }
}
cat("  reference ATAC cells: ", length(reference_cells), "\n\n", sep = "")

cat("[2/4] Building module3-style atac_z_block_v2.rds ...\n")
atac_z <- load_atac_z_block(reference_cells)
saveRDS(atac_z, file.path(OUT_DIR, "atac_z_block_v2.rds"))
cat("  saved atac_z_block_v2.rds: ", paste(dim(atac_z), collapse = " x "), "\n\n", sep = "")

fusion_rows <- list()
seed_rows <- list()

cat("[3/4] Building rna_z_block_v2_<trait>.rds from pseudo RNA profiles ...\n")
for (trait in TRAITS_TO_RUN) {
  cat("\n=== Trait ", trait, " RNA block ===\n", sep = "")
  pseudo_rna_obj <- readRDS(file.path(TRANSFER_DIR, paste0("granja_pseudo_rna_profile_", trait, ".rds")))
  res <- build_rna_z_block(pseudo_rna_obj, reference_cells)
  saveRDS(res$rna_z, file.path(OUT_DIR, paste0("rna_z_block_v2_", trait, ".rds")))

  fusion_rows[[trait]] <- data.frame(
    trait = trait,
    n_cells = nrow(res$rna_z),
    n_g_trait_used = nrow(pseudo_rna_obj$matrix),
    K_rna_picked = res$k_rna,
    cumvar_at_K = res$cumvar_at_k,
    rna_z_dim = paste(dim(res$rna_z), collapse = " x "),
    source = "module4 pseudo RNA profile",
    no_celltype_used = TRUE,
    stringsAsFactors = FALSE
  )
  cat("  saved rna_z_block_v2_", trait, ".rds: ",
      paste(dim(res$rna_z), collapse = " x "),
      " | K=", res$k_rna,
      " | cumvar@K=", round(res$cumvar_at_k, 3), "\n", sep = "")
  rm(pseudo_rna_obj, res)
  invisible(gc(verbose = FALSE))
}

cat("\n[4/4] Building seed_v2_<trait>.rds from pseudo zRNA top 1% ...\n")
for (trait in TRAITS_TO_RUN) {
  cat("\n=== Trait ", trait, " seed ===\n", sep = "")
  pseudo_z_obj <- readRDS(file.path(TRANSFER_DIR, paste0("granja_pseudo_zrna_", trait, ".rds")))
  seed_res <- make_seed(pseudo_z_obj, reference_cells, top_pct = SEED_TOP_PCT)
  saveRDS(seed_res$seed, file.path(OUT_DIR, paste0("seed_v2_", trait, ".rds")))

  seed_rows[[trait]] <- data.frame(
    trait = trait,
    n_cells = length(seed_res$seed),
    n_seed = seed_res$n_seed,
    seed_top_pct = SEED_TOP_PCT,
    thr_pseudo_z_rna = seed_res$threshold,
    score_min = min(seed_res$score),
    score_median = median(seed_res$score),
    score_max = max(seed_res$score),
    source = "module4 pseudo zRNA",
    no_celltype_used = TRUE,
    stringsAsFactors = FALSE
  )
  cat("  saved seed_v2_", trait, ".rds",
      " | n_seed=", seed_res$n_seed,
      " | threshold=", round(seed_res$threshold, 3), "\n", sep = "")
  rm(pseudo_z_obj, seed_res)
  invisible(gc(verbose = FALSE))
}

fusion_summary <- do.call(rbind, fusion_rows)
seed_summary <- do.call(rbind, seed_rows)
write.csv(fusion_summary, file.path(OUT_DIR, "fusion_summary.csv"), row.names = FALSE)
write.csv(seed_summary, file.path(OUT_DIR, "seed_summary_v2.csv"), row.names = FALSE)

manifest <- data.frame(
  key = c(
    "atac_z_block",
    paste0("rna_z_block_", TRAITS_TO_RUN),
    paste0("seed_", TRAITS_TO_RUN)
  ),
  path = c(
    file.path(OUT_DIR, "atac_z_block_v2.rds"),
    file.path(OUT_DIR, paste0("rna_z_block_v2_", TRAITS_TO_RUN, ".rds")),
    file.path(OUT_DIR, paste0("seed_v2_", TRAITS_TO_RUN, ".rds"))
  ),
  stringsAsFactors = FALSE
)
write.csv(manifest, file.path(OUT_DIR, "module3_v2_compat_manifest.csv"), row.names = FALSE)

cat("\n=== module3 v2-compatible inputs are ready ===\n")
cat("\nFusion summary:\n")
print(fusion_summary, row.names = FALSE)
cat("\nSeed summary:\n")
print(seed_summary, row.names = FALSE)
cat("\nOutputs are in: ", OUT_DIR, "\n", sep = "")
