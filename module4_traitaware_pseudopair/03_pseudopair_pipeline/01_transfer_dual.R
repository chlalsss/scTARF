# module4 / Phase 1.5.1: dual transfer from CITE RNA to Granja ATAC
# -----------------------------------------------------------------------------
# 目的:
#   基于 02_trait_mapping/01_align_granja_cite_by_trait.R 得到的 top-k RNA
#   邻居和权重, 把 CITE RNA 信息转移到每个 Granja ATAC 细胞上。
#
# 输出两路信息:
#   1. pseudo RNA profile:
#      每个 ATAC 细胞 = 它的 top-k CITE RNA 邻居 log-normalized expression
#      的加权平均。后面用于 ATAC+RNA 融合向量。
#
#   2. pseudo z_rna:
#      每个 ATAC 细胞 = 它的 top-k CITE RNA 邻居 scDRS norm_score 的
#      加权平均。后面用于 seed / RWR / TRS。
#
# 关键原则:
#   - 本脚本不做新的 cell 匹配, 只使用已经完成的 mapping。
#   - 本脚本不读取、不使用任何 cell_type 标签。
#   - ATAC 细胞仍然是主对象; RNA 只是为 ATAC 细胞补充辅助信息。
#
# 运行方式:
#   在 RStudio 里运行:
#     setwd("D:/RStudio/work/module4_traitaware_pseudopair")
#     source("03_pseudopair_pipeline/01_transfer_dual.R")
#
# 可选:
#   只跑一个 trait:
#     TRAITS_TO_RUN <- c("MONO")
#     source("03_pseudopair_pipeline/01_transfer_dual.R")
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
})

ROOT <- "D:/RStudio/work/module4_traitaware_pseudopair"
if (!exists("TRANSFER_MAP_DIR")) TRANSFER_MAP_DIR <- file.path(ROOT, "02_trait_mapping/output")
MAP_DIR <- TRANSFER_MAP_DIR
CITE_EXTRACT_DIR <- file.path(MAP_DIR, "cite_extract")
if (!exists("TRANSFER_OUT_DIR")) TRANSFER_OUT_DIR <- file.path(ROOT, "03_pseudopair_pipeline/output/transfer")
OUT_DIR <- TRANSFER_OUT_DIR
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

if (!exists("CITE_MTX")) CITE_MTX <- file.path(CITE_EXTRACT_DIR, "cite_trait_counts.mtx")
if (!exists("CITE_GENES")) CITE_GENES <- file.path(CITE_EXTRACT_DIR, "cite_trait_counts_genes.txt")
if (!exists("CITE_CELLS")) CITE_CELLS <- file.path(CITE_EXTRACT_DIR, "cite_trait_counts_cells.txt")
if (!exists("CITE_SCORE_RDS")) {
  CITE_SCORE_RDS <- file.path(ROOT, "01_data_prep_external_rna/output/cite_rna_scdrs_scores.rds")
}

if (!exists("TRAITS_ALL")) TRAITS_ALL <- c("MONO", "PLT", "LYM")
if (!exists("TRAITS_TO_RUN")) TRAITS_TO_RUN <- TRAITS_ALL
if (!exists("CPM_SCALE")) CPM_SCALE <- 1e4

stop_missing <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0L) {
    stop("Missing required file(s):\n  ", paste(missing, collapse = "\n  "),
         call. = FALSE)
  }
}

make_weight_matrix <- function(mapping, cite_cells) {
  n_atac <- length(mapping$atac_barcodes)
  k <- mapping$k
  neighbor_idx <- match(as.vector(mapping$rna_neighbor_barcode), cite_cells)
  if (anyNA(neighbor_idx)) {
    stop("Some RNA neighbor barcodes are not present in extracted CITE cell list.",
         call. = FALSE)
  }
  sparseMatrix(
    i = as.integer(neighbor_idx),
    j = rep(seq_len(n_atac), times = k),
    x = as.numeric(mapping$neighbor_weight),
    dims = c(length(cite_cells), n_atac),
    dimnames = list(cite_cells, mapping$atac_barcodes)
  )
}

log_normalize_counts <- function(counts_g_by_c, lib_size, scale_factor = 1e4) {
  lib <- lib_size[colnames(counts_g_by_c)]
  if (anyNA(lib) || any(!is.finite(lib))) {
    stop("Library size vector is missing values for some CITE cells.", call. = FALSE)
  }
  lib[lib < 1] <- 1
  norm <- Matrix::t(Matrix::t(counts_g_by_c) / lib) * scale_factor
  log1p(norm)
}

vector_stats <- function(x) {
  c(
    min = min(x, na.rm = TRUE),
    q05 = as.numeric(quantile(x, 0.05, na.rm = TRUE)),
    median = median(x, na.rm = TRUE),
    mean = mean(x, na.rm = TRUE),
    q95 = as.numeric(quantile(x, 0.95, na.rm = TRUE)),
    max = max(x, na.rm = TRUE)
  )
}

cat(">>> Phase 1.5.1 dual transfer: CITE RNA -> Granja ATAC\n\n")
cat("Traits: ", paste(TRAITS_TO_RUN, collapse = ", "), "\n\n", sep = "")

stop_missing(c(CITE_MTX, CITE_GENES, CITE_CELLS, CITE_SCORE_RDS))

cat("[1/4] Loading extracted CITE RNA count matrix ...\n")
cite_counts <- Matrix::readMM(CITE_MTX)
cite_counts <- as(cite_counts, "CsparseMatrix")
cite_genes <- readLines(CITE_GENES, warn = FALSE)
cite_cells <- readLines(CITE_CELLS, warn = FALSE)
rownames(cite_counts) <- cite_genes
colnames(cite_counts) <- cite_cells
cat("  CITE counts: ", paste(dim(cite_counts), collapse = " x "),
    " (gene x RNA cell)\n", sep = "")

cat("[2/4] Loading CITE scDRS scores and library sizes ...\n")
cite_scores <- readRDS(CITE_SCORE_RDS)
cite_scores <- as.data.frame(cite_scores)

# 注意: 这里刻意只使用 cell、scdrs_input_total_counts 和 *_norm_score。
# cell_type 列即使存在, 本脚本也不读取、不参与任何计算。
needed_score_cols <- paste0(TRAITS_TO_RUN, "_norm_score")
needed_cols <- c("cell", "scdrs_input_total_counts", needed_score_cols)
missing_cols <- setdiff(needed_cols, colnames(cite_scores))
if (length(missing_cols) > 0L) {
  stop("Missing required column(s) in CITE score RDS: ",
       paste(missing_cols, collapse = ", "), call. = FALSE)
}
cite_scores <- cite_scores[, needed_cols, drop = FALSE]

lib_size <- cite_scores$scdrs_input_total_counts[match(cite_cells, cite_scores$cell)]
names(lib_size) <- cite_cells
if (anyNA(lib_size) || any(!is.finite(lib_size))) {
  stop("Missing CITE library sizes. Check cite_rna_scdrs_scores.rds.", call. = FALSE)
}
lib_size[lib_size < 1] <- 1

summary_rows <- list()

cat("[3/4] Running transfer per trait ...\n")
for (trait in TRAITS_TO_RUN) {
  cat("\n=== Trait ", trait, " ===\n", sep = "")
  map_file <- file.path(MAP_DIR, paste0("granja_cite_mapping_", trait, ".rds"))
  stop_missing(map_file)
  mapping <- readRDS(map_file)

  genes <- mapping$g_trait_used
  missing_genes <- setdiff(genes, rownames(cite_counts))
  if (length(missing_genes) > 0L) {
    stop("CITE extracted matrix is missing mapping genes for ", trait, ": ",
         paste(head(missing_genes, 10), collapse = ", "), call. = FALSE)
  }

  cat("  building RNA-neighbor weight matrix ...\n")
  W <- make_weight_matrix(mapping, cite_cells)
  w_colsum <- Matrix::colSums(W)
  if (!all(abs(w_colsum - 1) < 1e-6)) {
    stop("Neighbor weights do not sum to 1 for all ATAC cells.", call. = FALSE)
  }

  cat("  log-normalizing CITE RNA counts on G_trait ...\n")
  rna_log <- log_normalize_counts(cite_counts[genes, , drop = FALSE],
                                  lib_size = lib_size,
                                  scale_factor = CPM_SCALE)

  cat("  transferring pseudo RNA profile: gene x ATAC ...\n")
  pseudo_rna <- rna_log %*% W
  pseudo_rna <- as(pseudo_rna, "CsparseMatrix")
  rownames(pseudo_rna) <- genes
  colnames(pseudo_rna) <- mapping$atac_barcodes

  cat("  transferring pseudo z_rna score ...\n")
  score_col <- paste0(trait, "_norm_score")
  score <- cite_scores[[score_col]][match(cite_cells, cite_scores$cell)]
  if (anyNA(score) || any(!is.finite(score))) {
    replacement <- median(score[is.finite(score)], na.rm = TRUE)
    warning("Some ", score_col, " values are missing; replacing with median = ",
            signif(replacement, 4))
    score[!is.finite(score) | is.na(score)] <- replacement
  }
  pseudo_z <- as.numeric(Matrix::crossprod(W, score))
  names(pseudo_z) <- mapping$atac_barcodes

  max_weight <- apply(mapping$neighbor_weight, 1, max)
  entropy <- -rowSums(mapping$neighbor_weight * log(pmax(mapping$neighbor_weight, 1e-12)))
  pseudo_z_dt <- data.table(
    barcode = mapping$atac_barcodes,
    trait = trait,
    pseudo_z_rna = pseudo_z,
    top1_distance = mapping$neighbor_distance[, 1],
    top1_cos_dist = mapping$neighbor_distance[, 1]^2 / 2,
    max_neighbor_weight = max_weight,
    neighbor_weight_entropy = entropy
  )

  pseudo_rna_obj <- list(
    trait = trait,
    trait_id = mapping$trait_id,
    matrix = pseudo_rna,
    orientation = "gene x ATAC_cell",
    genes = genes,
    atac_barcodes = mapping$atac_barcodes,
    source = "weighted average of top-k CITE RNA log-normalized expression",
    cpm_scale = CPM_SCALE,
    mapping_file = map_file,
    no_celltype_used = TRUE,
    created_at = as.character(Sys.time())
  )
  pseudo_z_obj <- list(
    trait = trait,
    trait_id = mapping$trait_id,
    scores = pseudo_z_dt,
    source = "weighted average of top-k CITE RNA scDRS norm_score",
    mapping_file = map_file,
    no_celltype_used = TRUE,
    created_at = as.character(Sys.time())
  )

  rna_out <- file.path(OUT_DIR, paste0("granja_pseudo_rna_profile_", trait, ".rds"))
  z_out <- file.path(OUT_DIR, paste0("granja_pseudo_zrna_", trait, ".rds"))
  z_csv <- file.path(OUT_DIR, paste0("granja_pseudo_zrna_", trait, ".csv"))
  saveRDS(pseudo_rna_obj, rna_out)
  saveRDS(pseudo_z_obj, z_out)
  fwrite(pseudo_z_dt, z_csv)

  z_stats <- vector_stats(pseudo_z)
  summary_rows[[trait]] <- data.frame(
    trait = trait,
    trait_id = mapping$trait_id,
    n_atac_cells = length(mapping$atac_barcodes),
    n_cite_cells = length(cite_cells),
    n_g_trait_used = length(genes),
    k = mapping$k,
    pseudo_rna_nnz = length(pseudo_rna@x),
    pseudo_rna_density = length(pseudo_rna@x) / max(prod(dim(pseudo_rna)), 1),
    pseudo_z_min = z_stats[["min"]],
    pseudo_z_q05 = z_stats[["q05"]],
    pseudo_z_median = z_stats[["median"]],
    pseudo_z_mean = z_stats[["mean"]],
    pseudo_z_q95 = z_stats[["q95"]],
    pseudo_z_max = z_stats[["max"]],
    weight_colsum_min = min(w_colsum),
    weight_colsum_max = max(w_colsum),
    no_celltype_used = TRUE,
    stringsAsFactors = FALSE
  )

  cat("  wrote: ", rna_out, "\n", sep = "")
  cat("  wrote: ", z_out, "\n", sep = "")
  cat("  wrote: ", z_csv, "\n", sep = "")

  rm(mapping, W, rna_log, pseudo_rna, pseudo_rna_obj,
     pseudo_z, pseudo_z_dt, pseudo_z_obj)
  invisible(gc(verbose = FALSE))
}

cat("\n[4/4] Writing transfer summary ...\n")
summary_df <- do.call(rbind, summary_rows)
summary_file <- file.path(OUT_DIR, "transfer_summary.csv")
write.csv(summary_df, summary_file, row.names = FALSE)

cat("\n=== Dual transfer finished ===\n")
print(summary_df, row.names = FALSE)
cat("\nOutputs are in: ", OUT_DIR, "\n", sep = "")
