# Top5 RNA-completion significance filtering for the final v3.3 non-paired path
# -----------------------------------------------------------------------------
# Purpose:
#   Starting from a per-trait ATAC->RNA top5 mapping, remove genes whose top5 RNA
#   average is only weakly nonzero because of smoothing. The gene universe is
#   mapping$g_trait_used, so each trait uses its own sparse mapping gene set.
#
# For each trait and ATAC cell:
#   real  = weighted average of its mapped top5 CITE RNA expression on G_trait
#   null  = B random top5 CITE RNA averages on the same G_trait
#   keep  = genes passing SIGNIFICANCE_METHOD, either empirical P or cell-level
#           BH-FDR, plus positive delta over the random top5 background median.
#
# Output is intentionally compatible with:
#   03_pseudopair_pipeline/02_prepare_module3_v2_inputs.R
#
# Run in RStudio through a driver such as RUN_V3_SPARSE_OTHER_TRAITS.R, or
# directly after defining MAP_DIR / TRANSFER_OUT_DIR / TRAITS_TO_RUN.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
})

if (!exists("ROOT")) ROOT <- Sys.getenv("SCTARF_NONPAIRED_ROOT", file.path("outputs", "nonpaired"))

if (!exists("MAP_DIR")) MAP_DIR <- file.path(ROOT, "02_trait_mapping/output")
if (!exists("TRANSFER_OUT_DIR")) {
  TRANSFER_OUT_DIR <- file.path(ROOT, "03_pseudopair_pipeline/output/v3_sparse_mapping_top5significant/transfer")
}
if (!exists("DIAG_DIR")) {
  DIAG_DIR <- file.path(ROOT, "02_trait_mapping/output/v3_sparse_mapping_top5significant/diagnostics")
}

CITE_EXTRACT_DIR <- file.path(MAP_DIR, "cite_extract")
if (!exists("CITE_MTX")) CITE_MTX <- file.path(CITE_EXTRACT_DIR, "cite_trait_counts.mtx")
if (!exists("CITE_GENES")) CITE_GENES <- file.path(CITE_EXTRACT_DIR, "cite_trait_counts_genes.txt")
if (!exists("CITE_CELLS")) CITE_CELLS <- file.path(CITE_EXTRACT_DIR, "cite_trait_counts_cells.txt")
if (!exists("CITE_SCORE_RDS")) {
  CITE_SCORE_RDS <- file.path(ROOT, "01_data_prep_external_rna/output/cite_rna_scdrs_scores.rds")
}

if (!exists("TRAITS_TO_RUN")) TRAITS_TO_RUN <- c("MONO")
if (!exists("EXPECTED_K_NEIGHBOR")) EXPECTED_K_NEIGHBOR <- 5L
if (!exists("B_NULL_TOP5")) B_NULL_TOP5 <- 200L
if (!exists("SIGNIFICANCE_METHOD")) SIGNIFICANCE_METHOD <- "empirical_p"
if (!exists("EMPIRICAL_P_THRESHOLD")) EMPIRICAL_P_THRESHOLD <- 0.05
if (!exists("FILTER_Q_THRESHOLD")) FILTER_Q_THRESHOLD <- 0.10
if (!exists("MIN_DELTA_OVER_NULL_MEDIAN")) MIN_DELTA_OVER_NULL_MEDIAN <- 0
if (!exists("MAX_RETAINED_GENES_PER_CELL")) MAX_RETAINED_GENES_PER_CELL <- Inf
if (!exists("FILTER_BATCH_SIZE")) FILTER_BATCH_SIZE <- 512L
if (!exists("CPM_SCALE")) CPM_SCALE <- 1e4
if (!exists("RNG_SEED")) RNG_SEED <- 20260518L
if (!exists("NULL_WEIGHT_MODE")) NULL_WEIGHT_MODE <- "sample_observed"

dir.create(TRANSFER_OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DIAG_DIR, recursive = TRUE, showWarnings = FALSE)

stop_missing <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0L) {
    stop("Missing required file(s):\n  ", paste(missing, collapse = "\n  "),
         call. = FALSE)
  }
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
    q05 = as.numeric(stats::quantile(x, 0.05, na.rm = TRUE)),
    median = median(x, na.rm = TRUE),
    mean = mean(x, na.rm = TRUE),
    q95 = as.numeric(stats::quantile(x, 0.95, na.rm = TRUE)),
    max = max(x, na.rm = TRUE)
  )
}

make_weight_matrix <- function(neighbor_idx, weight_mat, n_rna, cols) {
  k <- ncol(neighbor_idx)
  idx <- neighbor_idx[cols, , drop = FALSE]
  w <- weight_mat[cols, , drop = FALSE]
  sparseMatrix(
    i = as.integer(as.vector(t(idx))),
    j = rep(seq_along(cols), each = k),
    x = as.numeric(as.vector(t(w))),
    dims = c(n_rna, length(cols))
  )
}

make_null_weight_matrix <- function(mapping, n_rna, b_null, mode) {
  k <- mapping$k
  random_idx <- t(replicate(b_null, sample.int(n_rna, k, replace = FALSE)))
  if (identical(mode, "equal")) {
    random_w <- matrix(1 / k, nrow = b_null, ncol = k)
  } else if (identical(mode, "sample_observed")) {
    weight_rows <- sample.int(nrow(mapping$neighbor_weight), b_null, replace = TRUE)
    random_w <- mapping$neighbor_weight[weight_rows, , drop = FALSE]
    random_w <- random_w / rowSums(random_w)
  } else {
    stop("Unknown NULL_WEIGHT_MODE: ", mode, call. = FALSE)
  }

  sparseMatrix(
    i = as.integer(as.vector(t(random_idx))),
    j = rep(seq_len(b_null), each = k),
    x = as.numeric(as.vector(t(random_w))),
    dims = c(n_rna, b_null)
  )
}

compute_p_batch <- function(real_batch, null_sorted) {
  b_null <- ncol(null_sorted)
  p <- matrix(1, nrow = nrow(real_batch), ncol = ncol(real_batch))
  for (g in seq_len(nrow(real_batch))) {
    # Ties are rare for weighted log-normalized expression. This computes the
    # one-sided empirical tail P(null >= real) with the standard +1 correction.
    n_le <- findInterval(real_batch[g, ], null_sorted[g, ])
    p[g, ] <- (1 + b_null - n_le) / (1 + b_null)
  }
  p
}

filter_real_batch <- function(real_batch, null_sorted, null_median,
                              significance_method, p_threshold,
                              q_threshold, min_delta, max_retained) {
  p <- compute_p_batch(real_batch, null_sorted)
  positive_delta <- sweep(real_batch, 1L, null_median + min_delta, ">")

  if (identical(significance_method, "empirical_p")) {
    keep <- p <= p_threshold & positive_delta & real_batch > 0
  } else if (identical(significance_method, "cell_bh")) {
    keep <- matrix(FALSE, nrow = nrow(real_batch), ncol = ncol(real_batch))
    for (j in seq_len(ncol(real_batch))) {
      q <- stats::p.adjust(p[, j], method = "BH")
      keep[, j] <- q <= q_threshold & positive_delta[, j] & real_batch[, j] > 0
    }
  } else {
    stop("Unknown SIGNIFICANCE_METHOD: ", significance_method, call. = FALSE)
  }

  if (is.finite(max_retained)) {
    delta <- sweep(real_batch, 1L, null_median, "-")
    for (j in seq_len(ncol(keep))) {
      kept <- which(keep[, j])
      if (length(kept) > max_retained) {
        ord <- kept[order(p[kept, j], -delta[kept, j], decreasing = FALSE)]
        keep[kept, j] <- FALSE
        keep[head(ord, max_retained), j] <- TRUE
      }
    }
  }
  keep
}

build_top5_significant_matrix <- function(rna_log, mapping, cite_cells) {
  n_rna <- length(cite_cells)
  n_atac <- length(mapping$atac_barcodes)
  genes <- rownames(rna_log)

  neighbor_idx <- match(as.vector(mapping$rna_neighbor_barcode), cite_cells)
  if (anyNA(neighbor_idx)) {
    stop("Some mapped RNA neighbor barcodes are absent from the CITE cell list.",
         call. = FALSE)
  }
  neighbor_idx <- matrix(
    neighbor_idx,
    nrow = nrow(mapping$rna_neighbor_barcode),
    ncol = ncol(mapping$rna_neighbor_barcode)
  )

  cat("  building random top5 background: B=", B_NULL_TOP5,
      "; mode=", NULL_WEIGHT_MODE, " ...\n", sep = "")
  W_null <- make_null_weight_matrix(mapping, n_rna, B_NULL_TOP5, NULL_WEIGHT_MODE)
  null_mat <- as.matrix(rna_log %*% W_null)  # gene x random-top5 group
  rownames(null_mat) <- genes
  null_sorted <- t(apply(null_mat, 1L, sort))
  null_median <- apply(null_mat, 1L, median)

  null_gene_dt <- data.table(
    gene = genes,
    null_mean = rowMeans(null_mat),
    null_median = null_median,
    null_q95 = apply(null_mat, 1L, stats::quantile, probs = 0.95),
    null_max = apply(null_mat, 1L, max)
  )

  i_list <- list()
  j_list <- list()
  x_list <- list()
  retained_n <- integer(n_atac)
  real_sum <- numeric(n_atac)
  kept_sum <- numeric(n_atac)

  n_batch <- ceiling(n_atac / FILTER_BATCH_SIZE)
  t0 <- Sys.time()
  for (bi in seq_len(n_batch)) {
    cols <- ((bi - 1L) * FILTER_BATCH_SIZE + 1L):min(bi * FILTER_BATCH_SIZE, n_atac)
    W_real <- make_weight_matrix(neighbor_idx, mapping$neighbor_weight, n_rna, cols)
    real_batch <- as.matrix(rna_log %*% W_real)
    keep <- filter_real_batch(
      real_batch,
      null_sorted = null_sorted,
      null_median = null_median,
      significance_method = SIGNIFICANCE_METHOD,
      p_threshold = EMPIRICAL_P_THRESHOLD,
      q_threshold = FILTER_Q_THRESHOLD,
      min_delta = MIN_DELTA_OVER_NULL_MEDIAN,
      max_retained = MAX_RETAINED_GENES_PER_CELL
    )

    retained_n[cols] <- colSums(keep)
    real_sum[cols] <- colSums(real_batch)
    kept_batch <- real_batch * keep
    kept_sum[cols] <- colSums(kept_batch)

    ix <- which(keep, arr.ind = TRUE)
    if (nrow(ix) > 0L) {
      li <- length(i_list) + 1L
      i_list[[li]] <- ix[, 1]
      j_list[[li]] <- cols[ix[, 2]]
      x_list[[li]] <- real_batch[ix]
    }

    if (bi %% 10L == 0L || bi == n_batch) {
      elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      cat(sprintf(
        "    batch %d/%d finished (%d/%d ATAC cells, median retained=%s, %.1fs)\n",
        bi, n_batch, max(cols), n_atac,
        signif(stats::median(retained_n[seq_len(max(cols))]), 4),
        elapsed
      ))
      flush.console()
    }

    rm(W_real, real_batch, keep, kept_batch, ix)
    invisible(gc(verbose = FALSE))
  }

  if (length(i_list) == 0L) {
    filtered <- Matrix(0, nrow = length(genes), ncol = n_atac, sparse = TRUE)
    rownames(filtered) <- genes
    colnames(filtered) <- mapping$atac_barcodes
  } else {
    filtered <- sparseMatrix(
      i = as.integer(unlist(i_list, use.names = FALSE)),
      j = as.integer(unlist(j_list, use.names = FALSE)),
      x = as.numeric(unlist(x_list, use.names = FALSE)),
      dims = c(length(genes), n_atac),
      dimnames = list(genes, mapping$atac_barcodes)
    )
  }

  cell_dt <- data.table(
    barcode = mapping$atac_barcodes,
    retained_gene_count = retained_n,
    retained_gene_fraction = retained_n / length(genes),
    top5_real_expression_sum = real_sum,
    top5_significant_expression_sum = kept_sum,
    retained_expression_fraction = fifelse(real_sum > 0, kept_sum / real_sum, 0),
    no_celltype_used = TRUE
  )

  retained_value_max <- numeric(nrow(filtered))
  if (length(filtered@x) > 0L) {
    row_max <- tapply(filtered@x, filtered@i + 1L, max)
    retained_value_max[as.integer(names(row_max))] <- as.numeric(row_max)
  }

  gene_dt <- data.table(
    gene = genes,
    retained_cell_count = Matrix::rowSums(filtered != 0),
    retained_cell_rate = Matrix::rowSums(filtered != 0) / n_atac,
    retained_value_mean = Matrix::rowMeans(filtered),
    retained_value_max = retained_value_max
  )
  gene_dt <- merge(gene_dt, null_gene_dt, by = "gene", all.x = TRUE)

  list(
    matrix = as(filtered, "CsparseMatrix"),
    cell_diagnostics = cell_dt,
    gene_diagnostics = gene_dt
  )
}

cat(">>> Build top5-significant RNA-completion profiles\n\n")
cat("MAP_DIR          : ", MAP_DIR, "\n", sep = "")
cat("TRANSFER_OUT_DIR : ", TRANSFER_OUT_DIR, "\n", sep = "")
cat("DIAG_DIR         : ", DIAG_DIR, "\n", sep = "")
cat("TRAITS_TO_RUN    : ", paste(TRAITS_TO_RUN, collapse = ", "), "\n", sep = "")
cat("B_NULL_TOP5      : ", B_NULL_TOP5, "\n", sep = "")
cat("SIGNIFICANCE     : ", SIGNIFICANCE_METHOD, "\n", sep = "")
cat("EMPIRICAL_P      : ", EMPIRICAL_P_THRESHOLD, "\n", sep = "")
cat("FILTER_Q         : ", FILTER_Q_THRESHOLD, " (used only when SIGNIFICANCE_METHOD='cell_bh')\n", sep = "")
cat("MIN_DELTA        : ", MIN_DELTA_OVER_NULL_MEDIAN, "\n\n", sep = "")

set.seed(RNG_SEED)

stop_missing(c(CITE_MTX, CITE_GENES, CITE_CELLS, CITE_SCORE_RDS))
stop_missing(file.path(MAP_DIR, paste0("granja_cite_mapping_", TRAITS_TO_RUN, ".rds")))

cat("[1/3] Loading CITE counts and scDRS scores ...\n")
cite_counts <- Matrix::readMM(CITE_MTX)
cite_counts <- as(cite_counts, "CsparseMatrix")
cite_genes <- readLines(CITE_GENES, warn = FALSE)
cite_cells <- readLines(CITE_CELLS, warn = FALSE)
rownames(cite_counts) <- cite_genes
colnames(cite_counts) <- cite_cells

cite_scores <- as.data.frame(readRDS(CITE_SCORE_RDS))
score_cols <- paste0(TRAITS_TO_RUN, "_norm_score")
needed_cols <- c("cell", "scdrs_input_total_counts", score_cols)
missing_cols <- setdiff(needed_cols, colnames(cite_scores))
if (length(missing_cols) > 0L) {
  stop("CITE score RDS missing required column(s): ",
       paste(missing_cols, collapse = ", "), call. = FALSE)
}

lib_size <- cite_scores$scdrs_input_total_counts[match(cite_cells, cite_scores$cell)]
names(lib_size) <- cite_cells
if (anyNA(lib_size) || any(!is.finite(lib_size))) {
  stop("Missing CITE library sizes. Check CITE_SCORE_RDS.", call. = FALSE)
}
lib_size[lib_size < 1] <- 1

summary_rows <- list()

cat("[2/3] Filtering top5 RNA-completion profiles per trait ...\n")
for (trait in TRAITS_TO_RUN) {
  cat("\n=== Trait ", trait, " ===\n", sep = "")
  map_file <- file.path(MAP_DIR, paste0("granja_cite_mapping_", trait, ".rds"))
  mapping <- readRDS(map_file)
  if (!identical(as.integer(mapping$k), as.integer(EXPECTED_K_NEIGHBOR))) {
    stop("Expected mapping k=", EXPECTED_K_NEIGHBOR, " for ", trait,
         " but found k=", mapping$k, call. = FALSE)
  }

  genes <- mapping$g_trait_used
  missing_genes <- setdiff(genes, rownames(cite_counts))
  if (length(missing_genes) > 0L) {
    stop("CITE extracted matrix is missing mapping genes for ", trait, ": ",
         paste(head(missing_genes, 10), collapse = ", "), call. = FALSE)
  }

  cat("  log-normalizing CITE RNA counts on G_trait ...\n")
  rna_log <- log_normalize_counts(cite_counts[genes, , drop = FALSE],
                                  lib_size = lib_size,
                                  scale_factor = CPM_SCALE)

  res <- build_top5_significant_matrix(rna_log, mapping, cite_cells)
  pseudo_rna <- res$matrix

  score_col <- paste0(trait, "_norm_score")
  score <- cite_scores[[score_col]][match(cite_cells, cite_scores$cell)]
  if (anyNA(score) || any(!is.finite(score))) {
    replacement <- median(score[is.finite(score)], na.rm = TRUE)
    warning("Some ", score_col, " values are missing; replacing with median = ",
            signif(replacement, 4))
    score[!is.finite(score) | is.na(score)] <- replacement
  }

  neighbor_idx <- match(as.vector(mapping$rna_neighbor_barcode), cite_cells)
  neighbor_idx <- matrix(
    neighbor_idx,
    nrow = nrow(mapping$rna_neighbor_barcode),
    ncol = ncol(mapping$rna_neighbor_barcode)
  )
  W_all <- make_weight_matrix(
    neighbor_idx,
    mapping$neighbor_weight,
    n_rna = length(cite_cells),
    cols = seq_along(mapping$atac_barcodes)
  )
  pseudo_z <- as.numeric(Matrix::crossprod(W_all, score))
  names(pseudo_z) <- mapping$atac_barcodes

  max_weight <- apply(mapping$neighbor_weight, 1L, max)
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
    source = "top5 CITE RNA log-normalized expression filtered by random top5 background significance",
    cpm_scale = CPM_SCALE,
    mapping_file = map_file,
    filter = list(
      method = "top5_significant_random_background",
      significance_method = SIGNIFICANCE_METHOD,
      b_null_top5 = B_NULL_TOP5,
      empirical_p_threshold = EMPIRICAL_P_THRESHOLD,
      q_threshold = FILTER_Q_THRESHOLD,
      min_delta_over_null_median = MIN_DELTA_OVER_NULL_MEDIAN,
      max_retained_genes_per_cell = MAX_RETAINED_GENES_PER_CELL,
      null_weight_mode = NULL_WEIGHT_MODE,
      rng_seed = RNG_SEED
    ),
    no_celltype_used = TRUE,
    created_at = as.character(Sys.time())
  )
  pseudo_z_obj <- list(
    trait = trait,
    trait_id = mapping$trait_id,
    scores = pseudo_z_dt,
    source = "weighted average of top5 CITE RNA scDRS norm_score; unchanged by gene-level completion filter",
    mapping_file = map_file,
    no_celltype_used = TRUE,
    created_at = as.character(Sys.time())
  )

  rna_out <- file.path(TRANSFER_OUT_DIR, paste0("granja_pseudo_rna_profile_", trait, ".rds"))
  z_out <- file.path(TRANSFER_OUT_DIR, paste0("granja_pseudo_zrna_", trait, ".rds"))
  z_csv <- file.path(TRANSFER_OUT_DIR, paste0("granja_pseudo_zrna_", trait, ".csv"))
  saveRDS(pseudo_rna_obj, rna_out)
  saveRDS(pseudo_z_obj, z_out)
  fwrite(pseudo_z_dt, z_csv)

  cell_diag <- copy(res$cell_diagnostics)
  cell_diag[, trait := trait]
  setcolorder(cell_diag, c("trait", setdiff(names(cell_diag), "trait")))
  gene_diag <- copy(res$gene_diagnostics)
  gene_diag[, trait := trait]
  setcolorder(gene_diag, c("trait", setdiff(names(gene_diag), "trait")))

  fwrite(cell_diag, file.path(DIAG_DIR, paste0("top5_significant_cell_diagnostics_", trait, ".csv.gz")))
  fwrite(gene_diag, file.path(DIAG_DIR, paste0("top5_significant_gene_diagnostics_", trait, ".csv")))

  retained_stats <- vector_stats(cell_diag$retained_gene_count)
  summary_rows[[trait]] <- data.table(
    trait = trait,
    trait_id = mapping$trait_id,
    n_atac_cells = length(mapping$atac_barcodes),
    n_cite_cells = length(cite_cells),
    n_g_trait_used = length(genes),
    k = mapping$k,
    significance_method = SIGNIFICANCE_METHOD,
    b_null_top5 = B_NULL_TOP5,
    empirical_p_threshold = EMPIRICAL_P_THRESHOLD,
    filter_q_threshold = FILTER_Q_THRESHOLD,
    min_delta_over_null_median = MIN_DELTA_OVER_NULL_MEDIAN,
    max_retained_genes_per_cell = MAX_RETAINED_GENES_PER_CELL,
    null_weight_mode = NULL_WEIGHT_MODE,
    retained_gene_min = retained_stats[["min"]],
    retained_gene_q05 = retained_stats[["q05"]],
    retained_gene_median = retained_stats[["median"]],
    retained_gene_mean = retained_stats[["mean"]],
    retained_gene_q95 = retained_stats[["q95"]],
    retained_gene_max = retained_stats[["max"]],
    empty_cell_rate = mean(cell_diag$retained_gene_count == 0),
    pseudo_rna_nnz = length(pseudo_rna@x),
    pseudo_rna_density = length(pseudo_rna@x) / max(prod(dim(pseudo_rna)), 1),
    no_celltype_used = TRUE
  )

  cat("  retained genes/cell median=", retained_stats[["median"]],
      "; q05=", retained_stats[["q05"]],
      "; q95=", retained_stats[["q95"]],
      "; empty rate=", signif(mean(cell_diag$retained_gene_count == 0), 4),
      "\n", sep = "")
  cat("  wrote: ", rna_out, "\n", sep = "")
  cat("  wrote: ", z_out, "\n", sep = "")

  rm(mapping, rna_log, res, pseudo_rna, pseudo_rna_obj, pseudo_z_obj,
     pseudo_z_dt, W_all, cell_diag, gene_diag)
  invisible(gc(verbose = FALSE))
}

cat("\n[3/3] Writing top5-significant transfer summary ...\n")
summary_dt <- rbindlist(summary_rows, use.names = TRUE, fill = TRUE)
summary_file <- file.path(TRANSFER_OUT_DIR, "top5_significant_transfer_summary.csv")
fwrite(summary_dt, summary_file)
fwrite(summary_dt, file.path(DIAG_DIR, "top5_significant_transfer_summary.csv"))

cat("\n=== top5-significant RNA-completion profiles are ready ===\n")
print(summary_dt)
cat("\nTransfer outputs are in: ", TRANSFER_OUT_DIR, "\n", sep = "")
cat("Diagnostics are in: ", DIAG_DIR, "\n", sep = "")
