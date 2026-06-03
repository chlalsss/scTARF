# module4 v3.3: sparse trait-aware ATAC -> RNA mapping
# -----------------------------------------------------------------------------
# Purpose:
#   Move the teacher-feedback sparsification step upstream into the mapping
#   stage. Instead of finding RNA neighbors in the full dense G_trait vector,
#   build cell-specific sparse trait vectors for ATAC and RNA, then search top5
#   RNA references from those sparse vectors.
#
# Trait handling:
#   The script is generic. For each entry in TRAITS_TO_RUN it reads that trait's
#   gene set from GS_FILE, intersects it with ATAC gene activity and CITE RNA
#   genes, and writes granja_cite_mapping_<trait>.rds. MONO is only the default
#   standalone value when the script is sourced without an external driver.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
})

ROOT <- "D:/RStudio/work/module4_traitaware_pseudopair"

if (!exists("SPARSE_MAP_OUT_DIR")) {
  SPARSE_MAP_OUT_DIR <- file.path(ROOT, "02_trait_mapping/output/v3_sparse_mapping")
}
if (!exists("DIAG_DIR")) {
  DIAG_DIR <- file.path(SPARSE_MAP_OUT_DIR, "diagnostics")
}
if (!exists("CITE_EXTRACT_DIR")) {
  CITE_EXTRACT_DIR <- file.path(ROOT, "02_trait_mapping/output/cite_extract")
}

if (!exists("GA_PATH")) {
  GA_PATH <- "D:/SCAVENGE-reproducibility-main/blood_data/multiomics_output_v2/gene_activity_matrix.rds"
}
if (!exists("ATAC_META_PATH")) {
  ATAC_META_PATH <- file.path(ROOT, "04_evaluation/inputs/granja_cell_meta.rds")
}
if (!exists("GS_FILE")) {
  GS_FILE <- file.path(ROOT, "01_data_prep/03_rna_scdrs_python/granja_3traits.gs")
}
if (!exists("CITE_H5AD")) {
  CITE_H5AD <- file.path(
    ROOT,
    "01_data_prep_external_rna/03_cite_scdrs_python/input/cite_rna_raw_counts_gex_only.h5ad"
  )
}
if (!exists("CITE_SCORE_RDS")) {
  CITE_SCORE_RDS <- file.path(ROOT, "01_data_prep_external_rna/output/cite_rna_scdrs_scores.rds")
}
if (!exists("PY_HELPER")) {
  PY_HELPER <- file.path(ROOT, "02_trait_mapping/extract_cite_trait_counts.py")
}
if (!exists("PYTHON_EXE")) PYTHON_EXE <- "D:/python/python.exe"

if (!exists("TRAITS_TO_RUN")) TRAITS_TO_RUN <- c("MONO")
if (!exists("TRAITS_ALL")) {
  TRAITS_ALL <- c(
    MONO = "UKB_460K.blood_MONOCYTE_COUNT",
    PLT  = "UKB_460K.blood_PLATELET_COUNT",
    LYM  = "UKB_460K.blood_LYMPHOCYTE_COUNT"
  )
}

if (!exists("K_NEIGHBOR")) K_NEIGHBOR <- 5L
if (!exists("BATCH_SIZE")) BATCH_SIZE <- 512L
if (!exists("B_NULL_CELL")) B_NULL_CELL <- 1000L
if (!exists("RNA_EMPIRICAL_P_THRESHOLD")) RNA_EMPIRICAL_P_THRESHOLD <- 0.05
if (!exists("ATAC_EMPIRICAL_P_THRESHOLD")) ATAC_EMPIRICAL_P_THRESHOLD <- 0.05
if (!exists("MIN_DELTA_OVER_NULL_MEDIAN")) MIN_DELTA_OVER_NULL_MEDIAN <- 0
if (!exists("USE_GWAS_WEIGHT")) USE_GWAS_WEIGHT <- TRUE
if (!exists("CPM_SCALE")) CPM_SCALE <- 1e4
if (!exists("RNG_SEED")) RNG_SEED <- 20260519L
if (!exists("FORCE_REEXTRACT_CITE")) FORCE_REEXTRACT_CITE <- FALSE

EPS <- 1e-8

dir.create(SPARSE_MAP_OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DIAG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(CITE_EXTRACT_DIR, recursive = TRUE, showWarnings = FALSE)

TRAITS <- TRAITS_ALL[TRAITS_TO_RUN]
if (anyNA(TRAITS) || length(TRAITS) == 0L) {
  stop("TRAITS_TO_RUN contains unknown traits: ",
       paste(setdiff(TRAITS_TO_RUN, names(TRAITS_ALL)), collapse = ", "),
       call. = FALSE)
}

ATAC_TARGETS <- list(
  MONO = c("11_CD14.Mono.1", "12_CD14.Mono.2"),
  PLT  = c("02_Early.Eryth", "03_Late.Eryth", "04_Early.Baso"),
  LYM  = c("16_Pre.B", "17_B", "18_Plasma",
           "19_CD8.N", "23_CD8.EM", "24_CD8.CM",
           "20_CD4.N1", "21_CD4.N2", "22_CD4.M",
           "25_NK")
)

stop_missing <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0L) {
    stop("Missing required file(s):\n  ", paste(missing, collapse = "\n  "),
         call. = FALSE)
  }
}

parse_gs_file <- function(path) {
  lines <- readLines(path, warn = FALSE)
  out <- list()
  for (line in lines[-1]) {
    parts <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(parts) < 2) next
    pairs <- strsplit(strsplit(parts[2], ",", fixed = TRUE)[[1]], ":", fixed = TRUE)
    genes <- vapply(pairs, `[`, character(1), 1)
    weights <- as.numeric(vapply(
      pairs,
      function(x) if (length(x) >= 2) x[2] else "1",
      character(1)
    ))
    out[[parts[1]]] <- list(genes = genes, weights = weights)
  }
  out
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

gc_quiet <- function() invisible(gc(verbose = FALSE))

chunk_starts <- function(n, chunk_size) {
  if (n <= 0L) return(integer())
  seq.int(1L, n, by = chunk_size)
}

add_group_sums <- function(acc, rows, vals) {
  sums <- rowsum(vals, rows, reorder = FALSE)
  row_ids <- as.integer(rownames(sums))
  acc[row_ids] <- acc[row_ids] + as.numeric(sums[, 1L])
  acc
}

sparse_row_moments <- function(m, chunk_size = 1000000L) {
  n_genes <- nrow(m)
  n_cells <- ncol(m)
  row_sum <- numeric(n_genes)
  row_sum_sq <- numeric(n_genes)
  n_x <- length(m@x)

  for (start in chunk_starts(n_x, chunk_size)) {
    end <- min(start + chunk_size - 1L, n_x)
    ix <- start:end
    rows <- m@i[ix] + 1L
    vals <- as.numeric(m@x[ix])
    row_sum <- add_group_sums(row_sum, rows, vals)
    row_sum_sq <- add_group_sums(row_sum_sq, rows, vals * vals)
  }

  mu <- row_sum / n_cells
  mu2 <- row_sum_sq / n_cells
  list(mu = mu, sd = sqrt(pmax(mu2 - mu * mu, 0)))
}

nonzero_entries_by_row <- function(m, chunk_size = 1000000L) {
  n_genes <- nrow(m)
  out <- vector("list", n_genes)
  n_x <- length(m@x)

  for (start in chunk_starts(n_x, chunk_size)) {
    end <- min(start + chunk_size - 1L, n_x)
    ix <- start:end
    rows <- m@i[ix] + 1L
    pieces <- split(ix, rows)
    for (nm in names(pieces)) {
      g <- as.integer(nm)
      out[[g]] <- c(out[[g]], pieces[[nm]])
    }
  }

  out
}

compute_auc_rank <- function(score, label) {
  keep <- is.finite(score) & !is.na(label)
  score <- score[keep]
  label <- as.integer(label[keep])
  n_pos <- sum(label == 1L)
  n_neg <- sum(label == 0L)
  if (n_pos < 10L || n_neg < 10L) return(NA_real_)
  if (!is.finite(stats::var(score)) || stats::var(score) < 1e-12) return(NA_real_)
  r <- rank(score, ties.method = "average")
  (sum(r[label == 1L]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

top_metrics <- function(score, label, top_fraction) {
  keep <- is.finite(score) & !is.na(label)
  score <- score[keep]
  label <- as.integer(label[keep])
  n <- length(score)
  n_top <- max(1L, ceiling(n * top_fraction))
  idx <- order(score, decreasing = TRUE)[seq_len(n_top)]
  precision <- mean(label[idx] == 1L)
  target_rate <- mean(label == 1L)
  c(
    n_top = n_top,
    precision = precision,
    enrichment = precision / max(target_rate, 1e-12)
  )
}

extract_cite_counts_if_needed <- function(all_requested_genes) {
  requested_file <- file.path(CITE_EXTRACT_DIR, "requested_trait_genes.txt")
  out_mtx <- file.path(CITE_EXTRACT_DIR, "cite_trait_counts.mtx")
  out_genes <- file.path(CITE_EXTRACT_DIR, "cite_trait_counts_genes.txt")
  out_cells <- file.path(CITE_EXTRACT_DIR, "cite_trait_counts_cells.txt")
  out_summary <- file.path(CITE_EXTRACT_DIR, "cite_trait_counts_summary.csv")
  out_missing <- file.path(CITE_EXTRACT_DIR, "cite_trait_counts_missing_genes.txt")

  writeLines(sort(unique(all_requested_genes)), requested_file, useBytes = TRUE)

  outputs <- c(out_mtx, out_genes, out_cells, out_summary)
  if (!FORCE_REEXTRACT_CITE && all(file.exists(outputs))) {
    present_genes <- readLines(out_genes, warn = FALSE)
    missing_genes <- if (file.exists(out_missing)) readLines(out_missing, warn = FALSE) else character()
    covered_genes <- unique(c(present_genes, missing_genes))
    if (all(sort(unique(all_requested_genes)) %in% covered_genes)) {
      cat("[extract] Reusing existing CITE trait-count matrix.\n")
      return(list(mtx = out_mtx, genes = out_genes, cells = out_cells))
    }
    cat("[extract] Existing CITE extraction does not cover all requested genes; re-extracting.\n")
  }

  cat("[extract] Extracting CITE RNA counts from h5ad via Python ...\n")
  stop_missing(c(PYTHON_EXE, PY_HELPER, CITE_H5AD))
  args <- c(
    PY_HELPER,
    "--h5ad", CITE_H5AD,
    "--genes", requested_file,
    "--out-mtx", out_mtx,
    "--out-genes", out_genes,
    "--out-cells", out_cells,
    "--out-summary", out_summary
  )
  status <- system2(PYTHON_EXE, args = args)
  if (!identical(status, 0L)) {
    stop("Python CITE extraction failed with status ", status, call. = FALSE)
  }
  list(mtx = out_mtx, genes = out_genes, cells = out_cells)
}

build_sparse_trait_features <- function(mat_g_by_c, gene_weight, b_null,
                                        p_threshold, min_delta, use_weight,
                                        modality, trait) {
  genes <- rownames(mat_g_by_c)
  cells <- colnames(mat_g_by_c)
  n_genes <- length(genes)
  n_cells <- length(cells)

  cat("  [", modality, "] sparse empirical background filtering for ",
      n_genes, " genes x ", n_cells, " cells ...\n", sep = "")

  mat_g_by_c <- as(mat_g_by_c, "CsparseMatrix")
  gc_quiet()

  moments <- sparse_row_moments(mat_g_by_c)
  mu <- moments$mu
  sd <- moments$sd
  rm(moments)
  sd[sd < EPS] <- EPS

  b_eff <- min(as.integer(b_null), n_cells)
  bg_idx <- sample.int(n_cells, b_eff, replace = n_cells < b_eff)
  null_mat <- as.matrix(mat_g_by_c[, bg_idx, drop = FALSE])
  storage.mode(null_mat) <- "double"
  null_sorted <- t(apply(null_mat, 1L, sort))
  null_median <- apply(null_mat, 1L, median)

  i_list <- vector("list", n_genes)
  j_list <- vector("list", n_genes)
  x_list <- vector("list", n_genes)
  retained_n <- integer(n_cells)
  retained_sum <- numeric(n_cells)
  retained_cell_count <- integer(n_genes)
  retained_value_sum <- numeric(n_genes)
  retained_value_max <- rep(NA_real_, n_genes)

  cell_ptr <- mat_g_by_c@p
  nnz_by_row <- nonzero_entries_by_row(mat_g_by_c)
  gc_quiet()

  for (g in seq_len(n_genes)) {
    nz <- nnz_by_row[[g]]
    if (length(nz) == 0L) {
      retained_value_max[g] <- 0
      next
    }

    jj_all <- findInterval(nz - 1L, cell_ptr)
    obs_g <- as.numeric(mat_g_by_c@x[nz])
    positive <- obs_g > 0
    if (!any(positive)) {
      retained_value_max[g] <- 0
      next
    }
    jj_all <- jj_all[positive]
    obs_g <- obs_g[positive]

    n_lt <- findInterval(obs_g, null_sorted[g, ], left.open = TRUE)
    p_g <- (1 + b_eff - n_lt) / (1 + b_eff)
    z_g <- (obs_g - mu[g]) / sd[g]
    keep <- p_g <= p_threshold &
      obs_g > (null_median[g] + min_delta) &
      obs_g > 0 &
      z_g > 0

    if (any(keep)) {
      jj <- jj_all[keep]
      vals <- z_g[keep]
      if (isTRUE(use_weight)) vals <- vals * gene_weight[g]
      vals <- as.numeric(vals)
      i_list[[g]] <- rep.int(g, length(jj))
      j_list[[g]] <- jj
      x_list[[g]] <- vals
      retained_n[jj] <- retained_n[jj] + 1L
      retained_sum[jj] <- retained_sum[jj] + vals
      retained_cell_count[g] <- length(jj)
      retained_value_sum[g] <- sum(vals)
      retained_value_max[g] <- max(vals)
    } else {
      retained_value_max[g] <- 0
    }

    if (g %% 100L == 0L || g == n_genes) {
      cat("    [", modality, "] sparse feature genes ",
          g, "/", n_genes, " done\n", sep = "")
      flush.console()
    }
  }

  nnz <- sum(lengths(j_list))
  if (nnz == 0L) {
    feature <- Matrix(0, nrow = n_genes, ncol = n_cells, sparse = TRUE)
    dimnames(feature) <- list(genes, cells)
  } else {
    feature <- sparseMatrix(
      i = as.integer(unlist(i_list, use.names = FALSE)),
      j = as.integer(unlist(j_list, use.names = FALSE)),
      x = as.numeric(unlist(x_list, use.names = FALSE)),
      dims = c(n_genes, n_cells),
      dimnames = list(genes, cells)
    )
  }

  cell_dt <- data.table(
    trait = trait,
    modality = modality,
    barcode = cells,
    retained_gene_count = retained_n,
    retained_gene_fraction = retained_n / n_genes,
    retained_feature_sum = retained_sum,
    empty_sparse_vector = retained_n == 0L,
    no_celltype_used = TRUE
  )

  gene_dt <- data.table(
    trait = trait,
    modality = modality,
    gene = genes,
    retained_cell_count = retained_cell_count,
    retained_cell_rate = retained_cell_count / n_cells,
    retained_value_mean = retained_value_sum / pmax(retained_cell_count, 1L),
    retained_value_max = retained_value_max,
    null_median = null_median,
    null_q95 = apply(null_mat, 1L, stats::quantile, probs = 0.95),
    null_max = apply(null_mat, 1L, max)
  )

  rm(nnz_by_row, null_mat, null_sorted)
  gc_quiet()

  list(
    feature = feature,
    cell_diag = cell_dt,
    gene_diag = gene_dt,
    b_eff = b_eff
  )
}

normalize_columns <- function(m) {
  norms <- sqrt(Matrix::colSums(m * m))
  zero <- !is.finite(norms) | norms < EPS
  norms[zero] <- 1
  out <- Matrix::t(Matrix::t(m) / norms)
  attr(out, "zero_columns") <- zero
  out
}

batched_top_k_sparse_cosine <- function(query_g_by_c, ref_g_by_c, k, batch = 512L) {
  qn <- normalize_columns(query_g_by_c)
  rn <- normalize_columns(ref_g_by_c)
  query_zero <- attr(qn, "zero_columns")
  ref_zero <- attr(rn, "zero_columns")

  n_q <- ncol(qn)
  n_batch <- ceiling(n_q / batch)
  idx_all <- matrix(0L, n_q, k)
  dist_all <- matrix(0, n_q, k)
  sim1_all <- numeric(n_q)
  t0 <- Sys.time()

  for (bi in seq_len(n_batch)) {
    start <- (bi - 1L) * batch + 1L
    end <- min(bi * batch, n_q)
    cols <- start:end
    sim <- Matrix::crossprod(qn[, cols, drop = FALSE], rn)
    sim <- as.matrix(sim)
    sim[, ref_zero] <- -Inf
    nb <- nrow(sim)
    row_ix <- seq_len(nb)

    for (j in seq_len(k)) {
      mx <- max.col(sim, ties.method = "first")
      val <- sim[cbind(row_ix, mx)]
      val[!is.finite(val)] <- 0
      val <- pmin(pmax(val, 0), 1)
      idx_all[cols, j] <- mx
      dist_all[cols, j] <- sqrt(pmax(2 - 2 * val, 0))
      if (j == 1L) sim1_all[cols] <- val
      sim[cbind(row_ix, mx)] <- -Inf
    }

    if (bi %% 5L == 0L || bi == n_batch) {
      elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      cat(sprintf("       sparse cosine batch %d/%d finished (%d/%d ATAC cells, %.1fs)\n",
                  bi, n_batch, end, n_q, elapsed))
      flush.console()
    }
    rm(sim)
    invisible(gc(verbose = FALSE))
  }

  idx_all[query_zero, ] <- idx_all[query_zero, , drop = FALSE]

  list(
    nn.idx = idx_all,
    nn.dists = dist_all,
    top1_cosine = sim1_all,
    query_zero = query_zero,
    ref_zero = ref_zero
  )
}

softmax_per_cell <- function(d_mat) {
  tau <- apply(d_mat, 1, function(x) max(stats::median(x), EPS))
  w <- exp(-d_mat / tau)
  w <- w / rowSums(w)
  list(weight = w, tau = tau)
}

cat(">>> v3.3 sparse trait-aware mapping\n\n")
cat("Traits: ", paste(names(TRAITS), collapse = ", "), "\n", sep = "")
cat("Output: ", SPARSE_MAP_OUT_DIR, "\n", sep = "")
cat("B_NULL_CELL = ", B_NULL_CELL,
    "; RNA_P = ", RNA_EMPIRICAL_P_THRESHOLD,
    "; ATAC_P = ", ATAC_EMPIRICAL_P_THRESHOLD,
    "; K = ", K_NEIGHBOR, "\n\n", sep = "")

set.seed(RNG_SEED)
stop_missing(c(GA_PATH, ATAC_META_PATH, GS_FILE, CITE_SCORE_RDS))

cat("[1/6] Reading gene sets ...\n")
gs <- parse_gs_file(GS_FILE)
all_trait_genes <- unique(unlist(lapply(unname(TRAITS_ALL), function(id) gs[[id]]$genes)))

cat("[2/6] Preparing CITE RNA trait-gene count matrix ...\n")
cite_files <- extract_cite_counts_if_needed(all_trait_genes)

cat("[3/6] Loading Granja ATAC Gene Activity and ATAC cell list ...\n")
ga <- readRDS(GA_PATH)
ga <- as(ga, "CsparseMatrix")
atac_meta <- readRDS(ATAC_META_PATH)
if (inherits(atac_meta, "data.table")) atac_meta <- as.data.frame(atac_meta)
if (is.null(rownames(atac_meta)) ||
    identical(rownames(atac_meta), as.character(seq_len(nrow(atac_meta))))) {
  rownames(atac_meta) <- atac_meta$barcode
}
atac_cells <- rownames(atac_meta)[rownames(atac_meta) %in% colnames(ga)]
ga <- ga[, atac_cells, drop = FALSE]
cat("  ATAC Gene Activity: ", paste(dim(ga), collapse = " x "),
    " (gene x ATAC cell)\n", sep = "")

cat("[4/6] Loading extracted CITE RNA counts ...\n")
cite_counts <- Matrix::readMM(cite_files$mtx)
cite_counts <- as(cite_counts, "CsparseMatrix")
cite_genes <- readLines(cite_files$genes, warn = FALSE)
cite_cells <- readLines(cite_files$cells, warn = FALSE)
rownames(cite_counts) <- cite_genes
colnames(cite_counts) <- cite_cells
cat("  CITE counts: ", paste(dim(cite_counts), collapse = " x "),
    " (gene x RNA cell)\n", sep = "")

cat("[5/6] Loading CITE library sizes and scDRS scores ...\n")
cite_scores <- as.data.frame(readRDS(CITE_SCORE_RDS))
lib <- cite_scores$scdrs_input_total_counts[match(cite_cells, cite_scores$cell)]
if (any(!is.finite(lib)) || any(is.na(lib))) {
  warning("Some CITE library sizes missing; falling back to extracted-gene colSums for those cells.")
  fallback <- Matrix::colSums(cite_counts)
  bad <- !is.finite(lib) | is.na(lib)
  lib[bad] <- fallback[bad]
}
lib[lib < 1] <- 1
names(lib) <- cite_cells

cat("[6/6] Running sparse mapping per trait ...\n")
summary_rows <- list()
feature_summary_rows <- list()

for (short in names(TRAITS)) {
  trait_id <- TRAITS[[short]]
  cat("\n=== Trait ", short, " / ", trait_id, " ===\n", sep = "")

  g_full <- gs[[trait_id]]$genes
  w_full <- gs[[trait_id]]$weights
  g_use <- Reduce(intersect, list(g_full, rownames(ga), rownames(cite_counts)))
  if (length(g_use) < 50L) {
    stop("Too few shared genes for ", short, ": ", length(g_use), call. = FALSE)
  }
  w_g <- w_full[match(g_use, g_full)]
  cat("  shared G_trait genes: ", length(g_use), " / ", length(g_full), "\n", sep = "")

  cat("  building RNA sparse trait vectors ...\n")
  rna_sub <- cite_counts[g_use, , drop = FALSE]
  rna_norm <- Matrix::t(Matrix::t(rna_sub) / lib[colnames(rna_sub)]) * CPM_SCALE
  rm(rna_sub)
  gc_quiet()
  rna_log <- log1p(rna_norm)
  rm(rna_norm)
  gc_quiet()
  rna_sparse <- build_sparse_trait_features(
    rna_log,
    gene_weight = w_g,
    b_null = B_NULL_CELL,
    p_threshold = RNA_EMPIRICAL_P_THRESHOLD,
    min_delta = MIN_DELTA_OVER_NULL_MEDIAN,
    use_weight = USE_GWAS_WEIGHT,
    modality = "RNA",
    trait = short
  )
  rm(rna_log)
  gc_quiet()

  cat("  building ATAC sparse trait vectors ...\n")
  atac_sub <- ga[g_use, , drop = FALSE]
  atac_sparse <- build_sparse_trait_features(
    atac_sub,
    gene_weight = w_g,
    b_null = B_NULL_CELL,
    p_threshold = ATAC_EMPIRICAL_P_THRESHOLD,
    min_delta = MIN_DELTA_OVER_NULL_MEDIAN,
    use_weight = USE_GWAS_WEIGHT,
    modality = "ATAC",
    trait = short
  )
  rm(atac_sub)
  gc_quiet()

  cat("  searching sparse top-", K_NEIGHBOR, " CITE RNA neighbors ...\n", sep = "")
  t0 <- Sys.time()
  nn <- batched_top_k_sparse_cosine(
    query_g_by_c = atac_sparse$feature,
    ref_g_by_c = rna_sparse$feature,
    k = K_NEIGHBOR,
    batch = BATCH_SIZE
  )
  seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  sm <- softmax_per_cell(nn$nn.dists)

  neighbor_barcodes <- matrix(
    cite_cells[as.integer(nn$nn.idx)],
    nrow = nrow(nn$nn.idx),
    ncol = K_NEIGHBOR
  )

  score_col <- paste0(short, "_norm_score")
  if (score_col %in% names(cite_scores)) {
    score <- cite_scores[[score_col]][match(cite_cells, cite_scores$cell)]
    nn_score <- matrix(score[as.integer(nn$nn.idx)], nrow = nrow(nn$nn.idx), ncol = K_NEIGHBOR)
    pseudo_z <- rowSums(sm$weight * nn_score, na.rm = TRUE)
  } else {
    pseudo_z <- rep(NA_real_, length(atac_cells))
  }
  names(pseudo_z) <- atac_cells

  max_weight <- apply(sm$weight, 1L, max)
  entropy <- -rowSums(sm$weight * log(pmax(sm$weight, 1e-12)))
  query_nnz <- Matrix::colSums(atac_sparse$feature != 0)
  top1_cosine <- 1 - nn$nn.dists[, 1]^2 / 2

  cell_diag <- data.table(
    trait = short,
    barcode = atac_cells,
    atac_sparse_gene_count = as.integer(query_nnz),
    atac_sparse_empty = as.integer(query_nnz) == 0L,
    top1_distance = nn$nn.dists[, 1],
    top1_cosine = top1_cosine,
    max_neighbor_weight = max_weight,
    neighbor_weight_entropy = entropy,
    pseudo_z_rna = pseudo_z,
    no_celltype_used = TRUE
  )

  if (short %in% names(ATAC_TARGETS) && "cell_type" %in% names(atac_meta)) {
    target <- as.integer(atac_meta[atac_cells, "cell_type"] %in% ATAC_TARGETS[[short]])
    top5 <- top_metrics(pseudo_z, target, 0.05)
    posthoc <- data.table(
      trait = short,
      posthoc_target = paste(ATAC_TARGETS[[short]], collapse = "/"),
      n_atac_cells = length(target),
      n_target_cells = sum(target == 1L),
      target_rate = mean(target == 1L),
      pseudo_z_target_auroc = compute_auc_rank(pseudo_z, target),
      pseudo_z_top5pct_precision = unname(top5["precision"]),
      pseudo_z_top5pct_enrichment = unname(top5["enrichment"]),
      no_celltype_used_in_mapping = TRUE
    )
  } else {
    posthoc <- data.table(
      trait = short,
      posthoc_target = NA_character_,
      n_atac_cells = length(atac_cells),
      n_target_cells = NA_integer_,
      target_rate = NA_real_,
      pseudo_z_target_auroc = NA_real_,
      pseudo_z_top5pct_precision = NA_real_,
      pseudo_z_top5pct_enrichment = NA_real_,
      no_celltype_used_in_mapping = TRUE
    )
  }

  mapping <- list(
    trait = short,
    trait_id = trait_id,
    method = "sparse_trait_mapping_topk_softmax",
    reference = "NeurIPS2021_BMMC_CITE_GEX",
    k = K_NEIGHBOR,
    use_gwas_weight = USE_GWAS_WEIGHT,
    g_trait_used = g_use,
    w_g = w_g,
    atac_barcodes = atac_cells,
    rna_barcodes = cite_cells,
    rna_neighbor_idx = nn$nn.idx,
    rna_neighbor_barcode = neighbor_barcodes,
    neighbor_distance = nn$nn.dists,
    neighbor_weight = sm$weight,
    tau_per_cell = sm$tau,
    sparse_mapping = list(
      b_null_cell = B_NULL_CELL,
      rna_empirical_p_threshold = RNA_EMPIRICAL_P_THRESHOLD,
      atac_empirical_p_threshold = ATAC_EMPIRICAL_P_THRESHOLD,
      min_delta_over_null_median = MIN_DELTA_OVER_NULL_MEDIAN,
      rna_empty_cell_rate = mean(rna_sparse$cell_diag$empty_sparse_vector),
      atac_empty_cell_rate = mean(atac_sparse$cell_diag$empty_sparse_vector),
      rna_nnz = length(rna_sparse$feature@x),
      atac_nnz = length(atac_sparse$feature@x)
    ),
    input_paths = list(
      gene_activity = GA_PATH,
      cite_h5ad = CITE_H5AD,
      cite_extract_dir = CITE_EXTRACT_DIR,
      gene_set = GS_FILE,
      cite_scores = CITE_SCORE_RDS
    ),
    created_at = as.character(Sys.time())
  )

  out_file <- file.path(SPARSE_MAP_OUT_DIR, paste0("granja_cite_mapping_", short, ".rds"))
  saveRDS(mapping, out_file)

  fwrite(cell_diag, file.path(DIAG_DIR, paste0("sparse_mapping_cell_diagnostics_", short, ".csv.gz")))
  fwrite(posthoc, file.path(DIAG_DIR, paste0("sparse_mapping_posthoc_summary_", short, ".csv")))
  fwrite(rna_sparse$cell_diag, file.path(DIAG_DIR, paste0("sparse_feature_cell_diagnostics_RNA_", short, ".csv.gz")))
  fwrite(atac_sparse$cell_diag, file.path(DIAG_DIR, paste0("sparse_feature_cell_diagnostics_ATAC_", short, ".csv.gz")))
  fwrite(rna_sparse$gene_diag, file.path(DIAG_DIR, paste0("sparse_feature_gene_diagnostics_RNA_", short, ".csv")))
  fwrite(atac_sparse$gene_diag, file.path(DIAG_DIR, paste0("sparse_feature_gene_diagnostics_ATAC_", short, ".csv")))

  rna_stats <- vector_stats(rna_sparse$cell_diag$retained_gene_count)
  atac_stats <- vector_stats(atac_sparse$cell_diag$retained_gene_count)
  mapping_summary <- data.table(
    trait = short,
    trait_id = trait_id,
    n_atac_cells = length(atac_cells),
    n_cite_cells = length(cite_cells),
    n_g_trait_full = length(g_full),
    n_g_trait_used = length(g_use),
    k = K_NEIGHBOR,
    b_null_cell = B_NULL_CELL,
    rna_empirical_p_threshold = RNA_EMPIRICAL_P_THRESHOLD,
    atac_empirical_p_threshold = ATAC_EMPIRICAL_P_THRESHOLD,
    rna_retained_gene_median = rna_stats[["median"]],
    rna_retained_gene_q95 = rna_stats[["q95"]],
    rna_empty_cell_rate = mean(rna_sparse$cell_diag$empty_sparse_vector),
    atac_retained_gene_median = atac_stats[["median"]],
    atac_retained_gene_q95 = atac_stats[["q95"]],
    atac_empty_cell_rate = mean(atac_sparse$cell_diag$empty_sparse_vector),
    top1_distance_median = median(nn$nn.dists[, 1]),
    top1_cosine_median = median(top1_cosine),
    max_neighbor_weight_median = median(max_weight),
    neighbor_weight_entropy_median = median(entropy),
    pseudo_z_target_auroc_posthoc = posthoc$pseudo_z_target_auroc[1],
    pseudo_z_top5pct_precision_posthoc = posthoc$pseudo_z_top5pct_precision[1],
    seconds = seconds,
    no_celltype_used_in_mapping = TRUE
  )
  summary_rows[[short]] <- mapping_summary

  feature_summary_rows[[paste(short, "RNA", sep = "_")]] <- data.table(
    trait = short,
    modality = "RNA",
    retained_gene_median = rna_stats[["median"]],
    retained_gene_q95 = rna_stats[["q95"]],
    empty_cell_rate = mean(rna_sparse$cell_diag$empty_sparse_vector),
    nnz = length(rna_sparse$feature@x),
    density = length(rna_sparse$feature@x) / max(prod(dim(rna_sparse$feature)), 1)
  )
  feature_summary_rows[[paste(short, "ATAC", sep = "_")]] <- data.table(
    trait = short,
    modality = "ATAC",
    retained_gene_median = atac_stats[["median"]],
    retained_gene_q95 = atac_stats[["q95"]],
    empty_cell_rate = mean(atac_sparse$cell_diag$empty_sparse_vector),
    nnz = length(atac_sparse$feature@x),
    density = length(atac_sparse$feature@x) / max(prod(dim(atac_sparse$feature)), 1)
  )

  cat("  saved mapping: ", out_file, "\n", sep = "")
  print(mapping_summary)

  rm(rna_sparse, atac_sparse, nn, sm, mapping, cell_diag, posthoc)
  invisible(gc(verbose = FALSE))
}

summary_dt <- rbindlist(summary_rows, use.names = TRUE, fill = TRUE)
feature_summary_dt <- rbindlist(feature_summary_rows, use.names = TRUE, fill = TRUE)

fwrite(summary_dt, file.path(SPARSE_MAP_OUT_DIR, "sparse_trait_mapping_summary.csv"))
fwrite(summary_dt, file.path(DIAG_DIR, "sparse_trait_mapping_summary.csv"))
fwrite(feature_summary_dt, file.path(DIAG_DIR, "sparse_trait_feature_summary.csv"))

cat("\n=== Sparse trait mapping finished ===\n")
print(summary_dt)
cat("\nOutputs are in: ", SPARSE_MAP_OUT_DIR, "\n", sep = "")
cat("Diagnostics are in: ", DIAG_DIR, "\n", sep = "")
