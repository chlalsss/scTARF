# module4 / Phase 1.3: trait-aware alignment between Granja ATAC and CITE RNA
# -----------------------------------------------------------------------------
# 目的:
#   以 ATAC 细胞为主, CITE RNA 细胞为辅助.
#   对每个目标性状分别在该性状的 GWAS/SCDRS 基因集合 G_trait 上做相似细胞匹配.
#
# 关键原则:
#   1. 对齐只使用 Granja ATAC Gene Activity 和 CITE RNA expression.
#   2. 每个 trait 单独对齐, 不用全基因整体相似度.
#   3. 对齐过程完全不读、不用 cell_type 标签.
#   4. cell_type 只在 02_check_alignment_by_celltype.R 里用于事后检查.
#
# 运行方式:
#   在 RStudio 里把工作目录设到本项目根目录, 然后运行:
#     source("02_trait_mapping/01_align_granja_cite_by_trait.R")
#
# 可选:
#   只跑一个 trait:
#     TRAITS_TO_RUN <- c("MONO")
#     source("02_trait_mapping/01_align_granja_cite_by_trait.R")
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
})

ROOT <- "D:/RStudio/work/module4_traitaware_pseudopair"
if (!exists("ALIGN_OUT_DIR")) ALIGN_OUT_DIR <- file.path(ROOT, "02_trait_mapping/output")
OUT_DIR <- ALIGN_OUT_DIR
EXTRACT_DIR <- file.path(OUT_DIR, "cite_extract")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(EXTRACT_DIR, recursive = TRUE, showWarnings = FALSE)

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
if (!exists("TRAITS_TO_RUN")) TRAITS_TO_RUN <- c("MONO", "PLT", "LYM")
if (!exists("K_NEIGHBOR")) K_NEIGHBOR <- 5L
if (!exists("BATCH_SIZE")) BATCH_SIZE <- 512L
if (!exists("USE_GWAS_WEIGHT")) USE_GWAS_WEIGHT <- TRUE
if (!exists("FORCE_REEXTRACT_CITE")) FORCE_REEXTRACT_CITE <- FALSE

if (!exists("TRAITS_ALL")) {
  TRAITS_ALL <- c(
    MONO = "UKB_460K.blood_MONOCYTE_COUNT",
    PLT  = "UKB_460K.blood_PLATELET_COUNT",
    LYM  = "UKB_460K.blood_LYMPHOCYTE_COUNT"
  )
}
TRAITS <- TRAITS_ALL[TRAITS_TO_RUN]
if (anyNA(TRAITS) || length(TRAITS) == 0L) {
  stop("TRAITS_TO_RUN contains trait short names not present in TRAITS_ALL: ",
       paste(setdiff(TRAITS_TO_RUN, names(TRAITS_ALL)), collapse = ", "),
       call. = FALSE)
}
EPS <- 1e-8

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

build_trait_block <- function(mat_g_by_c, w_g, use_weight = TRUE) {
  mu <- Matrix::rowMeans(mat_g_by_c)
  m2 <- Matrix::rowMeans(mat_g_by_c * mat_g_by_c)
  sd <- sqrt(pmax(m2 - mu * mu, 0))
  sd[sd < EPS] <- EPS

  z <- as.matrix(mat_g_by_c)
  z <- (z - mu) / sd
  if (isTRUE(use_weight)) z <- z * w_g

  x <- t(z)
  norms <- sqrt(rowSums(x * x))
  norms[norms < EPS] <- EPS
  x / norms
}

batched_top_k_cosine <- function(x_query, x_ref, k, batch = 512L) {
  n_q <- nrow(x_query)
  n_batch <- ceiling(n_q / batch)
  x_ref_t <- t(x_ref)
  idx_all <- matrix(0L, n_q, k)
  dist_all <- matrix(0, n_q, k)
  t0 <- Sys.time()

  for (bi in seq_len(n_batch)) {
    start <- (bi - 1L) * batch + 1L
    end <- min(bi * batch, n_q)
    sim <- x_query[start:end, , drop = FALSE] %*% x_ref_t
    sim <- as.matrix(sim)
    nb <- nrow(sim)
    row_ix <- seq_len(nb)

    for (j in seq_len(k)) {
      mx <- max.col(sim, ties.method = "first")
      val <- sim[cbind(row_ix, mx)]
      val <- pmin(pmax(val, -1), 1)
      idx_all[start:end, j] <- mx
      dist_all[start:end, j] <- sqrt(pmax(2 - 2 * val, 0))
      sim[cbind(row_ix, mx)] <- -Inf
    }

    if (bi %% 5L == 0L || bi == n_batch) {
      elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      cat(sprintf("       batch %d/%d finished (%d/%d ATAC cells, %.1fs)\n",
                  bi, n_batch, end, n_q, elapsed))
      flush.console()
    }
    rm(sim)
    invisible(gc(verbose = FALSE))
  }

  list(nn.idx = idx_all, nn.dists = dist_all)
}

softmax_per_cell <- function(d_mat) {
  tau <- apply(d_mat, 1, function(x) max(stats::median(x), EPS))
  w <- exp(-d_mat / tau)
  w <- w / rowSums(w)
  list(weight = w, tau = tau)
}

stop_missing <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0L) {
    stop("Missing required file(s):\n  ", paste(missing, collapse = "\n  "), call. = FALSE)
  }
}

extract_cite_counts_if_needed <- function(all_requested_genes) {
  dir.create(EXTRACT_DIR, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(EXTRACT_DIR)) {
    stop("Could not create CITE extraction directory: ", EXTRACT_DIR, call. = FALSE)
  }

  requested_file <- file.path(EXTRACT_DIR, "requested_trait_genes.txt")
  out_mtx <- file.path(EXTRACT_DIR, "cite_trait_counts.mtx")
  out_genes <- file.path(EXTRACT_DIR, "cite_trait_counts_genes.txt")
  out_cells <- file.path(EXTRACT_DIR, "cite_trait_counts_cells.txt")
  out_summary <- file.path(EXTRACT_DIR, "cite_trait_counts_summary.csv")
  out_missing <- file.path(EXTRACT_DIR, "cite_trait_counts_missing_genes.txt")

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

cat(">>> Phase 1.3 Granja ATAC x CITE RNA trait-aware alignment\n\n")
cat("Traits: ", paste(names(TRAITS), collapse = ", "), "\n", sep = "")
cat("K_NEIGHBOR = ", K_NEIGHBOR, "; BATCH_SIZE = ", BATCH_SIZE,
    "; USE_GWAS_WEIGHT = ", USE_GWAS_WEIGHT, "\n\n", sep = "")

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
if (is.null(rownames(atac_meta)) || identical(rownames(atac_meta), as.character(seq_len(nrow(atac_meta))))) {
  rownames(atac_meta) <- atac_meta$barcode
}
atac_cells <- intersect(rownames(atac_meta), colnames(ga))
atac_cells <- rownames(atac_meta)[rownames(atac_meta) %in% atac_cells]
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

cat("[5/6] Loading CITE library sizes for RNA normalization ...\n")
cite_scores <- readRDS(CITE_SCORE_RDS)
cite_scores <- as.data.frame(cite_scores)
lib <- cite_scores$scdrs_input_total_counts[match(cite_cells, cite_scores$cell)]
if (any(!is.finite(lib)) || any(is.na(lib))) {
  warning("Some CITE library sizes missing; falling back to extracted-gene colSums for those cells.")
  fallback <- Matrix::colSums(cite_counts)
  bad <- !is.finite(lib) | is.na(lib)
  lib[bad] <- fallback[bad]
}
lib[lib < 1] <- 1
names(lib) <- cite_cells
rm(cite_scores)
invisible(gc(verbose = FALSE))

cat("[6/6] Running per-trait alignment ...\n")
summary_rows <- list()

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

  cat("  building ATAC trait-space block ...\n")
  x_atac <- build_trait_block(ga[g_use, , drop = FALSE], w_g, use_weight = USE_GWAS_WEIGHT)

  cat("  building CITE RNA trait-space block ...\n")
  rna_sub <- cite_counts[g_use, , drop = FALSE]
  rna_norm <- Matrix::t(Matrix::t(rna_sub) / lib[colnames(rna_sub)]) * 1e4
  rna_log <- log1p(rna_norm)
  x_rna <- build_trait_block(rna_log, w_g, use_weight = USE_GWAS_WEIGHT)

  cat("  searching top-", K_NEIGHBOR, " CITE RNA neighbors for each ATAC cell ...\n", sep = "")
  t0 <- Sys.time()
  nn <- batched_top_k_cosine(
    x_query = x_atac,
    x_ref = x_rna,
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

  mapping <- list(
    trait = short,
    trait_id = trait_id,
    method = "trait_gene_cosine_topk_softmax",
    reference = "NeurIPS2021_BMMC_CITE_GEX",
    k = K_NEIGHBOR,
    use_gwas_weight = USE_GWAS_WEIGHT,
    g_trait_used = g_use,
    w_g = w_g,
    atac_barcodes = colnames(ga),
    rna_barcodes = cite_cells,
    rna_neighbor_idx = nn$nn.idx,
    rna_neighbor_barcode = neighbor_barcodes,
    neighbor_distance = nn$nn.dists,
    neighbor_weight = sm$weight,
    tau_per_cell = sm$tau,
    input_paths = list(
      gene_activity = GA_PATH,
      cite_h5ad = CITE_H5AD,
      gene_set = GS_FILE,
      cite_scores = CITE_SCORE_RDS
    ),
    created_at = as.character(Sys.time())
  )

  out_file <- file.path(OUT_DIR, paste0("granja_cite_mapping_", short, ".rds"))
  saveRDS(mapping, out_file)
  cat("  saved: ", out_file, "\n", sep = "")

  d <- nn$nn.dists
  wsum <- rowSums(sm$weight)
  summary_rows[[short]] <- data.frame(
    trait = short,
    trait_id = trait_id,
    n_atac_cells = nrow(x_atac),
    n_cite_cells = nrow(x_rna),
    n_g_trait_full = length(g_full),
    n_g_trait_used = length(g_use),
    k = K_NEIGHBOR,
    use_gwas_weight = USE_GWAS_WEIGHT,
    euc_dist_median = median(d),
    cos_dist_median = median(d^2 / 2),
    cos_dist_q05 = as.numeric(quantile(d^2 / 2, 0.05)),
    cos_dist_q95 = as.numeric(quantile(d^2 / 2, 0.95)),
    tau_median = median(sm$tau),
    weight_rowsum_min = min(wsum),
    weight_rowsum_max = max(wsum),
    seconds = seconds,
    stringsAsFactors = FALSE
  )

  rm(x_atac, x_rna, rna_sub, rna_norm, rna_log, nn, sm, neighbor_barcodes, mapping)
  invisible(gc(verbose = FALSE))
}

summary_df <- do.call(rbind, summary_rows)
write.csv(summary_df, file.path(OUT_DIR, "granja_cite_mapping_summary.csv"), row.names = FALSE)

cat("\n=== Alignment finished ===\n")
print(summary_df, row.names = FALSE)
cat("\nOutputs are in: ", OUT_DIR, "\n", sep = "")
