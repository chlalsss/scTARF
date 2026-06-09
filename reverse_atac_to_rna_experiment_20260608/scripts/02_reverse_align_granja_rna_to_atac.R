# Reverse Step 1: trait-aware RNA -> ATAC neighbor mapping.
#
# Thin reverse counterpart of:
#   module4_traitaware_pseudopair/02_trait_mapping/01_align_granja_cite_by_trait.R
#
# Difference:
#   query = Granja RNA cells
#   ref   = Granja ATAC cells
#   each RNA cell gets top-k ATAC neighbors in the same trait-gene space.

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
})

ROOT <- "D:/RStudio/work"
MODULE4_ROOT <- file.path(ROOT, "module4_traitaware_pseudopair")
REV_ROOT <- file.path(ROOT, "reverse_atac_to_rna_experiment_20260608")

SCRIPT_ENV <- environment()
has_local <- function(x) exists(x, envir = SCRIPT_ENV, inherits = FALSE)

if (!has_local("TRAITS_TO_RUN")) TRAITS_TO_RUN <- c("PLT", "MCH", "RBC")
if (!has_local("TRAIT_IDS")) {
  TRAIT_IDS <- c(
    PLT = "UKB_460K.blood_PLATELET_COUNT",
    MCH = "UKB_460K.blood_MEAN_CORPUSCULAR_HEMOGLOBIN",
    RBC = "UKB_460K.blood_RED_COUNT"
  )
}
if (!has_local("K_NEIGHBOR")) K_NEIGHBOR <- 5L
if (!has_local("BATCH_SIZE")) BATCH_SIZE <- 512L
if (!has_local("USE_GWAS_WEIGHT")) USE_GWAS_WEIGHT <- TRUE

if (!has_local("REVERSE_TAG")) REVERSE_TAG <- "reverse_granja_rna_main_plt_mch_rbc"
if (!has_local("ALIGN_OUT_DIR")) {
  ALIGN_OUT_DIR <- file.path(REV_ROOT, "outputs", REVERSE_TAG, "mapping")
}
if (!has_local("RNA_EXTRACT_DIR")) RNA_EXTRACT_DIR <- file.path(ALIGN_OUT_DIR, "granja_rna_extract")

if (!has_local("GA_PATH")) {
  GA_PATH <- "D:/SCAVENGE-reproducibility-main/blood_data/multiomics_output_v2/gene_activity_matrix.rds"
}
if (!has_local("ATAC_META_PATH")) {
  ATAC_META_PATH <- file.path(MODULE4_ROOT, "04_evaluation/inputs/granja_cell_meta.rds")
}
if (!has_local("GRANJA_RNA_H5AD")) {
  GRANJA_RNA_H5AD <- file.path(
    MODULE4_ROOT,
    "01_data_prep/03_rna_scdrs_python/input/granja_rna_raw_counts.h5ad"
  )
}
if (!has_local("GRANJA_RNA_META_RDS")) {
  GRANJA_RNA_META_RDS <- file.path(MODULE4_ROOT, "01_data_prep/output/granja_rna_cell_meta.rds")
}
if (!has_local("GS_FILES")) {
  GS_FILES <- c(
    file.path(MODULE4_ROOT, "01_data_prep/03_rna_scdrs_python/granja_3traits.gs"),
    file.path(MODULE4_ROOT, "01_data_prep_external_rna/03_cite_scdrs_python/cite_trait_screening.gs")
  )
}
if (!has_local("PYTHON_EXE")) PYTHON_EXE <- "D:/python/python.exe"
if (!has_local("PY_HELPER")) {
  PY_HELPER <- file.path(MODULE4_ROOT, "02_trait_mapping/extract_cite_trait_counts.py")
}
if (!has_local("FORCE_REEXTRACT_RNA")) FORCE_REEXTRACT_RNA <- FALSE
if (!has_local("ALLOW_OUTSIDE_REV_ROOT")) ALLOW_OUTSIDE_REV_ROOT <- FALSE

EPS <- 1e-8

stop_missing <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0L) {
    stop("Missing required file(s):\n  ", paste(missing, collapse = "\n  "), call. = FALSE)
  }
}

assert_write_under_rev <- function(paths) {
  if (isTRUE(ALLOW_OUTSIDE_REV_ROOT)) return(invisible(TRUE))
  root <- normalizePath(REV_ROOT, winslash = "/", mustWork = TRUE)
  for (path in paths) {
    resolved <- normalizePath(path, winslash = "/", mustWork = FALSE)
    if (!startsWith(resolved, paste0(root, "/")) && !identical(resolved, root)) {
      stop("Refusing to write outside reverse workspace: ", path,
           "\nSet ALLOW_OUTSIDE_REV_ROOT <- TRUE only if this is intentional.",
           call. = FALSE)
    }
  }
  invisible(TRUE)
}

parse_gs_file <- function(path) {
  lines <- readLines(path, warn = FALSE)
  out <- list()
  for (line in lines[-1]) {
    parts <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(parts) < 2L) next
    pairs <- strsplit(strsplit(parts[2], ",", fixed = TRUE)[[1]], ":", fixed = TRUE)
    genes <- vapply(pairs, `[`, character(1), 1)
    weights <- as.numeric(vapply(
      pairs,
      function(x) if (length(x) >= 2L) x[2] else "1",
      character(1)
    ))
    out[[parts[1]]] <- list(genes = genes, weights = weights, gs_file = path)
  }
  out
}

parse_gs_files <- function(paths) {
  stop_missing(paths)
  out <- list()
  for (path in paths) {
    one <- parse_gs_file(path)
    for (nm in names(one)) {
      if (!nm %in% names(out)) out[[nm]] <- one[[nm]]
    }
  }
  out
}

build_trait_block <- function(mat_g_by_c, w_g, use_weight = TRUE) {
  mu <- Matrix::rowMeans(mat_g_by_c)
  m2 <- Matrix::rowMeans(mat_g_by_c * mat_g_by_c)
  sd <- sqrt(pmax(m2 - mu * mu, 0))
  sd[sd < EPS | !is.finite(sd)] <- EPS

  z <- as.matrix(mat_g_by_c)
  z <- (z - mu) / sd
  if (isTRUE(use_weight)) z <- z * w_g

  x <- t(z)
  norms <- sqrt(rowSums(x * x))
  norms[norms < EPS | !is.finite(norms)] <- EPS
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
    sim <- as.matrix(x_query[start:end, , drop = FALSE] %*% x_ref_t)
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
      cat(sprintf("       batch %d/%d finished (%d/%d RNA cells, %.1fs)\n",
                  bi, n_batch, end, n_q, elapsed))
      flush.console()
    }
    rm(sim)
    invisible(gc(verbose = FALSE))
  }

  list(nn.idx = idx_all, nn.dists = dist_all)
}

softmax_per_cell <- function(d_mat) {
  tau <- apply(d_mat, 1L, function(x) max(stats::median(x), EPS))
  w <- exp(-d_mat / tau)
  w <- w / rowSums(w)
  list(weight = w, tau = tau)
}

extract_rna_counts_if_needed <- function(all_requested_genes) {
  requested_file <- file.path(RNA_EXTRACT_DIR, "requested_trait_genes.txt")
  out_mtx <- file.path(RNA_EXTRACT_DIR, "cite_trait_counts.mtx")
  out_genes <- file.path(RNA_EXTRACT_DIR, "cite_trait_counts_genes.txt")
  out_cells <- file.path(RNA_EXTRACT_DIR, "cite_trait_counts_cells.txt")
  out_summary <- file.path(RNA_EXTRACT_DIR, "cite_trait_counts_summary.csv")
  out_missing <- file.path(RNA_EXTRACT_DIR, "cite_trait_counts_missing_genes.txt")

  writeLines(sort(unique(all_requested_genes)), requested_file, useBytes = TRUE)

  outputs <- c(out_mtx, out_genes, out_cells, out_summary)
  if (!FORCE_REEXTRACT_RNA && all(file.exists(outputs))) {
    present_genes <- readLines(out_genes, warn = FALSE)
    missing_genes <- if (file.exists(out_missing)) readLines(out_missing, warn = FALSE) else character()
    covered_genes <- unique(c(present_genes, missing_genes))
    if (all(sort(unique(all_requested_genes)) %in% covered_genes)) {
      cat("[extract] Reusing existing Granja RNA trait-count matrix.\n")
      return(list(mtx = out_mtx, genes = out_genes, cells = out_cells))
    }
  }

  cat("[extract] Extracting Granja RNA trait counts from h5ad via existing Python helper ...\n")
  stop_missing(c(PYTHON_EXE, PY_HELPER, GRANJA_RNA_H5AD))
  status <- system2(PYTHON_EXE, args = c(
    PY_HELPER,
    "--h5ad", GRANJA_RNA_H5AD,
    "--genes", requested_file,
    "--out-mtx", out_mtx,
    "--out-genes", out_genes,
    "--out-cells", out_cells,
    "--out-summary", out_summary
  ))
  if (!identical(status, 0L)) {
    stop("Python h5ad extraction failed with status ", status, call. = FALSE)
  }
  list(mtx = out_mtx, genes = out_genes, cells = out_cells)
}

load_meta_by_barcode <- function(path) {
  meta <- readRDS(path)
  if (inherits(meta, "data.table")) meta <- as.data.frame(meta)
  if (is.null(rownames(meta)) || identical(rownames(meta), as.character(seq_len(nrow(meta))))) {
    if ("barcode" %in% colnames(meta)) rownames(meta) <- as.character(meta$barcode)
  }
  meta
}

cat(">>> Reverse Step 1: Granja RNA -> Granja ATAC trait-aware mapping\n\n")
cat("Traits: ", paste(TRAITS_TO_RUN, collapse = ", "), "\n", sep = "")
cat("Output: ", ALIGN_OUT_DIR, "\n\n", sep = "")

stop_missing(c(GA_PATH, ATAC_META_PATH, GRANJA_RNA_H5AD, GRANJA_RNA_META_RDS, PYTHON_EXE, PY_HELPER))
assert_write_under_rev(c(ALIGN_OUT_DIR, RNA_EXTRACT_DIR))
dir.create(ALIGN_OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RNA_EXTRACT_DIR, recursive = TRUE, showWarnings = FALSE)

gs <- parse_gs_files(GS_FILES)
trait_ids <- TRAIT_IDS[TRAITS_TO_RUN]
if (anyNA(trait_ids)) {
  stop("TRAITS_TO_RUN contains names absent from TRAIT_IDS: ",
       paste(setdiff(TRAITS_TO_RUN, names(TRAIT_IDS)), collapse = ", "),
       call. = FALSE)
}
missing_gs <- trait_ids[!trait_ids %in% names(gs)]
if (length(missing_gs) > 0L) {
  stop("Missing gene-set rows for trait IDs:\n  ", paste(missing_gs, collapse = "\n  "),
       call. = FALSE)
}
all_trait_genes <- unique(unlist(lapply(unname(trait_ids), function(id) gs[[id]]$genes)))

cat("[1/5] Preparing Granja RNA trait-gene count matrix ...\n")
rna_files <- extract_rna_counts_if_needed(all_trait_genes)

cat("[2/5] Loading Granja ATAC Gene Activity ...\n")
ga <- readRDS(GA_PATH)
ga <- as(ga, "CsparseMatrix")
atac_meta <- load_meta_by_barcode(ATAC_META_PATH)
atac_cells <- rownames(atac_meta)[rownames(atac_meta) %in% colnames(ga)]
ga <- ga[, atac_cells, drop = FALSE]
cat("  GA: ", paste(dim(ga), collapse = " x "), " (gene x ATAC cell)\n", sep = "")

cat("[3/5] Loading Granja RNA extracted counts and metadata ...\n")
rna_counts <- Matrix::readMM(rna_files$mtx)
rna_counts <- as(rna_counts, "CsparseMatrix")
rna_genes <- readLines(rna_files$genes, warn = FALSE)
rna_cells <- readLines(rna_files$cells, warn = FALSE)
rownames(rna_counts) <- rna_genes
colnames(rna_counts) <- rna_cells

rna_meta <- load_meta_by_barcode(GRANJA_RNA_META_RDS)
rna_cells <- rna_cells[rna_cells %in% rownames(rna_meta)]
rna_counts <- rna_counts[, rna_cells, drop = FALSE]
lib_col <- intersect(c("scdrs_input_total_counts", "nUMI", "nUMI_pre"), colnames(rna_meta))[1]
if (is.na(lib_col)) {
  stop("Granja RNA metadata needs one of scdrs_input_total_counts/nUMI/nUMI_pre.", call. = FALSE)
}
rna_lib <- as.numeric(rna_meta[rna_cells, lib_col])
names(rna_lib) <- rna_cells
bad_lib <- !is.finite(rna_lib) | is.na(rna_lib) | rna_lib < 1
if (any(bad_lib)) rna_lib[bad_lib] <- Matrix::colSums(rna_counts)[bad_lib]
rna_lib[rna_lib < 1 | !is.finite(rna_lib)] <- 1
cat("  RNA counts: ", paste(dim(rna_counts), collapse = " x "), " (gene x RNA cell)\n", sep = "")

cat("[4/5] Running per-trait reverse alignment ...\n")
summary_rows <- list()

for (short in TRAITS_TO_RUN) {
  trait_id <- TRAIT_IDS[[short]]
  out_file <- file.path(ALIGN_OUT_DIR, paste0("granja_rna_to_atac_mapping_", short, ".rds"))
  if (file.exists(out_file) &&
      !isTRUE(get0("FORCE_REMAP", ifnotfound = FALSE, inherits = FALSE))) {
    cat("\n=== Trait ", short, " already mapped; skipping ===\n", sep = "")
    mapping <- readRDS(out_file)
    summary_rows[[short]] <- data.frame(
      trait = short,
      trait_id = trait_id,
      n_rna_cells = length(mapping$rna_barcodes),
      n_atac_cells = length(mapping$atac_barcodes),
      n_g_trait_used = length(mapping$g_trait_used),
      k = mapping$k,
      skipped_existing = TRUE,
      stringsAsFactors = FALSE
    )
    next
  }

  cat("\n=== Trait ", short, " / ", trait_id, " ===\n", sep = "")
  g_full <- gs[[trait_id]]$genes
  w_full <- gs[[trait_id]]$weights
  g_use <- Reduce(intersect, list(g_full, rownames(ga), rownames(rna_counts)))
  if (length(g_use) < 50L) {
    stop("Too few shared genes for ", short, ": ", length(g_use), call. = FALSE)
  }
  w_g <- w_full[match(g_use, g_full)]
  cat("  shared G_trait genes: ", length(g_use), " / ", length(g_full), "\n", sep = "")

  cat("  building RNA trait-space query block ...\n")
  rna_sub <- rna_counts[g_use, , drop = FALSE]
  rna_norm <- Matrix::t(Matrix::t(rna_sub) / rna_lib[colnames(rna_sub)]) * 1e4
  rna_log <- log1p(rna_norm)
  x_rna <- build_trait_block(rna_log, w_g, use_weight = USE_GWAS_WEIGHT)

  cat("  building ATAC Gene Activity trait-space reference block ...\n")
  x_atac <- build_trait_block(ga[g_use, , drop = FALSE], w_g, use_weight = USE_GWAS_WEIGHT)

  cat("  searching top-", K_NEIGHBOR, " ATAC neighbors for each RNA cell ...\n", sep = "")
  t0 <- Sys.time()
  nn <- batched_top_k_cosine(
    x_query = x_rna,
    x_ref = x_atac,
    k = K_NEIGHBOR,
    batch = BATCH_SIZE
  )
  seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  sm <- softmax_per_cell(nn$nn.dists)

  neighbor_barcodes <- matrix(
    colnames(ga)[as.integer(nn$nn.idx)],
    nrow = nrow(nn$nn.idx),
    ncol = K_NEIGHBOR
  )

  mapping <- list(
    trait = short,
    trait_id = trait_id,
    method = "reverse_trait_gene_cosine_topk_softmax",
    query = "Granja RNA",
    reference = "Granja ATAC Gene Activity",
    k = K_NEIGHBOR,
    use_gwas_weight = USE_GWAS_WEIGHT,
    g_trait_used = g_use,
    w_g = w_g,
    rna_barcodes = colnames(rna_counts),
    atac_barcodes = colnames(ga),
    atac_neighbor_idx = nn$nn.idx,
    atac_neighbor_barcode = neighbor_barcodes,
    neighbor_distance = nn$nn.dists,
    neighbor_weight = sm$weight,
    tau_per_cell = sm$tau,
    input_paths = list(
      gene_activity = GA_PATH,
      granja_rna_h5ad = GRANJA_RNA_H5AD,
      granja_rna_meta = GRANJA_RNA_META_RDS,
      gene_sets = GS_FILES
    ),
    no_celltype_used = TRUE,
    created_at = as.character(Sys.time())
  )
  saveRDS(mapping, out_file)
  cat("  saved: ", out_file, "\n", sep = "")

  d <- nn$nn.dists
  wsum <- rowSums(sm$weight)
  summary_rows[[short]] <- data.frame(
    trait = short,
    trait_id = trait_id,
    n_rna_cells = nrow(x_rna),
    n_atac_cells = nrow(x_atac),
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
    skipped_existing = FALSE,
    stringsAsFactors = FALSE
  )

  rm(x_rna, x_atac, rna_sub, rna_norm, rna_log, nn, sm, neighbor_barcodes, mapping)
  invisible(gc(verbose = FALSE))
}

cat("\n[5/5] Writing summary ...\n")
summary_df <- rbindlist(summary_rows, fill = TRUE)
fwrite(summary_df, file.path(ALIGN_OUT_DIR, "reverse_mapping_summary.csv"))

cat("\n=== Reverse alignment finished ===\n")
print(summary_df)
cat("\nOutputs are in: ", ALIGN_OUT_DIR, "\n", sep = "")
