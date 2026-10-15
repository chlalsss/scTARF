source("scripts/_paths.R")
# Reverse Step 3a: build Granja RNA main-space z block.
#
# This produces a module3-compatible "atac_z_block_v2.rds" override:
#   rows = Granja RNA cells
#   cols = RNA PCs, z-scored
#
# It reuses the existing h5ad extraction helper and the same sparse
# log-normalize + irlba PCA pattern used elsewhere in module4.

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(irlba)
  library(matrixStats)
})

if (!exists("ROOT")) ROOT <- Sys.getenv("SCTARF_WORK_ROOT", normalizePath(".", winslash = "/", mustWork = FALSE))
if (!exists("MODULE4_ROOT")) MODULE4_ROOT <- Sys.getenv("SCTARF_NONPAIRED_ROOT", file.path(ROOT, "outputs", "nonpaired"))
if (!exists("REV_ROOT")) REV_ROOT <- Sys.getenv("SCTARF_REVERSE_ROOT", file.path(ROOT, "outputs", "reverse"))
LOCAL_NONPAIRED_SCRIPT_DIR <- if (exists("SCTARF_NONPAIRED_SCRIPT_DIR")) SCTARF_NONPAIRED_SCRIPT_DIR else file.path(ROOT, "scripts", "03_nonpaired_analysis")

SCRIPT_ENV <- environment()
has_local <- function(x) exists(x, envir = SCRIPT_ENV, inherits = FALSE)

if (!has_local("REVERSE_TAG")) REVERSE_TAG <- "reverse_granja_rna_main_plt_mch_rbc"
if (!has_local("OUT_DIR")) OUT_DIR <- file.path(REV_ROOT, "outputs", REVERSE_TAG, "_inputs")
if (!has_local("RNA_MAIN_Z_RDS")) RNA_MAIN_Z_RDS <- file.path(OUT_DIR, "granja_rna_main_z_block_v2.rds")
if (!has_local("RNA_MAIN_EXTRACT_DIR")) RNA_MAIN_EXTRACT_DIR <- file.path(OUT_DIR, "granja_rna_main_extract")

if (!has_local("GRANJA_RNA_H5AD")) {
  GRANJA_RNA_H5AD <- file.path(
    MODULE4_ROOT,
    "01_data_prep/03_rna_scdrs_python/input/granja_rna_raw_counts.h5ad"
  )
}
if (!has_local("GRANJA_RNA_FEATURES_RDS")) {
  GRANJA_RNA_FEATURES_RDS <- file.path(MODULE4_ROOT, "01_data_prep/output/granja_rna_features.rds")
}
if (!has_local("GRANJA_RNA_META_RDS")) {
  GRANJA_RNA_META_RDS <- file.path(MODULE4_ROOT, "01_data_prep/output/granja_rna_cell_meta.rds")
}
if (!has_local("PYTHON_EXE")) PYTHON_EXE <- sctarf_python()
if (!has_local("PY_HELPER")) {
  PY_HELPER <- file.path(LOCAL_NONPAIRED_SCRIPT_DIR, "41_extract_cite_trait_counts.py")
}

if (!has_local("RNA_N_HVG")) RNA_N_HVG <- 3000L
if (!has_local("RNA_N_PCS")) RNA_N_PCS <- 30L
if (!has_local("CPM_SCALE")) CPM_SCALE <- 1e4
if (!has_local("FORCE_REBUILD_RNA_MAIN_Z")) FORCE_REBUILD_RNA_MAIN_Z <- FALSE
if (!has_local("FORCE_REEXTRACT_RNA_MAIN")) FORCE_REEXTRACT_RNA_MAIN <- FALSE
if (!has_local("ALLOW_OUTSIDE_REV_ROOT")) ALLOW_OUTSIDE_REV_ROOT <- FALSE

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

extract_feature_names <- function(x) {
  if (is.data.frame(x)) {
    for (nm in c("gene", "genes", "gene_name", "name", "feature", "features")) {
      if (nm %in% colnames(x)) return(as.character(x[[nm]]))
    }
  }
  if (is.vector(x)) return(as.character(x))
  stop("Could not extract gene names from GRANJA_RNA_FEATURES_RDS.", call. = FALSE)
}

load_meta_by_barcode <- function(path) {
  meta <- readRDS(path)
  if (inherits(meta, "data.table")) meta <- as.data.frame(meta)
  if (is.null(rownames(meta)) || identical(rownames(meta), as.character(seq_len(nrow(meta))))) {
    if ("barcode" %in% colnames(meta)) rownames(meta) <- as.character(meta$barcode)
  }
  meta
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

bad_rownames <- function(x) {
  is.null(rownames(x)) || anyNA(rownames(x)) || any(!nzchar(rownames(x)))
}

repair_existing_z_block <- function(path) {
  z <- readRDS(path)
  if (!bad_rownames(z)) return(TRUE)

  cells_file <- file.path(RNA_MAIN_EXTRACT_DIR, "cite_trait_counts_cells.txt")
  cells <- character()
  if (file.exists(cells_file)) {
    cells <- readLines(cells_file, warn = FALSE)
    if (file.exists(GRANJA_RNA_META_RDS)) {
      meta <- load_meta_by_barcode(GRANJA_RNA_META_RDS)
      cells <- cells[cells %in% rownames(meta)]
    }
  }

  if (length(cells) != nrow(z)) {
    transfer_dir <- file.path(REV_ROOT, "outputs", REVERSE_TAG, "transfer")
    profile_files <- list.files(
      transfer_dir,
      pattern = "^granja_pseudo_rna_profile_.*\\.rds$",
      full.names = TRUE
    )
    if (length(profile_files) > 0L) {
      obj <- readRDS(profile_files[[1]])
      cells <- obj$reference_cells
      if (is.null(cells)) cells <- obj$rna_barcodes
      if (is.null(cells)) cells <- colnames(obj$matrix)
    }
  }

  if (length(cells) != nrow(z)) return(FALSE)

  rownames(z) <- cells
  saveRDS(z, path)
  cat("Existing RNA main z-block had missing rownames; repaired from extracted RNA cell list.\n")
  TRUE
}

extract_all_rna_counts_if_needed <- function(genes) {
  requested_file <- file.path(RNA_MAIN_EXTRACT_DIR, "requested_main_genes.txt")
  out_mtx <- file.path(RNA_MAIN_EXTRACT_DIR, "cite_trait_counts.mtx")
  out_genes <- file.path(RNA_MAIN_EXTRACT_DIR, "cite_trait_counts_genes.txt")
  out_cells <- file.path(RNA_MAIN_EXTRACT_DIR, "cite_trait_counts_cells.txt")
  out_summary <- file.path(RNA_MAIN_EXTRACT_DIR, "cite_trait_counts_summary.csv")
  out_missing <- file.path(RNA_MAIN_EXTRACT_DIR, "cite_trait_counts_missing_genes.txt")
  writeLines(sort(unique(genes)), requested_file, useBytes = TRUE)

  outputs <- c(out_mtx, out_genes, out_cells, out_summary)
  if (!FORCE_REEXTRACT_RNA_MAIN && all(file.exists(outputs))) {
    present_genes <- readLines(out_genes, warn = FALSE)
    missing_genes <- if (file.exists(out_missing)) readLines(out_missing, warn = FALSE) else character()
    covered <- unique(c(present_genes, missing_genes))
    if (all(sort(unique(genes)) %in% covered)) {
      cat("[extract] Reusing existing Granja RNA main count matrix.\n")
      return(list(mtx = out_mtx, genes = out_genes, cells = out_cells))
    }
  }

  cat("[extract] Extracting Granja RNA all-gene count matrix via existing Python helper ...\n")
  stop_missing(c(PYTHON_EXE, PY_HELPER, GRANJA_RNA_H5AD))
  status <- system2(PYTHON_EXE, args = shQuote(c(
    PY_HELPER,
    "--h5ad", GRANJA_RNA_H5AD,
    "--genes", requested_file,
    "--out-mtx", out_mtx,
    "--out-genes", out_genes,
    "--out-cells", out_cells,
    "--out-summary", out_summary
  )))
  if (!identical(status, 0L)) {
    stop("Python h5ad extraction failed with status ", status, call. = FALSE)
  }
  list(mtx = out_mtx, genes = out_genes, cells = out_cells)
}

cat(">>> Reverse Step 3a: build Granja RNA main z-block\n\n")
cat("Output: ", RNA_MAIN_Z_RDS, "\n\n", sep = "")
assert_write_under_rev(c(OUT_DIR, RNA_MAIN_EXTRACT_DIR, RNA_MAIN_Z_RDS))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RNA_MAIN_EXTRACT_DIR, recursive = TRUE, showWarnings = FALSE)

if (file.exists(RNA_MAIN_Z_RDS) && !FORCE_REBUILD_RNA_MAIN_Z) {
  if (!repair_existing_z_block(RNA_MAIN_Z_RDS)) {
    stop(
      "Existing RNA main z-block is present but rownames are missing and cannot be repaired.\n",
      "Set FORCE_REBUILD_RNA_MAIN_Z <- TRUE and rerun RUN_STEPS <- c('rna_main_z').",
      call. = FALSE
    )
  }
  cat("Existing RNA main z-block found; skipping rebuild:\n  ", RNA_MAIN_Z_RDS, "\n", sep = "")
} else {
  stop_missing(c(GRANJA_RNA_H5AD, GRANJA_RNA_FEATURES_RDS, GRANJA_RNA_META_RDS, PYTHON_EXE, PY_HELPER))

  cat("[1/5] Reading Granja RNA gene list ...\n")
  genes <- extract_feature_names(readRDS(GRANJA_RNA_FEATURES_RDS))
  genes <- genes[!is.na(genes) & nzchar(genes)]
  cat("  genes requested: ", length(unique(genes)), "\n", sep = "")

  cat("[2/5] Extracting/loading RNA counts ...\n")
  rna_files <- extract_all_rna_counts_if_needed(genes)
  counts <- Matrix::readMM(rna_files$mtx)
  counts <- as(counts, "CsparseMatrix")
  rownames(counts) <- readLines(rna_files$genes, warn = FALSE)
  colnames(counts) <- readLines(rna_files$cells, warn = FALSE)

  meta <- load_meta_by_barcode(GRANJA_RNA_META_RDS)
  cells <- colnames(counts)[colnames(counts) %in% rownames(meta)]
  counts <- counts[, cells, drop = FALSE]
  lib_col <- intersect(c("scdrs_input_total_counts", "nUMI", "nUMI_pre"), colnames(meta))[1]
  if (is.na(lib_col)) {
    stop("Granja RNA metadata needs one of scdrs_input_total_counts/nUMI/nUMI_pre.", call. = FALSE)
  }
  lib <- as.numeric(meta[cells, lib_col])
  names(lib) <- cells
  bad <- !is.finite(lib) | is.na(lib) | lib < 1
  if (any(bad)) lib[bad] <- Matrix::colSums(counts)[bad]
  lib[lib < 1 | !is.finite(lib)] <- 1
  cat("  counts: ", paste(dim(counts), collapse = " x "), " (gene x RNA cell)\n", sep = "")

  cat("[3/5] Log-normalizing and selecting HVGs ...\n")
  norm <- Matrix::t(Matrix::t(counts) / lib[colnames(counts)]) * CPM_SCALE
  log_norm <- log1p(norm)
  gene_mean <- Matrix::rowMeans(log_norm)
  gene_m2 <- Matrix::rowMeans(log_norm * log_norm)
  gene_var <- pmax(gene_m2 - gene_mean * gene_mean, 0)
  gene_var[!is.finite(gene_var)] <- 0
  n_hvg <- min(RNA_N_HVG, length(gene_var))
  hvg <- names(sort(gene_var, decreasing = TRUE))[seq_len(n_hvg)]
  cat("  HVGs selected: ", length(hvg), "\n", sep = "")

  cat("[4/5] Running sparse PCA with irlba ...\n")
  x <- Matrix::t(log_norm[hvg, , drop = FALSE])
  pc_n <- min(RNA_N_PCS, nrow(x) - 1L, ncol(x) - 1L)
  if (pc_n < 2L) stop("Too few cells/features for PCA.", call. = FALSE)
  center <- Matrix::colMeans(x)
  m2 <- Matrix::colMeans(x * x)
  scale <- sqrt(pmax(m2 - center * center, 0))
  scale[!is.finite(scale) | scale < 1e-8] <- 1
  pca <- irlba::prcomp_irlba(x, n = pc_n, center = center, scale. = scale)
  z <- zscore_columns(pca$x)
  rownames(z) <- rownames(x)
  colnames(z) <- paste0("RNA_PC_", seq_len(ncol(z)))
  if (is.null(rownames(z)) || anyNA(rownames(z)) || any(!nzchar(rownames(z)))) {
    stop("RNA PCA z-block rownames are missing after PCA. Check extracted RNA cell names.",
         call. = FALSE)
  }

  cat("[5/5] Saving z-block and summary ...\n")
  saveRDS(z, RNA_MAIN_Z_RDS)
  summary <- data.frame(
    n_cells = nrow(z),
    n_genes_extracted = nrow(counts),
    n_hvg = length(hvg),
    n_pcs = ncol(z),
    source_h5ad = GRANJA_RNA_H5AD,
    output_rds = RNA_MAIN_Z_RDS,
    stringsAsFactors = FALSE
  )
  fwrite(summary, file.path(OUT_DIR, "granja_rna_main_z_block_summary.csv"))
  writeLines(hvg, file.path(OUT_DIR, "granja_rna_main_hvg_genes.txt"), useBytes = TRUE)

  cat("\n=== RNA main z-block ready ===\n")
  print(summary)
}
