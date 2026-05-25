# Pure-ATAC original SCAVENGE baseline robustness under raw ATAC read drop-out.
#
# Pipeline per p_drop x seed:
#   raw peak-by-cell ATAC counts
#     -> binomial read drop-out
#     -> TF-IDF + LSI
#     -> SCAVENGE mutual KNN
#     -> gchromVAR trait z-score on the degraded counts
#     -> SCAVENGE seed cells from ATAC z-score
#     -> SCAVENGE randomWalk_sparse + TRS post-processing
#     -> optional AUROC against configured cell-type targets
#
# Run:
#   setwd("D:/RStudio/work/scavenge_atac_baseline_robustness")
#   source("01_run_scavenge_atac_dropout.R")

this_file <- tryCatch(normalizePath(sys.frame(1)$ofile, winslash = "/"),
                      error = function(e) "")
candidate_dirs <- character()
if (nzchar(this_file)) candidate_dirs <- c(candidate_dirs, dirname(this_file))
if (exists("ROBUST_SCRIPT_DIR")) candidate_dirs <- c(candidate_dirs, ROBUST_SCRIPT_DIR)
candidate_dirs <- c(
  candidate_dirs,
  "D:/RStudio/work/scavenge_atac_baseline_robustness",
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
)
candidate_dirs <- unique(normalizePath(candidate_dirs, winslash = "/", mustWork = FALSE))
config_hits <- candidate_dirs[file.exists(file.path(candidate_dirs, "00_config.R"))]
if (length(config_hits) == 0L) {
  stop("Cannot locate 00_config.R. Candidate dirs:\n  ",
       paste(candidate_dirs, collapse = "\n  "), call. = FALSE)
}
ROBUST_SCRIPT_DIR <- config_hits[1L]
source(file.path(ROBUST_SCRIPT_DIR, "00_config.R"), local = environment())

required_pkgs <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) {
    stop("Missing R package(s): ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
}

required_pkgs(c(
  "Matrix", "SummarizedExperiment", "S4Vectors", "GenomicRanges",
  "chromVAR", "gchromVAR", "SCAVENGE", "data.table", "dplyr", "pROC"
))

suppressPackageStartupMessages({
  library(Matrix)
  library(SummarizedExperiment)
  library(data.table)
})

stop_missing <- function(paths) {
  paths <- paths[nzchar(paths)]
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0L) {
    stop("Missing required file(s):\n  ", paste(missing, collapse = "\n  "),
         call. = FALSE)
  }
}

format_p <- function(p_drop) sprintf("%.2f", round(as.numeric(p_drop), 2))

cache_seed_for <- function(p_drop, seed) {
  if (isTRUE(all.equal(as.numeric(p_drop), 0))) 0L else as.integer(seed)
}

run_key <- function(p_drop, seed) {
  paste0("p", format_p(p_drop), "_seed", as.integer(seed))
}

has_gc_bias <- function(SE) {
  rd_names <- names(SummarizedExperiment::rowData(SE))
  any(rd_names %in% c("GC", "bias", "gc_bias"))
}

load_genome_object <- function() {
  required_pkgs(GENOME_PACKAGE)
  getExportedValue(GENOME_PACKAGE, GENOME_OBJECT)
}

counts_assay <- function(SE) {
  assay_names <- SummarizedExperiment::assayNames(SE)
  if (!ATAC_ASSAY %in% assay_names) {
    stop("ATAC assay '", ATAC_ASSAY, "' not found. Available assays: ",
         paste(assay_names, collapse = ", "), call. = FALSE)
  }
  counts <- SummarizedExperiment::assay(SE, ATAC_ASSAY)
  if (!inherits(counts, "dgCMatrix")) counts <- Matrix::Matrix(counts, sparse = TRUE)
  counts
}

replace_counts_assay <- function(SE, counts) {
  SummarizedExperiment::assay(SE, ATAC_ASSAY) <- counts
  SE
}

ensure_gc_se <- function(SE) {
  if (has_gc_bias(SE)) {
    cat("  SE already has a GC/bias rowData column; using it as-is.\n")
    return(SE)
  }

  gc_file <- file.path(CACHE_DIR, "SE_atac_gcbiased.rds")
  if (SKIP_IF_EXISTS && file.exists(gc_file)) {
    cat("  loading cached GC-biased SE: ", gc_file, "\n", sep = "")
    cached <- tryCatch(readRDS(gc_file), error = function(e) NULL)
    if (!is.null(cached)) return(cached)
    cat("  cached GC-biased SE could not be read; recomputing.\n")
  }

  cat("  adding GC bias with ", GENOME_PACKAGE, "::", GENOME_OBJECT, " ...\n", sep = "")
  genome <- load_genome_object()
  SE <- chromVAR::addGCBias(object = SE, genome = genome)
  saveRDS(SE, gc_file)
  SE
}

load_meta <- function(SE) {
  if (nzchar(META_FILE) && file.exists(META_FILE)) {
    meta <- readRDS(META_FILE)
  } else {
    meta <- as.data.frame(SummarizedExperiment::colData(SE))
  }
  if (inherits(meta, "data.table")) meta <- as.data.frame(meta)

  if ((is.null(rownames(meta)) ||
       identical(rownames(meta), as.character(seq_len(nrow(meta))))) &&
      META_CELL_COL %in% colnames(meta)) {
    rownames(meta) <- as.character(meta[[META_CELL_COL]])
  }
  meta
}

infer_chr_prefix <- function(SE) {
  rr_chr <- as.character(GenomicRanges::seqnames(SummarizedExperiment::rowRanges(SE)))
  if (any(grepl("^chr", rr_chr))) "chr" else ""
}

normalize_chr <- function(chr, chr_prefix) {
  chr <- as.character(chr)
  chr <- sub("^CHR", "chr", chr, ignore.case = TRUE)
  chr <- sub("^chrM$", "chrMT", chr, ignore.case = TRUE)
  chr <- sub("^M$", "MT", chr, ignore.case = TRUE)
  out <- if (identical(chr_prefix, "chr")) {
    ifelse(grepl("^chr", chr), chr, paste0("chr", chr))
  } else {
    sub("^chr", "", chr)
  }
  out[is.na(chr)] <- NA_character_
  out
}

parse_variant_to_bed <- function(variant, chr_prefix) {
  variant <- as.character(variant)
  m <- regexec("^([^:]+):([0-9]+)_([^_]+)_([^_]+)$", variant)
  parts <- regmatches(variant, m)
  ok <- lengths(parts) == 5L

  chr <- rep(NA_character_, length(variant))
  pos <- rep(NA_integer_, length(variant))
  if (any(ok)) {
    chr[ok] <- vapply(parts[ok], `[`, character(1), 2)
    pos[ok] <- as.integer(vapply(parts[ok], `[`, character(1), 3))
  }

  data.table(
    ok = ok & !is.na(pos),
    chr = normalize_chr(chr, chr_prefix),
    start = pmax(pos - 1L, 0L),
    end = pos
  )
}

gwas_bed_for_trait <- function(gwas, trait, chr_prefix) {
  trait_vec <- as.character(gwas[[GWAS_TRAIT_COL]])
  dt <- gwas[trait_vec == trait]
  if (nrow(dt) == 0L) stop("No GWAS rows for trait: ", trait, call. = FALSE)

  if (all(c(GWAS_CHR_COL, GWAS_START_COL, GWAS_END_COL) %in% names(dt))) {
    bed_dt <- data.table(
      ok = TRUE,
      chr = normalize_chr(dt[[GWAS_CHR_COL]], chr_prefix),
      start = as.integer(dt[[GWAS_START_COL]]),
      end = as.integer(dt[[GWAS_END_COL]])
    )
  } else {
    variant_col <- GWAS_VARIANT_COL
    if (!variant_col %in% names(dt) && identical(variant_col, "var") &&
        "variant" %in% names(dt)) {
      variant_col <- "variant"
    }
    if (!variant_col %in% names(dt)) {
      stop("GWAS variant column not found: ", GWAS_VARIANT_COL, call. = FALSE)
    }
    bed_dt <- parse_variant_to_bed(dt[[variant_col]], chr_prefix = chr_prefix)
  }

  dt <- dt[bed_dt$ok]
  bed_dt <- bed_dt[bed_dt$ok]
  if (nrow(dt) == 0L) {
    stop("No parseable GWAS rows for trait: ", trait, call. = FALSE)
  }
  if (!GWAS_SCORE_COL %in% names(dt)) {
    stop("GWAS score column not found: ", GWAS_SCORE_COL, call. = FALSE)
  }

  bed <- data.table(
    chr = bed_dt$chr,
    start = as.integer(bed_dt$start),
    end = as.integer(bed_dt$end),
    name = trait,
    score = as.numeric(dt[[GWAS_SCORE_COL]])
  )
  bed[end <= start, end := start + 1L]
  bed
}

make_trait_imports <- function(SE, gwas) {
  chr_prefix <- infer_chr_prefix(SE)
  imports <- list()
  rows <- list()

  for (trait in TRAITS) {
    bed <- gwas_bed_for_trait(gwas, trait, chr_prefix)
    bed_file <- file.path(BED_DIR, paste0("trait_", make.names(trait), ".bed"))
    data.table::fwrite(bed, bed_file, sep = "\t", col.names = FALSE)

    imports[[trait]] <- gchromVAR::importBedScore(
      SummarizedExperiment::rowRanges(SE),
      bed_file,
      colidx = 5
    )
    rows[[trait]] <- data.table(trait = trait, n_gwas_rows = nrow(bed),
                                bed_file = bed_file)
    cat("  trait import ", trait, ": ", nrow(bed), " rows\n", sep = "")
  }
  data.table::fwrite(rbindlist(rows), file.path(OUTPUT_DIR, "trait_import_manifest.csv"))
  imports
}

sparsify_binomial <- function(counts, p_keep, seed) {
  if (p_keep >= 1.0) return(counts)
  set.seed(seed)
  out <- counts
  out@x <- as.numeric(stats::rbinom(
    n = length(counts@x),
    size = as.integer(counts@x),
    prob = p_keep
  ))
  Matrix::drop0(out)
}

get_background_peaks <- function(SE_sp, key, use_clean_bg) {
  bg_clean_file <- file.path(CACHE_DIR, "bg_peaks_clean.rds")
  if (use_clean_bg && file.exists(bg_clean_file)) {
    return(readRDS(bg_clean_file))
  }

  bg_file <- if (use_clean_bg) {
    bg_clean_file
  } else {
    file.path(CACHE_DIR, paste0("bg_peaks_", key, ".rds"))
  }

  if (SKIP_IF_EXISTS && file.exists(bg_file)) {
    return(readRDS(bg_file))
  }

  cat("  computing background peaks: ", basename(bg_file), "\n", sep = "")
  bg <- tryCatch(
    chromVAR::getBackgroundPeaks(object = SE_sp, niterations = BG_NITER),
    error = function(e) {
      rd <- SummarizedExperiment::rowData(SE_sp)
      bias <- NULL
      if ("GC" %in% names(rd)) bias <- rd$GC
      if (is.null(bias) && "bias" %in% names(rd)) bias <- rd$bias
      if (is.null(bias)) stop(e)
      chromVAR::getBackgroundPeaks(
        object = SE_sp,
        bias = bias,
        niterations = BG_NITER
      )
    }
  )
  saveRDS(bg, bg_file)
  bg
}

prepare_degraded_inputs <- function(SE_gc, counts_full, p_drop, seed, key) {
  p_keep <- 1 - p_drop
  counts_sp <- sparsify_binomial(counts_full, p_keep = p_keep, seed = seed)

  if (DROP_ALL_ZERO_PEAKS) {
    keep_peaks <- Matrix::rowSums(counts_sp) > 0
  } else {
    keep_peaks <- rep(TRUE, nrow(counts_sp))
  }
  if (!any(keep_peaks)) {
    stop("All peaks are zero after drop-out for ", key, call. = FALSE)
  }

  n_drop_peaks <- sum(!keep_peaks)
  if (n_drop_peaks > 0L) {
    cat("  dropping all-zero peaks: ", n_drop_peaks, " / ", length(keep_peaks), "\n", sep = "")
  }

  counts_sp <- counts_sp[keep_peaks, , drop = FALSE]
  SE_sp <- SE_gc[keep_peaks, ]
  SE_sp <- replace_counts_assay(SE_sp, counts_sp)

  bg_peaks <- get_background_peaks(
    SE_sp = SE_sp,
    key = key,
    use_clean_bg = n_drop_peaks == 0L
  )

  list(
    SE_sp = SE_sp,
    counts_sp = counts_sp,
    keep_peaks = keep_peaks,
    bg_peaks = bg_peaks,
    n_drop_peaks = n_drop_peaks
  )
}

orient_lsi <- function(lsi_mat, cell_ids) {
  lsi_mat <- as.matrix(lsi_mat)
  if (nrow(lsi_mat) == length(cell_ids)) {
    if (is.null(rownames(lsi_mat))) rownames(lsi_mat) <- cell_ids
    return(lsi_mat[cell_ids, , drop = FALSE])
  }
  if (ncol(lsi_mat) == length(cell_ids)) {
    lsi_mat <- t(lsi_mat)
    rownames(lsi_mat) <- cell_ids
    return(lsi_mat)
  }
  stop("LSI dimensions do not match cell count: ",
       paste(dim(lsi_mat), collapse = " x "), " vs ", length(cell_ids),
       call. = FALSE)
}

build_or_load_graph <- function(counts_sp, cell_ids, key) {
  lsi_file <- file.path(CACHE_DIR, paste0("atac_lsi_", key, ".rds"))
  knn_file <- file.path(CACHE_DIR, paste0("mutualknn_", key, ".rds"))

  if (SKIP_IF_EXISTS && file.exists(knn_file)) {
    mknn <- readRDS(knn_file)
    return(mknn)
  }

  if (SKIP_IF_EXISTS && file.exists(lsi_file)) {
    lsi_mat <- readRDS(lsi_file)
  } else {
    cat("  TF-IDF + LSI dims=", LSI_DIMS, " ...\n", sep = "")
    tfidf_mat <- SCAVENGE::tfidf(
      bmat = counts_sp,
      mat_binary = TRUE,
      TF = TRUE,
      log_TF = TRUE
    )
    lsi_mat <- SCAVENGE::do_lsi(mat = tfidf_mat, dims = LSI_DIMS)
    lsi_mat <- orient_lsi(lsi_mat, cell_ids)
    saveRDS(lsi_mat, lsi_file)
  }

  lsi_mat <- orient_lsi(lsi_mat, cell_ids)
  cat("  SCAVENGE mutual KNN k=", MUTUAL_K, " ...\n", sep = "")
  mknn <- SCAVENGE::getmutualknn(lsimat = lsi_mat, num_k = MUTUAL_K)
  if (is.null(rownames(mknn))) rownames(mknn) <- rownames(lsi_mat)
  if (is.null(colnames(mknn))) colnames(mknn) <- rownames(lsi_mat)
  saveRDS(mknn, knn_file)
  mknn
}

subset_trait_import <- function(trait_import, keep_peaks) {
  if (is.null(dim(trait_import))) {
    return(trait_import[keep_peaks])
  }
  if (nrow(trait_import) == length(keep_peaks)) {
    return(trait_import[keep_peaks, , drop = FALSE])
  }
  if (ncol(trait_import) == length(keep_peaks)) {
    return(trait_import[, keep_peaks, drop = FALSE])
  }
  stop("Cannot subset trait_import; dimensions are ",
       paste(dim(trait_import), collapse = " x "), call. = FALSE)
}

compute_or_load_z <- function(SE_sp, trait_import, bg_peaks, trait, key) {
  z_file <- file.path(CACHE_DIR, paste0("gchromvar_z_", key, "_", make.names(trait), ".rds"))
  if (SKIP_IF_EXISTS && file.exists(z_file)) return(readRDS(z_file))

  cat("  gchromVAR z-score: ", trait, "\n", sep = "")
  SE_DEV <- gchromVAR::computeWeightedDeviations(
    object = SE_sp,
    weights = trait_import,
    background_peaks = bg_peaks
  )
  z_score <- as.numeric(t(SummarizedExperiment::assays(SE_DEV)[["z"]]))
  if (length(z_score) != ncol(SE_sp)) {
    stop("Unexpected z-score length for ", trait, ": ", length(z_score),
         " vs cells ", ncol(SE_sp), call. = FALSE)
  }
  names(z_score) <- colnames(SE_sp)
  saveRDS(z_score, z_file)
  z_score
}

seed_cells_from_z <- function(z_graph) {
  idx <- SCAVENGE::seedindex(z_score = z_graph, percent_cut = SEED_PERCENT_CUT)
  if (is.logical(idx)) {
    names(z_graph)[idx]
  } else {
    names(z_graph)[as.integer(idx)]
  }
}

compute_auroc <- function(score, label) {
  keep <- is.finite(score) & !is.na(label)
  score <- score[keep]
  label <- as.integer(label[keep])
  if (length(score) < 20L) return(NA_real_)
  if (sum(label == 1L) < 10L || sum(label == 0L) < 10L) return(NA_real_)
  v <- stats::var(score)
  if (!is.finite(v) || v < 1e-12) return(NA_real_)
  roc_obj <- suppressMessages(pROC::roc(
    label,
    score,
    quiet = TRUE,
    direction = "<",
    levels = c(0L, 1L)
  ))
  as.numeric(pROC::auc(roc_obj))
}

run_trait_scavenge <- function(mknn, z_score, trait, display_key, meta) {
  graph_cells <- rownames(mknn)
  z_graph <- z_score[graph_cells]
  seed_cells <- seed_cells_from_z(z_graph)
  if (length(seed_cells) == 0L) {
    stop("No seed cells for ", trait, " / ", display_key, call. = FALSE)
  }

  scale_factor <- SCAVENGE::cal_scalefactor(
    z_score = z_score,
    percent_cut = SCALE_PERCENT_CUT
  )
  np_score <- SCAVENGE::randomWalk_sparse(
    intM = mknn,
    queryCells = seed_cells,
    gamma = GAMMA
  )
  if (is.null(names(np_score))) names(np_score) <- rownames(mknn)

  omit_idx <- np_score == 0
  np_use <- np_score[!omit_idx]
  keep_cells <- names(np_use)

  trs <- SCAVENGE::capOutlierQuantile(x = np_use, q_ceiling = CAP_Q)
  trs <- SCAVENGE::max_min_scale(trs)
  trs <- trs * scale_factor

  out_df <- data.frame(
    cell = keep_cells,
    z_score = as.numeric(z_score[keep_cells]),
    seed = keep_cells %in% seed_cells,
    np_score = as.numeric(np_use),
    TRS = as.numeric(trs),
    stringsAsFactors = FALSE
  )

  if (!is.null(meta) && nrow(meta) > 0L) {
    meta_out <- meta
    meta_out$cell <- rownames(meta_out)
    out_df <- dplyr::left_join(out_df, meta_out, by = "cell")
  }

  trs_file <- file.path(
    TRS_DIR,
    paste0("SCAVENGE_ATAC_", display_key, "_", make.names(trait), ".tsv.gz")
  )
  data.table::fwrite(out_df, trs_file, sep = "\t", compress = "gzip")

  auc <- NA_real_
  n_positive <- NA_integer_
  n_eval <- NA_integer_
  if (!is.null(meta) && META_CELLTYPE_COL %in% colnames(meta) && trait %in% names(TARGETS)) {
    common_eval <- intersect(out_df$cell, rownames(meta))
    label <- as.integer(as.character(meta[common_eval, META_CELLTYPE_COL]) %in% TARGETS[[trait]])
    score <- out_df$TRS[match(common_eval, out_df$cell)]
    auc <- compute_auroc(score, label)
    n_positive <- sum(label)
    n_eval <- length(label)
  }

  list(
    trs_file = trs_file,
    auroc = auc,
    n_positive = n_positive,
    n_eval = n_eval,
    n_seed = length(seed_cells),
    n_graph_cells = length(graph_cells),
    n_omit_zero = sum(omit_idx),
    scale_factor = scale_factor
  )
}

cat(">>> Pure-ATAC original SCAVENGE drop-out robustness\n\n")
cat("Run tag     : ", RUN_TAG, "\n", sep = "")
cat("Output dir  : ", OUTPUT_DIR, "\n", sep = "")
cat("Traits      : ", paste(TRAITS, collapse = ", "), "\n", sep = "")
cat("P_DROP_LIST : ", paste(P_DROP_LIST, collapse = ", "), "\n", sep = "")
cat("SEED_LIST   : ", paste(SEED_LIST, collapse = ", "), "\n\n", sep = "")

stop_missing(c(RAW_ATAC_SE_FILE, GWAS_FILE, META_FILE))

cat("[1/4] Loading inputs ...\n")
cat("  [1.1] read raw ATAC SE ...\n")
SE_raw <- readRDS(RAW_ATAC_SE_FILE)
cat("  [1.2] check GC bias ...\n")
SE_gc <- ensure_gc_se(SE_raw)
cat("  [1.3] extract counts assay ...\n")
counts_full <- counts_assay(SE_gc)
if (is.null(colnames(counts_full))) {
  colnames(counts_full) <- colnames(SE_gc)
}
if (is.null(colnames(SE_gc)) && !is.null(colnames(counts_full))) {
  colnames(SE_gc) <- colnames(counts_full)
}
cat("  [1.4] load metadata ...\n")
meta <- load_meta(SE_gc)

cat("  [1.5] load GWAS ...\n")
gwas <- data.table::fread(GWAS_FILE)
if (!GWAS_TRAIT_COL %in% names(gwas)) {
  stop("GWAS trait column not found: ", GWAS_TRAIT_COL, call. = FALSE)
}
gwas <- gwas[as.character(gwas[[GWAS_TRAIT_COL]]) %in% TRAITS]
if (nrow(gwas) == 0L) stop("No GWAS rows left after TRAITS filter.", call. = FALSE)

cat("  SE peak x cell: ", nrow(SE_gc), " x ", ncol(SE_gc), "\n", sep = "")
cat("  counts nnz    : ", length(counts_full@x), "\n", sep = "")
cat("  GWAS rows     : ", nrow(gwas), "\n", sep = "")
cat("  metadata rows : ", if (!is.null(meta)) nrow(meta) else 0L, "\n\n", sep = "")

cat("[2/4] Precomputing trait peak weights ...\n")
trait_imports <- make_trait_imports(SE_gc, gwas)

eval_rows <- list()
manifest_rows <- list()
row_i <- 0L
manifest_i <- 0L

cat("\n[3/4] Running drop-out grid ...\n")
for (p_drop in P_DROP_LIST) {
  p_drop <- round(as.numeric(p_drop), 2)
  for (seed in SEED_LIST) {
    seed <- as.integer(seed)
    cache_seed <- cache_seed_for(p_drop, seed)
    cache_key <- run_key(p_drop, cache_seed)
    display_key <- run_key(p_drop, seed)

    cat("\n", strrep("=", 72), "\n", sep = "")
    cat("Run ", display_key, " (cache ", cache_key, ")\n", sep = "")
    cat(strrep("=", 72), "\n", sep = "")

    degraded <- prepare_degraded_inputs(
      SE_gc = SE_gc,
      counts_full = counts_full,
      p_drop = p_drop,
      seed = cache_seed,
      key = cache_key
    )
    mknn <- build_or_load_graph(
      counts_sp = degraded$counts_sp,
      cell_ids = colnames(degraded$SE_sp),
      key = cache_key
    )

    manifest_i <- manifest_i + 1L
    manifest_rows[[manifest_i]] <- data.table(
      p_drop = p_drop,
      seed = seed,
      cache_seed = cache_seed,
      n_peaks = nrow(degraded$SE_sp),
      n_cells = ncol(degraded$SE_sp),
      n_nonzero_counts = length(degraded$counts_sp@x),
      n_dropped_all_zero_peaks = degraded$n_drop_peaks,
      knn_nnz = length(mknn@x)
    )

    for (trait in TRAITS) {
      trs_file <- file.path(
        TRS_DIR,
        paste0("SCAVENGE_ATAC_", display_key, "_", make.names(trait), ".tsv.gz")
      )

      if (SKIP_IF_EXISTS && file.exists(trs_file)) {
        cat("  skip existing TRS: ", basename(trs_file), "\n", sep = "")
        trs_dt <- data.table::fread(trs_file)
        auc <- NA_real_
        n_positive <- NA_integer_
        n_eval <- NA_integer_
        if (META_CELLTYPE_COL %in% colnames(meta) && trait %in% names(TARGETS)) {
          common_eval <- intersect(trs_dt$cell, rownames(meta))
          label <- as.integer(as.character(meta[common_eval, META_CELLTYPE_COL]) %in% TARGETS[[trait]])
          score <- trs_dt$TRS[match(common_eval, trs_dt$cell)]
          auc <- compute_auroc(score, label)
          n_positive <- sum(label)
          n_eval <- length(label)
        }
        res <- list(
          trs_file = trs_file,
          auroc = auc,
          n_positive = n_positive,
          n_eval = n_eval,
          n_seed = sum(trs_dt$seed),
          n_graph_cells = nrow(trs_dt),
          n_omit_zero = NA_integer_,
          scale_factor = NA_real_
        )
      } else {
        trait_import_current <- subset_trait_import(
          trait_imports[[trait]],
          degraded$keep_peaks
        )
        z_score <- compute_or_load_z(
          SE_sp = degraded$SE_sp,
          trait_import = trait_import_current,
          bg_peaks = degraded$bg_peaks,
          trait = trait,
          key = cache_key
        )
        res <- run_trait_scavenge(
          mknn = mknn,
          z_score = z_score,
          trait = trait,
          display_key = display_key,
          meta = meta
        )
      }

      row_i <- row_i + 1L
      eval_rows[[row_i]] <- data.table(
        method = "original_scavenge_atac",
        p_drop = p_drop,
        seed = seed,
        trait = trait,
        target = if (trait %in% names(TARGETS)) trait else NA_character_,
        auroc = res$auroc,
        n_positive = res$n_positive,
        n_eval = res$n_eval,
        n_seed = res$n_seed,
        n_graph_cells = res$n_graph_cells,
        n_omit_zero = res$n_omit_zero,
        scale_factor = res$scale_factor,
        trs_file = res$trs_file
      )

      cat(sprintf(
        "  %-10s AUROC=%s n_seed=%d\n",
        trait,
        ifelse(is.na(res$auroc), "NA", sprintf("%.4f", res$auroc)),
        as.integer(res$n_seed)
      ))
    }
  }
}

cat("\n[4/4] Writing evaluation tables ...\n")
eval_dt <- rbindlist(eval_rows, use.names = TRUE, fill = TRUE)
if (APPEND_EXISTING_EVAL && file.exists(EVAL_FILE)) {
  old <- data.table::fread(EVAL_FILE)
  eval_dt <- rbindlist(list(old, eval_dt), use.names = TRUE, fill = TRUE)
}
eval_dt <- unique(eval_dt, by = c("method", "p_drop", "seed", "trait", "target"))
setorder(eval_dt, trait, p_drop, seed)
data.table::fwrite(eval_dt, EVAL_FILE)

manifest_dt <- rbindlist(manifest_rows, use.names = TRUE, fill = TRUE)
if (file.exists(RUN_MANIFEST_FILE)) {
  old_manifest <- data.table::fread(RUN_MANIFEST_FILE)
  manifest_dt <- rbindlist(list(old_manifest, manifest_dt), use.names = TRUE, fill = TRUE)
}
manifest_dt <- unique(manifest_dt, by = c("p_drop", "seed"))
setorder(manifest_dt, p_drop, seed)
data.table::fwrite(manifest_dt, RUN_MANIFEST_FILE)

summary_dt <- eval_dt[
  method == "original_scavenge_atac" & trait %in% TRAITS & is.finite(auroc),
  .(
    mean_auroc = mean(auroc, na.rm = TRUE),
    sd_auroc = stats::sd(auroc, na.rm = TRUE),
    n_runs = .N,
    n_positive = n_positive[1],
    n_eval = n_eval[1]
  ),
  by = .(trait, target, p_drop)
]
summary_dt[, se_auroc := sd_auroc / sqrt(n_runs)]
summary_dt[, lo95 := mean_auroc - 1.96 * se_auroc]
summary_dt[, hi95 := mean_auroc + 1.96 * se_auroc]
setorder(summary_dt, trait, p_drop)
data.table::fwrite(summary_dt, SUMMARY_FILE)

cat("Wrote:\n")
cat("  ", EVAL_FILE, "\n", sep = "")
cat("  ", SUMMARY_FILE, "\n", sep = "")
cat("  ", RUN_MANIFEST_FILE, "\n", sep = "")
cat("\nDone. Run 02_plot_robustness.R to draw AUROC curves.\n")
