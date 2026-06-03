# module4 robustness: v3 sparse MONO full raw ATAC drop-out pipeline
# -----------------------------------------------------------------------------
# Purpose:
#   Run a paired-module3-style robustness test for the current module4 v3 route:
#
#     raw ATAC counts drop-out
#       -> recompute ATAC LSI z block
#       -> recompute ATAC Gene Activity
#       -> rerun sparse ATAC->RNA mapping
#       -> rerun top5-significant RNA completion
#       -> rerun RNA neighborhood-FDR evidence from the new mapping
#       -> rerun module3 KNN/RWR/AUROC
#
# Boundary:
#   This is intentionally heavier than 09_run_v3_current_code_atac_dropout.R.
#   That older script only perturbs the final ATAC z block. This script starts
#   from raw ATAC counts and reruns the current v3 pipeline for each
#   p_drop x seed condition.
#
# Run in RStudio:
#   setwd("D:/RStudio/work/module4_traitaware_pseudopair")
#   source("05_robustness/10_run_v3_sparse_mono_raw_atac_full_pipeline_dropout.R")
#
# Useful quick sanity:
#   ROBUST_P_DROP_LIST <- c(0)
#   ROBUST_SEED_LIST <- c(1)
#   source("05_robustness/10_run_v3_sparse_mono_raw_atac_full_pipeline_dropout.R")
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Matrix)
  library(SummarizedExperiment)
  library(SCAVENGE)
  library(data.table)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(ensembldb)
  library(EnsDb.Hsapiens.v75)
})

ROOT <- "D:/RStudio/work/module4_traitaware_pseudopair"

RUN_ENV <- environment()
has_local <- function(name) exists(name, envir = RUN_ENV, inherits = FALSE)

if (!has_local("ROBUST_OUTPUT_TAG")) {
  ROBUST_OUTPUT_TAG <- "v3_sparse_mono_raw_atac_full_pipeline_dropout"
}
if (!has_local("ROBUST_P_DROP_LIST")) {
  ROBUST_P_DROP_LIST <- round(seq(0.0, 0.9, by = 0.1), 2)
}
if (!has_local("ROBUST_SEED_LIST")) ROBUST_SEED_LIST <- c(1L, 2L, 3L)
if (!has_local("ROBUST_ATAC_W_LIST")) ROBUST_ATAC_W_LIST <- c(0.7, 1.0)
if (!has_local("ROBUST_TRAITS_TO_RUN")) ROBUST_TRAITS_TO_RUN <- c("MONO")
if (!has_local("ROBUST_RUN_LABEL")) ROBUST_RUN_LABEL <- paste(tolower(ROBUST_TRAITS_TO_RUN), collapse = "_")
if (!has_local("ROBUST_NBD_MODE")) ROBUST_NBD_MODE <- "recompute"
if (!has_local("ROBUST_SEED_MODE")) ROBUST_SEED_MODE <- "neighborhood_fdr"
if (!has_local("ROBUST_P0_USE_CLEAN_RAW_PRODUCTS")) ROBUST_P0_USE_CLEAN_RAW_PRODUCTS <- TRUE
if (!has_local("ROBUST_CACHE_DEGRADED_COUNTS")) ROBUST_CACHE_DEGRADED_COUNTS <- FALSE
if (!has_local("ROBUST_SPARSE_MAP_BATCH_SIZE")) ROBUST_SPARSE_MAP_BATCH_SIZE <- 128L
if (!has_local("ROBUST_SKIP_COMPLETED_CONDITIONS")) ROBUST_SKIP_COMPLETED_CONDITIONS <- TRUE
if (!has_local("ROBUST_AUROC_RUNNER")) {
  ROBUST_AUROC_RUNNER <- "04_evaluation/run_neighborhood_fdr_seed_auroc.R"
}
if (!has_local("ROBUST_AUROC_DIR_LABEL")) {
  ROBUST_AUROC_DIR_LABEL <- paste0(ROBUST_RUN_LABEL, "_neighborhood_fdr_ge0p80")
}
if (!has_local("ROBUST_METHOD_LABEL")) ROBUST_METHOD_LABEL <- "v3_raw_full_pipeline"

VALID_STEPS <- c(
  "raw_products",
  "sparse_mapping",
  "top5sig_transfer",
  "neighborhood_fdr",
  "module3_blocks",
  "knn_rwr",
  "neighborhood_seed",
  "neighborhood_rwr",
  "auroc",
  "summarize"
)
if (has_local("ROBUST_RESUME_STEPS")) {
  ROBUST_RUN_STEPS <- ROBUST_RESUME_STEPS
} else if (!has_local("ROBUST_RUN_STEPS")) {
  ROBUST_RUN_STEPS <- VALID_STEPS
}
bad_steps <- setdiff(ROBUST_RUN_STEPS, VALID_STEPS)
if (length(bad_steps) > 0L) {
  stop("Unknown ROBUST_RUN_STEPS: ", paste(bad_steps, collapse = ", "),
       "\nValid steps: ", paste(VALID_STEPS, collapse = ", "),
       call. = FALSE)
}
if (!ROBUST_NBD_MODE %in% c("recompute", "reuse_clean")) {
  stop("ROBUST_NBD_MODE must be 'recompute' or 'reuse_clean'.", call. = FALSE)
}
if (!ROBUST_SEED_MODE %in% c("neighborhood_fdr", "sparse_neighborhood_screen")) {
  stop("ROBUST_SEED_MODE must be 'neighborhood_fdr' or 'sparse_neighborhood_screen'.",
       call. = FALSE)
}

if (!has_local("RAW_ATAC_SE_FILE")) {
  RAW_ATAC_SE_FILE <- "D:/SCAVENGE-reproducibility-main/blood_data/scATAC-Healthy-Hematopoiesis-191120.rds"
}
if (!has_local("CLEAN_GA_PATH")) {
  CLEAN_GA_PATH <- "D:/SCAVENGE-reproducibility-main/blood_data/multiomics_output_v2/gene_activity_matrix.rds"
}
if (!has_local("CLEAN_COMPAT_DIR")) {
  CLEAN_COMPAT_DIR <- file.path(
    ROOT,
    "03_pseudopair_pipeline/output/v3_sparse_mapping_top5significant",
    "module3_mono_neighborhood_fdr_ge0p80"
  )
}
if (!has_local("CLEAN_ATAC_Z_PATH")) {
  CLEAN_ATAC_Z_PATH <- file.path(CLEAN_COMPAT_DIR, "atac_z_block_v2.rds")
}
if (!has_local("CLEAN_NBD_DIR")) {
  CLEAN_NBD_DIR <- file.path(ROOT, "04_evaluation/module4_v2_diagnostics/neighborhood_fdr")
}
if (!has_local("CLEAN_AUROC_CSV")) {
  CLEAN_AUROC_CSV <- file.path(
    ROOT,
    "04_evaluation/v3_sparse_mapping_top5significant/mono_neighborhood_fdr_ge0p80",
    "v3_sparse_mapping_top5significant_mono_neighborhood_fdr_ge0p80_auroc.csv"
  )
}

if (!has_local("BASE_CITE_EXTRACT_DIR")) {
  BASE_CITE_EXTRACT_DIR <- file.path(ROOT, "02_trait_mapping/output/cite_extract")
}
if (!has_local("ROBUST_GS_FILE")) {
  ROBUST_GS_FILE <- file.path(ROOT, "01_data_prep/03_rna_scdrs_python/granja_3traits.gs")
}
if (!has_local("ROBUST_CITE_H5AD")) ROBUST_CITE_H5AD <- NULL
if (!has_local("CITE_SCORE_RDS")) {
  CITE_SCORE_RDS <- file.path(ROOT, "01_data_prep_external_rna/output/cite_rna_scdrs_scores.rds")
}
if (!has_local("ROBUST_TRAITS_ALL")) {
  ROBUST_TRAITS_ALL <- c(
    MONO = "UKB_460K.blood_MONOCYTE_COUNT",
    PLT  = "UKB_460K.blood_PLATELET_COUNT",
    LYM  = "UKB_460K.blood_LYMPHOCYTE_COUNT",
    MCH  = "UKB_460K.blood_MEAN_CORPUSCULAR_HEMOGLOBIN",
    RDW  = "UKB_460K.blood_RBC_DISTRIB_WIDTH",
    RBC  = "UKB_460K.blood_RED_COUNT",
    EOS  = "UKB_460K.blood_EOSINOPHIL_COUNT",
    WBC  = "UKB_460K.blood_WHITE_COUNT"
  )
}
if (!has_local("CITE_FULL_SCORE_DIR")) {
  CITE_FULL_SCORE_DIR <- file.path(ROOT, "01_data_prep_external_rna/03_cite_scdrs_python/output")
}
if (!has_local("ATAC_META_FILE")) {
  ATAC_META_FILE <- file.path(ROOT, "04_evaluation/inputs/granja_cell_meta.rds")
}
if (!has_local("ROBUST_RNA_GROUP_FILE")) {
  ROBUST_RNA_GROUP_FILE <- file.path(
    ROOT,
    "04_evaluation/module4_v2_diagnostics/neighborhood_fdr/rna_unsupervised_group_membership.csv.gz"
  )
}
if (!has_local("ROBUST_RNA_GROUP_TEST_FILES")) {
  ROBUST_RNA_GROUP_TEST_FILES <- c(
    file.path(ROOT, "04_evaluation/module4_v2_diagnostics/neighborhood_fdr/neighborhood_fdr_group_tests.csv"),
    file.path(ROOT, "01_data_prep_external_rna/output/trait_screening/cite_trait_screening_neighborhood_fdr_group_tests.csv")
  )
}
if (!has_local("ROBUST_GROUP_TEST_TRAIT_MAP")) {
  ROBUST_GROUP_TEST_TRAIT_MAP <- c(EOS = "EOSINOPHIL", WBC = "WHITE")
}
if (!has_local("ROBUST_GROUP_FDR_THRESHOLD")) ROBUST_GROUP_FDR_THRESHOLD <- 0.10

ROBUST_DIR <- file.path(ROOT, "05_robustness/output", ROBUST_OUTPUT_TAG)
RAW_DIR <- file.path(ROBUST_DIR, "raw_products")
RUNS_DIR <- file.path(ROBUST_DIR, "runs")
dir.create(RAW_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RUNS_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_ALL <- file.path(ROBUST_DIR, paste0(ROBUST_OUTPUT_TAG, "_auroc_all.csv"))
OUT_PLOT <- file.path(ROBUST_DIR, paste0(ROBUST_OUTPUT_TAG, "_plot_data.csv"))
OUT_MANIFEST <- file.path(ROBUST_DIR, paste0(ROBUST_OUTPUT_TAG, "_run_manifest.csv"))
OUT_P0_COMPARE <- file.path(ROBUST_DIR, paste0(ROBUST_OUTPUT_TAG, "_p0_clean_compare.csv"))
OUT_MISSING_GRID <- file.path(ROBUST_DIR, paste0(ROBUST_OUTPUT_TAG, "_missing_grid.csv"))

stop_missing <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0L) {
    stop("Missing required file(s):\n  ", paste(missing, collapse = "\n  "),
         call. = FALSE)
  }
}

step_enabled <- function(x) x %in% ROBUST_RUN_STEPS

condition_auroc_csv <- function(p_drop, seed) {
  run_id <- sprintf("p%.2f_seed%d", p_drop, seed)
  file.path(
    RUNS_DIR,
    run_id,
    "04_evaluation",
    ROBUST_AUROC_DIR_LABEL,
    paste0(run_id, "_v3_sparse_", ROBUST_RUN_LABEL, "_full_pipeline_auroc.csv")
  )
}

condition_is_complete <- function(p_drop, seed, traits = ROBUST_TRAITS_TO_RUN) {
  csv <- condition_auroc_csv(p_drop, seed)
  if (!file.exists(csv)) return(FALSE)
  dt <- tryCatch(fread(csv), error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0L) return(FALSE)
  trait_col <- intersect(c("trait", "target"), names(dt))
  if (length(trait_col) == 0L) return(TRUE)
  seen <- unique(as.character(dt[[trait_col[1L]]]))
  all(traits %in% seen)
}

merge_existing_auroc <- function(path, old_dt) {
  if (is.null(old_dt) || nrow(old_dt) == 0L || !file.exists(path)) {
    return(invisible(FALSE))
  }
  new_dt <- tryCatch(fread(path), error = function(e) NULL)
  if (is.null(new_dt) || nrow(new_dt) == 0L) return(invisible(FALSE))
  merged <- rbindlist(list(old_dt, new_dt), use.names = TRUE, fill = TRUE)
  key_cols <- intersect(c("method", "target", "atac_w"), names(merged))
  if (length(key_cols) > 0L) {
    merged <- unique(merged, by = key_cols, fromLast = TRUE)
  }
  fwrite(merged, path)
  invisible(TRUE)
}

run_script <- function(script, vars = list()) {
  e <- new.env(parent = globalenv())
  for (nm in names(vars)) assign(nm, vars[[nm]], envir = e)
  cat("\n>>> source(", script, ")\n", sep = "")
  source(file.path(ROOT, script), local = e)
  invisible(e)
}

zscore_columns <- function(x) {
  x <- as.matrix(x)
  mu <- colMeans(x)
  sd <- apply(x, 2L, stats::sd)
  sd[!is.finite(sd) | sd == 0] <- 1
  x <- sweep(x, 2L, mu, "-")
  sweep(x, 2L, sd, "/")
}

sparsify_binomial <- function(counts, p_keep, seed) {
  if (p_keep >= 1.0) return(counts)
  set.seed(seed)
  out <- counts
  out@x <- as.numeric(stats::rbinom(
    length(counts@x),
    size = as.integer(counts@x),
    prob = p_keep
  ))
  Matrix::drop0(out)
}

get_counts <- function(se) {
  assay_names <- SummarizedExperiment::assayNames(se)
  if ("counts" %in% assay_names) {
    return(SummarizedExperiment::assay(se, "counts"))
  }
  SummarizedExperiment::assay(se, 1L)
}

build_gene_regions_hg19 <- function(peak_gr) {
  peak_has_chr <- any(grepl("^chr", as.character(GenomeInfoDb::seqlevels(peak_gr))))
  style <- if (peak_has_chr) "UCSC" else "NCBI"

  gene_coords <- ensembldb::genes(
    EnsDb.Hsapiens.v75,
    filter = ~ gene_biotype == "protein_coding",
    columns = c("gene_name", "gene_biotype")
  )
  GenomeInfoDb::seqlevelsStyle(gene_coords) <- style
  gene_coords <- GenomeInfoDb::keepStandardChromosomes(
    gene_coords,
    pruning.mode = "coarse"
  )
  gene_coords <- gene_coords[!duplicated(gene_coords$gene_name)]
  names(gene_coords) <- gene_coords$gene_name

  gene_regions <- gene_coords
  sp <- as.character(GenomicRanges::strand(gene_regions))
  strand_pos <- sp == "+"
  strand_neg <- sp == "-"
  strand_star <- sp == "*"
  GenomicRanges::start(gene_regions)[strand_pos | strand_star] <-
    pmax(GenomicRanges::start(gene_regions)[strand_pos | strand_star] - 2000L, 1L)
  GenomicRanges::end(gene_regions)[strand_neg] <-
    GenomicRanges::end(gene_regions)[strand_neg] + 2000L
  gene_regions
}

gene_peak_map_cache_path <- function() {
  file.path(RAW_DIR, "granja_hg19_promoter2kb_gene_peak_map.rds")
}

get_gene_peak_map <- function(peak_gr, n_peaks) {
  cache <- gene_peak_map_cache_path()
  if (file.exists(cache)) {
    obj <- readRDS(cache)
    if (as.integer(obj$n_peaks) == as.integer(n_peaks)) return(obj)
    cat("  cached gene-peak map has different peak count; rebuilding.\n")
  }

  cat("  building Granja hg19 promoter2kb gene-peak map ...\n")
  gene_regions <- build_gene_regions_hg19(peak_gr)
  ov <- GenomicRanges::findOverlaps(gene_regions, peak_gr)
  gene_peak_map <- Matrix::sparseMatrix(
    i = S4Vectors::queryHits(ov),
    j = S4Vectors::subjectHits(ov),
    x = 1,
    dims = c(length(gene_regions), n_peaks)
  )
  rownames(gene_peak_map) <- names(gene_regions)
  obj <- list(
    gene_peak_map = gene_peak_map,
    gene_names = names(gene_regions),
    n_peaks = n_peaks,
    created_at = as.character(Sys.time())
  )
  saveRDS(obj, cache)
  obj
}

build_gene_activity <- function(gene_peak_map, counts, out_file) {
  if (file.exists(out_file)) return(out_file)
  cat("  computing degraded Gene Activity ...\n")
  ga <- gene_peak_map %*% counts
  rownames(ga) <- rownames(gene_peak_map)
  colnames(ga) <- colnames(counts)
  keep <- Matrix::rowSums(ga) > 0
  ga <- ga[keep, , drop = FALSE]
  saveRDS(ga, out_file)
  out_file
}

build_atac_z <- function(counts, all_cells, out_file) {
  if (file.exists(out_file)) return(out_file)
  cat("  recomputing ATAC TF-IDF + LSI z block ...\n")
  tfidf_mat <- SCAVENGE::tfidf(
    bmat = counts,
    mat_binary = TRUE,
    TF = TRUE,
    log_TF = TRUE
  )
  lsi_mat <- SCAVENGE::do_lsi(mat = tfidf_mat, dims = 30)
  if (ncol(lsi_mat) != 30 && nrow(lsi_mat) == 30) lsi_mat <- t(lsi_mat)
  if (is.null(rownames(lsi_mat)) && nrow(lsi_mat) == length(all_cells)) {
    rownames(lsi_mat) <- all_cells
  }
  if (nrow(lsi_mat) != length(all_cells)) {
    stop("LSI cell dimension mismatch: ", nrow(lsi_mat),
         " vs ", length(all_cells), call. = FALSE)
  }
  rownames(lsi_mat) <- all_cells
  atac_z <- zscore_columns(lsi_mat)
  rownames(atac_z) <- all_cells
  colnames(atac_z) <- paste0("ATAC_LSI_", seq_len(ncol(atac_z)))
  saveRDS(atac_z, out_file)
  out_file
}

prepare_raw_products <- function(p_drop, seed, se, counts_full, gene_peak_obj) {
  prefix <- sprintf("p%.2f_seed%d", p_drop, seed)
  ga_file <- file.path(RAW_DIR, paste0("gene_activity_", prefix, ".rds"))
  z_file <- file.path(RAW_DIR, paste0("atac_z_block_", prefix, ".rds"))

  if (p_drop == 0 && isTRUE(ROBUST_P0_USE_CLEAN_RAW_PRODUCTS)) {
    return(list(
      ga_path = CLEAN_GA_PATH,
      atac_z_path = CLEAN_ATAC_Z_PATH,
      raw_product_mode = "clean_existing_products"
    ))
  }

  p_keep <- 1 - p_drop
  cat("  sparsifying raw ATAC counts: p_drop=", p_drop,
      " seed=", seed, "\n", sep = "")
  counts_sp_file <- file.path(RAW_DIR, paste0("counts_", prefix, ".rds"))
  if (isTRUE(ROBUST_CACHE_DEGRADED_COUNTS) && file.exists(counts_sp_file)) {
    counts_sp <- readRDS(counts_sp_file)
  } else {
    counts_sp <- sparsify_binomial(counts_full, p_keep = p_keep, seed = seed)
    if (isTRUE(ROBUST_CACHE_DEGRADED_COUNTS)) saveRDS(counts_sp, counts_sp_file)
  }

  build_gene_activity(gene_peak_obj$gene_peak_map, counts_sp, ga_file)
  build_atac_z(counts_sp, colnames(counts_full), z_file)

  list(
    ga_path = ga_file,
    atac_z_path = z_file,
    raw_product_mode = "recomputed_from_degraded_raw_counts"
  )
}

mapping_rng_seed <- function(p_drop, seed) {
  if (p_drop == 0) return(20260519L)
  20260519L + as.integer(seed) + as.integer(round(p_drop * 100))
}

completion_rng_seed <- function(p_drop, seed) {
  if (p_drop == 0) return(20260520L)
  20260520L + as.integer(seed) + as.integer(round(p_drop * 100))
}

run_one_condition <- function(p_drop, seed, products) {
  run_id <- sprintf("p%.2f_seed%d", p_drop, seed)
  run_root <- file.path(RUNS_DIR, run_id)
  sparse_map_dir <- file.path(run_root, "02_trait_mapping", "v3_sparse_mapping")
  sparse_map_diag <- file.path(sparse_map_dir, "diagnostics")
  transfer_dir <- file.path(run_root, "03_pseudopair_pipeline",
                            paste0("transfer_", ROBUST_RUN_LABEL))
  completion_diag <- file.path(run_root, "02_trait_mapping",
                               paste0(ROBUST_RUN_LABEL, "_completion_diagnostics"))
  seed_diag_dir <- file.path(run_root, "02_trait_mapping",
                             paste0(ROBUST_RUN_LABEL, "_neighborhood_seed_diagnostics"))
  nbd_dir <- if (identical(ROBUST_NBD_MODE, "recompute")) {
    file.path(run_root, "04_evaluation",
              paste0("neighborhood_fdr_", ROBUST_RUN_LABEL))
  } else {
    CLEAN_NBD_DIR
  }
  compat_dir <- file.path(run_root, "03_pseudopair_pipeline",
                          paste0("module3_", ROBUST_RUN_LABEL, "_top5significant"))
  nbd_compat_dir <- file.path(run_root, "03_pseudopair_pipeline",
                              paste0("module3_", ROBUST_RUN_LABEL, "_neighborhood_fdr_ge0p80"))
  auroc_dir <- file.path(run_root, "04_evaluation", ROBUST_AUROC_DIR_LABEL)
  auroc_csv <- paste0(run_id, "_v3_sparse_", ROBUST_RUN_LABEL,
                      "_full_pipeline_auroc.csv")

  if (step_enabled("sparse_mapping")) {
    run_script("02_trait_mapping/08_build_sparse_trait_mapping.R", list(
      SPARSE_MAP_OUT_DIR = sparse_map_dir,
      DIAG_DIR = sparse_map_diag,
      CITE_EXTRACT_DIR = BASE_CITE_EXTRACT_DIR,
      GA_PATH = products$ga_path,
      ATAC_META_PATH = ATAC_META_FILE,
      GS_FILE = ROBUST_GS_FILE,
      CITE_H5AD = ROBUST_CITE_H5AD,
      CITE_SCORE_RDS = CITE_SCORE_RDS,
      TRAITS_ALL = ROBUST_TRAITS_ALL,
      TRAITS_TO_RUN = ROBUST_TRAITS_TO_RUN,
      BATCH_SIZE = ROBUST_SPARSE_MAP_BATCH_SIZE,
      B_NULL_CELL = 1000L,
      RNA_EMPIRICAL_P_THRESHOLD = 0.05,
      ATAC_EMPIRICAL_P_THRESHOLD = 0.05,
      MIN_DELTA_OVER_NULL_MEDIAN = 0,
      RNG_SEED = mapping_rng_seed(p_drop, seed)
    ))
  }

  if (step_enabled("top5sig_transfer")) {
    run_script("02_trait_mapping/07_build_top5_significant_completion.R", list(
      MAP_DIR = sparse_map_dir,
      TRANSFER_OUT_DIR = transfer_dir,
      DIAG_DIR = completion_diag,
      CITE_MTX = file.path(BASE_CITE_EXTRACT_DIR, "cite_trait_counts.mtx"),
      CITE_GENES = file.path(BASE_CITE_EXTRACT_DIR, "cite_trait_counts_genes.txt"),
      CITE_CELLS = file.path(BASE_CITE_EXTRACT_DIR, "cite_trait_counts_cells.txt"),
      CITE_SCORE_RDS = CITE_SCORE_RDS,
      TRAITS_TO_RUN = ROBUST_TRAITS_TO_RUN,
      B_NULL_TOP5 = 200L,
      SIGNIFICANCE_METHOD = "empirical_p",
      EMPIRICAL_P_THRESHOLD = 0.05,
      FILTER_Q_THRESHOLD = 0.10,
      MIN_DELTA_OVER_NULL_MEDIAN = 0,
      MAX_RETAINED_GENES_PER_CELL = Inf,
      FILTER_BATCH_SIZE = 512L,
      RNG_SEED = completion_rng_seed(p_drop, seed),
      NULL_WEIGHT_MODE = "sample_observed"
    ))
  }

  if (step_enabled("neighborhood_fdr") && identical(ROBUST_NBD_MODE, "recompute")) {
    if (identical(ROBUST_SEED_MODE, "sparse_neighborhood_screen")) {
      seed_vars <- list(
        MAP_DIR = sparse_map_dir,
        OUT_DIR = seed_diag_dir,
        TRAITS_TO_SCREEN = ROBUST_TRAITS_TO_RUN,
        RNA_GROUP_FILE = ROBUST_RNA_GROUP_FILE,
        RNA_GROUP_TEST_FILES = ROBUST_RNA_GROUP_TEST_FILES,
        GROUP_TEST_TRAIT_MAP = ROBUST_GROUP_TEST_TRAIT_MAP,
        GROUP_FDR_THRESHOLD = ROBUST_GROUP_FDR_THRESHOLD
      )
      if (has_local("ROBUST_POSTHOC_TARGETS_LIST")) {
        seed_vars$POSTHOC_TARGETS_LIST <- ROBUST_POSTHOC_TARGETS_LIST
      }
      if (has_local("ROBUST_POSTHOC_TARGET_NAMES")) {
        seed_vars$POSTHOC_TARGET_NAMES <- ROBUST_POSTHOC_TARGET_NAMES
      }
      run_script("02_trait_mapping/09_screen_sparse_mapping_neighborhood_seeds.R", seed_vars)
    } else {
      run_script("04_evaluation/diagnose_scdrs_neighborhood_fdr.R", list(
        MAP_DIR = sparse_map_dir,
        CITE_MTX = file.path(BASE_CITE_EXTRACT_DIR, "cite_trait_counts.mtx"),
        CITE_GENES = file.path(BASE_CITE_EXTRACT_DIR, "cite_trait_counts_genes.txt"),
        CITE_CELLS = file.path(BASE_CITE_EXTRACT_DIR, "cite_trait_counts_cells.txt"),
        ATAC_META_FILE = ATAC_META_FILE,
        CITE_SCORE_RDS = CITE_SCORE_RDS,
        FULL_SCORE_DIR = CITE_FULL_SCORE_DIR,
        OUT_DIR = nbd_dir,
        TRAITS = ROBUST_TRAITS_TO_RUN,
        GROUP_FDR_THRESHOLDS = c(ROBUST_GROUP_FDR_THRESHOLD)
      ))
    }
  }

  if (step_enabled("module3_blocks")) {
    if (identical(ROBUST_SEED_MODE, "sparse_neighborhood_screen")) {
      run_script("03_pseudopair_pipeline/07_prepare_trait_screening_module3_inputs.R", list(
        MODULE3_TRANSFER_DIR = transfer_dir,
        MODULE3_COMPAT_OUT_DIR = nbd_compat_dir,
        TRAITS_TO_RUN = ROBUST_TRAITS_TO_RUN,
        ATAC_Z_BLOCK_OVERRIDE = products$atac_z_path,
        SEED_CELL_FILE = file.path(seed_diag_dir, "sparse_neighborhood_seed_cell_level.csv.gz"),
        SEED_SUMMARY_FILE = file.path(nbd_compat_dir, "seed_summary_neighborhood_fdr_ge0p80.csv")
      ))
    } else {
      run_script("03_pseudopair_pipeline/02_prepare_module3_v2_inputs.R", list(
        MODULE3_TRANSFER_DIR = transfer_dir,
        MODULE3_COMPAT_OUT_DIR = compat_dir,
        TRAITS_TO_RUN = ROBUST_TRAITS_TO_RUN,
        ATAC_Z_BLOCK_OVERRIDE = products$atac_z_path
      ))
    }
  }

  if (step_enabled("knn_rwr")) {
    run_script("03_pseudopair_pipeline/03_run_module3_v2_knn_rwr.R", list(
      MODULE3_COMPAT_OUT_DIR = if (identical(ROBUST_SEED_MODE, "sparse_neighborhood_screen")) {
        nbd_compat_dir
      } else {
        compat_dir
      },
      PATCHED_DIR_OVERRIDE = file.path(ROBUST_DIR, "_patched_module3_scripts", run_id),
      TRAITS_TO_RUN = ROBUST_TRAITS_TO_RUN,
      ATAC_W_LIST = ROBUST_ATAC_W_LIST
    ))
  }

  if (step_enabled("neighborhood_seed") && identical(ROBUST_SEED_MODE, "neighborhood_fdr")) {
    run_script("03_pseudopair_pipeline/04_prepare_neighborhood_fdr_seed_inputs.R", list(
      BASE_COMPAT_DIR = compat_dir,
      NBD_DIR = nbd_dir,
      OUT_DIR = nbd_compat_dir,
      NBD_TRAITS_TO_PREPARE = ROBUST_TRAITS_TO_RUN,
      NBD_ATAC_W_LIST = ROBUST_ATAC_W_LIST,
      NBD_GROUP_FDR_THRESHOLD = ROBUST_GROUP_FDR_THRESHOLD,
      NBD_SEED_CANDIDATE_COLUMN = "candidate_group_weight_ge0p80"
    ))
  }

  if (step_enabled("neighborhood_rwr") && identical(ROBUST_SEED_MODE, "neighborhood_fdr")) {
    run_script("03_pseudopair_pipeline/05_run_module3_v2_neighborhood_fdr_rwr.R", list(
      OUT_DIR = nbd_compat_dir,
      NBD_TRAITS_TO_RUN = ROBUST_TRAITS_TO_RUN,
      NBD_ATAC_W_LIST = ROBUST_ATAC_W_LIST
    ))
  }

  if (step_enabled("auroc")) {
    auroc_path <- file.path(auroc_dir, auroc_csv)
    old_auroc <- if (file.exists(auroc_path)) {
      tryCatch(fread(auroc_path), error = function(e) NULL)
    } else {
      NULL
    }
    auroc_vars <- list(
      V2_FILES_DIR_NEW = nbd_compat_dir,
      OUT_DIR_NEW = auroc_dir,
      OUT_CSV_NAME_NEW = auroc_csv,
      PATCHED_SCRIPT = file.path(auroc_dir, paste0("eval_auroc.", run_id, ".patched.R")),
      METHOD_LABEL_NEW = ROBUST_METHOD_LABEL
    )
    if (has_local("ROBUST_TARGET_REPLACEMENT")) {
      auroc_vars$TARGET_REPLACEMENT <- ROBUST_TARGET_REPLACEMENT
    }
    if (has_local("ROBUST_PRED_REPLACEMENT")) {
      auroc_vars$PRED_REPLACEMENT <- ROBUST_PRED_REPLACEMENT
    }
    run_script(ROBUST_AUROC_RUNNER, auroc_vars)
    merge_existing_auroc(auroc_path, old_auroc)
  }

  data.table(
    run_id = run_id,
    p_drop = p_drop,
    stress_seed = seed,
    run_root = run_root,
    ga_path = products$ga_path,
    atac_z_path = products$atac_z_path,
    raw_product_mode = products$raw_product_mode,
    nbd_mode = ROBUST_NBD_MODE,
    auroc_csv = file.path(auroc_dir, auroc_csv)
  )
}

summarize_runs <- function(manifest) {
  rows <- list()
  for (i in seq_len(nrow(manifest))) {
    f <- manifest$auroc_csv[i]
    if (!file.exists(f)) next
    dt <- fread(f)
    dt <- dt[method != "baseline_atac"]
    if (nrow(dt) == 0L) next
    dt[, `:=`(
      run_id = manifest$run_id[i],
      p_drop = manifest$p_drop[i],
      stress_seed = manifest$stress_seed[i]
    )]
    rows[[length(rows) + 1L]] <- dt
  }
  if (length(rows) == 0L) {
    warning("No AUROC files found to summarize.")
    return(invisible(NULL))
  }

  all_dt <- rbindlist(rows, use.names = TRUE, fill = TRUE)
  setcolorder(all_dt, c("run_id", "p_drop", "stress_seed"))
  setorder(all_dt, target, atac_w, p_drop, stress_seed)
  fwrite(all_dt, OUT_ALL)

  summary_dt <- all_dt[, .(
    auroc = mean(auroc, na.rm = TRUE),
    sd_auroc = stats::sd(auroc, na.rm = TRUE),
    n_runs = .N,
    n_positive = n_positive[1],
    n_total = n_total[1]
  ), by = .(target, atac_w, p_drop)]
  summary_dt[, se_auroc := sd_auroc / sqrt(n_runs)]
  summary_dt[, lo95 := auroc - 1.96 * se_auroc]
  summary_dt[, hi95 := auroc + 1.96 * se_auroc]
  summary_dt[, method := sprintf("w=%.1f", atac_w)]
  summary_dt[, trait := target]

  plot_dt <- summary_dt[, .(
    trait, target, p_drop, method, atac_w,
    auroc, lo95, hi95, sd_auroc, n_runs, n_positive, n_total
  )]
  setorder(plot_dt, trait, target, atac_w, p_drop)
  fwrite(plot_dt, OUT_PLOT)

  observed_grid <- unique(plot_dt[, .(
    target,
    method,
    p_drop = round(p_drop, 2)
  )])
  expected_grid <- data.table::CJ(
    target = sort(unique(plot_dt$target)),
    method = sort(unique(plot_dt$method)),
    p_drop = round(ROBUST_P_DROP_LIST, 2),
    unique = TRUE
  )
  missing_grid <- fsetdiff(expected_grid, observed_grid)
  if (nrow(missing_grid) > 0L) {
    fwrite(missing_grid, OUT_MISSING_GRID)
    warning("Missing target/method/p_drop combinations. See: ",
            OUT_MISSING_GRID, call. = FALSE)
  } else if (file.exists(OUT_MISSING_GRID)) {
    unlink(OUT_MISSING_GRID)
  }

  if (file.exists(CLEAN_AUROC_CSV)) {
    clean <- fread(CLEAN_AUROC_CSV)
    clean <- clean[!is.na(atac_w)]
    p0 <- plot_dt[p_drop == 0, .(target, atac_w, full_pipeline_p0_auroc = auroc)]
    cmp <- merge(
      p0,
      clean[, .(target, atac_w, current_clean_auroc = auroc)],
      by = c("target", "atac_w"),
      all.x = TRUE
    )
    cmp[, delta_vs_current_clean := full_pipeline_p0_auroc - current_clean_auroc]
    fwrite(cmp, OUT_P0_COMPARE)
    cat("\np_drop=0 comparison against current clean AUROC:\n")
    print(cmp)
  }

  cat("\nPlot-ready full-pipeline data:\n")
  print(plot_dt[, .(
    trait,
    target,
    p_drop,
    method,
    auroc = round(auroc, 4),
    lo95 = round(lo95, 4),
    hi95 = round(hi95, 4),
    n_runs
  )])

  cat("\nWrote:\n")
  cat("  ", OUT_ALL, "\n", sep = "")
  cat("  ", OUT_PLOT, "\n", sep = "")
  cat("  ", OUT_P0_COMPARE, "\n", sep = "")
  if (file.exists(OUT_MISSING_GRID)) cat("  ", OUT_MISSING_GRID, "\n", sep = "")
  invisible(plot_dt)
}

cat(">>> module4 v3 sparse raw ATAC full-pipeline drop-out robustness\n\n")
if (has_local("RUN_STEPS")) {
  cat("Ignoring RUN_STEPS from current R session.\n")
}
if (has_local("ROBUST_RUN_STEPS") && !has_local("ROBUST_RESUME_STEPS")) {
  cat("Using default full workflow; set ROBUST_RESUME_STEPS only for partial resume.\n")
}
cat("ROBUST_RUN_STEPS : ", paste(ROBUST_RUN_STEPS, collapse = ", "), "\n", sep = "")
cat("P_DROP_LIST      : ", paste(ROBUST_P_DROP_LIST, collapse = ", "), "\n", sep = "")
cat("SEED_LIST        : ", paste(ROBUST_SEED_LIST, collapse = ", "), "\n", sep = "")
cat("ATAC_W_LIST      : ", paste(ROBUST_ATAC_W_LIST, collapse = ", "), "\n", sep = "")
cat("NBD_MODE         : ", ROBUST_NBD_MODE, "\n", sep = "")
cat("SEED_MODE        : ", ROBUST_SEED_MODE, "\n", sep = "")
cat("RUN_LABEL        : ", ROBUST_RUN_LABEL, "\n", sep = "")
cat("AUROC_RUNNER     : ", ROBUST_AUROC_RUNNER, "\n", sep = "")
cat("MAP_BATCH_SIZE   : ", ROBUST_SPARSE_MAP_BATCH_SIZE, "\n", sep = "")
cat("Output           : ", ROBUST_DIR, "\n\n", sep = "")

needs_raw_input <- step_enabled("raw_products") &&
  (any(ROBUST_P_DROP_LIST > 0) || !isTRUE(ROBUST_P0_USE_CLEAN_RAW_PRODUCTS))

required_paths <- c(
  CLEAN_GA_PATH,
  CLEAN_ATAC_Z_PATH,
  ATAC_META_FILE,
  ROBUST_GS_FILE,
  CITE_SCORE_RDS,
  file.path(BASE_CITE_EXTRACT_DIR, "cite_trait_counts.mtx"),
  file.path(BASE_CITE_EXTRACT_DIR, "cite_trait_counts_genes.txt"),
  file.path(BASE_CITE_EXTRACT_DIR, "cite_trait_counts_cells.txt")
)
if (needs_raw_input) required_paths <- c(RAW_ATAC_SE_FILE, required_paths)
if (identical(ROBUST_SEED_MODE, "sparse_neighborhood_screen")) {
  required_paths <- c(required_paths, ROBUST_RNA_GROUP_FILE, ROBUST_RNA_GROUP_TEST_FILES)
}
stop_missing(required_paths)

se <- NULL
counts_full <- NULL
gene_peak_obj <- NULL
if (needs_raw_input) {
  cat("[raw] Loading raw ATAC SE and preparing reusable gene-peak map ...\n")
  se <- readRDS(RAW_ATAC_SE_FILE)
  counts_full <- get_counts(se)
  counts_full <- as(counts_full, "CsparseMatrix")
  gene_peak_obj <- get_gene_peak_map(
    peak_gr = SummarizedExperiment::rowRanges(se),
    n_peaks = nrow(counts_full)
  )
  rownames(gene_peak_obj$gene_peak_map) <- gene_peak_obj$gene_names
} else if (step_enabled("raw_products")) {
  cat("[raw] p_drop=0 clean-product run; raw ATAC SE loading skipped.\n")
}

manifest_rows <- list()
for (p_drop in ROBUST_P_DROP_LIST) {
  for (seed in ROBUST_SEED_LIST) {
    run_id <- sprintf("p%.2f_seed%d", p_drop, seed)
    cat("\n", strrep("=", 78), "\n", sep = "")
    cat("Run: ", run_id, "\n", sep = "")
    cat(strrep("=", 78), "\n", sep = "")

    if (isTRUE(ROBUST_SKIP_COMPLETED_CONDITIONS) &&
        condition_is_complete(p_drop, seed)) {
      cat("  [skip] existing AUROC output detected; reusing saved results.\n")
      next
    }

    if (step_enabled("raw_products")) {
      products <- prepare_raw_products(p_drop, seed, se, counts_full, gene_peak_obj)
    } else {
      raw_prefix <- sprintf("p%.2f_seed%d", p_drop, seed)
      products <- list(
        ga_path = if (p_drop == 0 && isTRUE(ROBUST_P0_USE_CLEAN_RAW_PRODUCTS)) {
          CLEAN_GA_PATH
        } else {
          file.path(RAW_DIR, paste0("gene_activity_", raw_prefix, ".rds"))
        },
        atac_z_path = if (p_drop == 0 && isTRUE(ROBUST_P0_USE_CLEAN_RAW_PRODUCTS)) {
          CLEAN_ATAC_Z_PATH
        } else {
          file.path(RAW_DIR, paste0("atac_z_block_", raw_prefix, ".rds"))
        },
        raw_product_mode = "reused_existing_raw_products"
      )
      stop_missing(c(products$ga_path, products$atac_z_path))
    }

    manifest_rows[[run_id]] <- run_one_condition(p_drop, seed, products)
    fwrite(rbindlist(manifest_rows, use.names = TRUE, fill = TRUE), OUT_MANIFEST)
  }
}

manifest <- rbindlist(manifest_rows, use.names = TRUE, fill = TRUE)
fwrite(manifest, OUT_MANIFEST)

if (step_enabled("summarize")) {
  summarize_runs(manifest)
}

cat("\n=== module4 v3 sparse raw ATAC full-pipeline dropout driver finished ===\n")
cat("Output: ", ROBUST_DIR, "\n", sep = "")
