# module4 robustness: dynamic ATAC-only baseline under raw read drop-out
# -----------------------------------------------------------------------------
# Why this script exists:
#   The previous "baseline" line in the robustness plot was a fixed reference
#   from the original non-degraded SCAVENGE output. That is not the right line
#   when the x-axis is ATAC read drop-out.
#
# What this script computes:
#   For each p_drop x stress seed x trait, it recomputes an ATAC-only baseline
#   on the same degraded ATAC data used by the robustness sweep:
#
#     degraded raw ATAC counts
#       -> gchromVAR trait z-score
#       -> ATAC-derived seed cells, top 5%
#       -> reuse existing ATAC-only graph for this p_drop/seed, w=1.0
#       -> RWR/TRS using the same small RWR helper as module4 robustness
#       -> AUROC by posthoc cell type
#
# Important:
#   Cell type is used only for the final AUROC check, not for seed or RWR.
#
# Run in RStudio:
#   setwd("D:/RStudio/work/module4_traitaware_pseudopair")
#   source("05_robustness/05_run_top3_raw_atac_dropout_dynamic_baseline.R")
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Matrix)
  library(SummarizedExperiment)
  library(gchromVAR)
  library(SCAVENGE)
  library(data.table)
  library(pROC)
})

if (!exists("ROOT")) {
  ROOT <- "D:/RStudio/work/module4_traitaware_pseudopair"
}
if (!exists("ROBUST_DIR")) {
  ROBUST_DIR <- file.path(ROOT, "05_robustness/output/top3_raw_atac_dropout")
}
if (!exists("BASELINE_DIR")) {
  BASELINE_DIR <- file.path(ROBUST_DIR, "dynamic_baseline")
}
dir.create(BASELINE_DIR, recursive = TRUE, showWarnings = FALSE)

if (!exists("RAW_ATAC_SE_FILE")) {
  RAW_ATAC_SE_FILE <- "D:/SCAVENGE-reproducibility-main/blood_data/scATAC-Healthy-Hematopoiesis-191120.rds"
}
if (!exists("BG_PEAKS_FILE")) {
  BG_PEAKS_FILE <- "D:/SCAVENGE-reproducibility-main/blood_data/scavenge_out_small_all_traits/bg_peaks_n200.rds"
}
if (!exists("GWAS_FILE")) {
  GWAS_FILE <- "D:/SCAVENGE-reproducibility-main/blood_data/CS.trans_all.PP01.strong_atac_overlap.tsv"
}
if (!exists("META_FILE")) {
  META_FILE <- file.path(ROOT, "04_evaluation/inputs/granja_cell_meta.rds")
}

if (!exists("TRAITS")) TRAITS <- c("MONO", "RDW", "RBC")
if (!exists("TARGETS")) {
  TARGETS <- list(
    MONO = c("11_CD14.Mono.1", "12_CD14.Mono.2"),
    RDW = c("02_Early.Eryth", "03_Late.Eryth", "04_Early.Baso"),
    RBC = c("02_Early.Eryth", "03_Late.Eryth", "04_Early.Baso")
  )
}

if (!exists("P_DROP_LIST")) P_DROP_LIST <- seq(0.0, 0.9, by = 0.1)
if (!exists("SEED_LIST")) SEED_LIST <- c(1L, 2L, 3L)

GAMMA <- 0.05
TOL <- 1e-6
MAX_ITER_RWR <- 200
CAP_Q <- 0.95
SEED_TOP_FRAC <- 0.05

stop_missing <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0L) {
    stop("Missing required file(s):\n  ", paste(missing, collapse = "\n  "),
         call. = FALSE)
  }
}

parse_var_to_bed <- function(variant, chr_prefix = "chr") {
  m <- regexec("^([^:]+):([0-9]+)_([^_]+)_([^_]+)$", variant)
  parts <- regmatches(variant, m)
  ok <- lengths(parts) == 5L
  chr <- rep(NA_character_, length(variant))
  pos <- rep(NA_integer_, length(variant))
  chr[ok] <- vapply(parts[ok], `[`, character(1), 2)
  pos[ok] <- as.integer(vapply(parts[ok], `[`, character(1), 3))
  if (nzchar(chr_prefix)) {
    chr <- ifelse(!is.na(chr) & !grepl("^chr", chr), paste0(chr_prefix, chr), chr)
  }
  data.table(
    ok = ok & !is.na(pos),
    chr = chr,
    start = pmax(pos - 1L, 0L),
    end = pos
  )
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

rwr <- function(A_norm, seed_vec, gamma = GAMMA, tol = TOL,
                max_iter = MAX_ITER_RWR) {
  if (sum(seed_vec) <= 0) stop("seed_vec has no positive entries.", call. = FALSE)
  s <- seed_vec / sum(seed_vec)
  x <- s
  for (iter in seq_len(max_iter)) {
    x_new <- gamma * s + (1 - gamma) * as.numeric(A_norm %*% x)
    delta <- sum(abs(x_new - x))
    x <- x_new
    if (delta < tol) return(x)
  }
  x
}

post_process <- function(trs, cap_q = CAP_Q) {
  cap <- as.numeric(stats::quantile(trs, cap_q, na.rm = TRUE))
  trs_c <- pmin(trs, cap)
  rng <- range(trs_c, na.rm = TRUE)
  if (is.finite(rng[2] - rng[1]) && rng[2] - rng[1] > 0) {
    (trs_c - rng[1]) / (rng[2] - rng[1])
  } else {
    rep(0, length(trs))
  }
}

compute_auroc <- function(predictor, labels) {
  keep <- is.finite(predictor) & !is.na(labels)
  predictor <- predictor[keep]
  labels <- labels[keep]
  if (length(predictor) < 20L) return(NA_real_)
  if (sum(labels == 1L) < 10L || sum(labels == 0L) < 10L) return(NA_real_)
  if (!is.finite(stats::var(predictor)) || stats::var(predictor) < 1e-12) {
    return(NA_real_)
  }
  roc_obj <- suppressMessages(pROC::roc(
    labels, predictor,
    quiet = TRUE,
    direction = "<",
    levels = c(0L, 1L)
  ))
  as.numeric(pROC::auc(roc_obj))
}

make_trait_imports <- function(SE, gwas) {
  rr_chr <- as.character(GenomicRanges::seqnames(SummarizedExperiment::rowRanges(SE)))
  chr_prefix <- if (any(grepl("^chr", rr_chr))) "chr" else ""

  imports <- list()
  bed_dir <- file.path(BASELINE_DIR, "trait_beds")
  dir.create(bed_dir, recursive = TRUE, showWarnings = FALSE)

  for (trait in TRAITS) {
    dt <- gwas[gwas$trait == trait]
    bed_dt <- parse_var_to_bed(dt$variant, chr_prefix = chr_prefix)
    dt <- dt[bed_dt$ok]
    bed_dt <- bed_dt[bed_dt$ok]
    if (nrow(dt) == 0L) stop("No GWAS variants after parsing for ", trait)

    bed <- data.table(
      chr = bed_dt$chr,
      start = as.integer(bed_dt$start),
      end = as.integer(bed_dt$end),
      name = trait,
      score = as.numeric(dt$PP)
    )
    bed[end <= start, end := start + 1L]
    bed_file <- file.path(bed_dir, paste0("trait_", trait, ".PP001.bed"))
    fwrite(bed, bed_file, sep = "\t", col.names = FALSE)

    imports[[trait]] <- gchromVAR::importBedScore(
      SummarizedExperiment::rowRanges(SE),
      bed_file,
      colidx = 5
    )
    cat("  trait import ", trait, ": ", nrow(bed), " variant rows\n", sep = "")
  }

  imports
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
  stop("Cannot subset trait_import with keep_peaks; unexpected dimensions: ",
       paste(dim(trait_import), collapse = " x "), call. = FALSE)
}

prepare_degraded_inputs <- function(SE_template, counts_full, p_keep, seed,
                                    p_drop, bg_peaks_full) {
  cat("  sparsify counts p_keep=", p_keep, " ...\n", sep = "")
  counts_sp <- sparsify_binomial(counts_full, p_keep, seed)

  keep_peaks <- Matrix::rowSums(counts_sp) > 0
  n_drop_peaks <- sum(!keep_peaks)
  if (n_drop_peaks > 0L) {
    cat("  filtering all-zero peaks after drop-out: ",
        n_drop_peaks, " / ", length(keep_peaks), "\n", sep = "")
  }

  counts_sp <- counts_sp[keep_peaks, , drop = FALSE]
  SE_sp <- SE_template[keep_peaks, ]
  SummarizedExperiment::assay(SE_sp, "counts") <- counts_sp

  if (n_drop_peaks == 0L) {
    bg_current <- bg_peaks_full
  } else {
    bg_file <- file.path(
      BASELINE_DIR,
      sprintf("bg_peaks_p%.2f_seed%d_nzpeaks.rds", p_drop, seed)
    )
    if (file.exists(bg_file)) {
      bg_current <- readRDS(bg_file)
    } else {
      bias <- SummarizedExperiment::rowData(SE_sp)$GC
      if (is.null(bias)) {
        stop("Cannot recompute background peaks: rowData(SE_sp)$GC is missing.",
             call. = FALSE)
      }
      cat("  recomputing background peaks for nonzero peak set ...\n")
      bg_current <- chromVAR::getBackgroundPeaks(
        SE_sp,
        bias = bias,
        niterations = ncol(bg_peaks_full)
      )
      saveRDS(bg_current, bg_file)
    }
  }

  list(SE_sp = SE_sp, keep_peaks = keep_peaks, bg_peaks = bg_current)
}

compute_z <- function(SE_sp, trait_import, bg_peaks_current) {
  SE_DEV <- gchromVAR::computeWeightedDeviations(
    object = SE_sp,
    weights = trait_import,
    background_peaks = bg_peaks_current
  )
  z_score <- t(SummarizedExperiment::assays(SE_DEV)[["z"]]) |> c()
  names(z_score) <- colnames(SE_sp)
  z_score
}

load_meta <- function(path) {
  meta <- readRDS(path)
  if (inherits(meta, "data.table")) meta <- as.data.frame(meta)
  if (is.null(rownames(meta)) ||
      identical(rownames(meta), as.character(seq_len(nrow(meta))))) {
    if ("barcode" %in% colnames(meta)) rownames(meta) <- meta$barcode
  }
  meta
}

cat(">>> Recompute dynamic ATAC-only baseline under raw read drop-out\n\n")
cat("P_DROP_LIST: ", paste(P_DROP_LIST, collapse = ", "), "\n", sep = "")
cat("SEED_LIST  : ", paste(SEED_LIST, collapse = ", "), "\n\n", sep = "")

stop_missing(c(
  RAW_ATAC_SE_FILE,
  BG_PEAKS_FILE,
  GWAS_FILE,
  META_FILE
))

cat("[1/5] Loading raw ATAC SE, bg peaks, GWAS, and metadata ...\n")
SE <- readRDS(RAW_ATAC_SE_FILE)
counts_full <- SummarizedExperiment::assay(SE, "counts")
bg_peaks_full <- readRDS(BG_PEAKS_FILE)
gwas <- fread(GWAS_FILE)
if ("var" %in% names(gwas) && !"variant" %in% names(gwas)) {
  setnames(gwas, "var", "variant")
}
gwas <- gwas[trait %in% TRAITS]
meta <- load_meta(META_FILE)

cat("  SE peak x cell: ", nrow(SE), " x ", ncol(SE), "\n", sep = "")
cat("  bg peaks dim  : ", paste(dim(bg_peaks_full), collapse = " x "), "\n", sep = "")
cat("  GWAS rows     : ", nrow(gwas), "\n\n", sep = "")

cat("[2/5] Precomputing trait peak weights ...\n")
trait_imports <- make_trait_imports(SE, gwas)

eval_rows <- list()
row_i <- 0L

cat("\n[3/5] Running dynamic baseline ...\n")
for (p_drop in P_DROP_LIST) {
  p_keep <- 1 - p_drop
  for (sd in SEED_LIST) {
    cat("\n=== p_drop=", p_drop, " seed=", sd, " ===\n", sep = "")

    z_prefix <- sprintf("p%.2f_seed%d", p_drop, sd)
    degraded <- NULL

    for (trait in TRAITS) {
      z_file <- file.path(BASELINE_DIR, sprintf("baseline_z_%s_%s.rds", z_prefix, trait))
      trs_file <- file.path(BASELINE_DIR, sprintf("baseline_TRS_%s_%s.tsv.gz", z_prefix, trait))
      knn_file <- file.path(ROBUST_DIR, sprintf("knn_p%.2f_seed%d_MONO_w1.0.rds", p_drop, sd))
      stop_missing(knn_file)

      if (file.exists(trs_file)) {
        cat("  skip exists: ", basename(trs_file), "\n", sep = "")
        trs_dt <- fread(trs_file)
      } else {
        if (file.exists(z_file)) {
          z_score <- readRDS(z_file)
        } else {
          if (is.null(degraded)) {
            degraded <- prepare_degraded_inputs(
              SE_template = SE,
              counts_full = counts_full,
              p_keep = p_keep,
              seed = sd,
              p_drop = p_drop,
              bg_peaks_full = bg_peaks_full
            )
          }
          cat("  compute gchromVAR z: ", trait, "\n", sep = "")
          trait_import_current <- subset_trait_import(
            trait_imports[[trait]],
            degraded$keep_peaks
          )
          z_score <- compute_z(
            degraded$SE_sp,
            trait_import_current,
            degraded$bg_peaks
          )
          saveRDS(z_score, z_file)
        }

        A_norm <- readRDS(knn_file)
        graph_cells <- rownames(A_norm)
        z_graph <- z_score[graph_cells]
        seed_idx <- SCAVENGE::seedindex(z_score = z_graph, percent_cut = SEED_TOP_FRAC)
        seed_cells <- names(z_graph)[seed_idx]
        seed_cells <- intersect(seed_cells, graph_cells)
        if (length(seed_cells) == 0L) {
          stop("No dynamic baseline seed cells for ", trait,
               " p_drop=", p_drop, " seed=", sd, call. = FALSE)
        }

        seed_vec <- numeric(length(graph_cells))
        names(seed_vec) <- graph_cells
        seed_vec[seed_cells] <- 1

        trs_raw <- rwr(A_norm, seed_vec)
        trs <- post_process(trs_raw)
        names(trs) <- graph_cells

        trs_dt <- data.table(
          barcode = names(trs),
          TRS = as.numeric(trs),
          z_score = as.numeric(z_graph[names(trs)]),
          seed = names(trs) %in% seed_cells
        )
        fwrite(trs_dt, trs_file, sep = "\t", compress = "gzip")
        cat("  wrote ", basename(trs_file), "; n_seed=", length(seed_cells), "\n", sep = "")
      }

      common_eval <- intersect(trs_dt$barcode, rownames(meta))
      lbl <- as.integer(as.character(meta[common_eval, "cell_type"]) %in% TARGETS[[trait]])
      score <- trs_dt$TRS[match(common_eval, trs_dt$barcode)]
      auc <- compute_auroc(score, lbl)

      row_i <- row_i + 1L
      eval_rows[[row_i]] <- data.table(
        p_drop = p_drop,
        seed = sd,
        trait = trait,
        target = trait,
        method = "dynamic_baseline_atac_seed",
        auroc = auc,
        n_positive = sum(lbl),
        n_total = length(lbl)
      )
      cat(sprintf("  AUROC baseline %-4s = %.4f\n", trait, auc))
    }
  }
}

cat("\n[4/5] Writing baseline tables ...\n")
eval_dt <- rbindlist(eval_rows)
eval_dt <- unique(eval_dt, by = c("p_drop", "seed", "trait", "target", "method"))
setorder(eval_dt, p_drop, seed, trait)

eval_file <- file.path(BASELINE_DIR, "dynamic_baseline_eval.csv")
summary_file <- file.path(BASELINE_DIR, "dynamic_baseline_summary.csv")
fwrite(eval_dt, eval_file)

summary_dt <- eval_dt[, .(
  baseline_mean_auroc = mean(auroc, na.rm = TRUE),
  baseline_sd_auroc = stats::sd(auroc, na.rm = TRUE),
  n_runs = .N,
  n_positive = n_positive[1],
  n_total = n_total[1]
), by = .(p_drop, trait, target)]
summary_dt[, baseline_se_auroc := baseline_sd_auroc / sqrt(n_runs)]
summary_dt[, baseline_lo95 := baseline_mean_auroc - 1.96 * baseline_se_auroc]
summary_dt[, baseline_hi95 := baseline_mean_auroc + 1.96 * baseline_se_auroc]
setorder(summary_dt, trait, p_drop)
fwrite(summary_dt, summary_file)

cat("\n[5/5] Dynamic baseline summary:\n")
print(summary_dt[, .(
  trait,
  p_drop,
  baseline_mean_auroc = round(baseline_mean_auroc, 4),
  n_runs
)])

cat("\nWrote:\n")
cat("  ", eval_file, "\n", sep = "")
cat("  ", summary_file, "\n", sep = "")
