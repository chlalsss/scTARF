# Run light-version SCAVENGE TRS with WNN-lite KNN graph.
#
# This script changes only the cell-cell graph:
#   original ATAC-only mutualKNN -> WNN-lite ATAC/RNA-reference mutualKNN.
# GWAS overlap, background peaks, z-score/seed calculation, and propagation
# parameters are kept aligned with the original blood_small.R workflow.

source("D:/RStudio/work/module2_wnn_lite_fusion/_bootstrap.R")

required_pkgs(c(
  "SCAVENGE", "chromVAR", "gchromVAR", "SummarizedExperiment",
  "data.table", "dplyr", "Matrix", "jsonlite"
))
suppressPackageStartupMessages({
  library(SCAVENGE)
  library(chromVAR)
  library(gchromVAR)
  library(SummarizedExperiment)
  library(data.table)
  library(dplyr)
  library(Matrix)
  library(jsonlite)
})

paths <- wnn_lite_paths()
ensure_wnn_lite_dirs(paths)

cutoff <- get_single_cutoff()
label <- cutoff_label(cutoff)
k <- env_integer("WNN_LITE_KNN_K", env_integer("MODULE2_KNN_K", 30))
gamma <- env_numeric("WNN_LITE_RW_GAMMA", env_numeric("MODULE2_RW_GAMMA", 0.05))
traits <- strsplit(Sys.getenv("MODULE2_TRAITS", "MONO,PLT,LYM"), ",", fixed = TRUE)[[1]]
traits <- trimws(traits)
traits <- traits[nzchar(traits)]
if (length(traits) == 0) stop("MODULE2_TRAITS is empty", call. = FALSE)

data_dir <- paths$data_dir
se_path <- Sys.getenv(
  "MODULE2_ATAC_SE",
  file.path(data_dir, "scATAC-Healthy-Hematopoiesis-191120.rds")
)
gwas_path <- Sys.getenv(
  "MODULE2_GWAS",
  file.path(data_dir, "CS.trans_all.PP01.strong_atac_overlap.tsv")
)
original_dir <- Sys.getenv(
  "MODULE2_ORIGINAL_SCAVENGE_DIR",
  file.path(data_dir, "scavenge_out_small_all_traits")
)
bg_peaks_path <- Sys.getenv(
  "MODULE2_BG_PEAKS",
  file.path(original_dir, "bg_peaks_n200.rds")
)
original_trs_path <- Sys.getenv(
  "MODULE2_ORIGINAL_TRS",
  file.path(original_dir, "TRS_matrix_small_all_traits_fullcells.rds")
)

default_knn_path <- file.path(
  paths$out_dir,
  sprintf("wnn_lite_mutualknn_cutoff_%s_k%d.rds", label, k)
)
wnn_knn_path <- Sys.getenv("WNN_LITE_KNN_PATH", Sys.getenv("MODULE2_KNN_PATH", default_knn_path))
diagnostics_path <- file.path(
  paths$out_dir,
  sprintf("wnn_lite_knn_diagnostics_cutoff_%s_k%d.csv", label, k)
)

default_run_label <- sprintf("trs_wnn_lite_cutoff_%s_k%d", label, k)
run_label <- Sys.getenv("WNN_LITE_RUN_LABEL", Sys.getenv("MODULE2_RUN_LABEL", default_run_label))
out_dir <- file.path(paths$out_dir, run_label)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("========== Module 2 WNN-lite: run SCAVENGE ==========\n")
cat("  traits        :", paste(traits, collapse = ", "), "\n")
cat("  WNN-lite KNN  :", wnn_knn_path, "\n")
cat("  diagnostics   :", diagnostics_path, "\n")
cat("  ATAC SE       :", se_path, "\n")
cat("  GWAS          :", gwas_path, "\n")
cat("  bg peaks      :", bg_peaks_path, "\n")
cat("  original TRS  :", original_trs_path, "\n")
cat("  gamma         :", gamma, "\n")
cat("  output dir    :", out_dir, "\n")

stopifnot(file.exists(se_path), file.exists(gwas_path), file.exists(bg_peaks_path), file.exists(wnn_knn_path))

if (file.exists(diagnostics_path)) {
  graph_diag <- fread(diagnostics_path)
  gate <- env_numeric("WNN_LITE_MAX_PCT_ISOLATED", 0.01)
  if ("pct_isolated" %in% colnames(graph_diag) && is.finite(graph_diag$pct_isolated[1])) {
    cat("  graph pct_isolated:", graph_diag$pct_isolated[1], "\n")
    if (graph_diag$pct_isolated[1] > gate) {
      warning(
        "WNN-lite graph pct_isolated exceeds gate (",
        graph_diag$pct_isolated[1], " > ", gate,
        "). TRS will still run after removing isolated cells."
      )
    }
  }
} else {
  graph_diag <- NULL
  warning("Missing WNN-lite diagnostics CSV: ", diagnostics_path)
}

parse_var_to_bed <- function(v) {
  v <- as.character(v)
  m <- regexec("^([0-9]+|X|Y|MT|M):(\\d+)_([^_]+)_([^_]+)$", v, ignore.case = TRUE)
  r <- regmatches(v, m)
  ok <- lengths(r) == 5
  chr <- start <- end <- rep(NA_character_, length(v))

  if (any(ok)) {
    raw_chr <- toupper(vapply(r[ok], `[[`, "", 2))
    raw_chr[raw_chr == "M"] <- "MT"
    chr[ok] <- paste0("chr", raw_chr)
    pos <- as.integer(vapply(r[ok], `[[`, "", 3))
    start[ok] <- as.character(pmax(pos - 1L, 0L))
    end[ok] <- as.character(pos)
  }

  data.table(
    chr = chr,
    start = as.integer(start),
    end = as.integer(end),
    ok = !is.na(chr) & !is.na(start) & !is.na(end)
  )
}

cat("\n========== Load inputs ==========\n")
SE_atac <- readRDS(se_path)
wnn_knn <- readRDS(wnn_knn_path)
bg_peaks <- readRDS(bg_peaks_path)
gwas <- data.table::fread(gwas_path)
original_trs <- if (file.exists(original_trs_path)) readRDS(original_trs_path) else NULL

col_var <- "var"
col_trait <- "trait"
col_pp <- "PP"
stopifnot(all(c(col_var, col_trait, col_pp) %in% names(gwas)))

if (!all(c("chr", "start", "end") %in% names(gwas))) {
  parsed <- parse_var_to_bed(gwas[[col_var]])
  bad_n <- sum(!parsed$ok)
  if (bad_n > 0) stop("Failed to parse ", bad_n, " GWAS variants from var column", call. = FALSE)
  gwas[, `:=`(
    chr = parsed$chr,
    start = parsed$start,
    end = parsed$end
  )]
  gwas[end <= start, end := start + 1L]
}

if (nrow(wnn_knn) != ncol(wnn_knn)) {
  stop("WNN-lite KNN must be a square matrix", call. = FALSE)
}
if (is.null(rownames(wnn_knn)) || is.null(colnames(wnn_knn))) {
  stop("WNN-lite KNN must have rownames and colnames", call. = FALSE)
}

se_cells <- colnames(SE_atac)
knn_cells <- rownames(wnn_knn)
original_cells <- if (!is.null(original_trs) && "cell" %in% colnames(original_trs)) {
  as.character(original_trs$cell)
} else {
  character()
}

matched_seed_wnn <- intersect(se_cells, knn_cells)
matched_seed_original <- intersect(se_cells, original_cells)
matched_all <- Reduce(intersect, list(se_cells, knn_cells, if (length(original_cells) > 0) original_cells else se_cells))

alignment_summary <- data.table(
  n_seed = length(se_cells),
  n_wnn_lite = length(knn_cells),
  n_original_trs = length(original_cells),
  n_matched_seed_wnn_lite = length(matched_seed_wnn),
  n_matched_seed_original = length(matched_seed_original),
  n_matched_all = length(matched_all),
  seed_wnn_lite_same_order = identical(se_cells, knn_cells)
)
print(alignment_summary)
fwrite(alignment_summary, file.path(out_dir, "cell_alignment_summary.csv"))

if (length(matched_seed_wnn) == 0) {
  stop("No overlapping cells between ATAC SE and WNN-lite KNN", call. = FALSE)
}
if (length(matched_seed_wnn) < length(se_cells) || length(matched_seed_wnn) < length(knn_cells)) {
  warning("ATAC SE and WNN-lite KNN do not have identical cell sets; subsetting to intersection.")
}

cell_ids <- se_cells[se_cells %in% matched_seed_wnn]
SE_use <- SE_atac[, cell_ids]
wnn_knn <- wnn_knn[cell_ids, cell_ids]
stopifnot(identical(colnames(SE_use), rownames(wnn_knn)))

degree <- Matrix::rowSums(wnn_knn)
cat("  WNN-lite KNN cells:", nrow(wnn_knn), "\n")
cat("  WNN-lite isolated :", sum(degree == 0), "(", round(mean(degree == 0) * 100, 3), "% )\n")

all_cell_ids <- cell_ids
graph_cell_ids <- names(degree)[degree > 0]
isolated_cells <- names(degree)[degree == 0]
if (length(graph_cell_ids) == 0) {
  stop("WNN-lite KNN has no non-isolated cells", call. = FALSE)
}
if (length(isolated_cells) > 0) {
  cat("  removing isolated cells before randomWalk_sparse:", length(isolated_cells), "\n")
  fwrite(
    data.table(cell = isolated_cells),
    file.path(out_dir, "wnnLite_isolated_cells_removed_before_random_walk.csv")
  )
}
wnn_knn_graph <- wnn_knn[graph_cell_ids, graph_cell_ids]
graph_degree <- Matrix::rowSums(wnn_knn_graph)
stopifnot(all(graph_degree > 0))

cat("\n========== Run traits ==========\n")
trs_list <- vector("list", length(traits))
names(trs_list) <- traits
trait_summary <- list()

for (tr in traits) {
  cat("\n[Trait]", tr, "\n")
  dt <- gwas[get(col_trait) == tr]
  if (nrow(dt) == 0) {
    warning("No GWAS rows for trait: ", tr)
    next
  }

  bed <- dt[, .(
    chr = as.character(chr),
    start = as.integer(start),
    end = as.integer(end),
    name = tr,
    score = as.numeric(get(col_pp))
  )]
  bed[end <= start, end := start + 1L]

  bed_file <- file.path(out_dir, paste0("trait_", make.names(tr), ".PP001.bed"))
  data.table::fwrite(bed, bed_file, sep = "\t", col.names = FALSE)

  trait_import <- importBedScore(rowRanges(SE_use), bed_file, colidx = 5)
  SE_DEV <- computeWeightedDeviations(
    object = SE_use,
    weights = trait_import,
    background_peaks = bg_peaks
  )

  z_score <- t(SummarizedExperiment::assays(SE_DEV)[["z"]]) |> c()
  names(z_score) <- colnames(SE_use)

  seed_idx <- seedindex(z_score = z_score, percent_cut = 0.05)
  scale_factor <- cal_scalefactor(z_score = z_score, percent_cut = 0.01)
  seed_cells_total <- rownames(wnn_knn)[seed_idx]
  seed_cells <- intersect(seed_cells_total, rownames(wnn_knn_graph))
  seed_isolated <- setdiff(seed_cells_total, seed_cells)
  if (length(seed_cells) == 0) {
    warning("No non-isolated seed cells for trait: ", tr)
    next
  }

  np_score <- randomWalk_sparse(
    intM = wnn_knn_graph,
    queryCells = seed_cells,
    gamma = gamma
  )

  omit_idx <- np_score == 0
  np_score_use <- np_score[!omit_idx]
  keep_cells <- names(np_score_use)
  TRS <- capOutlierQuantile(x = np_score_use, q_ceiling = 0.95) |> max_min_scale()
  TRS <- TRS * scale_factor

  out_df <- data.frame(
    cell = keep_cells,
    z_score = z_score[keep_cells],
    seed = keep_cells %in% seed_cells,
    np_score = as.numeric(np_score_use),
    TRS = as.numeric(TRS),
    stringsAsFactors = FALSE
  )

  meta <- as.data.frame(colData(SE_use))
  meta$cell <- rownames(meta)
  out_df <- dplyr::left_join(out_df, meta, by = "cell")

  out_file <- file.path(out_dir, paste0("SCAVENGE_wnnLite_", make.names(tr), ".tsv.gz"))
  data.table::fwrite(out_df, out_file, sep = "\t")

  trs_list[[tr]] <- out_df[, c("cell", "TRS")]
  colnames(trs_list[[tr]]) <- c("cell", make.names(tr))

  trait_summary[[tr]] <- data.table(
    trait = tr,
    n_gwas_rows = nrow(dt),
    n_seed_total = length(seed_cells_total),
    n_seed_in_graph = length(seed_cells),
    n_seed_isolated = length(seed_isolated),
    n_graph_cells = length(graph_cell_ids),
    n_isolated_cells = length(isolated_cells),
    n_omit = sum(omit_idx),
    pct_omit = mean(omit_idx),
    scale_factor = scale_factor,
    mean_np_score = mean(np_score_use),
    mean_TRS = mean(TRS),
    median_TRS = median(TRS),
    output_file = out_file
  )

  cat("  seed:", length(seed_cells), "omit:", sum(omit_idx), "\n")
  cat("  saved:", out_file, "\n")
}

trs_list <- trs_list[!vapply(trs_list, is.null, logical(1))]
if (length(trs_list) == 0) stop("No trait results were generated", call. = FALSE)

trs_mat <- Reduce(function(x, y) dplyr::full_join(x, y, by = "cell"), trs_list)
trs_full <- data.frame(cell = all_cell_ids, stringsAsFactors = FALSE) |>
  dplyr::left_join(trs_mat, by = "cell")
trait_cols <- setdiff(colnames(trs_full), "cell")
for (cc in trait_cols) trs_full[[cc]][is.na(trs_full[[cc]])] <- 0

trs_rds <- file.path(out_dir, "TRS_matrix_wnnLite_selected_traits_fullcells.rds")
trs_tsv <- file.path(out_dir, "TRS_matrix_wnnLite_selected_traits_fullcells.tsv.gz")
saveRDS(trs_full, trs_rds, compress = FALSE)
data.table::fwrite(trs_full, trs_tsv, sep = "\t")

trait_summary_dt <- rbindlist(trait_summary, use.names = TRUE, fill = TRUE)
summary_out <- file.path(out_dir, "wnnLite_trait_run_summary.csv")
fwrite(trait_summary_dt, summary_out)

manifest <- list(
  script = "03_run_scavenge_wnn_lite_knn.R",
  mode = "wnn_lite_knn_only",
  traits = traits,
  cutoff = cutoff,
  k = k,
  gamma = gamma,
  wnn_lite_knn = wnn_knn_path,
  wnn_lite_diagnostics = diagnostics_path,
  atac_se = se_path,
  gwas = gwas_path,
  bg_peaks = bg_peaks_path,
  original_trs = original_trs_path,
  cell_alignment_summary = file.path(out_dir, "cell_alignment_summary.csv"),
  trait_summary = summary_out,
  trs_rds = trs_rds,
  trs_tsv = trs_tsv,
  created_at = as.character(Sys.time())
)
jsonlite::write_json(
  manifest,
  file.path(out_dir, "wnnLite_run_manifest.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

cat("\n========== Done ==========\n")
cat("  TRS RDS :", trs_rds, "\n")
cat("  summary :", summary_out, "\n")

