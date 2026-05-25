# Run light-version SCAVENGE TRS with fused KNN graph.
#
# This script changes only the cell-cell graph:
#   original ATAC-only mutualKNN -> fused ATAC+RNA mutualKNN.
# GWAS overlap, background peaks, z-score/seed calculation, and propagation
# parameters are kept aligned with the original blood_small.R workflow.

source("D:/RStudio/work/module2_adaptive_fusion/_bootstrap.R")

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

paths <- module2_paths()
ensure_module2_dirs(paths)

label <- get_active_label()
cutoff <- if (label == "top1" || startsWith(label, "top1disease")) NA_real_ else get_single_cutoff()
npcs_raw <- Sys.getenv("MODULE2_FUSED_NPCS", "auto")
npcs <- if (tolower(npcs_raw) == "auto") {
  read_auto_value(
    file.path(paths$out_dir, sprintf("fused_embedding_auto_npcs_cutoff_%s.txt", label)),
    "fused PCA dimension"
  )
} else {
  as.integer(npcs_raw)
}
k_raw <- Sys.getenv("MODULE2_KNN_K", "auto")
k <- if (tolower(k_raw) == "auto") {
  read_auto_value(
    file.path(paths$out_dir, sprintf("fused_knn_auto_k_cutoff_%s_pca%d.txt", label, npcs)),
    "KNN k"
  )
} else {
  as.integer(k_raw)
}
gamma <- as.numeric(Sys.getenv("MODULE2_RW_GAMMA", "0.05"))
traits <- strsplit(Sys.getenv("MODULE2_TRAITS", "MONO,PLT,LYM"), ",", fixed = TRUE)[[1]]
traits <- trimws(traits)

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
  sprintf("fused_mutualknn_cutoff_%s_pca%d_k%d.rds", label, npcs, k)
)
fused_knn_path <- Sys.getenv("MODULE2_KNN_PATH", default_knn_path)

default_run_label <- sprintf("trs_fused_knn_cutoff_%s_pca%d_k%d", label, npcs, k)
run_label <- Sys.getenv("MODULE2_RUN_LABEL", default_run_label)
out_dir <- file.path(paths$out_dir, run_label)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("========== Module 2: run fused-KNN SCAVENGE ==========\n")
cat("  traits        :", paste(traits, collapse = ", "), "\n")
cat("  fused KNN     :", fused_knn_path, "\n")
cat("  ATAC SE       :", se_path, "\n")
cat("  GWAS          :", gwas_path, "\n")
cat("  bg peaks      :", bg_peaks_path, "\n")
cat("  original TRS  :", original_trs_path, "\n")
cat("  output dir    :", out_dir, "\n")

stopifnot(file.exists(se_path), file.exists(gwas_path), file.exists(bg_peaks_path), file.exists(fused_knn_path))

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
fused_knn <- readRDS(fused_knn_path)
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
  if (bad_n > 0) stop("Failed to parse ", bad_n, " GWAS variants from var column")
  gwas[, `:=`(
    chr = parsed$chr,
    start = parsed$start,
    end = parsed$end
  )]
  gwas[end <= start, end := start + 1L]
}

se_cells <- colnames(SE_atac)
knn_cells <- rownames(fused_knn)
original_cells <- if (!is.null(original_trs) && "cell" %in% colnames(original_trs)) {
  as.character(original_trs$cell)
} else {
  character()
}

matched_seed_fused <- intersect(se_cells, knn_cells)
matched_seed_original <- intersect(se_cells, original_cells)
matched_all <- Reduce(intersect, list(se_cells, knn_cells, if (length(original_cells) > 0) original_cells else se_cells))

alignment_summary <- data.table(
  n_seed = length(se_cells),
  n_fused = length(knn_cells),
  n_original_trs = length(original_cells),
  n_matched_seed_fused = length(matched_seed_fused),
  n_matched_seed_original = length(matched_seed_original),
  n_matched_all = length(matched_all),
  seed_fused_same_order = identical(se_cells, knn_cells)
)
print(alignment_summary)
fwrite(alignment_summary, file.path(out_dir, "cell_alignment_summary.csv"))

if (length(matched_seed_fused) == 0) {
  stop("No overlapping cells between ATAC SE and fused KNN")
}
if (length(matched_seed_fused) < length(se_cells) || length(matched_seed_fused) < length(knn_cells)) {
  warning("ATAC SE and fused KNN do not have identical cell sets; subsetting to intersection.")
}

cell_ids <- se_cells[se_cells %in% matched_seed_fused]
SE_use <- SE_atac[, cell_ids]
fused_knn <- fused_knn[cell_ids, cell_ids]
stopifnot(identical(colnames(SE_use), rownames(fused_knn)))

degree <- Matrix::rowSums(fused_knn)
cat("  fused KNN cells:", nrow(fused_knn), "\n")
cat("  fused KNN isolated:", sum(degree == 0), "(", round(mean(degree == 0) * 100, 3), "% )\n")

all_cell_ids <- cell_ids
graph_cell_ids <- names(degree)[degree > 0]
isolated_cells <- names(degree)[degree == 0]
if (length(graph_cell_ids) == 0) {
  stop("Fused KNN has no non-isolated cells")
}
if (length(isolated_cells) > 0) {
  cat("  removing isolated cells before randomWalk_sparse:", length(isolated_cells), "\n")
}
fused_knn_graph <- fused_knn[graph_cell_ids, graph_cell_ids]
graph_degree <- Matrix::rowSums(fused_knn_graph)
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
  seed_cells_total <- rownames(fused_knn)[seed_idx]
  seed_cells <- intersect(seed_cells_total, rownames(fused_knn_graph))
  seed_isolated <- setdiff(seed_cells_total, seed_cells)
  if (length(seed_cells) == 0) {
    warning("No non-isolated seed cells for trait: ", tr)
    next
  }

  np_score <- randomWalk_sparse(
    intM = fused_knn_graph,
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

  out_file <- file.path(out_dir, paste0("SCAVENGE_fusedKNN_", make.names(tr), ".tsv.gz"))
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
if (length(trs_list) == 0) stop("No trait results were generated")

trs_mat <- Reduce(function(x, y) dplyr::full_join(x, y, by = "cell"), trs_list)
trs_full <- data.frame(cell = all_cell_ids, stringsAsFactors = FALSE) |>
  dplyr::left_join(trs_mat, by = "cell")
trait_cols <- setdiff(colnames(trs_full), "cell")
for (cc in trait_cols) trs_full[[cc]][is.na(trs_full[[cc]])] <- 0

trs_rds <- file.path(out_dir, "TRS_matrix_fusedKNN_selected_traits_fullcells.rds")
trs_tsv <- file.path(out_dir, "TRS_matrix_fusedKNN_selected_traits_fullcells.tsv.gz")
saveRDS(trs_full, trs_rds, compress = FALSE)
data.table::fwrite(trs_full, trs_tsv, sep = "\t")

trait_summary_dt <- rbindlist(trait_summary, use.names = TRUE, fill = TRUE)
summary_out <- file.path(out_dir, "fusedKNN_trait_run_summary.csv")
fwrite(trait_summary_dt, summary_out)

manifest <- list(
  script = "07_run_scavenge_fused_knn.R",
  mode = "fused_knn_only",
  traits = traits,
  cutoff = cutoff,
  npcs = npcs,
  k = k,
  gamma = gamma,
  fused_knn = fused_knn_path,
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
  file.path(out_dir, "fusedKNN_run_manifest.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

cat("\n========== Done ==========\n")
cat("  TRS RDS :", trs_rds, "\n")
cat("  summary :", summary_out, "\n")

