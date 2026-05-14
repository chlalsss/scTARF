# 03_scavenge_baseline/03_run_baseline_traits.R
#
# Forked from module2_wnn_lite_fusion/03_run_scavenge_wnn_lite_knn.R lines
# 200-330 (the trait-loop block; cleaner than blood_small.R because it handles
# isolated cells before randomWalk_sparse). Replaces wnn_knn with the standard
# mutualKNN graph from 02_bg_peaks_and_lsi.R — this is the ATAC-only baseline
# whose TRS we'll compare new dual-modality methods against.
#
# Inputs:
#   <out_dir>/03_scavenge_baseline/SE_atac_gcbiased.rds
#   <out_dir>/03_scavenge_baseline/bg_peaks_n200.rds
#   <out_dir>/03_scavenge_baseline/mutualknn_k30.rds
#   paths$gwas_tsv  CS.trans_all.PP01.strong_atac_overlap.tsv (var/trait/PP)
#
# Outputs (per trait):
#   <out_dir>/03_scavenge_baseline/SCAVENGE_baseline_<trait>.tsv.gz
#                                  (cell, z_score, seed, np_score, TRS, +meta)
#   <out_dir>/03_scavenge_baseline/TRS_matrix_baseline.rds
#                                  cell x trait wide matrix
#   <out_dir>/03_scavenge_baseline/baseline_trait_summary.csv

source(file.path(dirname(sys.frame(1)$ofile %||% "."), "..", "_bootstrap.R"))
required_pkgs(c(
  "data.table", "SummarizedExperiment", "chromVAR", "gchromVAR", "SCAVENGE",
  "Matrix", "dplyr"
))

paths <- m3_paths()
out_dir <- file.path(paths$out_dir, "03_scavenge_baseline")

# --- 1. inputs ------------------------------------------------------------

SE       <- readRDS(file.path(out_dir, "SE_atac_gcbiased.rds"))
bg_peaks <- readRDS(file.path(out_dir, "bg_peaks_n200.rds"))
mknn     <- readRDS(file.path(out_dir, "mutualknn_k30.rds"))
message("Inputs: SE ", nrow(SE), "x", ncol(SE),
        " | bg_peaks ", nrow(bg_peaks), "x", ncol(bg_peaks),
        " | mutualKNN ", nrow(mknn), "x", ncol(mknn))

if (!file.exists(paths$gwas_tsv)) stop("Missing GWAS TSV: ", paths$gwas_tsv)
gwas <- data.table::fread(paths$gwas_tsv)
if ("var" %in% colnames(gwas) && !"variant" %in% colnames(gwas)) {
  data.table::setnames(gwas, "var", "variant")
}
keep_traits <- c("MONO", "PLT", "LYM")
gwas <- gwas[trait %in% keep_traits]
message("GWAS rows after MONO/PLT/LYM filter: ", nrow(gwas))

# Variant -> bed coords
bed_dt <- parse_var_to_bed(gwas$variant)
gwas <- gwas[bed_dt$ok]
bed_dt <- bed_dt[bed_dt$ok]
gwas[, `:=`(chr = bed_dt$chr, start = bed_dt$start, end = bed_dt$end)]
gwas[end <= start, end := start + 1L]

col_pp    <- "PP"
col_trait <- "trait"
gamma     <- 0.05

# --- 2. KNN bookkeeping ---------------------------------------------------

degree <- Matrix::rowSums(mknn)
graph_cells    <- names(degree)[degree > 0]
isolated_cells <- names(degree)[degree == 0]
if (length(graph_cells) == 0) stop("All cells isolated in mutualKNN.")
message("KNN: connected cells = ", length(graph_cells),
        " | isolated cells = ", length(isolated_cells))
if (length(isolated_cells) > 0) {
  data.table::fwrite(
    data.table::data.table(cell = isolated_cells),
    file.path(out_dir, "baseline_isolated_cells.csv")
  )
}
mknn_g <- mknn[graph_cells, graph_cells]
stopifnot(all(Matrix::rowSums(mknn_g) > 0))

# --- 3. trait loop --------------------------------------------------------

trs_list <- vector("list", length(keep_traits))
names(trs_list) <- keep_traits
trait_summary <- list()

for (tr in keep_traits) {
  message("\n[Trait] ", tr)
  dt <- gwas[get(col_trait) == tr]
  if (nrow(dt) == 0) {
    warning("No GWAS rows for trait: ", tr); next
  }

  bed <- dt[, .(
    chr   = as.character(chr),
    start = as.integer(start),
    end   = as.integer(end),
    name  = tr,
    score = as.numeric(get(col_pp))
  )]
  bed[end <= start, end := start + 1L]

  bed_file <- file.path(out_dir, paste0("trait_", make.names(tr), ".PP01.bed"))
  data.table::fwrite(bed, bed_file, sep = "\t", col.names = FALSE)

  trait_import <- gchromVAR::importBedScore(
    SummarizedExperiment::rowRanges(SE), bed_file, colidx = 5
  )
  SE_DEV <- gchromVAR::computeWeightedDeviations(
    object = SE,
    weights = trait_import,
    background_peaks = bg_peaks
  )

  z_score <- t(SummarizedExperiment::assays(SE_DEV)[["z"]]) |> c()
  names(z_score) <- colnames(SE)

  seed_idx     <- SCAVENGE::seedindex(z_score = z_score, percent_cut = 0.05)
  scale_factor <- SCAVENGE::cal_scalefactor(z_score = z_score, percent_cut = 0.01)
  seed_cells_total <- rownames(mknn)[seed_idx]
  seed_cells       <- intersect(seed_cells_total, rownames(mknn_g))
  seed_isolated    <- setdiff(seed_cells_total, seed_cells)
  if (length(seed_cells) == 0) {
    warning("No non-isolated seed cells for trait: ", tr); next
  }

  np_score <- SCAVENGE::randomWalk_sparse(
    intM = mknn_g,
    queryCells = seed_cells,
    gamma = gamma
  )

  omit_idx <- np_score == 0
  np_use   <- np_score[!omit_idx]
  keep     <- names(np_use)
  TRS <- SCAVENGE::capOutlierQuantile(x = np_use, q_ceiling = 0.95) |>
    SCAVENGE::max_min_scale()
  TRS <- TRS * scale_factor

  out_df <- data.frame(
    cell     = keep,
    z_score  = z_score[keep],
    seed     = keep %in% seed_cells,
    np_score = as.numeric(np_use),
    TRS      = as.numeric(TRS),
    stringsAsFactors = FALSE
  )

  meta <- as.data.frame(SummarizedExperiment::colData(SE))
  meta$cell <- rownames(meta)
  out_df <- dplyr::left_join(out_df, meta, by = "cell")

  out_file <- file.path(out_dir, paste0("SCAVENGE_baseline_", make.names(tr), ".tsv.gz"))
  data.table::fwrite(out_df, out_file, sep = "\t")
  message("  wrote ", out_file)

  trs_list[[tr]] <- out_df[, c("cell", "TRS")]
  colnames(trs_list[[tr]]) <- c("cell", make.names(tr))

  trait_summary[[tr]] <- data.table::data.table(
    trait              = tr,
    n_gwas_rows        = nrow(dt),
    n_seed_total       = length(seed_cells_total),
    n_seed_in_graph    = length(seed_cells),
    n_seed_isolated    = length(seed_isolated),
    n_graph_cells      = length(graph_cells),
    n_omit_zero        = sum(omit_idx),
    scale_factor       = scale_factor,
    z_max              = max(z_score, na.rm = TRUE),
    z_p95              = stats::quantile(z_score, 0.95, na.rm = TRUE)
  )
}

# --- 4. merged outputs ---------------------------------------------------

trs_list <- trs_list[!vapply(trs_list, is.null, logical(1))]
if (length(trs_list) > 0) {
  trs_mat <- Reduce(function(x, y) dplyr::full_join(x, y, by = "cell"), trs_list)
  saveRDS(trs_mat, file.path(out_dir, "TRS_matrix_baseline.rds"))
  data.table::fwrite(
    trs_mat, file.path(out_dir, "TRS_matrix_baseline.tsv.gz"), sep = "\t"
  )
}

if (length(trait_summary) > 0) {
  data.table::fwrite(
    data.table::rbindlist(trait_summary),
    file.path(out_dir, "baseline_trait_summary.csv")
  )
}

message("\nBaseline done. Outputs in: ", out_dir)
