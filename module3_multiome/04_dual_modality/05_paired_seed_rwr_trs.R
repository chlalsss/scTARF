# 04_dual_modality/05_paired_seed_rwr_trs.R
#
# §2.6 配对一致性 seed selection + §2.7 带权 RWR + TRS 后处理.
#
# §2.6 seed:
#   primary: seed = (z_atac top 5%) ∩ (z_rna top 5%)
#   fallback (intersect < 200):
#     rank_atac = rank(z_atac)/N; rank_rna = rank(z_rna)/N
#     combined = rank_atac + rank_rna; seed = top 5% by combined
#   合规性: z_atac 来自 g-chromVAR (variant→peak), z_rna 来自 SCDRS
#   (variant→peak→gene→expr); 都不读 cell_type, 双重 GWAS 证据 ≠ 双重答案.
#
# §2.7 RWR + TRS:
#   randomWalk_sparse on §2.5's row-normalized weighted KNN, γ = 0.05
#   capOutlierQuantile(q=0.95) → max_min_scale → × scale_factor
#
# Inputs:
#   <out>/03_scavenge_baseline/SE_atac_gcbiased.rds       SE for chromVAR
#   <out>/03_scavenge_baseline/bg_peaks_n200.rds          bg peaks
#   <out>/03_scavenge_baseline/trait_<trait>.PP01.bed     reused PP-weighted bed
#   paths$gwas_tsv                                         backup if bed missing
#   <out>/04_dual_modality/zrna_<trait>.rds                §2.3 cell-level z_rna
#   <out>/04_dual_modality/knn_graph_<trait>.rds           §2.5 weighted KNN
#
# Outputs:
#   <out>/04_dual_modality/zatac_<trait>.rds               cell-level z_atac (named)
#   <out>/04_dual_modality/seed_paired_<trait>.rds         seed cell barcodes (chr vec)
#   <out>/04_dual_modality/TRS_dual_<trait>.tsv.gz         per-cell long table
#   <out>/04_dual_modality/TRS_dual_matrix.rds             cell × trait wide matrix
#   <out>/04_dual_modality/paired_seed_summary.csv         per-trait seed bookkeeping
#
# Red line: NO cell_type / cell_meta in this script.

source(file.path(dirname(sys.frame(1)$ofile %||% "."), "..", "_bootstrap.R"))
required_pkgs(c(
  "data.table", "Matrix",
  "SummarizedExperiment", "chromVAR", "gchromVAR", "SCAVENGE", "dplyr"
))

paths   <- m3_paths()
out_dir <- file.path(paths$out_dir, "04_dual_modality")
base_dir <- file.path(paths$out_dir, "03_scavenge_baseline")

# --- locked params --------------------------------------------------------
TRAITS         <- c("MONO", "PLT", "LYM")
SEED_TOP_FRAC  <- 0.05
# 2026-05-04: lowered from 200 to 50 after BMMC empirics. With G_trait
# constructed by GWAS-only (no cell-type prior), z_atac top 5% and z_rna
# top 5% overlap at chance (~178 cells). The rank-sum fallback dilutes
# seeds toward "generally-active" cell types (T cells, monocytes) and
# wipes out trait specificity. Using the small but trait-faithful
# intersection directly gives RWR a high-confidence seed core that
# subsequent KNN propagation expands.
INTERSECT_FLOOR<- 50L
GAMMA          <- 0.05
TRS_Q_CEILING  <- 0.95
SCALE_PCT      <- 0.01

# --- 1. shared chromVAR inputs -------------------------------------------

SE       <- readRDS(file.path(base_dir, "SE_atac_gcbiased.rds"))
bg_peaks <- readRDS(file.path(base_dir, "bg_peaks_n200.rds"))
message("SE: ", nrow(SE), " peaks x ", ncol(SE), " cells")

# --- 2. per-trait pipeline -----------------------------------------------

trs_long_list <- list()
trs_wide_list <- list()
sum_rows      <- list()

for (tr in TRAITS) {
  message("\n[Trait] ", tr)

  # 2a. recompute z_atac via g-chromVAR (don't reuse baseline TSV — that one
  # drops cells with np_score==0 and we want all cells for seed selection).
  bed_file <- file.path(base_dir, paste0("trait_", tr, ".PP01.bed"))
  if (!file.exists(bed_file)) stop("Missing trait bed: ", bed_file)
  trait_imp <- gchromVAR::importBedScore(
    SummarizedExperiment::rowRanges(SE), bed_file, colidx = 5
  )
  SE_DEV <- gchromVAR::computeWeightedDeviations(
    object = SE, weights = trait_imp, background_peaks = bg_peaks
  )
  z_atac <- t(SummarizedExperiment::assays(SE_DEV)[["z"]]) |> c()
  names(z_atac) <- colnames(SE)
  saveRDS(z_atac, file.path(out_dir, paste0("zatac_", tr, ".rds")))
  message(sprintf("  z_atac: n=%d, q95=%.3f, max=%.3f",
                  length(z_atac),
                  as.numeric(stats::quantile(z_atac, 0.95, na.rm = TRUE)),
                  max(z_atac, na.rm = TRUE)))

  # 2b. load z_rna and align to common cell set
  z_rna <- readRDS(file.path(out_dir, paste0("zrna_", tr, ".rds")))
  common <- intersect(names(z_atac), names(z_rna))
  if (length(common) < 0.99 * length(z_atac)) {
    warning("Only ", length(common), " / ", length(z_atac),
            " cells overlap between z_atac and z_rna for ", tr)
  }
  z_atac_c <- z_atac[common]
  z_rna_c  <- z_rna[common]

  # 2c. paired seed: primary intersection
  q_atac <- stats::quantile(z_atac_c, 1 - SEED_TOP_FRAC, na.rm = TRUE)
  q_rna  <- stats::quantile(z_rna_c,  1 - SEED_TOP_FRAC, na.rm = TRUE)
  hi_atac <- common[!is.na(z_atac_c) & z_atac_c >= q_atac]
  hi_rna  <- common[!is.na(z_rna_c)  & z_rna_c  >= q_rna]
  seed_inter <- intersect(hi_atac, hi_rna)

  fallback <- length(seed_inter) < INTERSECT_FLOOR
  if (fallback) {
    message(sprintf(
      "  intersection only %d (< %d); falling back to rank-sum top 5%%",
      length(seed_inter), INTERSECT_FLOOR
    ))
    n_c <- length(common)
    r_atac <- rank(z_atac_c, ties.method = "average") / n_c
    r_rna  <- rank(z_rna_c,  ties.method = "average") / n_c
    combined <- r_atac + r_rna
    n_seed <- ceiling(SEED_TOP_FRAC * n_c)
    ord <- order(combined, decreasing = TRUE)
    seed_cells <- common[ord[seq_len(n_seed)]]
  } else {
    seed_cells <- seed_inter
  }
  saveRDS(seed_cells, file.path(out_dir, paste0("seed_paired_", tr, ".rds")))
  message(sprintf(
    "  hi_atac=%d, hi_rna=%d, intersect=%d, fallback=%s, seed_used=%d",
    length(hi_atac), length(hi_rna), length(seed_inter),
    fallback, length(seed_cells)
  ))

  # 2d. load weighted KNN, subset to non-isolated, intersect with seeds
  g <- readRDS(file.path(out_dir, paste0("knn_graph_", tr, ".rds")))
  deg <- Matrix::rowSums(g)
  graph_cells <- rownames(g)[deg > 0]
  g_sub <- g[graph_cells, graph_cells]
  seed_in_graph <- intersect(seed_cells, graph_cells)
  seed_dropped  <- setdiff(seed_cells, graph_cells)
  message(sprintf(
    "  graph cells=%d (isolated=%d) | seeds in graph=%d (dropped=%d)",
    length(graph_cells), length(deg) - length(graph_cells),
    length(seed_in_graph), length(seed_dropped)
  ))
  if (length(seed_in_graph) == 0) {
    warning("No seeds left after graph filter; skip trait"); next
  }

  # 2e. scale factor uses z_atac on full graph cells (matching baseline)
  z_atac_g <- z_atac[graph_cells]
  scale_factor <- SCAVENGE::cal_scalefactor(z_score = z_atac_g,
                                            percent_cut = SCALE_PCT)

  # 2f. RWR on weighted graph
  np_score <- SCAVENGE::randomWalk_sparse(
    intM = g_sub, queryCells = seed_in_graph, gamma = GAMMA
  )

  # 2g. TRS post-processing (matches baseline 03_run_baseline_traits.R)
  omit_idx <- np_score == 0
  np_use   <- np_score[!omit_idx]
  keep     <- names(np_use)
  TRS <- SCAVENGE::capOutlierQuantile(x = np_use, q_ceiling = TRS_Q_CEILING) |>
    SCAVENGE::max_min_scale()
  TRS <- TRS * scale_factor

  # 2h. long table out
  out_df <- data.frame(
    cell    = keep,
    z_atac  = z_atac[keep],
    z_rna   = z_rna[keep],
    seed    = keep %in% seed_in_graph,
    np_score= as.numeric(np_use),
    TRS     = as.numeric(TRS),
    stringsAsFactors = FALSE
  )
  out_file <- file.path(out_dir,
                        paste0("TRS_dual_", make.names(tr), ".tsv.gz"))
  data.table::fwrite(out_df, out_file, sep = "\t")
  message("  wrote ", out_file, " (rows=", nrow(out_df), ")")

  trs_long_list[[tr]] <- out_df
  trs_wide_list[[tr]] <- data.frame(cell = keep, x = as.numeric(TRS))
  colnames(trs_wide_list[[tr]])[2] <- tr

  sum_rows[[tr]] <- data.table::data.table(
    trait            = tr,
    n_cells_total    = length(z_atac),
    n_cells_common   = length(common),
    hi_atac          = length(hi_atac),
    hi_rna           = length(hi_rna),
    intersect_size   = length(seed_inter),
    fallback_used    = fallback,
    n_seed_used      = length(seed_cells),
    n_seed_in_graph  = length(seed_in_graph),
    n_graph_cells    = length(graph_cells),
    n_TRS_cells      = nrow(out_df),
    scale_factor     = scale_factor,
    TRS_q95          = as.numeric(stats::quantile(out_df$TRS, 0.95)),
    TRS_q99          = as.numeric(stats::quantile(out_df$TRS, 0.99)),
    TRS_max          = max(out_df$TRS)
  )
}

# --- 3. wide TRS matrix (cell × trait) -----------------------------------

trs_wide_list <- trs_wide_list[!vapply(trs_wide_list, is.null, logical(1))]
if (length(trs_wide_list) > 0) {
  trs_mat <- Reduce(function(x, y) dplyr::full_join(x, y, by = "cell"),
                    trs_wide_list)
  saveRDS(trs_mat, file.path(out_dir, "TRS_dual_matrix.rds"))
  data.table::fwrite(trs_mat,
                     file.path(out_dir, "TRS_dual_matrix.tsv.gz"), sep = "\t")
}

summary_dt <- data.table::rbindlist(sum_rows)
data.table::fwrite(summary_dt, file.path(out_dir, "paired_seed_summary.csv"))
print(summary_dt)

message("\nDual-modality TRS built. Outputs in: ", out_dir)
