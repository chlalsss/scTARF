# 04_dual_modality/test_atacseed_dualgraph.R
#
# Ablation test: ATAC-only seed (z_atac top 5%, baseline-style) on the dual
# weighted KNN graph (§2.5 product). Isolates which side hurts AUROC:
#   baseline                 : z_atac seed  + ATAC-only k=30 binary KNN  → AUROC_base
#   dual paired (current)    : intersect    + dual k=90 weighted KNN     → AUROC_paired
#   dual ATAC seed (this test): z_atac seed  + dual k=90 weighted KNN     → AUROC_test
#
#   AUROC_test vs base : graph contribution alone (positive = dual graph helps)
#   AUROC_test vs paired: seed contribution alone (positive = ATAC seed beats paired)
#
# Inputs:
#   zatac_<trait>.rds      from §2.7 (we already paid the chromVAR cost)
#   knn_graph_<trait>.rds  from §2.5 (the dual weighted KNN)
#   BMMC_cell_meta.rds     for AUROC eval only
#
# Outputs:
#   TRS_atacseed_dualgraph_<trait>.tsv.gz     long table
#   atacseed_dualgraph_summary.csv            per trait AUROC + 3-way compare

source(file.path(dirname(sys.frame(1)$ofile %||% "."), "..", "_bootstrap.R"))
required_pkgs(c("data.table", "Matrix", "SCAVENGE"))

paths   <- m3_paths()
out_dir <- file.path(paths$out_dir, "04_dual_modality")
base_dir <- file.path(paths$out_dir, "03_scavenge_baseline")

TRAITS         <- c("MONO", "PLT", "LYM")
SEED_TOP_FRAC  <- 0.05
GAMMA          <- 0.05
TRS_Q_CEILING  <- 0.95
SCALE_PCT      <- 0.01

auroc_mwu <- function(scores, labels) {
  ok <- !is.na(scores) & !is.na(labels)
  s  <- scores[ok]; l <- as.logical(labels[ok])
  n_pos <- sum(l); n_neg <- sum(!l)
  if (n_pos == 0 || n_neg == 0) return(NA_real_)
  r <- rank(s)
  (sum(r[l]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

# meta only loaded for end-of-script AUROC eval
meta <- readRDS(file.path(paths$out_dir, "01_data_qc", "BMMC_cell_meta.rds"))
meta_dt <- data.table::as.data.table(meta)

# baseline / paired AUROC reference (from prior runs)
auroc_eval_csv <- file.path(out_dir, "evaluation_auroc.csv")
if (file.exists(auroc_eval_csv)) {
  prev <- data.table::fread(auroc_eval_csv)
} else {
  prev <- NULL
}

rows <- list()

for (tr in TRAITS) {
  message("\n[Trait] ", tr)
  z_atac <- readRDS(file.path(out_dir, paste0("zatac_", tr, ".rds")))
  q_atac <- stats::quantile(z_atac, 1 - SEED_TOP_FRAC, na.rm = TRUE)
  seed_cells <- names(z_atac)[!is.na(z_atac) & z_atac >= q_atac]
  message("  z_atac seed: ", length(seed_cells), " cells (top ",
          100 * SEED_TOP_FRAC, "%)")

  g <- readRDS(file.path(out_dir, paste0("knn_graph_", tr, ".rds")))
  deg <- Matrix::rowSums(g)
  graph_cells <- rownames(g)[deg > 0]
  g_sub <- g[graph_cells, graph_cells]
  seed_in_graph <- intersect(seed_cells, graph_cells)
  message("  graph cells=", length(graph_cells),
          " | seeds in graph=", length(seed_in_graph))

  scale_factor <- SCAVENGE::cal_scalefactor(
    z_score = z_atac[graph_cells], percent_cut = SCALE_PCT
  )

  np_score <- SCAVENGE::randomWalk_sparse(
    intM = g_sub, queryCells = seed_in_graph, gamma = GAMMA
  )

  omit_idx <- np_score == 0
  np_use   <- np_score[!omit_idx]
  keep     <- names(np_use)
  TRS <- SCAVENGE::capOutlierQuantile(np_use, q_ceiling = TRS_Q_CEILING) |>
    SCAVENGE::max_min_scale()
  TRS <- TRS * scale_factor

  out_df <- data.frame(
    cell    = keep,
    z_atac  = z_atac[keep],
    TRS     = as.numeric(TRS),
    stringsAsFactors = FALSE
  )
  out_file <- file.path(out_dir,
                        paste0("TRS_atacseed_dualgraph_", tr, ".tsv.gz"))
  data.table::fwrite(out_df, out_file, sep = "\t")

  # AUROC against doc-primary target (and LYM_TNK as bonus)
  ct_lookup <- meta_dt[, .(cell = barcode, cell_type)]
  out_dt <- merge(data.table::as.data.table(out_df), ct_lookup, by = "cell")

  primary_set <- bmmc_target_types()[[tr]]
  is_primary <- out_dt$cell_type %in% primary_set
  auroc_primary <- auroc_mwu(out_dt$TRS, is_primary)

  auroc_secondary <- NA_real_
  sec_label <- ""
  if (tr == "LYM") {
    sec_set <- bmmc_target_types_secondary()$LYM_TNK
    auroc_secondary <- auroc_mwu(out_dt$TRS, out_dt$cell_type %in% sec_set)
    sec_label <- "LYM_TNK"
  }

  rows[[tr]] <- data.table::data.table(
    trait        = tr,
    n_seed       = length(seed_in_graph),
    n_TRS_cells  = nrow(out_df),
    AUROC_test   = auroc_primary,
    AUROC_secondary = auroc_secondary,
    secondary_label = sec_label
  )
  message(sprintf(
    "  AUROC primary (%s) = %.3f%s",
    tr, auroc_primary,
    if (!is.na(auroc_secondary))
      sprintf(" | %s = %.3f", sec_label, auroc_secondary) else ""
  ))
}

summary_dt <- data.table::rbindlist(rows, fill = TRUE)

# attach baseline + paired numbers for 3-way compare
if (!is.null(prev)) {
  base_map <- prev[method == "baseline" & target_label %in% TRAITS,
                   .(trait = target_label, AUROC_base = AUROC)]
  paired_map <- prev[method == "dual" & target_label %in% TRAITS,
                     .(trait = target_label, AUROC_paired = AUROC)]
  summary_dt <- merge(summary_dt, base_map, by = "trait", all.x = TRUE)
  summary_dt <- merge(summary_dt, paired_map, by = "trait", all.x = TRUE)
  summary_dt[, delta_test_vs_base   := AUROC_test - AUROC_base]
  summary_dt[, delta_test_vs_paired := AUROC_test - AUROC_paired]
  data.table::setcolorder(summary_dt,
    c("trait", "n_seed", "n_TRS_cells",
      "AUROC_base", "AUROC_paired", "AUROC_test",
      "delta_test_vs_base", "delta_test_vs_paired",
      "AUROC_secondary", "secondary_label"))
}

data.table::fwrite(summary_dt,
                   file.path(out_dir, "atacseed_dualgraph_summary.csv"))
print(summary_dt)

message("\nDone. Output: ", out_dir)
