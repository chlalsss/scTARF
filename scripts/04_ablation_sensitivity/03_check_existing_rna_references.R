# Label-free RNA reference suitability scoring.
# -----------------------------------------------------------------------------
# This script ranks candidate RNA references before the downstream sparse
# ATAC-RNA fusion workflow.
#
# It uses two checks only:
#   1. trait gene coverage in RNA and ATAC gene-activity feature spaces
#   2. intrinsic RNA scDRS signal strength from existing score files
#
# It does not use annotation labels, mapping output, pseudo pairs, KNN, RWR, or
# AUROC. The output is a trait-reference suitability score. For each trait, the
# reference with the highest finite score is selected.
#
# Run from the repository root with optional SCTARF_* path overrides.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

source(file.path("scripts", "_paths.R"))
ROOT <- sctarf_nonpaired_root()
OUT_DIR <- file.path(ROOT, "01_data_prep_rna_reference_screening/output")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

TRAITS <- data.table(
  trait = c("MONO", "PLT", "LYM", "MCH", "RDW", "RBC", "EOS", "WBC"),
  trait_id = c(
    "UKB_460K.blood_MONOCYTE_COUNT",
    "UKB_460K.blood_PLATELET_COUNT",
    "UKB_460K.blood_LYMPHOCYTE_COUNT",
    "UKB_460K.blood_MEAN_CORPUSCULAR_HEMOGLOBIN",
    "UKB_460K.blood_RBC_DISTRIB_WIDTH",
    "UKB_460K.blood_RED_COUNT",
    "UKB_460K.blood_EOSINOPHIL_COUNT",
    "UKB_460K.blood_WHITE_COUNT"
  )
)

GENE_SET_FILES <- c(
  file.path(ROOT, "01_data_prep/03_rna_scdrs_python/granja_3traits.gs"),
  file.path(ROOT, "01_data_prep_external_rna/03_cite_scdrs_python/cite_3traits.gs"),
  file.path(ROOT, "01_data_prep_external_rna/03_cite_scdrs_python/cite_trait_screening.gs")
)

ATAC_GENE_TXT <- file.path(ROOT, "01_data_prep_external_rna/output/granja_ga_genes.txt")
ATAC_GENE_ACTIVITY_RDS <- file.path(sctarf_blood_root(), "multiomics_output_v2", "gene_activity_matrix.rds")

RNA_REFERENCES <- list(
  granja_rna = list(
    label = "Granja healthy hematopoiesis RNA",
    gene_source = "rds_feature_table",
    gene_file = file.path(ROOT, "01_data_prep/output/granja_rna_features.rds"),
    score_dirs = c(
      file.path(ROOT, "01_data_prep/03_rna_scdrs_python/output"),
      file.path(ROOT, "01_data_prep/03_rna_scdrs_python/output_missing_traits")
    )
  ),
  cite_rna = list(
    label = "NeurIPS 2021 BMMC CITE GEX",
    gene_source = "h5ad_var",
    gene_file = file.path(
      ROOT,
      "01_data_prep_external_rna/03_cite_scdrs_python/input/cite_rna_raw_counts_gex_only.h5ad"
    ),
    score_dirs = c(
      file.path(ROOT, "01_data_prep_external_rna/03_cite_scdrs_python/output"),
      file.path(ROOT, "01_data_prep_external_rna/03_cite_scdrs_python/output_trait_screening")
    )
  )
)

SUMMARY_FILE <- file.path(OUT_DIR, "rna_suitability_summary.csv")
METRIC_FILE <- file.path(OUT_DIR, "rna_suitability_metrics_long.csv")
RANKED_FILE <- file.path(OUT_DIR, "rna_suitability_ranked_references.csv")
COMPARISON_FILE <- file.path(OUT_DIR, "rna_suitability_trait_reference_comparison.csv")
BEST_FILE <- file.path(OUT_DIR, "rna_suitability_best_reference.csv")
REPORT_FILE <- file.path(OUT_DIR, "rna_suitability_decision_report.md")

clamp01 <- function(x) {
  pmin(1, pmax(0, x))
}

gene_key <- function(x) {
  toupper(trimws(as.character(x)))
}

stop_missing <- function(paths, label) {
  paths <- paths[!is.na(paths)]
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0L) {
    stop(label, " missing file(s):\n  ", paste(missing, collapse = "\n  "),
         call. = FALSE)
  }
}

read_text_genes <- function(path) {
  stop_missing(path, "gene list")
  x <- fread(path, header = FALSE, col.names = "gene")$gene
  unique(gene_key(x[nzchar(x)]))
}

read_rds_feature_genes <- function(path) {
  stop_missing(path, "RNA feature table")
  x <- readRDS(path)
  if (is.data.frame(x)) {
    candidates <- c("gene", "gene_name", "symbol", "feature", "features")
    hit <- intersect(candidates, names(x))[1]
    if (!is.na(hit)) return(unique(gene_key(x[[hit]])))
  }
  rn <- rownames(x)
  if (!is.null(rn)) return(unique(gene_key(rn)))
  stop("Could not find genes in RDS feature file: ", path, call. = FALSE)
}

read_rds_row_genes <- function(path) {
  stop_missing(path, "ATAC gene activity matrix")
  x <- readRDS(path)
  rn <- rownames(x)
  if (is.null(rn)) {
    stop("ATAC gene activity object has no rownames: ", path, call. = FALSE)
  }
  unique(gene_key(rn))
}

read_h5ad_var_genes <- function(path) {
  stop_missing(path, "RNA h5ad")
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("Package hdf5r is required to read h5ad genes: ", path, call. = FALSE)
  }
  h5 <- hdf5r::H5File$new(path, mode = "r")
  on.exit(h5$close_all(), add = TRUE)
  if (!"var" %in% names(h5)) {
    stop("h5ad file has no var group: ", path, call. = FALSE)
  }
  var <- h5[["var"]]
  candidates <- c("gene_symbol", "_index", "index", "gene", "gene_name", "gene_id")
  for (nm in candidates) {
    if (nm %in% names(var)) {
      vals <- as.character(var[[nm]]$read())
      vals <- vals[!is.na(vals) & nzchar(vals)]
      if (length(vals) > 0L) return(unique(gene_key(vals)))
    }
  }
  stop("Could not find a usable gene field in h5ad var group: ", path,
       call. = FALSE)
}

read_reference_genes <- function(cfg) {
  if (identical(cfg$gene_source, "rds_feature_table")) {
    return(read_rds_feature_genes(cfg$gene_file))
  }
  if (identical(cfg$gene_source, "h5ad_var")) {
    return(read_h5ad_var_genes(cfg$gene_file))
  }
  stop("Unknown gene_source: ", cfg$gene_source, call. = FALSE)
}

parse_geneset_field <- function(trait_id, geneset) {
  parts <- strsplit(geneset, ",", fixed = TRUE)[[1]]
  parts <- parts[nzchar(parts)]
  genes <- sub(":.*$", "", parts)
  has_weight <- grepl(":", parts, fixed = TRUE)
  weights <- rep(1, length(parts))
  if (any(has_weight)) {
    parsed <- suppressWarnings(as.numeric(sub("^.*:", "", parts[has_weight])))
    parsed[!is.finite(parsed)] <- 1
    weights[has_weight] <- parsed
  }
  data.table(
    trait_id = trait_id,
    gene = genes,
    gene_key = gene_key(genes),
    weight_abs = abs(weights)
  )
}

read_gene_sets <- function(paths) {
  stop_missing(paths, "gene-set")
  rows <- list()
  idx <- 1L
  for (path in paths) {
    gs <- fread(path, sep = "\t")
    if (!all(c("TRAIT", "GENESET") %in% names(gs))) {
      stop("Gene-set file must contain TRAIT and GENESET columns: ", path,
           call. = FALSE)
    }
    gs <- gs[TRAIT %in% TRAITS$trait_id]
    if (nrow(gs) == 0L) next
    for (i in seq_len(nrow(gs))) {
      rows[[idx]] <- parse_geneset_field(gs$TRAIT[i], gs$GENESET[i])
      idx <- idx + 1L
    }
  }
  out <- rbindlist(rows, use.names = TRUE, fill = TRUE)
  out <- out[nzchar(gene_key)]
  out[, weight_abs := fifelse(is.finite(weight_abs) & weight_abs > 0, weight_abs, 1)]
  out <- out[, .(
    gene = gene[which.max(weight_abs)],
    weight_abs = max(weight_abs)
  ), by = .(trait_id, gene_key)]
  merge(out, TRAITS, by = "trait_id", all.x = TRUE, sort = FALSE)
}

score_gene_coverage <- function(gs_trait, rna_genes, atac_genes) {
  genes <- unique(gs_trait, by = "gene_key")
  genes[, in_rna := gene_key %in% rna_genes]
  genes[, in_atac := gene_key %in% atac_genes]
  genes[, in_shared := in_rna & in_atac]
  total_weight <- sum(genes$weight_abs)
  if (!is.finite(total_weight) || total_weight <= 0) total_weight <- nrow(genes)

  n_trait <- nrow(genes)
  n_in_rna <- sum(genes$in_rna)
  n_in_atac <- sum(genes$in_atac)
  n_shared <- sum(genes$in_shared)
  rna_weighted <- sum(genes[in_rna == TRUE, weight_abs]) / total_weight
  atac_weighted <- sum(genes[in_atac == TRUE, weight_abs]) / total_weight
  shared_weighted <- sum(genes[in_shared == TRUE, weight_abs]) / total_weight
  unweighted_shared <- n_shared / max(n_trait, 1L)
  score <- 100 * (
    0.45 * rna_weighted +
      0.45 * shared_weighted +
      0.10 * atac_weighted
  )

  data.table(
    n_trait_genes = n_trait,
    n_trait_genes_in_rna = n_in_rna,
    n_trait_genes_in_atac = n_in_atac,
    n_trait_genes_shared = n_shared,
    rna_trait_gene_coverage = n_in_rna / max(n_trait, 1L),
    atac_trait_gene_coverage = n_in_atac / max(n_trait, 1L),
    shared_trait_gene_coverage = unweighted_shared,
    rna_weighted_gene_coverage = rna_weighted,
    atac_weighted_gene_coverage = atac_weighted,
    shared_weighted_gene_coverage = shared_weighted,
    gene_coverage_score = round(score, 2)
  )
}

find_score_file <- function(score_dirs, trait_id) {
  candidates <- file.path(score_dirs, paste0(trait_id, ".score.gz"))
  hit <- candidates[file.exists(candidates)][1]
  if (length(hit) == 0L || is.na(hit)) return(NA_character_)
  hit
}

read_trait_score <- function(cfg, trait_id) {
  path <- find_score_file(cfg$score_dirs, trait_id)
  if (is.na(path)) return(NULL)
  dt <- fread(path)
  cell_col <- intersect(c("cell", "cell_id", "barcode", "V1"), names(dt))[1]
  score_col <- intersect(c("norm_score", "zscore", "score"), names(dt))[1]
  p_col <- intersect(c("pval", "mc_pval", "pval_ctrl", "pval_mc"), names(dt))[1]
  if (is.na(cell_col) || is.na(score_col)) {
    stop("Score file lacks required columns: ", path, call. = FALSE)
  }
  out <- data.table(
    cell = as.character(dt[[cell_col]]),
    norm_score = as.numeric(dt[[score_col]])
  )
  if (!is.na(p_col)) {
    out[, pval := as.numeric(dt[[p_col]])]
  } else {
    out[, pval := NA_real_]
  }
  out[!is.na(cell) & nzchar(cell)]
}

score_rna_signal <- function(score_dt) {
  if (is.null(score_dt)) {
    return(data.table(
      score_file_found = FALSE,
      n_scored_cells = 0L,
      score_mean = NA_real_,
      score_sd = NA_real_,
      score_iqr = NA_real_,
      score_q95 = NA_real_,
      score_q99 = NA_real_,
      top1pct_mean_score = NA_real_,
      background_mean_score = NA_real_,
      top1pct_delta_score = NA_real_,
      pval_available = FALSE,
      frac_pval_le_0p01 = NA_real_,
      frac_pval_le_0p05 = NA_real_,
      pval_0p01_enrichment = NA_real_,
      n_qval_le_0p10_positive = NA_integer_,
      qval_le_0p10_positive_fraction = NA_real_,
      rna_signal_score = NA_real_
    ))
  }

  dt <- score_dt[is.finite(norm_score)]
  n_scored <- nrow(dt)
  if (n_scored == 0L) {
    return(data.table(
      score_file_found = TRUE,
      n_scored_cells = 0L,
      score_mean = NA_real_,
      score_sd = NA_real_,
      score_iqr = NA_real_,
      score_q95 = NA_real_,
      score_q99 = NA_real_,
      top1pct_mean_score = NA_real_,
      background_mean_score = NA_real_,
      top1pct_delta_score = NA_real_,
      pval_available = FALSE,
      frac_pval_le_0p01 = NA_real_,
      frac_pval_le_0p05 = NA_real_,
      pval_0p01_enrichment = NA_real_,
      n_qval_le_0p10_positive = NA_integer_,
      qval_le_0p10_positive_fraction = NA_real_,
      rna_signal_score = NA_real_
    ))
  }

  ord <- order(dt$norm_score, decreasing = TRUE)
  top_n <- max(1L, ceiling(0.01 * n_scored))
  top_idx <- ord[seq_len(top_n)]
  top_mean <- mean(dt$norm_score[top_idx])
  bg_mean <- mean(dt$norm_score[-top_idx])
  delta <- top_mean - bg_mean

  has_pval <- any(is.finite(dt$pval))
  frac_p01 <- NA_real_
  frac_p05 <- NA_real_
  p01_enrichment <- NA_real_
  n_q10_pos <- NA_integer_
  q10_pos_frac <- NA_real_
  if (has_pval) {
    p <- dt$pval
    frac_p01 <- mean(p <= 0.01, na.rm = TRUE)
    frac_p05 <- mean(p <= 0.05, na.rm = TRUE)
    p01_enrichment <- frac_p01 / 0.01
    q <- p.adjust(p, method = "BH")
    n_q10_pos <- sum(q <= 0.10 & dt$norm_score > 0, na.rm = TRUE)
    q10_pos_frac <- n_q10_pos / n_scored
  }

  score_sd <- stats::sd(dt$norm_score)
  p_component <- if (has_pval) clamp01((p01_enrichment - 1) / 2) else 0
  q_component <- if (has_pval) {
    clamp01(n_q10_pos / max(20, ceiling(0.001 * n_scored)))
  } else {
    0
  }
  sd_component <- clamp01((score_sd - 1.0) / 1.5)
  signal_score <- if (has_pval) {
    100 * (0.50 * p_component + 0.30 * q_component + 0.20 * sd_component)
  } else {
    100 * (0.65 * sd_component + 0.35 * clamp01((delta - 2.5) / 2.5))
  }

  data.table(
    score_file_found = TRUE,
    n_scored_cells = n_scored,
    score_mean = mean(dt$norm_score),
    score_sd = score_sd,
    score_iqr = stats::IQR(dt$norm_score),
    score_q95 = as.numeric(stats::quantile(dt$norm_score, 0.95)),
    score_q99 = as.numeric(stats::quantile(dt$norm_score, 0.99)),
    top1pct_mean_score = top_mean,
    background_mean_score = bg_mean,
    top1pct_delta_score = delta,
    pval_available = has_pval,
    frac_pval_le_0p01 = frac_p01,
    frac_pval_le_0p05 = frac_p05,
    pval_0p01_enrichment = p01_enrichment,
    n_qval_le_0p10_positive = n_q10_pos,
    qval_le_0p10_positive_fraction = q10_pos_frac,
    rna_signal_score = round(signal_score, 2)
  )
}

suitability_score <- function(coverage_score, signal_score) {
  if (!is.finite(signal_score)) return(NA_real_)
  round(0.45 * coverage_score + 0.55 * signal_score, 2)
}

selection_note_one <- function(selected, has_score, score,
                               gap_to_best, margin_next,
                               coverage_score, signal_score) {
  score_txt <- if (is.finite(score)) sprintf("%.2f", score) else "NA"
  cov_txt <- if (is.finite(coverage_score)) sprintf("%.2f", coverage_score) else "NA"
  sig_txt <- if (is.finite(signal_score)) sprintf("%.2f", signal_score) else "NA"

  if (isTRUE(selected)) {
    if (is.finite(margin_next)) {
      return(sprintf(
        "Selected: highest RNA suitability score (%s); margin over next reference = %.2f. Coverage score = %s; RNA signal score = %s.",
        score_txt, margin_next, cov_txt, sig_txt
      ))
    }
    return(sprintf(
      "Selected: highest RNA suitability score (%s). Coverage score = %s; RNA signal score = %s.",
      score_txt, cov_txt, sig_txt
    ))
  }

  if (!isTRUE(has_score)) {
    return(sprintf(
      "Not selected: no comparable RNA suitability score. Coverage score = %s; RNA signal score = %s.",
      cov_txt, sig_txt
    ))
  }

  gap_txt <- if (is.finite(gap_to_best)) sprintf("%.2f", gap_to_best) else "NA"
  sprintf(
    "Not selected: lower than the selected reference by %s points. Score = %s; coverage score = %s; RNA signal score = %s.",
    gap_txt, score_txt, cov_txt, sig_txt
  )
}

add_selection_ranking <- function(summary_dt) {
  ranked <- copy(summary_dt)
  ranked[, has_suitability_score := is.finite(rna_suitability_score)]
  ranked[, score_for_order := fifelse(is.finite(rna_suitability_score),
                                      rna_suitability_score, -Inf)]

  setorder(ranked, trait, -has_suitability_score, -score_for_order, reference)
  ranked[, suitability_rank_for_trait := seq_len(.N), by = trait]
  ranked[, best_score := {
    vals <- rna_suitability_score[has_suitability_score == TRUE]
    if (length(vals) == 0L) NA_real_ else max(vals)
  }, by = trait]
  ranked[, next_score := {
    vals <- sort(rna_suitability_score[has_suitability_score == TRUE],
                 decreasing = TRUE)
    if (length(vals) < 2L) NA_real_ else vals[2]
  }, by = trait]
  ranked[, score_gap_to_best := fifelse(
    is.finite(best_score) & is.finite(rna_suitability_score),
    round(best_score - rna_suitability_score, 2),
    NA_real_
  )]
  ranked[, selected_by_score := has_suitability_score == TRUE &
           suitability_rank_for_trait == 1L]
  ranked[, score_margin_over_next := fifelse(
    selected_by_score == TRUE & is.finite(next_score),
    round(rna_suitability_score - next_score, 2),
    NA_real_
  )]
  ranked[, selection_note := mapply(
    selection_note_one,
    selected = selected_by_score,
    has_score = has_suitability_score,
    score = rna_suitability_score,
    gap_to_best = score_gap_to_best,
    margin_next = score_margin_over_next,
    coverage_score = gene_coverage_score,
    signal_score = rna_signal_score,
    USE.NAMES = FALSE
  )]
  ranked[, score_for_order := NULL]
  ranked[]
}

make_ref_wide <- function(ranked_dt, ref_name) {
  cols <- c(
    "rna_suitability_score",
    "gene_coverage_score",
    "rna_signal_score",
    "suitability_rank_for_trait",
    "selected_by_score"
  )
  out <- ranked_dt[reference == ref_name, c("trait", cols), with = FALSE]
  setnames(out, cols, paste0(ref_name, "_", cols))
  out
}

selection_summary_one <- function(recommended_reference, selected_score,
                                  granja_score, cite_score) {
  fmt <- function(x) {
    if (is.finite(x)) sprintf("%.2f", x) else "NA"
  }

  if (identical(recommended_reference, "none")) {
    return(sprintf(
      "No RNA reference selected because no comparable suitability score is available. Granja score=%s; CITE score=%s.",
      fmt(granja_score), fmt(cite_score)
    ))
  }

  if (identical(recommended_reference, "granja_rna")) {
    if (is.finite(cite_score)) {
      return(sprintf(
        "Selected granja_rna: score %.2f vs cite_rna %.2f; margin %.2f.",
        selected_score, cite_score, selected_score - cite_score
      ))
    }
    return(sprintf(
      "Selected granja_rna: score %.2f. CITE has no comparable suitability score.",
      selected_score
    ))
  }

  if (identical(recommended_reference, "cite_rna")) {
    if (is.finite(granja_score)) {
      return(sprintf(
        "Selected cite_rna: score %.2f vs granja_rna %.2f; margin %.2f.",
        selected_score, granja_score, selected_score - granja_score
      ))
    }
    return(sprintf(
      "Selected cite_rna: score %.2f. Granja has no comparable suitability score.",
      selected_score
    ))
  }

  sprintf("Selected %s with suitability score %s.", recommended_reference,
          fmt(selected_score))
}

make_trait_reference_comparison <- function(ranked_dt, best_dt) {
  ref_tables <- lapply(names(RNA_REFERENCES), function(ref_name) {
    make_ref_wide(ranked_dt, ref_name)
  })
  wide <- Reduce(function(x, y) merge(x, y, by = "trait", all = TRUE), ref_tables)

  selected <- best_dt[, .(
    trait,
    recommended_reference,
    selected_rna_suitability_score = fifelse(recommended_reference == "none",
                                             NA_real_, rna_suitability_score),
    score_margin_over_next
  )]
  out <- merge(selected, wide, by = "trait", all.x = TRUE, sort = FALSE)

  if (all(c("granja_rna_rna_suitability_score",
            "cite_rna_rna_suitability_score") %in% names(out))) {
    out[, granja_minus_cite_suitability_score := fifelse(
      is.finite(granja_rna_rna_suitability_score) &
        is.finite(cite_rna_rna_suitability_score),
      round(granja_rna_rna_suitability_score - cite_rna_rna_suitability_score, 2),
      NA_real_
    )]
    out[, selection_summary := mapply(
      selection_summary_one,
      recommended_reference = recommended_reference,
      selected_score = selected_rna_suitability_score,
      granja_score = granja_rna_rna_suitability_score,
      cite_score = cite_rna_rna_suitability_score,
      USE.NAMES = FALSE
    )]
  }

  front <- c(
    "trait",
    "recommended_reference",
    "selected_rna_suitability_score",
    "score_margin_over_next",
    "granja_minus_cite_suitability_score",
    "granja_rna_rna_suitability_score",
    "granja_rna_gene_coverage_score",
    "granja_rna_rna_signal_score",
    "cite_rna_rna_suitability_score",
    "cite_rna_gene_coverage_score",
    "cite_rna_rna_signal_score",
    "selection_summary"
  )
  setcolorder(out, c(intersect(front, names(out)), setdiff(names(out), front)))
  out[order(trait)]
}

metric_long_table <- function(summary_dt) {
  metric_cols <- c(
    "gene_coverage_score",
    "rna_signal_score",
    "rna_suitability_score",
    "shared_weighted_gene_coverage",
    "pval_0p01_enrichment",
    "frac_pval_le_0p05",
    "qval_le_0p10_positive_fraction"
  )
  melt(
    summary_dt,
    id.vars = c("reference", "reference_label", "trait", "trait_id"),
    measure.vars = metric_cols,
    variable.name = "metric",
    value.name = "value"
  )
}

md_table <- function(dt, cols) {
  if (nrow(dt) == 0L) return("")
  x <- copy(dt[, ..cols])
  for (nm in names(x)) x[[nm]] <- as.character(x[[nm]])
  header <- paste0("| ", paste(cols, collapse = " | "), " |")
  sep <- paste0("| ", paste(rep("---", length(cols)), collapse = " | "), " |")
  body <- apply(x, 1, function(row) paste0("| ", paste(row, collapse = " | "), " |"))
  paste(c(header, sep, body), collapse = "\n")
}

write_report <- function(summary_dt, ranked_dt, comparison_dt, best_dt) {
  report_summary <- copy(summary_dt)
  report_ranked <- copy(ranked_dt)
  report_comparison <- copy(comparison_dt)
  num_cols <- intersect(c(
    "rna_suitability_score",
    "gene_coverage_score",
    "rna_signal_score",
    "shared_weighted_gene_coverage",
    "pval_0p01_enrichment",
    "frac_pval_le_0p05"
  ), names(report_summary))
  for (nm in num_cols) report_summary[[nm]] <- round(report_summary[[nm]], 3)
  for (nm in intersect(num_cols, names(report_ranked))) {
    report_ranked[[nm]] <- round(report_ranked[[nm]], 3)
  }
  for (nm in names(report_comparison)) {
    if (is.numeric(report_comparison[[nm]])) {
      report_comparison[[nm]] <- round(report_comparison[[nm]], 3)
    }
  }

  best_report <- copy(best_dt)
  for (nm in intersect(num_cols, names(best_report))) {
    best_report[[nm]] <- round(best_report[[nm]], 3)
  }

  lines <- c(
    "# Label-Free RNA Suitability Score Selection",
    "",
    paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "",
    "## Method",
    "",
    "This is a score-based RNA reference selection tool. It uses no annotation labels and no mapping output.",
    "",
    "Each trait-reference pair receives two component scores.",
    "",
    "The gene coverage score measures how well the trait gene set is represented in the RNA reference and the ATAC gene-activity feature space.",
    "",
    "The RNA signal score measures whether existing RNA scDRS scores show a non-uniform positive tail using cell-level p-value enrichment, BH-adjusted positive counts, and score dispersion.",
    "",
    "Scores:",
    "",
    "- gene_coverage_score = 100 * (0.45 * RNA weighted gene coverage + 0.45 * shared RNA-ATAC weighted gene coverage + 0.10 * ATAC weighted gene coverage).",
    "- RNA signal p-value component = clamp((p<=0.01 enrichment - 1) / 2, 0, 1).",
    "- RNA signal q-value component = clamp(number of positive q<=0.10 cells / max(20, 0.1% of scored cells), 0, 1).",
    "- RNA signal dispersion component = clamp((score SD - 1.0) / 1.5, 0, 1).",
    "- rna_signal_score = 100 * (0.50 * p-value component + 0.30 * q-value component + 0.20 * dispersion component).",
    "- rna_suitability_score = 0.45 * gene_coverage_score + 0.55 * rna_signal_score.",
    "",
    "Selection rule: for each trait, select the candidate RNA reference with the highest finite rna_suitability_score. Candidates without an RNA scDRS score are kept in the tables but are not comparable for selection.",
    "",
    "## Recommended Reference Per Trait",
    "",
    md_table(
      best_report,
      c("trait", "recommended_reference", "selected_rna_suitability_score",
        "score_margin_over_next", "gene_coverage_score", "rna_signal_score")
    ),
    "",
    "## Trait-Reference Suitability Comparison",
    "",
    md_table(
      report_comparison,
      c("trait", "recommended_reference", "selected_rna_suitability_score",
        "score_margin_over_next", "granja_rna_rna_suitability_score",
        "cite_rna_rna_suitability_score", "granja_minus_cite_suitability_score",
        "selection_summary")
    ),
    "",
    "## Ranked Reference Comparison",
    "",
    md_table(
      report_ranked[order(trait, suitability_rank_for_trait)],
      c("trait", "suitability_rank_for_trait", "reference", "selected_by_score",
        "rna_suitability_score", "score_gap_to_best",
        "gene_coverage_score", "rna_signal_score", "selection_note")
    ),
    "",
    "## Full Summary",
    "",
    md_table(
      report_summary[order(trait, reference)],
      c("trait", "reference", "rna_suitability_score",
        "gene_coverage_score", "rna_signal_score",
        "shared_weighted_gene_coverage", "pval_0p01_enrichment",
        "frac_pval_le_0p05")
    ),
    "",
    "## Output Files",
    "",
    paste0("- ", SUMMARY_FILE),
    paste0("- ", METRIC_FILE),
    paste0("- ", RANKED_FILE),
    paste0("- ", COMPARISON_FILE),
    paste0("- ", BEST_FILE)
  )
  writeLines(lines, REPORT_FILE, useBytes = TRUE)
}

cat(">>> Label-free RNA suitability score selection\n\n")
cat("Output directory: ", OUT_DIR, "\n\n", sep = "")

gene_sets <- read_gene_sets(GENE_SET_FILES)
if (file.exists(ATAC_GENE_TXT)) {
  atac_genes <- read_text_genes(ATAC_GENE_TXT)
} else {
  atac_genes <- read_rds_row_genes(ATAC_GENE_ACTIVITY_RDS)
}

summary_rows <- list()
row_idx <- 1L

for (ref_name in names(RNA_REFERENCES)) {
  cfg <- RNA_REFERENCES[[ref_name]]
  cat("Reference: ", ref_name, " - ", cfg$label, "\n", sep = "")
  rna_genes <- read_reference_genes(cfg)
  cat("  RNA genes: ", length(rna_genes), "\n", sep = "")

  for (i in seq_len(nrow(TRAITS))) {
    trait_short <- TRAITS$trait[i]
    trait_id <- TRAITS$trait_id[i]
    gs_trait <- gene_sets[trait == trait_short]
    if (nrow(gs_trait) == 0L) {
      warning("No gene set found for trait: ", trait_short)
      next
    }

    cov_dt <- score_gene_coverage(gs_trait, rna_genes, atac_genes)
    score_dt <- read_trait_score(cfg, trait_id)
    sig_dt <- score_rna_signal(score_dt)
    total_score <- suitability_score(cov_dt$gene_coverage_score, sig_dt$rna_signal_score)

    summary_rows[[row_idx]] <- cbind(
      data.table(
        reference = ref_name,
        reference_label = cfg$label,
        trait = trait_short,
        trait_id = trait_id,
        rna_suitability_score = total_score
      ),
      cov_dt,
      sig_dt
    )
    row_idx <- row_idx + 1L

    cat("  ", trait_short,
        ": score=", ifelse(is.na(total_score), "NA", total_score),
        " | coverage_score=", cov_dt$gene_coverage_score,
        " | rna_signal_score=", ifelse(is.na(sig_dt$rna_signal_score),
                                       "NA", sig_dt$rna_signal_score),
        "\n", sep = "")
  }
  cat("\n")
}

summary_dt <- rbindlist(summary_rows, use.names = TRUE, fill = TRUE)
ranked_dt <- add_selection_ranking(summary_dt)
setorder(summary_dt, trait, -rna_suitability_score, reference)

best_dt <- ranked_dt[suitability_rank_for_trait == 1L]
best_dt[, recommended_reference := fifelse(
  selected_by_score == TRUE,
  reference,
  "none"
)]
best_dt[, selected_rna_suitability_score := fifelse(
  recommended_reference == "none",
  NA_real_,
  rna_suitability_score
)]
setcolorder(best_dt, c(
  "trait", "recommended_reference", "selected_rna_suitability_score",
  "score_margin_over_next", "reference",
  setdiff(names(best_dt), c("trait", "recommended_reference",
                            "selected_rna_suitability_score",
                            "score_margin_over_next",
                            "reference"))
))
setorder(best_dt, trait)

comparison_dt <- make_trait_reference_comparison(ranked_dt, best_dt)
metrics_long <- metric_long_table(summary_dt)

fwrite(summary_dt, SUMMARY_FILE)
fwrite(metrics_long, METRIC_FILE)
fwrite(ranked_dt, RANKED_FILE)
fwrite(comparison_dt, COMPARISON_FILE)
fwrite(best_dt, BEST_FILE)
write_report(summary_dt, ranked_dt, comparison_dt, best_dt)

cat("Recommended reference per trait by RNA suitability score:\n")
print(best_dt[, .(
  trait,
  recommended_reference,
  selected_rna_suitability_score,
  score_margin_over_next,
  gene_coverage_score,
  rna_signal_score
)])

cat("\nWrote:\n")
cat("  ", SUMMARY_FILE, "\n", sep = "")
cat("  ", METRIC_FILE, "\n", sep = "")
cat("  ", RANKED_FILE, "\n", sep = "")
cat("  ", COMPARISON_FILE, "\n", sep = "")
cat("  ", BEST_FILE, "\n", sep = "")
cat("  ", REPORT_FILE, "\n", sep = "")

cat("\nDone. No downstream fusion or mapping was run.\n")
