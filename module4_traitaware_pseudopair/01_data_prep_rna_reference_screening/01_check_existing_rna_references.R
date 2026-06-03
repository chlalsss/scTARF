# Check existing RNA references for trait suitability.
# -----------------------------------------------------------------------------
# This script does not run scDRS, mapping, pseudo RNA completion, KNN, RWR, or
# AUROC. It only reads existing RNA metadata and existing scDRS scores, then
# reports which trait-reference pairs are already usable or need more work.
#
# Run in RStudio:
#   setwd("D:/RStudio/work/module4_traitaware_pseudopair")
#   source("01_data_prep_rna_reference_screening/01_check_existing_rna_references.R")
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

ROOT <- "D:/RStudio/work/module4_traitaware_pseudopair"
OUT_DIR <- file.path(ROOT, "01_data_prep_rna_reference_screening/output")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

TRAITS_TO_CHECK <- c("MONO", "PLT", "LYM", "MCH", "RDW", "RBC", "EOS", "WBC")

REFERENCE_CONFIG <- list(
  granja_rna = list(
    label = "Granja scRNA Healthy Hematopoiesis",
    meta_type = "rds",
    meta_file = file.path(ROOT, "01_data_prep/output/granja_rna_cell_meta.rds"),
    score_files = c(
      file.path(ROOT, "01_data_prep/output/granja_rna_scdrs_scores.rds")
    ),
    cell_col_candidates = c("cell", "barcode"),
    celltype_col_candidates = c("cell_type", "BioClassification")
  ),
  cite_rna = list(
    label = "NeurIPS 2021 BMMC CITE GEX",
    meta_type = "csv",
    meta_file = file.path(
      ROOT,
      "01_data_prep_external_rna/03_cite_scdrs_python/cite_scdrs_cell_meta.csv.gz"
    ),
    score_files = c(
      file.path(ROOT, "01_data_prep_external_rna/output/cite_rna_scdrs_scores.rds"),
      file.path(
        ROOT,
        "01_data_prep_external_rna/output/trait_screening/cite_rna_scdrs_scores_trait_screening.rds"
      ),
      file.path(
        ROOT,
        "01_data_prep_external_rna/output/trait_screening/cite_rna_scdrs_scores_trait_screening_remaining.rds"
      )
    ),
    cell_col_candidates = c("cell", "barcode"),
    celltype_col_candidates = c("cell_type")
  )
)

stop_missing <- function(paths, label) {
  paths <- paths[!is.na(paths)]
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0L) {
    stop(label, " missing required file(s):\n  ",
         paste(missing, collapse = "\n  "), call. = FALSE)
  }
}

pick_col <- function(nms, candidates, required = TRUE, label = "table") {
  hit <- intersect(candidates, nms)[1]
  if (is.na(hit) && required) {
    stop("Could not find required column in ", label, ". Tried: ",
         paste(candidates, collapse = ", "), call. = FALSE)
  }
  hit
}

as_dt_preserve_rownames <- function(x) {
  if (inherits(x, "data.table")) {
    dt <- copy(x)
  } else {
    df <- as.data.frame(x, stringsAsFactors = FALSE)
    rn <- rownames(df)
    if (!is.null(rn) &&
        !identical(rn, as.character(seq_len(nrow(df)))) &&
        !"cell" %in% colnames(df) &&
        !"barcode" %in% colnames(df)) {
      df$cell <- rn
    }
    dt <- as.data.table(df)
  }
  dt
}

read_reference_meta <- function(cfg, ref_name) {
  stop_missing(cfg$meta_file, paste0(ref_name, " metadata"))

  if (identical(cfg$meta_type, "rds")) {
    dt <- as_dt_preserve_rownames(readRDS(cfg$meta_file))
  } else if (identical(cfg$meta_type, "csv")) {
    dt <- fread(cfg$meta_file)
  } else {
    stop("Unknown metadata type for ", ref_name, ": ", cfg$meta_type,
         call. = FALSE)
  }

  cell_col <- pick_col(names(dt), cfg$cell_col_candidates,
                       label = paste(ref_name, "metadata"))
  celltype_col <- pick_col(names(dt), cfg$celltype_col_candidates,
                           label = paste(ref_name, "metadata"))

  out <- dt[, .(
    cell = as.character(get(cell_col)),
    cell_type = as.character(get(celltype_col))
  )]
  out <- out[!is.na(cell) & cell != ""]
  unique(out, by = "cell")
}

read_one_score_file <- function(path) {
  if (!file.exists(path)) return(NULL)
  dt <- as_dt_preserve_rownames(readRDS(path))
  cell_col <- pick_col(names(dt), c("cell", "barcode"), label = path)
  setnames(dt, cell_col, "cell")
  dt[, cell := as.character(cell)]

  score_cols <- grep("_norm_score$", names(dt), value = TRUE)
  aux_cols <- intersect(c("scdrs_input_total_counts", "scdrs_input_n_genes"),
                        names(dt))
  keep <- unique(c("cell", score_cols, aux_cols))
  unique(dt[, ..keep], by = "cell")
}

merge_reference_scores <- function(paths) {
  score_list <- Filter(Negate(is.null), lapply(paths, read_one_score_file))
  if (length(score_list) == 0L) {
    return(data.table(cell = character()))
  }

  out <- score_list[[1]]
  if (length(score_list) > 1L) {
    for (i in 2:length(score_list)) {
      add <- score_list[[i]]
      add_cols <- setdiff(names(add), names(out))
      out <- merge(out, add[, c("cell", add_cols), with = FALSE],
                   by = "cell", all = TRUE, sort = FALSE)
    }
  }
  out
}

is_expected_for_trait <- function(trait, cell_type) {
  x <- tolower(cell_type)

  is_mono <- grepl("mono|cd14|cd16\\.mono|cd16\\+ mono", x)
  is_b <- grepl("(^|[^a-z])b($|[^a-z])|pre\\.b|plasma|b1 b|naive cd20", x)
  is_t <- grepl("cd4|cd8|mait|t reg|gdt|dnt|t naive|t activated|tigit", x)
  is_nk_ilc <- grepl("(^|[^a-z])nk($|[^a-z])|nk cd|(^|[^a-z])ilc($|[^a-z])|ilc1", x)
  is_lymph <- is_b | is_t | is_nk_ilc | grepl("lymph", x)
  is_ery <- grepl(
    "eryth|reticulocyte|normoblast|proeryth|early\\.eryth|late\\.eryth|baso|mk/e prog",
    x
  )
  is_mk <- grepl("megak|platelet|mk/e prog|\\bmk\\b", x)
  is_gran <- grepl("eosin|neut|gmp\\.neut|g/m prog|(^|[^a-z])gmp($|[^a-z])|granul|baso|early\\.baso", x)
  is_dc <- grepl("\\bdc\\b|cdc|pdc", x)
  is_prog <- grepl("hsc|cmp|lmpp|(^|[^a-z])gmp($|[^a-z])|clp|prog", x)

  if (identical(trait, "MONO")) return(is_mono)
  if (identical(trait, "LYM")) return(is_lymph)
  if (trait %in% c("MCH", "RDW", "RBC")) return(is_ery)
  if (identical(trait, "PLT")) return(is_mk)
  if (identical(trait, "EOS")) return(is_gran)
  if (identical(trait, "WBC")) {
    return((is_mono | is_lymph | is_gran | is_dc | is_prog) & !is_ery)
  }

  rep(FALSE, length(cell_type))
}

expected_target_label <- function(trait) {
  if (identical(trait, "MONO")) return("monocyte")
  if (identical(trait, "LYM")) return("lymphoid")
  if (trait %in% c("MCH", "RDW", "RBC")) return("erythroid_or_red_proxy")
  if (identical(trait, "PLT")) return("megakaryocyte_or_platelet")
  if (identical(trait, "EOS")) return("eosinophil_or_granulocyte")
  if (identical(trait, "WBC")) return("broad_non_erythroid_leukocyte")
  "unknown"
}

auc_rank <- function(score, label) {
  keep <- is.finite(score) & !is.na(label)
  score <- score[keep]
  label <- as.integer(label[keep])
  n_pos <- sum(label == 1L)
  n_neg <- sum(label == 0L)
  if (n_pos < 10L || n_neg < 10L) return(NA_real_)
  if (length(unique(score)) < 2L) return(NA_real_)
  r <- rank(score, ties.method = "average")
  (sum(r[label == 1L]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

top_label <- function(dt, group_col, n = 10L) {
  if (nrow(dt) == 0L) return("")
  tab <- dt[, .N, by = group_col][order(-N)]
  paste(paste0(tab[[group_col]][seq_len(min(n, nrow(tab)))], ":",
               tab$N[seq_len(min(n, nrow(tab)))]),
        collapse = "; ")
}

call_status <- function(trait, has_score, n_expected, n_total,
                        top_expected_frac, enrichment, auc) {
  if (n_total == 0L) return("FAIL_EMPTY_REFERENCE")
  if (n_expected < 10L) return("FAIL_NO_TARGET_CELLS")
  if (!has_score) return("NEEDS_SCDRS")
  if (identical(trait, "WBC")) {
    if (is.finite(top_expected_frac) && top_expected_frac >= 0.60) {
      return("BROAD_TARGET_REVIEW")
    }
    return("FAIL_WEAK_WBC_COMPONENT")
  }
  if (is.finite(top_expected_frac) &&
      is.finite(enrichment) &&
      top_expected_frac >= 0.60 &&
      enrichment >= 1.50) {
    return("PASS_STRONG")
  }
  if (is.finite(top_expected_frac) &&
      is.finite(enrichment) &&
      top_expected_frac >= 0.35 &&
      enrichment >= 1.10) {
    return("PASS_WEAK")
  }
  if (is.finite(auc) && auc < 0.50) return("FAIL_ANTI_TARGET")
  "FAIL_WRONG_SOURCE"
}

cat(">>> Check existing RNA references for trait suitability\n\n")

summary_rows <- list()
celltype_rows <- list()
top_rows <- list()

for (ref_name in names(REFERENCE_CONFIG)) {
  cfg <- REFERENCE_CONFIG[[ref_name]]
  cat("=== Reference: ", ref_name, " ===\n", sep = "")
  cat("Label: ", cfg$label, "\n", sep = "")

  meta <- read_reference_meta(cfg, ref_name)
  scores <- merge_reference_scores(cfg$score_files)
  ref_dt <- merge(meta, scores, by = "cell", all.x = TRUE, sort = FALSE)
  cat("Cells: ", nrow(ref_dt), "\n", sep = "")
  cat("Score columns: ",
      paste(grep("_norm_score$", names(ref_dt), value = TRUE), collapse = ", "),
      "\n\n", sep = "")

  for (trait in TRAITS_TO_CHECK) {
    score_col <- paste0(trait, "_norm_score")
    has_score <- score_col %in% names(ref_dt)
    expected <- is_expected_for_trait(trait, ref_dt$cell_type)
    expected[is.na(expected)] <- FALSE
    n_expected <- sum(expected)
    expected_rate <- n_expected / max(nrow(ref_dt), 1L)

    top_expected_frac <- NA_real_
    top_n <- NA_integer_
    top_celltypes <- ""
    auc <- NA_real_
    enrichment <- NA_real_
    score_median <- NA_real_
    score_max <- NA_real_

    if (has_score) {
      score <- as.numeric(ref_dt[[score_col]])
      score_median <- median(score, na.rm = TRUE)
      score_max <- max(score, na.rm = TRUE)
      auc <- auc_rank(score, expected)

      keep <- is.finite(score)
      if (any(keep)) {
        tmp_score <- data.table(
          cell = ref_dt$cell[keep],
          cell_type = ref_dt$cell_type[keep],
          score = score[keep],
          expected_target = expected[keep]
        )

        top_n <- max(1L, ceiling(nrow(tmp_score) * 0.01))
        ord <- order(tmp_score$score, decreasing = TRUE)
        top_cells <- tmp_score[ord[seq_len(min(top_n, length(ord)))]]
        top_expected_frac <- mean(top_cells$expected_target)
        enrichment <- top_expected_frac / max(expected_rate, 1e-12)
        top_celltypes <- top_label(top_cells, "cell_type", n = 12L)

        ct <- tmp_score[, .(
          n_cells = .N,
          expected_target_fraction = mean(expected_target),
          n_expected_target = sum(expected_target),
          mean_score = mean(score, na.rm = TRUE),
          median_score = median(score, na.rm = TRUE),
          q90_score = as.numeric(stats::quantile(score, 0.90, na.rm = TRUE)),
          n_top1pct = sum(cell %in% top_cells$cell)
        ), by = cell_type][order(-mean_score)]
        ct[, `:=`(
          reference = ref_name,
          reference_label = cfg$label,
          trait = trait
        )]
        setcolorder(ct, c("reference", "reference_label", "trait", "cell_type"))
        celltype_rows[[paste(ref_name, trait, sep = "_")]] <- ct

        top_dt <- copy(top_cells)
        top_dt[, `:=`(
          reference = ref_name,
          reference_label = cfg$label,
          trait = trait
        )]
        setcolorder(top_dt, c("reference", "reference_label", "trait"))
        top_rows[[paste(ref_name, trait, sep = "_")]] <- top_dt
      }
    }

    status <- call_status(
      trait = trait,
      has_score = has_score,
      n_expected = n_expected,
      n_total = nrow(ref_dt),
      top_expected_frac = top_expected_frac,
      enrichment = enrichment,
      auc = auc
    )

    note <- ""
    if (identical(trait, "PLT")) {
      note <- "Requires MK/platelet-like RNA source; red-lineage proxy alone is not enough."
    } else if (identical(trait, "EOS")) {
      note <- "Requires eosinophil/granulocyte-like RNA source."
    } else if (identical(trait, "WBC")) {
      note <- "Broad target; consider component-wise mono/lymph/granulocyte references."
    }

    summary_rows[[paste(ref_name, trait, sep = "_")]] <- data.table(
      reference = ref_name,
      reference_label = cfg$label,
      trait = trait,
      expected_target = expected_target_label(trait),
      status = status,
      has_existing_scdrs_score = has_score,
      n_cells = nrow(ref_dt),
      n_expected_target_cells = n_expected,
      expected_target_cell_fraction = expected_rate,
      score_median = score_median,
      score_max = score_max,
      target_auroc = auc,
      top1pct_n_cells = top_n,
      top1pct_expected_target_fraction = top_expected_frac,
      top1pct_expected_target_enrichment = enrichment,
      top1pct_celltypes = top_celltypes,
      note = note
    )
  }
}

summary_dt <- rbindlist(summary_rows, use.names = TRUE, fill = TRUE)
celltype_dt <- rbindlist(celltype_rows, use.names = TRUE, fill = TRUE)
top_dt <- rbindlist(top_rows, use.names = TRUE, fill = TRUE)

setorder(summary_dt, trait, reference)

summary_file <- file.path(OUT_DIR, "existing_rna_reference_trait_summary.csv")
celltype_file <- file.path(OUT_DIR, "existing_rna_reference_celltype_scores.csv")
top_file <- file.path(OUT_DIR, "existing_rna_reference_top1pct_cells.csv")

fwrite(summary_dt, summary_file)
fwrite(celltype_dt, celltype_file)
fwrite(top_dt, top_file)

recommended <- summary_dt[
  status %in% c("PASS_STRONG", "PASS_WEAK", "BROAD_TARGET_REVIEW", "NEEDS_SCDRS"),
  .SD[order(
    factor(status, levels = c("PASS_STRONG", "PASS_WEAK", "BROAD_TARGET_REVIEW", "NEEDS_SCDRS")),
    -top1pct_expected_target_fraction,
    -expected_target_cell_fraction
  )[1]],
  by = trait
]
recommended_file <- file.path(OUT_DIR, "recommended_existing_reference_map.csv")
fwrite(recommended, recommended_file)

cat("\n=== Existing RNA reference trait summary ===\n")
print(summary_dt[, .(
  trait,
  reference,
  status,
  has_existing_scdrs_score,
  n_expected_target_cells,
  expected_target_cell_fraction = signif(expected_target_cell_fraction, 3),
  target_auroc = signif(target_auroc, 3),
  top1pct_expected_target_fraction = signif(top1pct_expected_target_fraction, 3),
  top1pct_expected_target_enrichment = signif(top1pct_expected_target_enrichment, 3),
  top1pct_celltypes
)])

cat("\n=== Recommended existing reference map ===\n")
print(recommended[, .(
  trait,
  reference,
  status,
  has_existing_scdrs_score,
  n_expected_target_cells,
  top1pct_expected_target_fraction = signif(top1pct_expected_target_fraction, 3),
  note
)])

cat("\nWrote:\n")
cat("  ", summary_file, "\n", sep = "")
cat("  ", celltype_file, "\n", sep = "")
cat("  ", top_file, "\n", sep = "")
cat("  ", recommended_file, "\n", sep = "")

cat("\nDone. This was a reference check only; no downstream mapping was run.\n")
