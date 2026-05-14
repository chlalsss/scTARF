# module4 trait screening: summarize CITE RNA scDRS candidate traits
# -----------------------------------------------------------------------------
# Purpose:
#   Summarize RNA-side scDRS results for candidate blood traits before deciding
#   which traits deserve the full module4 pipeline.
#
# This script uses cell_type only for post hoc interpretation. It does not select
# cells for mapping or downstream analysis.
#
# Run in RStudio after scDRS finishes:
#   setwd("D:/RStudio/work/module4_traitaware_pseudopair")
#   source("01_data_prep_external_rna/08_summarize_cite_scdrs_trait_screening.R")
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

ROOT <- "D:/RStudio/work/module4_traitaware_pseudopair"
GS_FILE <- file.path(ROOT, "01_data_prep_external_rna/03_cite_scdrs_python/cite_trait_screening.gs")
SCORE_DIR <- file.path(ROOT, "01_data_prep_external_rna/03_cite_scdrs_python/output_trait_screening")
CITE_SCORE_RDS <- file.path(ROOT, "01_data_prep_external_rna/output/cite_rna_scdrs_scores.rds")
OUT_DIR <- file.path(ROOT, "01_data_prep_external_rna/output/trait_screening")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

if (!exists("TOP_PCT")) TOP_PCT <- 0.01
if (!exists("FDR_THRESHOLD")) FDR_THRESHOLD <- 0.10

score_col <- function(dt) {
  candidates <- c("norm_score", "zscore", "score")
  hit <- intersect(candidates, colnames(dt))[1]
  if (is.na(hit)) stop("No score column found in: ", paste(colnames(dt), collapse = ", "))
  hit
}

pval_col <- function(dt) {
  candidates <- c("pval", "pval_ctrl", "pval_mc")
  hit <- intersect(candidates, colnames(dt))[1]
  if (is.na(hit)) return(NA_character_)
  hit
}

cell_col <- function(dt) {
  candidates <- c("cell", "cell_id", "barcode", "V1")
  hit <- intersect(candidates, colnames(dt))[1]
  if (is.na(hit)) hit <- colnames(dt)[1]
  hit
}

short_name <- function(trait_id) {
  x <- sub("^UKB_460K\\.blood_", "", trait_id)
  x <- gsub("_COUNT$", "", x)
  x <- gsub("MEAN_CORPUSCULAR_HEMOGLOBIN", "MCH", x)
  x <- gsub("RBC_DISTRIB_WIDTH", "RDW", x)
  x <- gsub("RED", "RBC", x)
  x
}

top_table_label <- function(x, n = 10L) {
  x <- x[!is.na(x)]
  tab <- sort(table(x), decreasing = TRUE)
  if (length(tab) == 0L) return("")
  paste(paste0(names(head(tab, n)), ":", as.integer(head(tab, n))), collapse = "; ")
}

lineage_group <- function(cell_type) {
  out <- rep("other", length(cell_type))
  erythroid <- c("MK/E prog", "Proerythroblast", "Erythroblast", "Normoblast", "Reticulocyte")
  mono <- c("CD14+ Mono", "CD16+ Mono")
  myeloid_other <- c("G/M prog", "cDC1", "cDC2", "pDC")
  lymph <- c(
    "CD4+ T activated", "CD4+ T naive", "CD4+ T activated integrinB7+",
    "CD4+ T CD314+ CD45RA+", "CD8+ T naive",
    "CD8+ T naive CD127+ CD26- CD101-", "CD8+ T CD49f+",
    "CD8+ T CD57+ CD45RA+", "CD8+ T CD57+ CD45RO+",
    "CD8+ T CD69+ CD45RA+", "CD8+ T CD69+ CD45RO+",
    "CD8+ T TIGIT+ CD45RA+", "CD8+ T TIGIT+ CD45RO+",
    "T reg", "T prog cycling", "MAIT", "gdT CD158b+", "gdT TCRVD2+",
    "dnT", "NK", "NK CD158e1+", "ILC", "ILC1", "Lymph prog",
    "Naive CD20+ B IGKC+", "Naive CD20+ B IGKC-",
    "B1 B IGKC+", "B1 B IGKC-", "Transitional B",
    "Plasma cell IGKC+", "Plasma cell IGKC-",
    "Plasmablast IGKC+", "Plasmablast IGKC-"
  )
  out[cell_type %in% erythroid] <- "erythroid_mk"
  out[cell_type %in% mono] <- "mono"
  out[cell_type %in% myeloid_other] <- "myeloid_other"
  out[cell_type %in% lymph] <- "lymph"
  out
}

expected_lineage <- function(short) {
  if (short %in% c("RBC", "MCH", "RDW")) return("erythroid_mk")
  if (short %in% c("EOSINOPHIL", "WHITE")) return("leukocyte_non_erythroid")
  NA_character_
}

expected_hit <- function(lineage, expected) {
  if (is.na(expected)) return(rep(NA, length(lineage)))
  if (expected == "leukocyte_non_erythroid") {
    return(lineage %in% c("mono", "myeloid_other", "lymph"))
  }
  lineage == expected
}

safe_mean <- function(x) {
  if (length(x) == 0L) return(NA_real_)
  mean(x, na.rm = TRUE)
}

cat(">>> Summarize CITE RNA scDRS trait screening\n\n")
cat("Score dir: ", SCORE_DIR, "\n", sep = "")
cat("TOP_PCT = ", TOP_PCT, "; FDR_THRESHOLD = ", FDR_THRESHOLD, "\n\n", sep = "")

missing <- c(GS_FILE, CITE_SCORE_RDS)[!file.exists(c(GS_FILE, CITE_SCORE_RDS))]
if (length(missing) > 0L) {
  stop("Missing required file(s):\n  ", paste(missing, collapse = "\n  "), call. = FALSE)
}

gs <- fread(GS_FILE)
traits <- gs$TRAIT
cite_meta <- as.data.table(readRDS(CITE_SCORE_RDS))
if (!all(c("cell", "cell_type") %in% names(cite_meta))) {
  stop("CITE metadata must contain cell and cell_type.", call. = FALSE)
}
cite_meta <- cite_meta[, .(cell, cell_type)]
cite_meta[, lineage := lineage_group(cell_type)]

summary_rows <- list()
celltype_rows <- list()
lineage_rows <- list()

for (trait_id in traits) {
  short <- short_name(trait_id)
  score_file <- file.path(SCORE_DIR, paste0(trait_id, ".score.gz"))
  if (!file.exists(score_file)) {
    warning("Missing score file for ", trait_id, ": ", score_file)
    next
  }

  cat("=== ", short, " / ", trait_id, " ===\n", sep = "")
  dt <- fread(score_file)
  cc <- cell_col(dt)
  sc <- score_col(dt)
  pc <- pval_col(dt)

  score_dt <- data.table(
    cell = as.character(dt[[cc]]),
    norm_score = as.numeric(dt[[sc]])
  )
  if (!is.na(pc)) {
    score_dt[, pval := as.numeric(dt[[pc]])]
    score_dt[, qval := p.adjust(pval, method = "BH")]
  } else {
    score_dt[, `:=`(pval = NA_real_, qval = NA_real_)]
  }

  score_dt <- merge(score_dt, cite_meta, by = "cell", all.x = TRUE, sort = FALSE)
  if (anyNA(score_dt$cell_type)) {
    warning(sum(is.na(score_dt$cell_type)), " cells missing CITE cell_type metadata for ", trait_id)
  }

  top_thr <- as.numeric(quantile(score_dt$norm_score, 1 - TOP_PCT, na.rm = TRUE))
  score_dt[, is_top := norm_score >= top_thr & is.finite(norm_score)]
  score_dt[, is_fdr_sig := qval <= FDR_THRESHOLD & norm_score > 0]

  expected <- expected_lineage(short)
  score_dt[, is_expected_lineage_posthoc := expected_hit(lineage, expected)]

  top_dt <- score_dt[is_top == TRUE]
  sig_dt <- score_dt[is_fdr_sig == TRUE]

  summary_rows[[trait_id]] <- data.table(
    trait = short,
    trait_id = trait_id,
    n_cells = nrow(score_dt),
    score_min = min(score_dt$norm_score, na.rm = TRUE),
    score_median = median(score_dt$norm_score, na.rm = TRUE),
    score_max = max(score_dt$norm_score, na.rm = TRUE),
    min_pval = min(score_dt$pval, na.rm = TRUE),
    min_qval = min(score_dt$qval, na.rm = TRUE),
    fdr_threshold = FDR_THRESHOLD,
    n_fdr_sig_positive = nrow(sig_dt),
    top_pct = TOP_PCT,
    top_threshold = top_thr,
    n_top = nrow(top_dt),
    expected_lineage_posthoc = expected,
    top_expected_lineage_fraction_posthoc = safe_mean(top_dt$is_expected_lineage_posthoc),
    fdr_sig_expected_lineage_fraction_posthoc = safe_mean(sig_dt$is_expected_lineage_posthoc),
    top_celltypes_posthoc = top_table_label(top_dt$cell_type, 12L),
    fdr_sig_celltypes_posthoc = top_table_label(sig_dt$cell_type, 12L),
    no_celltype_used_in_scdrs = TRUE
  )

  by_ct <- score_dt[, .(
    n_cells = .N,
    mean_score = mean(norm_score, na.rm = TRUE),
    median_score = median(norm_score, na.rm = TRUE),
    q90_score = as.numeric(quantile(norm_score, 0.90, na.rm = TRUE)),
    n_top = sum(is_top, na.rm = TRUE),
    n_fdr_sig_positive = sum(is_fdr_sig, na.rm = TRUE)
  ), by = .(cell_type, lineage)][order(-mean_score)]
  by_ct[, `:=`(trait = short, trait_id = trait_id)]
  setcolorder(by_ct, c("trait", "trait_id", "cell_type", "lineage"))
  celltype_rows[[trait_id]] <- by_ct

  by_lin <- score_dt[, .(
    n_cells = .N,
    mean_score = mean(norm_score, na.rm = TRUE),
    median_score = median(norm_score, na.rm = TRUE),
    q90_score = as.numeric(quantile(norm_score, 0.90, na.rm = TRUE)),
    n_top = sum(is_top, na.rm = TRUE),
    n_fdr_sig_positive = sum(is_fdr_sig, na.rm = TRUE)
  ), by = lineage][order(-mean_score)]
  by_lin[, `:=`(trait = short, trait_id = trait_id)]
  setcolorder(by_lin, c("trait", "trait_id", "lineage"))
  lineage_rows[[trait_id]] <- by_lin

  cat("  min qval=", signif(summary_rows[[trait_id]]$min_qval, 4),
      "; q<=", FDR_THRESHOLD, " & score>0 cells=",
      summary_rows[[trait_id]]$n_fdr_sig_positive,
      "; top expected lineage fraction=",
      signif(summary_rows[[trait_id]]$top_expected_lineage_fraction_posthoc, 4),
      "\n", sep = "")
}

summary_dt <- rbindlist(summary_rows, use.names = TRUE, fill = TRUE)
celltype_dt <- rbindlist(celltype_rows, use.names = TRUE, fill = TRUE)
lineage_dt <- rbindlist(lineage_rows, use.names = TRUE, fill = TRUE)

fwrite(summary_dt, file.path(OUT_DIR, "cite_trait_screening_scdrs_summary.csv"))
fwrite(celltype_dt, file.path(OUT_DIR, "cite_trait_screening_by_celltype_posthoc.csv"))
fwrite(lineage_dt, file.path(OUT_DIR, "cite_trait_screening_by_lineage_posthoc.csv"))

cat("\n=== CITE RNA trait-screening summary ===\n")
print(summary_dt[, .(
  trait, min_qval, n_fdr_sig_positive,
  top_expected_lineage_fraction_posthoc,
  fdr_sig_expected_lineage_fraction_posthoc,
  top_celltypes_posthoc
)])
cat("\nOutputs are in: ", OUT_DIR, "\n", sep = "")
