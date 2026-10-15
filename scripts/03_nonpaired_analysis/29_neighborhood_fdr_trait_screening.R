# module4 trait screening: RNA neighborhood-level FDR for candidate traits
# -----------------------------------------------------------------------------
# Purpose:
#   Use the already computed scDRS full_score files for the candidate traits and
#   test whether existing unsupervised CITE RNA neighborhoods are trait-associated.
#
# This is still RNA-side screening only. It does not run ATAC mapping, transfer,
# RWR, or AUROC.
#
# Method boundary:
#   - RNA neighborhoods were built without celltype.
#   - FDR testing uses scDRS real scores and matched-control scores.
#   - celltype/lineage labels are used only for post hoc interpretation.
#
# Run from the repository root with optional SCTARF_* path overrides.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

source(file.path("scripts", "_paths.R"))
ROOT <- sctarf_nonpaired_root()
GS_FILE <- file.path(ROOT, "01_data_prep_external_rna/03_cite_scdrs_python/cite_trait_screening.gs")
FULL_SCORE_DIR <- file.path(ROOT, "01_data_prep_external_rna/03_cite_scdrs_python/output_trait_screening")
GROUP_FILE <- file.path(
  ROOT,
  "04_evaluation/module4_v2_diagnostics/neighborhood_fdr/rna_unsupervised_group_membership.csv.gz"
)
OUT_DIR <- file.path(ROOT, "01_data_prep_external_rna/output/trait_screening")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

if (!exists("GROUP_FDR_THRESHOLD")) GROUP_FDR_THRESHOLD <- 0.10

short_name <- function(trait_id) {
  x <- sub("^UKB_460K\\.blood_", "", trait_id)
  x <- gsub("_COUNT$", "", x)
  x <- gsub("MEAN_CORPUSCULAR_HEMOGLOBIN", "MCH", x)
  x <- gsub("RBC_DISTRIB_WIDTH", "RDW", x)
  x <- gsub("RED", "RBC", x)
  x
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

top_table_label <- function(x, n = 10L) {
  x <- x[!is.na(x)]
  tab <- sort(table(x), decreasing = TRUE)
  if (length(tab) == 0L) return("")
  paste(paste0(names(head(tab, n)), ":", as.integer(head(tab, n))), collapse = "; ")
}

read_full_score_needed <- function(path) {
  dt <- fread(path)
  ctrl_cols <- grep("^ctrl_norm_score_", names(dt), value = TRUE)
  required <- c("cell", "norm_score", ctrl_cols)
  if (!all(c("cell", "norm_score") %in% names(dt)) || length(ctrl_cols) == 0L) {
    stop("full_score file must contain cell, norm_score, and ctrl_norm_score_* columns: ",
         path, call. = FALSE)
  }
  dt[, ..required]
}

test_groups <- function(score_dt, group_dt) {
  dt <- merge(
    group_dt[, .(cell, group_id, cell_type_posthoc, lineage_posthoc)],
    score_dt,
    by = "cell",
    all.x = FALSE,
    all.y = FALSE,
    sort = FALSE
  )
  ctrl_cols <- grep("^ctrl_norm_score_", names(dt), value = TRUE)
  if (nrow(dt) == 0L) stop("No overlapping cells between group file and full_score.", call. = FALSE)

  group_ids <- sort(unique(dt$group_id))
  rows <- vector("list", length(group_ids))

  for (i in seq_along(group_ids)) {
    gid <- group_ids[i]
    sub <- dt[group_id == gid]
    obs <- mean(sub$norm_score, na.rm = TRUE)
    null <- colMeans(as.matrix(sub[, ..ctrl_cols]), na.rm = TRUE)
    p_mc <- (1 + sum(null >= obs, na.rm = TRUE)) / (1 + length(null))

    rows[[i]] <- data.table(
      group_id = gid,
      n_scored_cells = nrow(sub),
      observed_mean_score = obs,
      null_mean = mean(null, na.rm = TRUE),
      null_sd = sd(null, na.rm = TRUE),
      p_mc = p_mc,
      z_vs_null = if (sd(null, na.rm = TRUE) > 0) {
        (obs - mean(null, na.rm = TRUE)) / sd(null, na.rm = TRUE)
      } else {
        NA_real_
      },
      top_celltypes_posthoc = top_table_label(sub$cell_type_posthoc, 8L),
      top_lineages_posthoc = top_table_label(sub$lineage_posthoc, 8L)
    )
  }

  out <- rbindlist(rows, use.names = TRUE, fill = TRUE)
  out[, qval := p.adjust(p_mc, method = "BH")]
  setorder(out, qval, -observed_mean_score)
  out
}

cat(">>> RNA neighborhood-FDR trait screening\n\n")
cat("Group file: ", GROUP_FILE, "\n", sep = "")
cat("Full score dir: ", FULL_SCORE_DIR, "\n", sep = "")
cat("GROUP_FDR_THRESHOLD = ", GROUP_FDR_THRESHOLD, "\n\n", sep = "")

missing <- c(GS_FILE, GROUP_FILE)[!file.exists(c(GS_FILE, GROUP_FILE))]
if (length(missing) > 0L) {
  stop("Missing required file(s):\n  ", paste(missing, collapse = "\n  "), call. = FALSE)
}

gs <- fread(GS_FILE)
traits <- gs$TRAIT
group_dt <- fread(GROUP_FILE)
if (!all(c("cell", "group_id", "cell_type_posthoc") %in% names(group_dt))) {
  stop("Group file must contain cell, group_id, and cell_type_posthoc.", call. = FALSE)
}
group_dt[, lineage_posthoc := lineage_group(cell_type_posthoc)]

summary_rows <- list()
group_rows <- list()

for (trait_id in traits) {
  short <- short_name(trait_id)
  full_file <- file.path(FULL_SCORE_DIR, paste0(trait_id, ".full_score.gz"))
  if (!file.exists(full_file)) {
    warning("Missing full_score for ", trait_id, ": ", full_file)
    next
  }

  cat("=== ", short, " / ", trait_id, " ===\n", sep = "")
  score_dt <- read_full_score_needed(full_file)
  group_test <- test_groups(score_dt, group_dt)
  group_test[, `:=`(trait = short, trait_id = trait_id)]
  setcolorder(group_test, c("trait", "trait_id", "group_id"))
  group_rows[[trait_id]] <- group_test

  sig <- group_test[qval <= GROUP_FDR_THRESHOLD & observed_mean_score > 0]
  sig_cells <- group_dt[group_id %in% sig$group_id]
  expected <- expected_lineage(short)
  expected_vec <- expected_hit(sig_cells$lineage_posthoc, expected)

  summary_rows[[trait_id]] <- data.table(
    trait = short,
    trait_id = trait_id,
    n_groups_tested = nrow(group_test),
    min_group_qval = min(group_test$qval, na.rm = TRUE),
    group_fdr_threshold = GROUP_FDR_THRESHOLD,
    n_sig_groups = nrow(sig),
    n_sig_group_cells = nrow(sig_cells),
    expected_lineage_posthoc = expected,
    sig_group_expected_lineage_fraction_posthoc = mean(expected_vec, na.rm = TRUE),
    sig_group_top_celltypes_posthoc = top_table_label(sig_cells$cell_type_posthoc, 12L),
    sig_group_top_lineages_posthoc = top_table_label(sig_cells$lineage_posthoc, 8L),
    no_celltype_used_for_grouping_or_fdr = TRUE
  )

  cat("  groups tested=", nrow(group_test),
      "; min group qval=", signif(min(group_test$qval, na.rm = TRUE), 4),
      "; q<=", GROUP_FDR_THRESHOLD, " groups=", nrow(sig),
      "; expected lineage fraction=",
      signif(summary_rows[[trait_id]]$sig_group_expected_lineage_fraction_posthoc, 4),
      "\n", sep = "")
}

summary_dt <- rbindlist(summary_rows, use.names = TRUE, fill = TRUE)
group_dt_out <- rbindlist(group_rows, use.names = TRUE, fill = TRUE)

fwrite(summary_dt, file.path(OUT_DIR, "cite_trait_screening_neighborhood_fdr_summary.csv"))
fwrite(group_dt_out, file.path(OUT_DIR, "cite_trait_screening_neighborhood_fdr_group_tests.csv"))

cat("\n=== RNA neighborhood-FDR trait-screening summary ===\n")
print(summary_dt[, .(
  trait,
  min_group_qval,
  n_sig_groups,
  n_sig_group_cells,
  expected_lineage_posthoc,
  sig_group_expected_lineage_fraction_posthoc,
  sig_group_top_celltypes_posthoc
)])
cat("\nOutputs are in: ", OUT_DIR, "\n", sep = "")
