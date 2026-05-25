# 01: Build per-cell gamma vectors for the two ablation versions.
#
# Version A (discrete):
#   gamma_i = gamma_high  if cell i ∈ progenitor / transition celltype set
#   gamma_i = gamma_low   otherwise
#
# Version B (continuous, rna_reliability based):
#   r_i        = empirical rank percentile of rna_reliability_i (in [0, 1])
#   gamma_i    = gamma_low + (gamma_high - gamma_low) * (1 - r_i)
#
# Both versions produce one row per cell with a gamma value in (0, 1).
#
# Outputs:
#   gamma_versionA.csv   (cell, BioClassification, is_transition, gamma)
#   gamma_versionB.csv   (cell, rna_reliability, rank_pct, gamma)
#   gamma_summary.csv    (version, n_cells, mean, sd, q10, q50, q90)
#   gamma_assignment_manifest.json
#
# Parameters:
#   MODULE2_GAMMA_HIGH       (default 0.15)
#   MODULE2_GAMMA_LOW        (default 0.03)

source("D:/RStudio/work/module2_adaptive_gamma/_bootstrap.R")

required_pkgs(c("data.table", "Matrix", "SummarizedExperiment", "jsonlite"))
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(SummarizedExperiment)
  library(jsonlite)
})

paths <- adaptive_gamma_paths()
ensure_adaptive_gamma_dirs(paths)

gamma_high <- env_numeric("MODULE2_GAMMA_HIGH", 0.15)
gamma_low <- env_numeric("MODULE2_GAMMA_LOW", 0.03)
if (gamma_high <= gamma_low) stop("gamma_high must exceed gamma_low", call. = FALSE)
if (gamma_high >= 1 || gamma_low <= 0) stop("gamma values must be in (0, 1)", call. = FALSE)

cat("========== Module 2 adaptive gamma: 01 build gamma vectors ==========\n")
cat("  gamma_high :", gamma_high, "(progenitor / low-reliability)\n")
cat("  gamma_low  :", gamma_low,  "(mature / high-reliability)\n")
cat("  out_dir    :", paths$out_dir, "\n")

stopifnot(file.exists(paths$atac_se))
stopifnot(file.exists(paths$wnn_lite_weights))

# Anchor cell list to the original LSI mutualKNN graph rownames so downstream
# RWR sees a consistent gamma vector.
stopifnot(file.exists(paths$original_knn))
original_knn <- readRDS(paths$original_knn)
cells_in_graph <- rownames(original_knn)
n_cells <- length(cells_in_graph)
cat("  cells in graph:", n_cells, "\n")

# Version A -------------------------------------------------------------

cat("\n[1/3] Version A (discrete, by BioClassification)\n")
SE_atac <- readRDS(paths$atac_se)
meta <- as.data.table(as.data.frame(colData(SE_atac)))
meta[, cell := colnames(SE_atac)]
if (!"BioClassification" %in% names(meta)) {
  stop("ATAC metadata must have BioClassification column", call. = FALSE)
}
meta_a <- meta[, .(cell, BioClassification)]
meta_a <- meta_a[match(cells_in_graph, cell)]
if (any(is.na(meta_a$cell))) {
  stop(sum(is.na(meta_a$cell)), " cells in graph are missing from ATAC SE", call. = FALSE)
}
meta_a[, is_transition := is_transition_celltype(BioClassification)]
meta_a[, gamma := ifelse(is_transition, gamma_high, gamma_low)]

cat("  progenitor / transition cells:", sum(meta_a$is_transition),
    "(", round(mean(meta_a$is_transition) * 100, 2), "%)\n")
cat("  gamma summary (Version A):\n")
print(summary(meta_a$gamma))

versionA_out <- file.path(paths$out_dir, "gamma_versionA.csv")
fwrite(meta_a[, .(cell, BioClassification, is_transition, gamma)], versionA_out)

# Version B -------------------------------------------------------------

cat("\n[2/3] Version B (continuous, by rna_reliability)\n")
weights_dt <- fread(paths$wnn_lite_weights)
if (!all(c("atac_cell", "rna_reliability") %in% names(weights_dt))) {
  stop("wnn_lite_weights must contain atac_cell and rna_reliability", call. = FALSE)
}
rel_map <- setNames(as.numeric(weights_dt$rna_reliability), weights_dt$atac_cell)
rel_for_graph <- rel_map[cells_in_graph]
n_missing_rel <- sum(is.na(rel_for_graph))
if (n_missing_rel > 0) {
  warning(n_missing_rel, " cells missing rna_reliability; ",
          "imputing with median for gamma calculation.")
  rel_for_graph[is.na(rel_for_graph)] <- median(rel_for_graph, na.rm = TRUE)
}

# Empirical rank percentile in [0, 1]; ties get average rank.
rank_pct <- rank(rel_for_graph, ties.method = "average") / length(rel_for_graph)
gamma_B <- gamma_low + (gamma_high - gamma_low) * (1 - rank_pct)

versionB_dt <- data.table(
  cell = cells_in_graph,
  rna_reliability = as.numeric(rel_for_graph),
  rank_pct = as.numeric(rank_pct),
  gamma = as.numeric(gamma_B)
)
versionB_out <- file.path(paths$out_dir, "gamma_versionB.csv")
fwrite(versionB_dt, versionB_out)

cat("  rna_reliability summary:\n")
print(summary(versionB_dt$rna_reliability))
cat("  gamma summary (Version B):\n")
print(summary(versionB_dt$gamma))

# Cross-check: how does Version A overlap with Version B's high-gamma tail?
cat("\n[3/3] Cross-version diagnostic\n")
cross <- merge(
  meta_a[, .(cell, is_transition, gamma_A = gamma)],
  versionB_dt[, .(cell, rna_reliability, gamma_B = gamma)],
  by = "cell"
)
cor_AB <- cor(cross$gamma_A, cross$gamma_B, method = "spearman")
cat("  Spearman corr(gamma_A, gamma_B):", round(cor_AB, 3), "\n")

cross_summary_path <- file.path(paths$out_dir, "gamma_cross_version.csv")
fwrite(cross, cross_summary_path)

# Summary stats ---------------------------------------------------------

summary_dt <- rbindlist(list(
  data.table(
    version = "A_discrete",
    n_cells = nrow(meta_a),
    n_high_gamma = sum(meta_a$is_transition),
    pct_high_gamma = mean(meta_a$is_transition),
    mean_gamma = mean(meta_a$gamma),
    sd_gamma = sd(meta_a$gamma),
    q10 = quantile(meta_a$gamma, 0.1),
    q50 = quantile(meta_a$gamma, 0.5),
    q90 = quantile(meta_a$gamma, 0.9)
  ),
  data.table(
    version = "B_continuous",
    n_cells = nrow(versionB_dt),
    n_high_gamma = NA_integer_,
    pct_high_gamma = NA_real_,
    mean_gamma = mean(versionB_dt$gamma),
    sd_gamma = sd(versionB_dt$gamma),
    q10 = quantile(versionB_dt$gamma, 0.1),
    q50 = quantile(versionB_dt$gamma, 0.5),
    q90 = quantile(versionB_dt$gamma, 0.9)
  )
), use.names = TRUE, fill = TRUE)
summary_path <- file.path(paths$out_dir, "gamma_summary.csv")
fwrite(summary_dt, summary_path)
cat("\n  gamma_summary.csv:\n")
print(summary_dt)

# Manifest --------------------------------------------------------------

manifest <- list(
  script = "01_build_gamma_vectors.R",
  parameters = list(
    gamma_high = gamma_high,
    gamma_low = gamma_low
  ),
  inputs = list(
    atac_se = paths$atac_se,
    wnn_lite_weights = paths$wnn_lite_weights,
    original_knn = paths$original_knn
  ),
  outputs = list(
    versionA = versionA_out,
    versionB = versionB_out,
    cross = cross_summary_path,
    summary = summary_path
  ),
  n_cells = n_cells,
  n_progenitor_cells_versionA = sum(meta_a$is_transition),
  spearman_AB = cor_AB,
  notes = c(
    "Version A flags HSC/CMP/LMPP/CLP/GMP/Pre/Early/Unk celltypes (consistent with module2_wnn_lite_gated).",
    "Version B uses empirical rank percentile of rna_reliability (low reliability => high gamma).",
    "Both versions produce gamma in [gamma_low, gamma_high]."
  ),
  created_at = as.character(Sys.time())
)
manifest_out <- file.path(paths$out_dir, "gamma_assignment_manifest.json")
write_json(manifest, manifest_out, pretty = TRUE, auto_unbox = TRUE)

cat("\n========== Done ==========\n")
cat("  manifest:", manifest_out, "\n")
