# Project trait-level gene risk onto ATAC cells through Gene Activity.
#
# This script does not run TRS. It builds per-cell gene_risk_cell_score and
# checks whether it is independent enough from the original ATAC z_score.

source("D:/RStudio/work/module2_gwas_gene_risk/_bootstrap.R")

required_pkgs(c("data.table", "Matrix", "ggplot2", "jsonlite"))
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(ggplot2)
  library(jsonlite)
})

paths <- gwas_gene_paths()
ensure_gwas_gene_dirs(paths)

traits <- strsplit(Sys.getenv("MODULE2_TRAITS", "MONO,PLT,LYM"), ",", fixed = TRUE)[[1]]
traits <- trimws(traits)
traits <- traits[nzchar(traits)]
if (length(traits) == 0) stop("MODULE2_TRAITS is empty", call. = FALSE)

activity_transform <- Sys.getenv("MODULE2_GENE_ACTIVITY_TRANSFORM", "log1p")
if (!activity_transform %in% c("log1p", "raw")) {
  stop("MODULE2_GENE_ACTIVITY_TRANSFORM must be log1p or raw", call. = FALSE)
}

gene_risk_path <- Sys.getenv(
  "MODULE2_TRAIT_GENE_RISK",
  file.path(paths$out_dir, "trait_gene_risk", "trait_gene_risk_long.csv")
)
seed_signal_path <- Sys.getenv(
  "MODULE2_SEED_SIGNAL_BY_CELL",
  file.path(paths$out_dir, "seed_signal_diagnostics", "seed_signal_by_cell.csv.gz")
)
out_dir <- file.path(paths$out_dir, "gene_risk_cell_projection")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("========== Module 2 GWAS gene-risk: project gene risk to cells ==========\n")
cat("  traits        :", paste(traits, collapse = ", "), "\n")
cat("  gene risk     :", gene_risk_path, "\n")
cat("  gene activity :", paths$gene_activity, "\n")
cat("  seed signal   :", seed_signal_path, "\n")
cat("  transform     :", activity_transform, "\n")
cat("  output dir    :", out_dir, "\n")

stopifnot(file.exists(gene_risk_path), file.exists(paths$gene_activity), file.exists(seed_signal_path))

rank01 <- function(x) {
  x <- as.numeric(x)
  ok <- is.finite(x)
  out <- rep(NA_real_, length(x))
  if (!any(ok)) {
    out[] <- NA_real_
    return(out)
  }
  r <- rank(x[ok], ties.method = "average", na.last = "keep")
  out[ok] <- (r - 0.5) / length(r)
  out
}

transform_activity <- function(mat, mode) {
  if (mode == "raw") return(mat)
  out <- mat
  out@x <- log1p(out@x)
  out
}

weighted_gene_activity_score <- function(ga, genes, weights, transform) {
  keep <- genes %in% rownames(ga) & is.finite(weights) & weights > 0
  genes <- genes[keep]
  weights <- as.numeric(weights[keep])
  if (length(genes) == 0) {
    return(list(score = rep(NA_real_, ncol(ga)), genes_used = character(), weights_used = numeric()))
  }

  ga_sub <- ga[genes, , drop = FALSE]
  ga_sub <- transform_activity(ga_sub, transform)
  score <- as.numeric(Matrix::crossprod(weights, ga_sub)) / sum(weights)
  names(score) <- colnames(ga)
  list(score = score, genes_used = genes, weights_used = weights)
}

cat("\n========== Load inputs ==========\n")
gene_risk <- data.table::fread(gene_risk_path)
seed_signal <- data.table::fread(seed_signal_path)
gene_activity <- readRDS(paths$gene_activity)

if (!inherits(gene_activity, "dgCMatrix")) {
  gene_activity <- as(gene_activity, "dgCMatrix")
}
if (is.null(rownames(gene_activity)) || is.null(colnames(gene_activity))) {
  stop("Gene Activity matrix must have rownames and colnames", call. = FALSE)
}
if (!all(c("trait", "gene_id", "gene_risk_score") %in% colnames(gene_risk))) {
  stop("Gene risk table must contain trait, gene_id, gene_risk_score", call. = FALSE)
}
if (!all(c("trait", "cell", "z_score", "BioClassification", "target_type") %in% colnames(seed_signal))) {
  stop("Seed signal table must contain trait, cell, z_score, BioClassification, target_type", call. = FALSE)
}

common_cells <- intersect(colnames(gene_activity), unique(seed_signal$cell))
if (length(common_cells) == 0) {
  stop("No common cells between Gene Activity and seed signal table", call. = FALSE)
}
cell_ids <- colnames(gene_activity)[colnames(gene_activity) %in% common_cells]
gene_activity <- gene_activity[, cell_ids, drop = FALSE]

cat("  Gene Activity:", nrow(gene_activity), "genes x", ncol(gene_activity), "cells\n")
cat("  common cells :", length(cell_ids), "\n")

cell_score_list <- list()
coverage_list <- list()
correlation_list <- list()
top_enrichment_list <- list()

for (tr in traits) {
  cat("\n[Trait]", tr, "\n")
  gr <- gene_risk[trait == tr]
  if (nrow(gr) == 0) {
    warning("No gene risk rows for trait: ", tr)
    next
  }
  ss <- seed_signal[trait == tr & cell %in% cell_ids]
  if (nrow(ss) == 0) {
    warning("No seed signal rows for trait: ", tr)
    next
  }
  setkey(ss, cell)
  ss <- ss[cell_ids]
  stopifnot(identical(ss$cell, cell_ids))

  projection <- weighted_gene_activity_score(
    ga = gene_activity,
    genes = gr$gene_id,
    weights = gr$gene_risk_score,
    transform = activity_transform
  )
  score <- projection$score
  score_rank01 <- rank01(score)
  z_rank01 <- rank01(ss$z_score)

  dt <- data.table(
    trait = tr,
    cell = cell_ids,
    gene_risk_cell_score = as.numeric(score),
    gene_risk_cell_rank01 = score_rank01,
    atac_z_score = as.numeric(ss$z_score),
    atac_z_rank01 = z_rank01,
    seed = as.logical(ss$seed),
    BioClassification = ss$BioClassification,
    target_type = as.logical(ss$target_type),
    target_definition = ss$target_definition
  )
  cell_score_list[[tr]] <- dt

  n_genes_total <- nrow(gr)
  n_genes_used <- length(projection$genes_used)
  positive_cells <- is.finite(score) & score > 0
  coverage_list[[tr]] <- data.table(
    trait = tr,
    n_gene_risk_genes = n_genes_total,
    n_genes_used_in_gene_activity = n_genes_used,
    pct_genes_used = n_genes_used / n_genes_total,
    n_cells = length(score),
    n_cells_score_positive = sum(positive_cells),
    pct_cells_score_positive = mean(positive_cells),
    mean_gene_risk_cell_score = mean(score, na.rm = TRUE),
    median_gene_risk_cell_score = median(score, na.rm = TRUE),
    max_gene_risk_cell_score = max(score, na.rm = TRUE),
    activity_transform = activity_transform
  )

  correlation_list[[tr]] <- data.table(
    trait = tr,
    pearson_score_vs_atac_z = suppressWarnings(cor(dt$gene_risk_cell_score, dt$atac_z_score, method = "pearson", use = "pairwise.complete.obs")),
    spearman_score_vs_atac_z = suppressWarnings(cor(dt$gene_risk_cell_score, dt$atac_z_score, method = "spearman", use = "pairwise.complete.obs")),
    pearson_rank_vs_atac_rank = suppressWarnings(cor(dt$gene_risk_cell_rank01, dt$atac_z_rank01, method = "pearson", use = "pairwise.complete.obs")),
    spearman_rank_vs_atac_rank = suppressWarnings(cor(dt$gene_risk_cell_rank01, dt$atac_z_rank01, method = "spearman", use = "pairwise.complete.obs")),
    high_collinearity_flag = suppressWarnings(abs(cor(dt$gene_risk_cell_rank01, dt$atac_z_rank01, method = "spearman", use = "pairwise.complete.obs")) > 0.7)
  )

  for (score_col in c("gene_risk_cell_score", "atac_z_score")) {
    score_vec <- dt[[score_col]]
    for (top_frac in c(0.05, 0.10)) {
      n_top <- max(1L, ceiling(nrow(dt) * top_frac))
      top_idx <- order(score_vec, decreasing = TRUE, na.last = NA)[seq_len(n_top)]
      top_target <- mean(dt$target_type[top_idx], na.rm = TRUE)
      all_target <- mean(dt$target_type, na.rm = TRUE)
      top_enrichment_list[[paste(tr, score_col, top_frac, sep = "_")]] <- data.table(
        trait = tr,
        score = score_col,
        top_frac = top_frac,
        n_top = n_top,
        target_fraction_top = top_target,
        target_fraction_all = all_target,
        enrichment_ratio = if (is.finite(all_target) && all_target > 0) top_target / all_target else NA_real_,
        target_definition = target_definition(tr)
      )
    }
  }

  cat("  genes used:", n_genes_used, "/", n_genes_total, "\n")
  cat("  positive cells:", sum(positive_cells), "/", length(score), "\n")
}

if (length(cell_score_list) == 0) {
  stop("No gene-risk cell projections were generated", call. = FALSE)
}

cell_scores <- rbindlist(cell_score_list, use.names = TRUE, fill = TRUE)
coverage <- rbindlist(coverage_list, use.names = TRUE, fill = TRUE)
correlations <- rbindlist(correlation_list, use.names = TRUE, fill = TRUE)
top_enrichment <- rbindlist(top_enrichment_list, use.names = TRUE, fill = TRUE)

cell_scores_out <- file.path(out_dir, "gene_risk_cell_scores.csv.gz")
coverage_out <- file.path(out_dir, "gene_risk_cell_score_coverage.csv")
correlation_out <- file.path(out_dir, "gene_risk_vs_atac_z_correlation.csv")
top_enrichment_out <- file.path(out_dir, "gene_risk_cell_top_enrichment.csv")
data.table::fwrite(cell_scores, cell_scores_out)
data.table::fwrite(coverage, coverage_out)
data.table::fwrite(correlations, correlation_out)
data.table::fwrite(top_enrichment, top_enrichment_out)

plot_dt <- correlations[, .(
  trait,
  pearson_rank_vs_atac_rank,
  spearman_rank_vs_atac_rank,
  high_collinearity_flag
)]
plot_long <- melt(
  plot_dt,
  id.vars = c("trait", "high_collinearity_flag"),
  variable.name = "correlation_type",
  value.name = "correlation"
)
p_corr <- ggplot(plot_long, aes(x = trait, y = correlation, fill = correlation_type)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) +
  geom_hline(yintercept = 0.7, linetype = "dashed", color = "firebrick") +
  labs(
    title = "Gene-risk cell score vs original ATAC z-score",
    subtitle = "Dashed line marks the high-collinearity threshold (0.7)",
    x = NULL,
    y = "Correlation"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top")
fig_corr <- file.path(out_dir, "Fig_gene_risk_vs_atac_z_correlation.png")
ggsave(fig_corr, p_corr, width = 8, height = 4.8, dpi = 300)

p_top <- ggplot(top_enrichment, aes(x = score, y = target_fraction_top, fill = score)) +
  geom_col(width = 0.65) +
  geom_hline(aes(yintercept = target_fraction_all), linetype = "dashed", color = "grey35") +
  facet_grid(trait ~ top_frac, labeller = label_both) +
  labs(
    title = "Target-cell fraction in top-scoring cells",
    subtitle = "Dashed line is trait target-cell fraction among all ATAC cells",
    x = NULL,
    y = "Target-cell fraction"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")
fig_top <- file.path(out_dir, "Fig_gene_risk_cell_top_enrichment.png")
ggsave(fig_top, p_top, width = 8, height = 6, dpi = 300)

manifest <- list(
  script = "04_project_gene_risk_to_cells.R",
  traits = traits,
  gene_risk = gene_risk_path,
  gene_activity = paths$gene_activity,
  seed_signal = seed_signal_path,
  activity_transform = activity_transform,
  cell_scores = cell_scores_out,
  coverage = coverage_out,
  correlations = correlation_out,
  top_enrichment = top_enrichment_out,
  correlation_figure = fig_corr,
  top_enrichment_figure = fig_top,
  created_at = as.character(Sys.time())
)
jsonlite::write_json(
  manifest,
  file.path(out_dir, "gene_risk_cell_projection_manifest.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

cat("\n========== Done ==========\n")
cat("  cell scores    :", cell_scores_out, "\n")
cat("  coverage       :", coverage_out, "\n")
cat("  correlations   :", correlation_out, "\n")
cat("  top enrichment :", top_enrichment_out, "\n")
