# Run alpha scan for gene-risk seed prior on a fixed cell-cell graph.
#
# This isolates seed-prior effects by defaulting to the original LSI mutualKNN.
# It does not use WNN-lite unless MODULE2_SEED_PRIOR_KNN is explicitly set.

source("D:/RStudio/work/module2_gwas_gene_risk/_bootstrap.R")

required_pkgs(c("SCAVENGE", "data.table", "Matrix", "ggplot2", "jsonlite"))
suppressPackageStartupMessages({
  library(SCAVENGE)
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

alpha_grid <- strsplit(Sys.getenv("MODULE2_ALPHA_GRID", "0,0.1,0.25,0.5,0.75,1"), ",", fixed = TRUE)[[1]]
alpha_grid <- as.numeric(trimws(alpha_grid))
if (any(!is.finite(alpha_grid)) || any(alpha_grid < 0 | alpha_grid > 1)) {
  stop("MODULE2_ALPHA_GRID must contain numeric values in [0, 1]", call. = FALSE)
}
alpha_grid <- sort(unique(alpha_grid))

seed_percent <- env_numeric("MODULE2_SEED_PERCENT", 0.05)
gamma <- env_numeric("MODULE2_SEED_PRIOR_RW_GAMMA", env_numeric("MODULE2_RW_GAMMA", 0.05))
graph_path <- Sys.getenv("MODULE2_SEED_PRIOR_KNN", paths$original_knn)
graph_label <- Sys.getenv("MODULE2_SEED_PRIOR_GRAPH_LABEL", "original_lsi")
cell_score_path <- Sys.getenv(
  "MODULE2_GENE_RISK_CELL_SCORES",
  file.path(paths$out_dir, "gene_risk_cell_projection", "gene_risk_cell_scores.csv.gz")
)
out_dir <- file.path(paths$out_dir, "seed_prior_alpha_scan", graph_label)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("========== Module 2 GWAS gene-risk: seed-prior alpha scan ==========\n")
cat("  traits      :", paste(traits, collapse = ", "), "\n")
cat("  alpha grid  :", paste(alpha_grid, collapse = ", "), "\n")
cat("  seed percent:", seed_percent, "\n")
cat("  gamma       :", gamma, "\n")
cat("  graph       :", graph_path, "\n")
cat("  graph label :", graph_label, "\n")
cat("  cell scores :", cell_score_path, "\n")
cat("  output dir  :", out_dir, "\n")

stopifnot(file.exists(graph_path), file.exists(cell_score_path))

alpha_label <- function(alpha) {
  gsub("\\.", "p", sprintf("%.2f", alpha))
}

top_n_idx <- function(x, frac) {
  n_top <- max(1L, ceiling(length(x) * frac))
  order(x, decreasing = TRUE, na.last = NA)[seq_len(n_top)]
}

safe_scale_factor <- function(z_score) {
  out <- tryCatch(
    SCAVENGE::cal_scalefactor(z_score = z_score, percent_cut = 0.01),
    error = function(e) NA_real_
  )
  if (!is.finite(out)) out <- 1
  out
}

cat("\n========== Load inputs ==========\n")
cell_scores <- data.table::fread(cell_score_path)
graph <- readRDS(graph_path)
graph <- as(graph, "dgCMatrix")

required_cols <- c(
  "trait", "cell", "gene_risk_cell_rank01", "atac_z_rank01",
  "gene_risk_cell_score", "atac_z_score", "BioClassification", "target_type"
)
if (!all(required_cols %in% colnames(cell_scores))) {
  stop("Cell score table must contain: ", paste(required_cols, collapse = ", "), call. = FALSE)
}
if (nrow(graph) != ncol(graph) || is.null(rownames(graph)) || is.null(colnames(graph))) {
  stop("KNN graph must be a square matrix with rownames and colnames", call. = FALSE)
}
if (!identical(rownames(graph), colnames(graph))) {
  stop("KNN graph rownames and colnames must be identical", call. = FALSE)
}

common_cells <- intersect(rownames(graph), unique(cell_scores$cell))
if (length(common_cells) == 0) {
  stop("No common cells between graph and gene-risk cell score table", call. = FALSE)
}
cell_ids <- rownames(graph)[rownames(graph) %in% common_cells]
graph <- graph[cell_ids, cell_ids]
degree <- Matrix::rowSums(graph)
graph_cell_ids <- names(degree)[degree > 0]
isolated_cells <- names(degree)[degree == 0]
graph_use <- graph[graph_cell_ids, graph_cell_ids]

cat("  graph cells      :", nrow(graph), "\n")
cat("  graph nonisolated:", length(graph_cell_ids), "\n")
cat("  graph isolated   :", length(isolated_cells), "\n")

run_summary_list <- list()
seed_enrichment_list <- list()
top_enrichment_list <- list()
cell_result_list <- list()

for (tr in traits) {
  cat("\n[Trait]", tr, "\n")
  dt <- cell_scores[trait == tr & cell %in% cell_ids]
  if (nrow(dt) == 0) {
    warning("No cell score rows for trait: ", tr)
    next
  }
  setkey(dt, cell)
  dt <- dt[cell_ids]
  stopifnot(identical(dt$cell, cell_ids))

  target_vec <- as.logical(dt$target_type)
  all_target <- mean(target_vec, na.rm = TRUE)
  z_score <- setNames(as.numeric(dt$atac_z_score), dt$cell)
  scale_factor <- safe_scale_factor(z_score)

  for (alpha in alpha_grid) {
    label <- alpha_label(alpha)
    mixed_seed_score <- (1 - alpha) * as.numeric(dt$atac_z_rank01) +
      alpha * as.numeric(dt$gene_risk_cell_rank01)
    names(mixed_seed_score) <- dt$cell

    seed_idx <- top_n_idx(mixed_seed_score, seed_percent)
    seed_cells_total <- names(mixed_seed_score)[seed_idx]
    seed_cells <- intersect(seed_cells_total, graph_cell_ids)
    seed_isolated <- setdiff(seed_cells_total, seed_cells)
    if (length(seed_cells) == 0) {
      warning("No non-isolated seed cells for trait ", tr, " alpha ", alpha)
      next
    }

    np_score <- SCAVENGE::randomWalk_sparse(
      intM = graph_use,
      queryCells = seed_cells,
      gamma = gamma
    )
    np_full <- setNames(rep(0, length(cell_ids)), cell_ids)
    np_full[names(np_score)] <- as.numeric(np_score)

    keep <- is.finite(np_score) & np_score > 0
    trs_full <- setNames(rep(0, length(cell_ids)), cell_ids)
    if (any(keep)) {
      trs <- SCAVENGE::capOutlierQuantile(x = np_score[keep], q_ceiling = 0.95)
      trs <- SCAVENGE::max_min_scale(trs) * scale_factor
      trs_full[names(trs)] <- as.numeric(trs)
    }

    seed_logical <- dt$cell %in% seed_cells_total
    seed_target <- mean(target_vec[seed_logical], na.rm = TRUE)
    seed_enrichment_list[[paste(tr, label, sep = "_")]] <- data.table(
      trait = tr,
      alpha = alpha,
      graph = graph_label,
      n_seed_total = length(seed_cells_total),
      n_seed_in_graph = length(seed_cells),
      n_seed_isolated = length(seed_isolated),
      target_fraction_seed = seed_target,
      target_fraction_all = all_target,
      seed_enrichment_ratio = if (is.finite(all_target) && all_target > 0) seed_target / all_target else NA_real_,
      target_definition = target_definition(tr)
    )

    for (score_col in c("TRS", "mixed_seed_score", "atac_z_score", "gene_risk_cell_score")) {
      score_vec <- switch(
        score_col,
        TRS = trs_full,
        mixed_seed_score = mixed_seed_score,
        atac_z_score = setNames(dt$atac_z_score, dt$cell),
        gene_risk_cell_score = setNames(dt$gene_risk_cell_score, dt$cell)
      )
      for (top_frac in c(0.05, 0.10)) {
        top_idx <- top_n_idx(score_vec, top_frac)
        top_target <- mean(target_vec[top_idx], na.rm = TRUE)
        top_enrichment_list[[paste(tr, label, score_col, top_frac, sep = "_")]] <- data.table(
          trait = tr,
          alpha = alpha,
          graph = graph_label,
          score = score_col,
          top_frac = top_frac,
          n_top = length(top_idx),
          target_fraction_top = top_target,
          target_fraction_all = all_target,
          enrichment_ratio = if (is.finite(all_target) && all_target > 0) top_target / all_target else NA_real_,
          target_definition = target_definition(tr)
        )
      }
    }

    run_summary_list[[paste(tr, label, sep = "_")]] <- data.table(
      trait = tr,
      alpha = alpha,
      graph = graph_label,
      n_cells = length(cell_ids),
      n_graph_cells = length(graph_cell_ids),
      n_isolated_cells = length(isolated_cells),
      n_seed_total = length(seed_cells_total),
      n_seed_in_graph = length(seed_cells),
      n_seed_isolated = length(seed_isolated),
      n_np_positive = sum(np_full > 0),
      pct_np_positive = mean(np_full > 0),
      scale_factor = scale_factor,
      mean_TRS = mean(trs_full),
      median_TRS = median(trs_full),
      max_TRS = max(trs_full)
    )

    cell_result_list[[paste(tr, label, sep = "_")]] <- data.table(
      trait = tr,
      alpha = alpha,
      graph = graph_label,
      cell = dt$cell,
      mixed_seed_score = as.numeric(mixed_seed_score),
      seed = seed_logical,
      np_score = as.numeric(np_full),
      TRS = as.numeric(trs_full),
      atac_z_score = as.numeric(dt$atac_z_score),
      gene_risk_cell_score = as.numeric(dt$gene_risk_cell_score),
      BioClassification = dt$BioClassification,
      target_type = target_vec
    )

    cat("  alpha", alpha, "seed target", round(seed_target, 4), "TRS positive", sum(np_full > 0), "\n")
  }
}

if (length(run_summary_list) == 0) {
  stop("No alpha scan results were generated", call. = FALSE)
}

run_summary <- rbindlist(run_summary_list, use.names = TRUE, fill = TRUE)
seed_enrichment <- rbindlist(seed_enrichment_list, use.names = TRUE, fill = TRUE)
top_enrichment <- rbindlist(top_enrichment_list, use.names = TRUE, fill = TRUE)
cell_results <- rbindlist(cell_result_list, use.names = TRUE, fill = TRUE)

top_baseline <- top_enrichment[score == "TRS" & alpha == 0, .(
  trait,
  top_frac,
  baseline_target_fraction_top = target_fraction_top
)]
top_enrichment <- merge(top_enrichment, top_baseline, by = c("trait", "top_frac"), all.x = TRUE)
top_enrichment[, delta_vs_alpha0_pp := 100 * (target_fraction_top - baseline_target_fraction_top)]

seed_baseline <- seed_enrichment[alpha == 0, .(
  trait,
  baseline_target_fraction_seed = target_fraction_seed
)]
seed_enrichment <- merge(seed_enrichment, seed_baseline, by = "trait", all.x = TRUE)
seed_enrichment[, delta_seed_vs_alpha0_pp := 100 * (target_fraction_seed - baseline_target_fraction_seed)]

run_summary_out <- file.path(out_dir, "seed_prior_alpha_scan_run_summary.csv")
seed_enrichment_out <- file.path(out_dir, "seed_prior_alpha_scan_seed_enrichment.csv")
top_enrichment_out <- file.path(out_dir, "seed_prior_alpha_scan_top_enrichment.csv")
cell_results_out <- file.path(out_dir, "seed_prior_alpha_scan_cell_results.csv.gz")
data.table::fwrite(run_summary, run_summary_out)
data.table::fwrite(seed_enrichment, seed_enrichment_out)
data.table::fwrite(top_enrichment, top_enrichment_out)
data.table::fwrite(cell_results, cell_results_out)

plot_dt <- top_enrichment[score == "TRS" & top_frac == 0.05]
p_top <- ggplot(plot_dt, aes(x = factor(alpha), y = target_fraction_top, group = trait, color = trait)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_wrap(~ trait, scales = "free_y") +
  labs(
    title = "Top 5% TRS target-cell fraction across seed-prior alpha",
    subtitle = sprintf("Graph: %s; alpha=0 is ATAC-z-only seed prior", graph_label),
    x = "alpha for gene-risk seed prior",
    y = "Target-cell fraction in top 5% TRS"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")
fig_top <- file.path(out_dir, "Fig_seed_prior_alpha_scan_top5_TRS.png")
ggsave(fig_top, p_top, width = 8, height = 5, dpi = 300)

seed_plot_dt <- seed_enrichment
p_seed <- ggplot(seed_plot_dt, aes(x = factor(alpha), y = target_fraction_seed, group = trait, color = trait)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_wrap(~ trait, scales = "free_y") +
  labs(
    title = "Seed target-cell fraction across alpha",
    subtitle = sprintf("Graph: %s; seed selection uses mixed rank score", graph_label),
    x = "alpha for gene-risk seed prior",
    y = "Target-cell fraction among seed cells"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")
fig_seed <- file.path(out_dir, "Fig_seed_prior_alpha_scan_seed_fraction.png")
ggsave(fig_seed, p_seed, width = 8, height = 5, dpi = 300)

manifest <- list(
  script = "05_run_seed_prior_alpha_scan.R",
  traits = traits,
  alpha_grid = alpha_grid,
  seed_percent = seed_percent,
  gamma = gamma,
  graph = graph_path,
  graph_label = graph_label,
  cell_scores = cell_score_path,
  run_summary = run_summary_out,
  seed_enrichment = seed_enrichment_out,
  top_enrichment = top_enrichment_out,
  cell_results = cell_results_out,
  top5_figure = fig_top,
  seed_figure = fig_seed,
  created_at = as.character(Sys.time())
)
jsonlite::write_json(
  manifest,
  file.path(out_dir, "seed_prior_alpha_scan_manifest.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

cat("\n========== Done ==========\n")
cat("  run summary    :", run_summary_out, "\n")
cat("  seed enrichment:", seed_enrichment_out, "\n")
cat("  top enrichment :", top_enrichment_out, "\n")
cat("  cell results   :", cell_results_out, "\n")
