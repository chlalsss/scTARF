source("D:/RStudio/work/module2_wnn_lite_gated/_bootstrap.R")

required_pkgs(c("data.table", "jsonlite"))
suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

paths <- gated_paths("compare_gated_vs_ungated")
ensure_gated_dirs(paths)

cutoff <- get_single_cutoff()
k <- env_integer("WNN_LITE_KNN_K", env_integer("MODULE2_KNN_K", 30))
label <- cutoff_label(cutoff)

cell_paths <- gated_paths("celltype_gated")
rel_paths <- gated_paths("reliability_gated")

out_dir <- paths$out_dir
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

files <- list(
  ungated_diag = file.path(paths$ungated_dir, sprintf("wnn_lite_knn_diagnostics_cutoff_%s_k%d.csv", label, k)),
  cell_diag = file.path(cell_paths$out_dir, sprintf("wnn_lite_celltype_gated_knn_diagnostics_k%d.csv", k)),
  rel_diag = file.path(rel_paths$out_dir, sprintf("wnn_lite_reliability_gated_knn_diagnostics_k%d.csv", k)),
  ungated_same = file.path(paths$ungated_dir, sprintf("wnn_lite_same_type_by_celltype_cutoff_%s_k%d.csv", label, k)),
  cell_same = file.path(cell_paths$out_dir, sprintf("wnn_lite_celltype_gated_same_type_by_celltype_k%d.csv", k)),
  rel_same = file.path(rel_paths$out_dir, sprintf("wnn_lite_reliability_gated_same_type_by_celltype_k%d.csv", k)),
  cell_gated_types = file.path(cell_paths$out_dir, "gated_celltypes.csv"),
  cell_hop = file.path(cell_paths$out_dir, "seed_hop_diagnostics", "celltype_gated_seed_graph_hop_delta.csv"),
  rel_hop = file.path(rel_paths$out_dir, "seed_hop_diagnostics", "reliability_gated_seed_graph_hop_delta.csv")
)

missing <- names(files)[!file.exists(unlist(files))]
if (length(missing) > 0) {
  stop("Missing required comparison inputs: ", paste(missing, collapse = ", "), call. = FALSE)
}

read_diag <- function(path, graph) {
  dt <- fread(path)
  dt[, graph := graph]
  keep <- c(
    "graph", "n_cells", "n_gated", "pct_gated", "n_mutual_pairs",
    "median_degree", "mean_degree", "n_isolated", "pct_isolated",
    "mean_rna_weight", "mean_selected_rna_component_frac",
    "mean_selected_rna_component_frac_gated",
    "mean_selected_rna_component_frac_non_gated"
  )
  dt[, intersect(keep, names(dt)), with = FALSE]
}

graph_diag <- rbindlist(list(
  read_diag(files$ungated_diag, "ungated_wnn_lite"),
  read_diag(files$cell_diag, "celltype_gated"),
  read_diag(files$rel_diag, "reliability_gated")
), fill = TRUE)

graph_diag[, isolated_pct_display := sprintf("%.2f%%", pct_isolated * 100)]
graph_diag[, selected_rna_component_display := sprintf("%.2f%%", mean_selected_rna_component_frac * 100)]
if ("mean_selected_rna_component_frac_gated" %in% names(graph_diag)) {
  graph_diag[, gated_rna_component_display := fifelse(
    is.na(mean_selected_rna_component_frac_gated),
    NA_character_,
    sprintf("%.2f%%", mean_selected_rna_component_frac_gated * 100)
  )]
}

graph_diag_out <- file.path(out_dir, "graph_diagnostics_summary.csv")
fwrite(graph_diag, graph_diag_out)

base_same <- fread(files$ungated_same)
cell_same <- fread(files$cell_same)
rel_same <- fread(files$rel_same)
gated_types <- fread(files$cell_gated_types)

base_same <- base_same[, .(
  cell_type,
  n_cells_base = n_cells,
  base_degree = mean_degree,
  base_mutual = mean_mutual_same_type_frac
)]
cell_same <- cell_same[, .(
  cell_type,
  n_cells,
  n_gated,
  pct_gated,
  cell_degree = mean_degree,
  cell_mutual = mean_mutual_same_type_frac
)]
rel_same <- rel_same[, .(
  cell_type,
  rel_n_gated = n_gated,
  rel_pct_gated = pct_gated,
  rel_degree = mean_degree,
  rel_mutual = mean_mutual_same_type_frac
)]

same_join <- Reduce(function(x, y) merge(x, y, by = "cell_type", all = TRUE), list(base_same, cell_same, rel_same))
same_join[, is_celltype_gated_type := cell_type %in% gated_types$cell_type]
same_join[, cell_mutual_delta_pp := (cell_mutual - base_mutual) * 100]
same_join[, rel_mutual_delta_pp := (rel_mutual - base_mutual) * 100]
same_join[, cell_degree_delta := cell_degree - base_degree]
same_join[, rel_degree_delta := rel_degree - base_degree]

same_type_detail_out <- file.path(out_dir, "same_type_by_celltype_deltas.csv")
fwrite(same_join, same_type_detail_out)

weighted_same <- same_join[is_celltype_gated_type == TRUE, .(
  graph = c("ungated_wnn_lite", "celltype_gated", "reliability_gated"),
  weighted_mutual_same_type_fraction = c(
    sum(base_mutual * n_cells_base, na.rm = TRUE) / sum(n_cells_base, na.rm = TRUE),
    sum(cell_mutual * n_cells_base, na.rm = TRUE) / sum(n_cells_base, na.rm = TRUE),
    sum(rel_mutual * n_cells_base, na.rm = TRUE) / sum(n_cells_base, na.rm = TRUE)
  ),
  weighted_mean_degree = c(
    sum(base_degree * n_cells_base, na.rm = TRUE) / sum(n_cells_base, na.rm = TRUE),
    sum(cell_degree * n_cells_base, na.rm = TRUE) / sum(n_cells_base, na.rm = TRUE),
    sum(rel_degree * n_cells_base, na.rm = TRUE) / sum(n_cells_base, na.rm = TRUE)
  )
)]
baseline_same <- weighted_same[graph == "ungated_wnn_lite", weighted_mutual_same_type_fraction]
weighted_same[, delta_vs_ungated_pp := (weighted_mutual_same_type_fraction - baseline_same) * 100]
weighted_same[, acceptance_target_delta_pp := 5]
weighted_same[, passes_same_type_criterion := graph != "ungated_wnn_lite" & delta_vs_ungated_pp >= acceptance_target_delta_pp]
same_type_summary_out <- file.path(out_dir, "same_type_weighted_summary.csv")
fwrite(weighted_same, same_type_summary_out)

read_hop <- function(path, graph) {
  dt <- fread(path)
  dt[, gated_graph := graph]
  dt
}
hop_delta <- rbindlist(list(
  read_hop(files$cell_hop, "celltype_gated"),
  read_hop(files$rel_hop, "reliability_gated")
), fill = TRUE)
hop_delta[, within_0p5pp_of_original_or_better := gated_minus_original_pp >= -0.5]
hop_delta_out <- file.path(out_dir, "seed_hop_delta_summary.csv")
fwrite(hop_delta, hop_delta_out)

hop2_summary <- hop_delta[neighborhood == "hop2_exact", .(
  n_traits = .N,
  n_traits_within_0p5pp = sum(within_0p5pp_of_original_or_better, na.rm = TRUE),
  passes_hop2_criterion = sum(within_0p5pp_of_original_or_better, na.rm = TRUE) >= 2,
  mean_gated_minus_original_pp = mean(gated_minus_original_pp, na.rm = TRUE),
  mean_gated_minus_ungated_pp = mean(gated_minus_ungated_pp, na.rm = TRUE)
), by = gated_graph]
hop2_summary_out <- file.path(out_dir, "hop2_acceptance_summary.csv")
fwrite(hop2_summary, hop2_summary_out)

acceptance <- data.table(
  gated_graph = c("celltype_gated", "reliability_gated"),
  sanity_bound_pass = c(
    graph_diag[graph == "celltype_gated", pct_gated <= 0.50],
    graph_diag[graph == "reliability_gated", pct_gated <= 0.40]
  ),
  same_type_pass = c(
    weighted_same[graph == "celltype_gated", passes_same_type_criterion],
    weighted_same[graph == "reliability_gated", passes_same_type_criterion]
  ),
  hop2_pass = c(
    hop2_summary[gated_graph == "celltype_gated", passes_hop2_criterion],
    hop2_summary[gated_graph == "reliability_gated", passes_hop2_criterion]
  )
)
acceptance[, overall_pre_trs_pass := sanity_bound_pass & same_type_pass & hop2_pass]
acceptance[, recommendation := fifelse(
  overall_pre_trs_pass,
  "eligible_for_full_TRS",
  "do_not_prioritize_full_TRS; record as negative/pre-TRS failed"
)]
acceptance_out <- file.path(out_dir, "pre_trs_acceptance_summary.csv")
fwrite(acceptance, acceptance_out)

manifest <- list(
  script = "06_compare_all_graphs.R",
  output_dir = out_dir,
  graph_diagnostics_summary = graph_diag_out,
  same_type_by_celltype_deltas = same_type_detail_out,
  same_type_weighted_summary = same_type_summary_out,
  seed_hop_delta_summary = hop_delta_out,
  hop2_acceptance_summary = hop2_summary_out,
  pre_trs_acceptance_summary = acceptance_out,
  note = "This comparison intentionally summarizes pre-TRS graph and hop diagnostics. TRS top5 is omitted until 04_run_scavenge_*_gated.R is run.",
  created_at = as.character(Sys.time())
)
jsonlite::write_json(
  manifest,
  file.path(out_dir, "compare_gated_vs_ungated_manifest.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

cat("========== WNN-lite gated: compare all pre-TRS diagnostics ==========\n")
cat("\n[Graph diagnostics]\n")
print(graph_diag)
cat("\n[Same-type weighted summary]\n")
print(weighted_same)
cat("\n[Hop2 acceptance summary]\n")
print(hop2_summary)
cat("\n[Pre-TRS acceptance summary]\n")
print(acceptance)
cat("\n========== Done ==========\n")
cat("  graph diagnostics:", graph_diag_out, "\n")
cat("  same type summary:", same_type_summary_out, "\n")
cat("  hop2 summary     :", hop2_summary_out, "\n")
cat("  acceptance       :", acceptance_out, "\n")
