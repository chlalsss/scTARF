# Compare graph-level quality across KNN strategies.
#
# This script checks whether KNN graphs are technically usable for SCAVENGE:
#   - square sparse graph with matched ATAC cell names
#   - degree distribution
#   - isolated cells
#   - same-celltype neighbor fraction
#   - edge overlap between strategies
#
# It does not run TRS.

source("D:/RStudio/work/module2_wnn_lite_fusion/_bootstrap.R")

required_pkgs(c("data.table", "Matrix", "SummarizedExperiment", "ggplot2", "jsonlite"))
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(SummarizedExperiment)
  library(ggplot2)
  library(jsonlite)
})

paths <- wnn_lite_paths()
ensure_wnn_lite_dirs(paths)

cutoff <- get_single_cutoff()
label <- cutoff_label(cutoff)
k <- env_integer("WNN_LITE_KNN_K", env_integer("MODULE2_KNN_K", 30))
npcs <- env_integer("MODULE2_FUSED_NPCS", 50)

se_path <- Sys.getenv(
  "MODULE2_ATAC_SE",
  file.path(paths$data_dir, "scATAC-Healthy-Hematopoiesis-191120.rds")
)
original_knn_path <- Sys.getenv(
  "MODULE2_ORIGINAL_KNN",
  file.path(paths$data_dir, "scavenge_out_small_all_traits", "mutualknn_k30.rds")
)
fused_knn_path <- Sys.getenv(
  "MODULE2_FUSED_KNN",
  file.path(
    paths$data_dir,
    "module2_adaptive_fusion",
    sprintf("fused_mutualknn_cutoff_%s_pca%d_k%d.rds", label, npcs, k)
  )
)
wnn_knn_path <- Sys.getenv(
  "WNN_LITE_KNN_PATH",
  file.path(paths$out_dir, sprintf("wnn_lite_mutualknn_cutoff_%s_k%d.rds", label, k))
)

out_dir <- file.path(paths$out_dir, sprintf("knn_strategy_comparison_cutoff_%s_k%d", label, k))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("========== Module 2: compare all KNN strategies ==========\n")
cat("  original:", original_knn_path, "\n")
cat("  fused   :", fused_knn_path, "\n")
cat("  WNN-lite:", wnn_knn_path, "\n")
cat("  ATAC SE :", se_path, "\n")
cat("  output  :", out_dir, "\n")

stopifnot(file.exists(se_path), file.exists(original_knn_path), file.exists(fused_knn_path), file.exists(wnn_knn_path))

SE_atac <- readRDS(se_path)
meta <- as.data.table(as.data.frame(colData(SE_atac)))
meta[, cell := colnames(SE_atac)]
if (!"BioClassification" %in% colnames(meta)) {
  stop("ATAC metadata must contain BioClassification", call. = FALSE)
}
type_map <- setNames(as.character(meta$BioClassification), meta$cell)
se_cells <- colnames(SE_atac)

graphs <- list(
  original_lsi = readRDS(original_knn_path),
  concat_pca_fused = readRDS(fused_knn_path),
  wnn_lite = readRDS(wnn_knn_path)
)

validate_graph <- function(adj, name) {
  if (nrow(adj) != ncol(adj)) stop(name, " is not square", call. = FALSE)
  if (is.null(rownames(adj)) || is.null(colnames(adj))) {
    stop(name, " must have rownames and colnames", call. = FALSE)
  }
  if (!identical(rownames(adj), colnames(adj))) {
    stop(name, " rownames and colnames are not identical", call. = FALSE)
  }
  common <- intersect(se_cells, rownames(adj))
  data.table(
    method = name,
    n_graph_cells = nrow(adj),
    n_se_cells = length(se_cells),
    n_common_cells = length(common),
    same_order_as_se = identical(se_cells, rownames(adj))
  )
}

extract_edges <- function(adj, cells) {
  upper <- as(Matrix::triu(adj, k = 1), "TsparseMatrix")
  paste(cells[upper@i + 1L], cells[upper@j + 1L], sep = "|")
}

summary_graph <- function(adj, name) {
  degree <- Matrix::rowSums(adj)
  cells <- rownames(adj)
  adj_t <- as(adj, "TsparseMatrix")
  i <- adj_t@i + 1L
  j <- adj_t@j + 1L
  keep <- i != j
  i <- i[keep]
  j <- j[keep]
  same <- type_map[cells[i]] == type_map[cells[j]]
  per_cell_same <- tapply(same, i, mean)

  data.table(
    method = name,
    n_cells = nrow(adj),
    n_edges = length(extract_edges(adj, cells)),
    min_degree = min(degree),
    q10_degree = as.numeric(stats::quantile(degree, 0.10)),
    q25_degree = as.numeric(stats::quantile(degree, 0.25)),
    median_degree = as.numeric(stats::median(degree)),
    mean_degree = mean(degree),
    q75_degree = as.numeric(stats::quantile(degree, 0.75)),
    q90_degree = as.numeric(stats::quantile(degree, 0.90)),
    max_degree = max(degree),
    n_isolated = sum(degree == 0),
    pct_isolated = mean(degree == 0),
    mean_same_type_frac = mean(per_cell_same, na.rm = TRUE),
    median_same_type_frac = median(per_cell_same, na.rm = TRUE),
    usable_for_scavenge = length(cells) == length(se_cells) &&
      length(intersect(cells, se_cells)) == length(se_cells) &&
      mean(degree == 0) <= 0.01
  )
}

edge_overlap <- function(a, b, name_a, name_b) {
  common <- intersect(rownames(a), rownames(b))
  a <- a[common, common]
  b <- b[common, common]
  ea <- extract_edges(a, common)
  eb <- extract_edges(b, common)
  n_inter <- length(intersect(ea, eb))
  n_union <- length(ea) + length(eb) - n_inter
  data.table(
    graph_a = name_a,
    graph_b = name_b,
    n_common_cells = length(common),
    n_edges_a = length(ea),
    n_edges_b = length(eb),
    n_edges_common = n_inter,
    jaccard = if (n_union > 0) n_inter / n_union else NA_real_,
    pct_a_overlap_b = if (length(ea) > 0) n_inter / length(ea) else NA_real_,
    pct_b_overlap_a = if (length(eb) > 0) n_inter / length(eb) else NA_real_
  )
}

alignment_dt <- rbindlist(Map(validate_graph, graphs, names(graphs)))
summary_dt <- rbindlist(Map(summary_graph, graphs, names(graphs)))

pairs <- combn(names(graphs), 2, simplify = FALSE)
overlap_dt <- rbindlist(lapply(pairs, function(pair) {
  edge_overlap(graphs[[pair[1]]], graphs[[pair[2]]], pair[1], pair[2])
}))

alignment_out <- file.path(out_dir, "knn_strategy_alignment_summary.csv")
summary_out <- file.path(out_dir, "knn_strategy_graph_summary.csv")
overlap_out <- file.path(out_dir, "knn_strategy_edge_overlap.csv")
fwrite(alignment_dt, alignment_out)
fwrite(summary_dt, summary_out)
fwrite(overlap_dt, overlap_out)

degree_dt <- rbindlist(lapply(names(graphs), function(name) {
  data.table(method = name, degree = as.numeric(Matrix::rowSums(graphs[[name]])))
}))
degree_out <- file.path(out_dir, "knn_strategy_degree_by_cell.csv")
fwrite(degree_dt, degree_out)

p_degree <- ggplot(degree_dt, aes(x = degree, fill = method)) +
  geom_histogram(binwidth = 1, color = "white", linewidth = 0.1) +
  facet_wrap(~ method, ncol = 1) +
  labs(
    title = "KNN strategy degree distributions",
    x = "MutualKNN degree",
    y = "Cell count"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")
fig_degree <- file.path(out_dir, "Fig_knn_strategy_degree_distribution.png")
ggsave(fig_degree, p_degree, width = 7, height = 8, dpi = 300)

jsonlite::write_json(
  list(
    script = "05_compare_all_knn_strategies.R",
    cutoff = cutoff,
    k = k,
    npcs = npcs,
    original_knn = original_knn_path,
    fused_knn = fused_knn_path,
    wnn_lite_knn = wnn_knn_path,
    alignment_csv = alignment_out,
    summary_csv = summary_out,
    edge_overlap_csv = overlap_out,
    degree_csv = degree_out,
    degree_figure = fig_degree,
    created_at = as.character(Sys.time())
  ),
  file.path(out_dir, "knn_strategy_comparison_manifest.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

cat("\n[Alignment]\n")
print(alignment_dt)
cat("\n[Graph summary]\n")
print(summary_dt)
cat("\n[Edge overlap]\n")
print(overlap_dt)
cat("\n========== Done ==========\n")
cat("  summary:", summary_out, "\n")
cat("  overlap:", overlap_out, "\n")

