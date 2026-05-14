# Common seed-hop diagnostics for gated WNN-lite graphs.
#
# Entry scripts must define:
#   gated_version
#   gated_label
#   gated_knn_file

if (!exists("wnn_gated_code_dir", inherits = TRUE)) {
  source("D:/RStudio/work/module2_wnn_lite_gated/_bootstrap.R")
}

required_pkgs(c("data.table", "Matrix", "ggplot2", "jsonlite"))
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(ggplot2)
  library(jsonlite)
})

if (!exists("gated_version", inherits = TRUE)) {
  stop("Entry script must define gated_version", call. = FALSE)
}
if (!exists("gated_label", inherits = TRUE)) {
  stop("Entry script must define gated_label", call. = FALSE)
}
if (!exists("gated_knn_file", inherits = TRUE)) {
  stop("Entry script must define gated_knn_file", call. = FALSE)
}

paths <- gated_paths(gated_version)
ensure_gated_dirs(paths)

cutoff <- get_single_cutoff()
k <- env_integer("WNN_LITE_KNN_K", env_integer("MODULE2_KNN_K", 30))
traits <- strsplit(Sys.getenv("MODULE2_TRAITS", "MONO,PLT,LYM"), ",", fixed = TRUE)[[1]]
traits <- trimws(traits)
traits <- traits[nzchar(traits)]
if (length(traits) == 0) stop("MODULE2_TRAITS is empty", call. = FALSE)

ungated_knn_file <- file.path(
  paths$ungated_dir,
  sprintf("wnn_lite_mutualknn_cutoff_%s_k%d.rds", cutoff_label(cutoff), k)
)
diag_dir <- file.path(paths$out_dir, "seed_hop_diagnostics")
dir.create(diag_dir, showWarnings = FALSE, recursive = TRUE)

stopifnot(
  file.exists(paths$seed_signal_by_cell),
  file.exists(paths$original_knn),
  file.exists(ungated_knn_file),
  file.exists(gated_knn_file)
)

cat("========== WNN-lite gated: seed-hop diagnostics ==========\n")
cat("  version       :", gated_version, "\n")
cat("  traits        :", paste(traits, collapse = ", "), "\n")
cat("  seed signal   :", paths$seed_signal_by_cell, "\n")
cat("  original KNN  :", paths$original_knn, "\n")
cat("  ungated KNN   :", ungated_knn_file, "\n")
cat("  gated KNN     :", gated_knn_file, "\n")
cat("  output dir    :", diag_dir, "\n")

seed_signal <- fread(paths$seed_signal_by_cell)
needed_seed_cols <- c("trait", "cell", "seed", "target_type", "BioClassification")
if (!all(needed_seed_cols %in% names(seed_signal))) {
  stop("seed_signal_by_cell must contain: ", paste(needed_seed_cols, collapse = ", "), call. = FALSE)
}
seed_signal[, seed := as.logical(seed)]
seed_signal[, target_type := as.logical(target_type)]

graphs <- list(
  original_lsi = readRDS(paths$original_knn),
  ungated_wnn_lite = readRDS(ungated_knn_file),
  gated_wnn_lite = readRDS(gated_knn_file)
)

validate_graph <- function(adj, graph_name) {
  if (nrow(adj) != ncol(adj)) stop(graph_name, " is not square", call. = FALSE)
  if (is.null(rownames(adj)) || is.null(colnames(adj))) {
    stop(graph_name, " must have rownames and colnames", call. = FALSE)
  }
  if (!identical(rownames(adj), colnames(adj))) {
    stop(graph_name, " rownames and colnames are not identical", call. = FALSE)
  }
  data.table(
    graph = graph_name,
    n_cells = nrow(adj),
    n_edges = length(adj@x) / 2,
    n_isolated = sum(Matrix::rowSums(adj) == 0),
    pct_isolated = mean(Matrix::rowSums(adj) == 0)
  )
}

alignment <- rbindlist(Map(validate_graph, graphs, names(graphs)))
fwrite(alignment, file.path(diag_dir, sprintf("%s_graph_alignment_summary.csv", gated_label)))
print(alignment)

common_cells <- Reduce(intersect, c(lapply(graphs, rownames), list(unique(seed_signal$cell))))
if (length(common_cells) == 0) {
  stop("No common cells across seed signal and graphs", call. = FALSE)
}

cell_order <- rownames(graphs$gated_wnn_lite)
cell_order <- cell_order[cell_order %in% common_cells]
if (length(cell_order) < length(common_cells)) {
  warning("Using gated graph cell order for common cells: ", length(cell_order), " / ", length(common_cells))
}

graphs <- lapply(graphs, function(adj) as_binary_graph(adj[cell_order, cell_order]))

hop_enrichment <- function(adj, seed_cells, target_vec, graph_name, trait) {
  cells <- rownames(adj)
  target_vec <- as.logical(target_vec[cells])
  seed_cells <- intersect(seed_cells, cells)
  seed_idx <- match(seed_cells, cells)
  seed_logical <- rep(FALSE, length(cells))
  seed_logical[seed_idx] <- TRUE

  degree <- Matrix::rowSums(adj)
  seed_in_nonisolated <- seed_logical & degree > 0

  if (length(seed_idx) == 0) {
    groups <- list(
      seed_self = seed_logical,
      hop1_exact = rep(FALSE, length(cells)),
      hop2_exact = rep(FALSE, length(cells)),
      within1_including_seed = seed_logical,
      within2_including_seed = seed_logical
    )
  } else {
    hop1 <- Matrix::rowSums(adj[, seed_idx, drop = FALSE]) > 0
    hop1_exact <- hop1 & !seed_logical

    if (any(hop1_exact)) {
      hop1_idx <- which(hop1_exact)
      hop2 <- Matrix::rowSums(adj[, hop1_idx, drop = FALSE]) > 0
    } else {
      hop2 <- rep(FALSE, length(cells))
    }
    hop2_exact <- hop2 & !seed_logical & !hop1_exact

    groups <- list(
      seed_self = seed_logical,
      hop1_exact = hop1_exact,
      hop2_exact = hop2_exact,
      within1_including_seed = seed_logical | hop1_exact,
      within2_including_seed = seed_logical | hop1_exact | hop2_exact
    )
  }

  all_target <- mean(target_vec, na.rm = TRUE)
  rbindlist(lapply(names(groups), function(group_name) {
    idx <- groups[[group_name]]
    frac <- if (any(idx)) mean(target_vec[idx], na.rm = TRUE) else NA_real_
    data.table(
      trait = trait,
      graph = graph_name,
      neighborhood = group_name,
      n_cells = sum(idx),
      n_target_cells = if (any(idx)) sum(target_vec[idx], na.rm = TRUE) else 0L,
      target_fraction = frac,
      target_fraction_all = all_target,
      enrichment_ratio = if (is.finite(all_target) && all_target > 0) frac / all_target else NA_real_,
      n_seed_total = length(seed_cells),
      n_seed_nonisolated = sum(seed_in_nonisolated),
      n_seed_isolated = sum(seed_logical & degree == 0),
      target_definition = target_definition(trait)
    )
  }))
}

hop_list <- list()
for (tr in traits) {
  dt <- seed_signal[trait == tr & cell %in% cell_order]
  if (nrow(dt) == 0) {
    warning("No seed signal rows for trait: ", tr)
    next
  }
  dt <- dt[match(cell_order, cell)]
  if (!identical(dt$cell, cell_order)) {
    stop("Seed signal alignment failed for trait: ", tr, call. = FALSE)
  }
  seed_cells <- dt[seed == TRUE, cell]
  target_vec <- setNames(dt$target_type, dt$cell)
  for (graph_name in names(graphs)) {
    hop_list[[paste(tr, graph_name, sep = "_")]] <- hop_enrichment(
      adj = graphs[[graph_name]],
      seed_cells = seed_cells,
      target_vec = target_vec,
      graph_name = graph_name,
      trait = tr
    )
  }
}

hop_dt <- rbindlist(hop_list, use.names = TRUE, fill = TRUE)
hop_out <- file.path(diag_dir, sprintf("%s_seed_graph_hop_enrichment.csv", gated_label))
fwrite(hop_dt, hop_out)

delta_dt <- dcast(
  hop_dt,
  trait + neighborhood ~ graph,
  value.var = "target_fraction"
)
if (all(c("original_lsi", "ungated_wnn_lite", "gated_wnn_lite") %in% names(delta_dt))) {
  delta_dt[, gated_minus_original_pp := (gated_wnn_lite - original_lsi) * 100]
  delta_dt[, gated_minus_ungated_pp := (gated_wnn_lite - ungated_wnn_lite) * 100]
  delta_dt[, ungated_minus_original_pp := (ungated_wnn_lite - original_lsi) * 100]
}
delta_out <- file.path(diag_dir, sprintf("%s_seed_graph_hop_delta.csv", gated_label))
fwrite(delta_dt, delta_out)

plot_dt <- hop_dt[
  neighborhood %in% c("seed_self", "hop1_exact", "hop2_exact", "within1_including_seed", "within2_including_seed")
]
plot_dt[, neighborhood := factor(
  neighborhood,
  levels = c("seed_self", "hop1_exact", "hop2_exact", "within1_including_seed", "within2_including_seed")
)]

p <- ggplot(plot_dt, aes(x = neighborhood, y = target_fraction, fill = graph)) +
  geom_col(position = position_dodge(width = 0.76), width = 0.68) +
  geom_hline(aes(yintercept = target_fraction_all), linetype = "dashed", color = "grey35") +
  facet_wrap(~ trait, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c(
    original_lsi = "#8C8C8C",
    ungated_wnn_lite = "#4C78A8",
    gated_wnn_lite = "#F58518"
  )) +
  labs(
    title = sprintf("Seed-hop target enrichment: %s", gated_label),
    subtitle = "Dashed line is target-cell fraction among all ATAC cells",
    x = NULL,
    y = "Target-cell fraction",
    fill = "Graph"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )
fig_out <- file.path(diag_dir, sprintf("Fig_%s_seed_hop_target_enrichment.png", gated_label))
ggsave(fig_out, p, width = 10, height = 8, dpi = 300)

hop2_summary <- delta_dt[neighborhood == "hop2_exact"]
cat("\n[hop2_exact delta]\n")
print(hop2_summary)

manifest <- list(
  script = "_diagnostics_gated_common.R",
  gated_version = gated_version,
  gated_label = gated_label,
  traits = traits,
  seed_signal = paths$seed_signal_by_cell,
  original_knn = paths$original_knn,
  ungated_knn = ungated_knn_file,
  gated_knn = gated_knn_file,
  hop_enrichment_csv = hop_out,
  hop_delta_csv = delta_out,
  figure = fig_out,
  created_at = as.character(Sys.time())
)
jsonlite::write_json(
  manifest,
  file.path(diag_dir, sprintf("%s_seed_hop_diagnostics_manifest.json", gated_label)),
  pretty = TRUE,
  auto_unbox = TRUE
)

cat("\n========== Done ==========\n")
cat("  hop enrichment:", hop_out, "\n")
cat("  hop delta     :", delta_out, "\n")
cat("  figure        :", fig_out, "\n")
