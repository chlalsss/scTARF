# Diagnose whether SCAVENGE seeds or graph neighborhoods explain TRS changes.
#
# This script does not modify TRS. It recomputes the original per-trait
# chromVAR/gchromVAR z_score and seed cells, then compares target-cell
# enrichment at seed, 1-hop, and 2-hop levels on the original LSI graph and
# the WNN-lite graph.

source("D:/RStudio/work/module2_gwas_gene_risk/_bootstrap.R")

required_pkgs(c(
  "SCAVENGE", "chromVAR", "gchromVAR", "SummarizedExperiment",
  "data.table", "dplyr", "Matrix", "ggplot2", "jsonlite"
))
suppressPackageStartupMessages({
  library(SCAVENGE)
  library(chromVAR)
  library(gchromVAR)
  library(SummarizedExperiment)
  library(data.table)
  library(dplyr)
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

seed_percent <- env_numeric("MODULE2_SEED_PERCENT", 0.05)
if (seed_percent <= 0 || seed_percent >= 1) {
  stop("MODULE2_SEED_PERCENT must be in (0, 1)", call. = FALSE)
}

out_dir <- file.path(paths$out_dir, "seed_signal_diagnostics")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "trait_bed"), showWarnings = FALSE, recursive = TRUE)

cat("========== Module 2 GWAS gene-risk: seed signal diagnostics ==========\n")
cat("  traits       :", paste(traits, collapse = ", "), "\n")
cat("  seed percent :", seed_percent, "\n")
cat("  ATAC SE      :", paths$atac_se, "\n")
cat("  GWAS         :", paths$gwas, "\n")
cat("  bg peaks     :", paths$bg_peaks, "\n")
cat("  original KNN :", paths$original_knn, "\n")
cat("  WNN-lite KNN :", paths$wnn_lite_knn, "\n")
cat("  output dir   :", out_dir, "\n")

stopifnot(
  file.exists(paths$atac_se),
  file.exists(paths$gwas),
  file.exists(paths$bg_peaks),
  file.exists(paths$original_knn),
  file.exists(paths$wnn_lite_knn)
)

cat("\n========== Load inputs ==========\n")
SE_atac <- readRDS(paths$atac_se)
bg_peaks <- readRDS(paths$bg_peaks)
gwas <- data.table::fread(paths$gwas)
original_knn <- readRDS(paths$original_knn)
wnn_lite_knn <- readRDS(paths$wnn_lite_knn)

meta <- as.data.table(as.data.frame(SummarizedExperiment::colData(SE_atac)))
meta[, cell := colnames(SE_atac)]
if (!"BioClassification" %in% colnames(meta)) {
  stop("ATAC metadata must contain BioClassification", call. = FALSE)
}
meta <- meta[, .(cell, BioClassification)]

col_var <- "var"
col_trait <- "trait"
col_pp <- "PP"
if (!all(c(col_var, col_trait, col_pp) %in% names(gwas))) {
  stop("GWAS table must contain var, trait, and PP columns", call. = FALSE)
}

if (!all(c("chr", "start", "end") %in% names(gwas))) {
  parsed <- parse_var_to_bed(gwas[[col_var]])
  bad_n <- sum(!parsed$ok)
  if (bad_n > 0) {
    stop("Failed to parse ", bad_n, " GWAS variants from var column", call. = FALSE)
  }
  gwas[, `:=`(
    chr = parsed$chr,
    start = parsed$start,
    end = parsed$end
  )]
  gwas[end <= start, end := start + 1L]
}

validate_graph <- function(adj, graph_name, se_cells) {
  if (nrow(adj) != ncol(adj)) stop(graph_name, " is not square", call. = FALSE)
  if (is.null(rownames(adj)) || is.null(colnames(adj))) {
    stop(graph_name, " must have rownames and colnames", call. = FALSE)
  }
  if (!identical(rownames(adj), colnames(adj))) {
    stop(graph_name, " rownames and colnames are not identical", call. = FALSE)
  }

  common <- intersect(se_cells, rownames(adj))
  data.table(
    graph = graph_name,
    n_se_cells = length(se_cells),
    n_graph_cells = nrow(adj),
    n_common_cells = length(common),
    same_order_as_se = identical(se_cells, rownames(adj)),
    n_missing_from_graph = length(setdiff(se_cells, rownames(adj))),
    n_extra_in_graph = length(setdiff(rownames(adj), se_cells))
  )
}

se_cells <- colnames(SE_atac)
alignment <- rbindlist(list(
  validate_graph(original_knn, "original_lsi", se_cells),
  validate_graph(wnn_lite_knn, "wnn_lite", se_cells)
))
print(alignment)
fwrite(alignment, file.path(out_dir, "seed_signal_alignment_summary.csv"))

common_cells <- Reduce(intersect, list(se_cells, rownames(original_knn), rownames(wnn_lite_knn)))
if (length(common_cells) == 0) {
  stop("No common cells across ATAC SE, original KNN, and WNN-lite KNN", call. = FALSE)
}
if (length(common_cells) < length(se_cells)) {
  warning("Subsetting ATAC SE and graphs to ", length(common_cells), " common cells.")
}

cell_ids <- se_cells[se_cells %in% common_cells]
SE_use <- SE_atac[, cell_ids]
meta_use <- meta[match(cell_ids, meta$cell)]
stopifnot(identical(meta_use$cell, cell_ids))

graphs <- list(
  original_lsi = original_knn[cell_ids, cell_ids],
  wnn_lite = wnn_lite_knn[cell_ids, cell_ids]
)
graphs <- lapply(graphs, as_binary_graph)

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

seed_self_list <- list()
hop_list <- list()
cell_signal_list <- list()

cat("\n========== Run trait diagnostics ==========\n")
for (tr in traits) {
  cat("\n[Trait]", tr, "\n")
  dt <- gwas[get(col_trait) == tr]
  if (nrow(dt) == 0) {
    warning("No GWAS rows for trait: ", tr)
    next
  }

  bed <- dt[, .(
    chr = as.character(chr),
    start = as.integer(start),
    end = as.integer(end),
    name = tr,
    score = as.numeric(get(col_pp))
  )]
  bed[end <= start, end := start + 1L]

  bed_file <- file.path(out_dir, "trait_bed", paste0("trait_", make.names(tr), ".PP001.bed"))
  data.table::fwrite(bed, bed_file, sep = "\t", col.names = FALSE)

  trait_import <- importBedScore(SummarizedExperiment::rowRanges(SE_use), bed_file, colidx = 5)
  SE_DEV <- computeWeightedDeviations(
    object = SE_use,
    weights = trait_import,
    background_peaks = bg_peaks
  )

  z_score <- t(SummarizedExperiment::assays(SE_DEV)[["z"]]) |> c()
  names(z_score) <- colnames(SE_use)

  seed_idx <- seedindex(z_score = z_score, percent_cut = seed_percent)
  seed_cells <- names(z_score)[seed_idx]
  target_vec <- setNames(is_target_type(tr, meta_use$BioClassification), meta_use$cell)

  cell_dt <- data.table(
    trait = tr,
    cell = names(z_score),
    z_score = as.numeric(z_score),
    seed = names(z_score) %in% seed_cells
  )
  cell_dt <- merge(cell_dt, meta_use, by = "cell", all.x = TRUE)
  cell_dt[, `:=`(
    target_type = is_target_type(tr, BioClassification),
    target_definition = target_definition(tr)
  )]
  cell_signal_list[[tr]] <- cell_dt

  seed_target <- cell_dt[seed == TRUE]
  all_target <- mean(cell_dt$target_type, na.rm = TRUE)
  seed_frac <- mean(seed_target$target_type, na.rm = TRUE)
  seed_self_list[[tr]] <- data.table(
    trait = tr,
    n_gwas_rows = nrow(dt),
    n_cells = nrow(cell_dt),
    n_seed = nrow(seed_target),
    seed_percent = seed_percent,
    target_fraction_seed = seed_frac,
    target_fraction_all = all_target,
    enrichment_ratio = if (is.finite(all_target) && all_target > 0) seed_frac / all_target else NA_real_,
    median_z_seed = median(seed_target$z_score),
    median_z_all = median(cell_dt$z_score),
    target_definition = target_definition(tr)
  )

  for (graph_name in names(graphs)) {
    hop_list[[paste(tr, graph_name, sep = "_")]] <- hop_enrichment(
      adj = graphs[[graph_name]],
      seed_cells = seed_cells,
      target_vec = target_vec,
      graph_name = graph_name,
      trait = tr
    )
  }

  cat("  GWAS rows:", nrow(dt), "seed:", length(seed_cells), "\n")
}

if (length(seed_self_list) == 0) {
  stop("No trait diagnostics were generated", call. = FALSE)
}

seed_self_dt <- rbindlist(seed_self_list, use.names = TRUE, fill = TRUE)
hop_dt <- rbindlist(hop_list, use.names = TRUE, fill = TRUE)
cell_signal_dt <- rbindlist(cell_signal_list, use.names = TRUE, fill = TRUE)

seed_self_out <- file.path(out_dir, "seed_self_enrichment.csv")
hop_out <- file.path(out_dir, "seed_graph_hop_enrichment.csv")
cell_signal_out <- file.path(out_dir, "seed_signal_by_cell.csv.gz")
fwrite(seed_self_dt, seed_self_out)
fwrite(hop_dt, hop_out)
fwrite(cell_signal_dt, cell_signal_out)

plot_dt <- hop_dt[
  neighborhood %in% c("seed_self", "hop1_exact", "hop2_exact", "within1_including_seed", "within2_including_seed")
]
plot_dt[, neighborhood := factor(
  neighborhood,
  levels = c("seed_self", "hop1_exact", "hop2_exact", "within1_including_seed", "within2_including_seed")
)]

p <- ggplot(plot_dt, aes(x = neighborhood, y = target_fraction, fill = graph)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.68) +
  geom_hline(aes(yintercept = target_fraction_all), linetype = "dashed", color = "grey35") +
  facet_wrap(~ trait, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c(original_lsi = "#8C8C8C", wnn_lite = "#4C78A8")) +
  labs(
    title = "Seed and graph-neighborhood target-cell enrichment",
    subtitle = "Dashed line is target-cell fraction among all ATAC cells",
    x = NULL,
    y = "Target-cell fraction"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )
fig_out <- file.path(out_dir, "Fig_seed_hop_target_enrichment.png")
ggsave(fig_out, p, width = 10, height = 8, dpi = 300)

manifest <- list(
  script = "01_seed_signal_diagnostics.R",
  traits = traits,
  seed_percent = seed_percent,
  atac_se = paths$atac_se,
  gwas = paths$gwas,
  bg_peaks = paths$bg_peaks,
  original_knn = paths$original_knn,
  wnn_lite_knn = paths$wnn_lite_knn,
  output_dir = out_dir,
  alignment_csv = file.path(out_dir, "seed_signal_alignment_summary.csv"),
  seed_self_csv = seed_self_out,
  hop_enrichment_csv = hop_out,
  seed_signal_by_cell_csv = cell_signal_out,
  figure = fig_out,
  target_definitions = setNames(lapply(traits, target_definition), traits),
  created_at = as.character(Sys.time())
)
jsonlite::write_json(
  manifest,
  file.path(out_dir, "seed_signal_diagnostics_manifest.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

cat("\n========== Done ==========\n")
cat("  seed self :", seed_self_out, "\n")
cat("  graph hops:", hop_out, "\n")
cat("  per cell  :", cell_signal_out, "\n")
cat("  figure    :", fig_out, "\n")
