# Define Version A cell-type-aware RNA gating set.

source("D:/RStudio/work/module2_wnn_lite_gated/_bootstrap.R")

required_pkgs(c("data.table", "SummarizedExperiment", "jsonlite"))
suppressPackageStartupMessages({
  library(data.table)
  library(SummarizedExperiment)
  library(jsonlite)
})

paths <- gated_paths("celltype_gated")
ensure_gated_dirs(paths)

cat("========== WNN-lite gated: define celltype gating ==========\n")
cat("  ATAC SE      :", paths$atac_se, "\n")
cat("  same type CSV:", paths$base_same_type_by_celltype, "\n")
cat("  seed signal  :", paths$seed_signal_by_cell, "\n")
cat("  output dir   :", paths$out_dir, "\n")

meta <- load_atac_metadata(paths$atac_se)
same_type <- if (file.exists(paths$base_same_type_by_celltype)) {
  data.table::fread(paths$base_same_type_by_celltype)
} else {
  NULL
}
seed_signal <- if (file.exists(paths$seed_signal_by_cell)) {
  data.table::fread(paths$seed_signal_by_cell)
} else {
  NULL
}

celltype_counts <- meta[, .(n_cells = .N), by = BioClassification][order(BioClassification)]
celltype_counts[, gated := is_transition_celltype(BioClassification)]
celltype_counts[, reason := fifelse(
  gated,
  "biological_prior: progenitor_or_transition_state",
  "not_gated"
)]
setnames(celltype_counts, "BioClassification", "cell_type")

if (!is.null(same_type)) {
  join_col <- intersect(c("BioClassification", "cell_type", "atac_cell_type"), names(same_type))[1]
  if (!is.na(join_col)) {
    same_type_small <- copy(same_type)
    setnames(same_type_small, join_col, "cell_type")
    celltype_counts <- merge(celltype_counts, same_type_small, by = "cell_type", all.x = TRUE)
  }
}

gated_celltypes <- celltype_counts[gated == TRUE]
gated_cells <- meta[BioClassification %in% gated_celltypes$cell_type]
gated_cells[, gating_reason := "celltype_prior"]

gated_fraction <- nrow(gated_cells) / nrow(meta)
summary_dt <- data.table(
  version = "celltype_gated",
  n_cells_total = nrow(meta),
  n_cells_gated = nrow(gated_cells),
  pct_cells_gated = gated_fraction,
  n_celltypes_total = nrow(celltype_counts),
  n_celltypes_gated = nrow(gated_celltypes),
  sanity_bound = 0.50,
  passes_sanity_bound = gated_fraction <= 0.50
)

seed_overlap <- data.table()
if (!is.null(seed_signal)) {
  seed_overlap <- seed_signal[, .(
    n_seed = sum(seed),
    n_seed_gated = sum(seed & cell %in% gated_cells$cell),
    pct_seed_gated = if (sum(seed) > 0) sum(seed & cell %in% gated_cells$cell) / sum(seed) else NA_real_,
    n_cells = .N,
    n_cells_gated = sum(cell %in% gated_cells$cell),
    pct_cells_gated = mean(cell %in% gated_cells$cell)
  ), by = trait][order(trait)]
  seed_overlap[, flag_seed_gated_over_10pct := pct_seed_gated > 0.10]
}

gated_celltypes_out <- file.path(paths$out_dir, "gated_celltypes.csv")
gated_cells_out <- file.path(paths$out_dir, "gated_cells.csv")
summary_out <- file.path(paths$out_dir, "gating_summary.csv")
seed_overlap_out <- file.path(paths$out_dir, "gated_seed_overlap.csv")
data.table::fwrite(gated_celltypes, gated_celltypes_out)
data.table::fwrite(gated_cells, gated_cells_out)
data.table::fwrite(summary_dt, summary_out)
if (nrow(seed_overlap) > 0) data.table::fwrite(seed_overlap, seed_overlap_out)

jsonlite::write_json(
  list(
    script = "celltype_gated/01_define_celltype_gating.R",
    rule = "BioClassification contains HSC|CMP|LMPP|CLP|GMP|Pre|Early|Unk",
    atac_se = paths$atac_se,
    same_type_by_celltype = paths$base_same_type_by_celltype,
    seed_signal = paths$seed_signal_by_cell,
    gated_celltypes = gated_celltypes_out,
    gated_cells = gated_cells_out,
    summary = summary_out,
    seed_overlap = if (nrow(seed_overlap) > 0) seed_overlap_out else NULL,
    created_at = as.character(Sys.time())
  ),
  file.path(paths$out_dir, "celltype_gating_manifest.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

cat("\n[Summary]\n")
print(summary_dt)
if (nrow(seed_overlap) > 0) {
  cat("\n[Seed overlap]\n")
  print(seed_overlap)
}
cat("\n========== Done ==========\n")
cat("  gated celltypes:", gated_celltypes_out, "\n")
cat("  gated cells    :", gated_cells_out, "\n")
