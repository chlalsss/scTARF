# Build Version A celltype-gated WNN-lite weights.

source("D:/RStudio/work/module2_wnn_lite_gated/_bootstrap.R")

required_pkgs(c("data.table", "jsonlite"))
suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

paths <- gated_paths("celltype_gated")
ensure_gated_dirs(paths)

gated_cells_path <- file.path(paths$out_dir, "gated_cells.csv")
out_weights <- file.path(paths$out_dir, "wnn_lite_celltype_gated_weights.csv")

cat("========== WNN-lite gated: build celltype-gated weights ==========\n")
cat("  base weights:", paths$base_weights, "\n")
cat("  gated cells :", gated_cells_path, "\n")
cat("  output      :", out_weights, "\n")

stopifnot(file.exists(paths$base_weights), file.exists(gated_cells_path))

weights <- data.table::fread(paths$base_weights)
gated_cells <- data.table::fread(gated_cells_path)

cell_col <- intersect(c("cell", "atac_cell"), names(weights))[1]
if (is.na(cell_col)) stop("Base weights must contain cell or atac_cell column", call. = FALSE)
if (!all(c("rna_weight", "atac_weight") %in% names(weights))) {
  stop("Base weights must contain rna_weight and atac_weight", call. = FALSE)
}

weights[, gated := get(cell_col) %in% gated_cells$cell]
weights[, gating_version := "celltype_gated"]
weights[, gating_reason := fifelse(gated, "celltype_prior", "not_gated")]
weights[gated == TRUE, `:=`(
  rna_weight_original = rna_weight,
  atac_weight_original = atac_weight,
  rna_weight = 0,
  atac_weight = 1
)]
weights[gated == FALSE, `:=`(
  rna_weight_original = rna_weight,
  atac_weight_original = atac_weight
)]

summary_dt <- weights[, .(
  n_cells = .N,
  n_gated = sum(gated),
  pct_gated = mean(gated),
  mean_rna_weight_original = mean(rna_weight_original, na.rm = TRUE),
  mean_rna_weight_gated = mean(rna_weight, na.rm = TRUE),
  sanity_bound = 0.50,
  passes_sanity_bound = mean(gated) <= 0.50
)]

data.table::fwrite(weights, out_weights)
summary_out <- file.path(paths$out_dir, "gated_weight_summary.csv")
data.table::fwrite(summary_dt, summary_out)

jsonlite::write_json(
  list(
    script = "celltype_gated/02_build_gated_weights.R",
    base_weights = paths$base_weights,
    gated_cells = gated_cells_path,
    output_weights = out_weights,
    summary = summary_out,
    created_at = as.character(Sys.time())
  ),
  file.path(paths$out_dir, "celltype_gated_weight_manifest.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

cat("\n[Summary]\n")
print(summary_dt)
cat("\n========== Done ==========\n")
cat("  weights:", out_weights, "\n")
