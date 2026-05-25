# Define Version B reliability-aware RNA gating set.

source("D:/RStudio/work/module2_wnn_lite_gated/_bootstrap.R")

required_pkgs(c("data.table", "jsonlite"))
suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

paths <- gated_paths("reliability_gated")
ensure_gated_dirs(paths)

threshold <- env_numeric("MODULE2_RELIABILITY_GATING_THRESHOLD", 0.108)

cat("========== WNN-lite gated: define reliability gating ==========\n")
cat("  base weights:", paths$base_weights, "\n")
cat("  threshold   :", threshold, "\n")
cat("  output dir  :", paths$out_dir, "\n")

stopifnot(file.exists(paths$base_weights))

weights <- data.table::fread(paths$base_weights)
cell_col <- intersect(c("cell", "atac_cell"), names(weights))[1]
if (is.na(cell_col)) stop("Base weights must contain cell or atac_cell column", call. = FALSE)

reliability_col <- intersect(c("rna_reliability", "reliability", "rna_quality_score"), names(weights))[1]
if (is.na(reliability_col)) {
  stop(
    "Base weights must contain an RNA reliability column. Looked for: ",
    "rna_reliability, reliability, rna_quality_score",
    call. = FALSE
  )
}

weights[, gated := get(reliability_col) < threshold]
gated_cells <- weights[gated == TRUE, .(
  cell = get(cell_col),
  rna_reliability = get(reliability_col),
  gating_reason = sprintf("rna_reliability_lt_%.3f", threshold)
)]

summary_dt <- weights[, .(
  n_cells = .N,
  n_gated = sum(gated),
  pct_gated = mean(gated),
  threshold = threshold,
  reliability_col = reliability_col,
  q10 = as.numeric(stats::quantile(get(reliability_col), 0.10, na.rm = TRUE)),
  q25 = as.numeric(stats::quantile(get(reliability_col), 0.25, na.rm = TRUE)),
  median = as.numeric(stats::median(get(reliability_col), na.rm = TRUE)),
  sanity_bound = 0.40,
  passes_sanity_bound = mean(gated) <= 0.40
)]

gated_cells_out <- file.path(paths$out_dir, "gated_cells.csv")
summary_out <- file.path(paths$out_dir, "gating_summary.csv")
data.table::fwrite(gated_cells, gated_cells_out)
data.table::fwrite(summary_dt, summary_out)

jsonlite::write_json(
  list(
    script = "reliability_gated/01_define_reliability_gating.R",
    base_weights = paths$base_weights,
    reliability_col = reliability_col,
    threshold = threshold,
    gated_cells = gated_cells_out,
    summary = summary_out,
    created_at = as.character(Sys.time())
  ),
  file.path(paths$out_dir, "reliability_gating_manifest.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

cat("\n[Summary]\n")
print(summary_dt)
cat("\n========== Done ==========\n")
cat("  gated cells:", gated_cells_out, "\n")
