# 01: Build cell-level RNA-inferred type distribution for every ATAC cell.
#
# Inputs:
#   - atac_topk_rna_neighbors.csv/.rds (atac_cell, rank, rna_cell, distance)
#   - rna_cells.csv (cell, cell_type)
#   - atac_cells.csv (cell, cell_type)            # for sanity-check only
#
# Outputs (in MODULE2_SEED_FILTER_OUT_DIR):
#   - cell_rna_neighbor_dominant.csv      # one row per ATAC cell
#   - cell_trait_rna_target_fraction.csv  # one row per (atac_cell, trait)
#   - rna_inferred_type_manifest.json
#
# Parameters (env vars):
#   MODULE2_SEED_FILTER_K        : top-K neighbors used for voting (default 30)
#   MODULE2_TRAITS               : comma-separated traits (default "MONO,PLT,LYM")

source("D:/RStudio/work/module2_seed_rna_filter/_bootstrap.R")

required_pkgs(c("data.table", "jsonlite"))
suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

paths <- seed_filter_paths()
ensure_seed_filter_dirs(paths)

k_neighbors <- env_integer("MODULE2_SEED_FILTER_K", 30L)
traits <- parse_csv_env("MODULE2_TRAITS", "MONO,PLT,LYM")
if (length(traits) == 0) stop("MODULE2_TRAITS is empty", call. = FALSE)

cat("========== Module 2 seed RNA filter: 01 build RNA-inferred type ==========\n")
cat("  k_neighbors :", k_neighbors, "\n")
cat("  traits      :", paste(traits, collapse = ", "), "\n")
cat("  topk csv    :", paths$topk_neighbors_csv, "\n")
cat("  rna_cells   :", paths$rna_cells, "\n")
cat("  out_dir     :", paths$out_dir, "\n")

stopifnot(
  file.exists(paths$rna_cells),
  file.exists(paths$atac_cells)
)

# Load --------------------------------------------------------------------

cat("\n[1/4] Load RNA cell types\n")
rna_cells <- fread(paths$rna_cells)
stopifnot(all(c("cell", "cell_type") %in% names(rna_cells)))
setnames(rna_cells, c("cell", "cell_type"), c("rna_cell", "rna_cell_type"))
cat("  RNA cells:", nrow(rna_cells), "\n")

cat("\n[2/4] Load top-K RNA neighbors\n")
neighbors <- load_topk_neighbors(paths, k_max = k_neighbors)
n_atac_in_topk <- length(unique(neighbors$atac_cell))
cat("  topk rows :", nrow(neighbors), "\n")
cat("  unique ATAC cells in topk:", n_atac_in_topk, "\n")

# Join in rna cell types
neighbors <- merge(neighbors, rna_cells, by = "rna_cell", all.x = TRUE)
n_missing <- sum(is.na(neighbors$rna_cell_type))
if (n_missing > 0) {
  warning("Missing rna_cell_type for ", n_missing, " neighbor rows; dropping them.")
  neighbors <- neighbors[!is.na(rna_cell_type)]
}

# Per-cell dominant type --------------------------------------------------

cat("\n[3/4] Compute per-cell dominant RNA-inferred type\n")
dominant_dt <- neighbors[
  ,
  {
    n_used <- .N
    tab <- sort(table(rna_cell_type), decreasing = TRUE)
    n_unique <- length(tab)
    top_names <- names(tab)
    top_fracs <- as.numeric(tab) / n_used
    list(
      n_neighbors = n_used,
      n_unique_types = n_unique,
      top1_type = top_names[1],
      top1_fraction = top_fracs[1],
      top2_type = if (n_unique >= 2) top_names[2] else NA_character_,
      top2_fraction = if (n_unique >= 2) top_fracs[2] else NA_real_,
      top3_type = if (n_unique >= 3) top_names[3] else NA_character_,
      top3_fraction = if (n_unique >= 3) top_fracs[3] else NA_real_
    )
  },
  by = atac_cell
]

dominant_out <- file.path(paths$out_dir, "cell_rna_neighbor_dominant.csv")
fwrite(dominant_dt, dominant_out)
cat("  wrote:", dominant_out, "\n")

# Per (cell, trait) target fraction ---------------------------------------

cat("\n[4/4] Compute per (cell, trait) RNA-target fraction\n")
trait_fraction_list <- vector("list", length(traits))
for (i in seq_along(traits)) {
  tr <- traits[i]
  cat("  trait:", tr, "(", target_definition(tr), ")\n")
  neighbors[, is_target_tmp := is_target_type(tr, rna_cell_type)]
  trait_dt <- neighbors[
    ,
    list(
      trait = tr,
      n_neighbors = .N,
      n_target_neighbors = sum(is_target_tmp),
      target_fraction = sum(is_target_tmp) / .N
    ),
    by = atac_cell
  ]
  trait_fraction_list[[i]] <- trait_dt
}
neighbors[, is_target_tmp := NULL]
trait_fraction_dt <- rbindlist(trait_fraction_list, use.names = TRUE)

trait_out <- file.path(paths$out_dir, "cell_trait_rna_target_fraction.csv")
fwrite(trait_fraction_dt, trait_out)
cat("  wrote:", trait_out, "\n")

# Manifest ----------------------------------------------------------------

manifest <- list(
  script = "01_build_rna_inferred_type.R",
  k_neighbors = k_neighbors,
  traits = traits,
  target_definitions = setNames(lapply(traits, target_definition), traits),
  inputs = list(
    rna_cells = paths$rna_cells,
    topk_neighbors = paths$topk_neighbors_csv
  ),
  outputs = list(
    cell_dominant = dominant_out,
    cell_trait_target_fraction = trait_out
  ),
  n_atac_cells_with_neighbors = n_atac_in_topk,
  n_neighbors_dropped_missing_rna_type = n_missing,
  created_at = as.character(Sys.time())
)
manifest_out <- file.path(paths$out_dir, "rna_inferred_type_manifest.json")
write_json(manifest, manifest_out, pretty = TRUE, auto_unbox = TRUE)

cat("\n========== Done ==========\n")
cat("  manifest:", manifest_out, "\n")
