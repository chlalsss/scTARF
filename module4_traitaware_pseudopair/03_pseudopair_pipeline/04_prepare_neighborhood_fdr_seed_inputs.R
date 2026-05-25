# module4: prepare module3 v2 inputs with neighborhood-FDR seeds
# -----------------------------------------------------------------------------
# Purpose:
#   Keep the existing module3 v2 fused KNN/RWR pipeline, but replace the old
#   pseudo-zRNA top-1% seeds with seeds derived from scDRS-significant RNA
#   neighborhoods.
#
# Main idea:
#   1. scDRS first identifies trait-associated RNA neighborhoods by group FDR.
#   2. For each ATAC cell, we already know how much of its mapped RNA-neighbor
#      weight falls inside those significant RNA neighborhoods.
#   3. ATAC cells with strong support from significant RNA neighborhoods become
#      RWR seeds.
#
# Default seed rule:
#   candidate_group_weight_ge0p80:
#     keep ATAC cells whose mapped RNA-neighbor weight to significant
#     trait-associated RNA neighborhoods is at least 0.80.
#
# Method boundary:
#   - Cell type is NOT used to build the RNA neighborhoods.
#   - Cell type is NOT used for scDRS FDR.
#   - Cell type is NOT used to select ATAC seeds here.
#   - PLT is excluded by default because post hoc diagnostics showed the current
#     neighborhood-FDR signal maps almost entirely to myeloid RNA sources and
#     gives near-zero PLT seed precision. That is a diagnostic decision, not a
#     selection rule inside the method.
#
# Run in RStudio:
#   setwd("D:/RStudio/work/module4_traitaware_pseudopair")
#   source("03_pseudopair_pipeline/04_prepare_neighborhood_fdr_seed_inputs.R")
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

ROOT <- "D:/RStudio/work/module4_traitaware_pseudopair"
if (!exists("BASE_COMPAT_DIR")) {
  BASE_COMPAT_DIR <- file.path(ROOT, "03_pseudopair_pipeline/output/module3_v2_compat")
}
if (!exists("NBD_DIR")) {
  NBD_DIR <- file.path(ROOT, "04_evaluation/module4_v2_diagnostics/neighborhood_fdr")
}
if (!exists("OUT_DIR")) {
  OUT_DIR <- file.path(ROOT, "03_pseudopair_pipeline/output/module3_v2_neighborhood_fdr_ge0p80")
}

# Use neighborhood-specific variable names so old variables left in the RStudio
# global environment by earlier scripts cannot accidentally pull PLT back in.
if (!exists("NBD_TRAITS_TO_PREPARE")) NBD_TRAITS_TO_PREPARE <- c("MONO", "LYM")
if (!exists("NBD_ATAC_W_LIST")) NBD_ATAC_W_LIST <- c(0.7, 1.0)
if (!exists("NBD_GROUP_FDR_THRESHOLD")) NBD_GROUP_FDR_THRESHOLD <- 0.10
if (!exists("NBD_SEED_CANDIDATE_COLUMN")) {
  NBD_SEED_CANDIDATE_COLUMN <- "candidate_group_weight_ge0p80"
}

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

stop_missing <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0L) {
    stop("Missing required file(s):\n  ", paste(missing, collapse = "\n  "),
         call. = FALSE)
  }
}

copy_required <- function(from, to) {
  if (!file.exists(from)) stop("Missing file to copy: ", from, call. = FALSE)
  ok <- file.copy(from, to, overwrite = TRUE)
  if (!ok) stop("Failed to copy:\n  from: ", from, "\n  to: ", to, call. = FALSE)
}

cat(">>> Prepare module3 v2 inputs with neighborhood-FDR seeds\n\n")
cat("Base compat dir : ", BASE_COMPAT_DIR, "\n", sep = "")
cat("Neighborhood dir: ", NBD_DIR, "\n", sep = "")
cat("Output dir      : ", OUT_DIR, "\n", sep = "")
cat("Traits          : ", paste(NBD_TRAITS_TO_PREPARE, collapse = ", "), "\n", sep = "")
cat("Group FDR       : ", NBD_GROUP_FDR_THRESHOLD, "\n", sep = "")
cat("Seed rule       : ", NBD_SEED_CANDIDATE_COLUMN, "\n\n", sep = "")

cell_file <- file.path(NBD_DIR, "neighborhood_fdr_cell_level.csv.gz")
summary_file <- file.path(NBD_DIR, "neighborhood_fdr_seed_summary.csv")
stop_missing(c(cell_file, summary_file, file.path(BASE_COMPAT_DIR, "atac_z_block_v2.rds")))

needed_cols <- c(
  "trait",
  "group_fdr_threshold",
  "barcode",
  "group_sig_source_weight",
  NBD_SEED_CANDIDATE_COLUMN,
  "no_celltype_used_to_select"
)

cat("[1/4] Loading neighborhood-FDR ATAC-level evidence ...\n")
cell_dt <- fread(cell_file, select = needed_cols)
if (!all(needed_cols %in% names(cell_dt))) {
  stop("Neighborhood cell-level file does not contain required columns.",
       call. = FALSE)
}
cell_dt <- cell_dt[group_fdr_threshold == NBD_GROUP_FDR_THRESHOLD]
cell_dt <- cell_dt[trait %in% NBD_TRAITS_TO_PREPARE]
if (nrow(cell_dt) == 0L) {
  stop("No rows found for requested traits and FDR threshold.", call. = FALSE)
}
if (!all(cell_dt$no_celltype_used_to_select %in% TRUE)) {
  stop("Input says cell type may have been used for selection; aborting.",
       call. = FALSE)
}

cat("[2/4] Copying reusable module3 v2 blocks and KNN files ...\n")
copy_required(
  file.path(BASE_COMPAT_DIR, "atac_z_block_v2.rds"),
  file.path(OUT_DIR, "atac_z_block_v2.rds")
)

for (trait in NBD_TRAITS_TO_PREPARE) {
  copy_required(
    file.path(BASE_COMPAT_DIR, paste0("rna_z_block_v2_", trait, ".rds")),
    file.path(OUT_DIR, paste0("rna_z_block_v2_", trait, ".rds"))
  )
  for (w in NBD_ATAC_W_LIST) {
    copy_required(
      file.path(BASE_COMPAT_DIR, sprintf("knn_v2_%s_w%.1f.rds", trait, w)),
      file.path(OUT_DIR, sprintf("knn_v2_%s_w%.1f.rds", trait, w))
    )
  }
}

cat("[3/4] Writing seed_v2_<trait>.rds from neighborhood-FDR evidence ...\n")
seed_rows <- list()
for (trait_name in NBD_TRAITS_TO_PREPARE) {
  dt <- cell_dt[trait == trait_name]
  if (anyDuplicated(dt$barcode)) {
    stop("Duplicated ATAC barcodes in neighborhood cell-level file for ", trait_name,
         call. = FALSE)
  }

  candidate <- as.logical(dt[[NBD_SEED_CANDIDATE_COLUMN]])
  candidate[is.na(candidate)] <- FALSE

  score <- dt$group_sig_source_weight
  score[!is.finite(score) | is.na(score)] <- 0
  score[!candidate] <- 0
  names(score) <- dt$barcode

  if (sum(score > 0) == 0L) {
    stop("No positive neighborhood-FDR seeds for trait ", trait_name, call. = FALSE)
  }

  saveRDS(score, file.path(OUT_DIR, paste0("seed_v2_", trait_name, ".rds")))

  seed_rows[[trait_name]] <- data.frame(
    trait = trait_name,
    n_cells = length(score),
    n_seed_positive = sum(score > 0),
    group_fdr_threshold = NBD_GROUP_FDR_THRESHOLD,
    seed_candidate_column = NBD_SEED_CANDIDATE_COLUMN,
    seed_score_source = "group_sig_source_weight among selected ATAC cells",
    seed_score_min_positive = min(score[score > 0]),
    seed_score_median_positive = median(score[score > 0]),
    seed_score_max = max(score),
    no_celltype_used = TRUE,
    stringsAsFactors = FALSE
  )

  cat("  ", trait_name, ": n_seed_positive=", sum(score > 0),
      "; median positive score=", signif(median(score[score > 0]), 4),
      "; max=", signif(max(score), 4), "\n", sep = "")
}

seed_summary <- rbindlist(seed_rows, use.names = TRUE, fill = TRUE)
fwrite(seed_summary, file.path(OUT_DIR, "seed_summary_neighborhood_fdr.csv"))

cat("[4/4] Writing manifest ...\n")
manifest <- data.frame(
  key = c(
    "atac_z_block",
    paste0("rna_z_block_", NBD_TRAITS_TO_PREPARE),
    paste0("seed_", NBD_TRAITS_TO_PREPARE),
    as.vector(outer(
      paste0("knn_", NBD_TRAITS_TO_PREPARE),
      sprintf("_w%.1f", NBD_ATAC_W_LIST),
      paste0
    ))
  ),
  path = c(
    file.path(OUT_DIR, "atac_z_block_v2.rds"),
    file.path(OUT_DIR, paste0("rna_z_block_v2_", NBD_TRAITS_TO_PREPARE, ".rds")),
    file.path(OUT_DIR, paste0("seed_v2_", NBD_TRAITS_TO_PREPARE, ".rds")),
    as.vector(sapply(NBD_TRAITS_TO_PREPARE, function(trait) {
      file.path(OUT_DIR, sprintf("knn_v2_%s_w%.1f.rds", trait, NBD_ATAC_W_LIST))
    }))
  ),
  stringsAsFactors = FALSE
)
fwrite(manifest, file.path(OUT_DIR, "module3_v2_neighborhood_fdr_manifest.csv"))

cat("\n=== Neighborhood-FDR seed inputs are ready ===\n")
print(seed_summary)
cat("\nOutputs are in: ", OUT_DIR, "\n", sep = "")
