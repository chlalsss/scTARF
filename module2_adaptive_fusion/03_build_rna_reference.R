# Build weighted RNA reference latent for each ATAC cell.

source("D:/RStudio/work/module2_adaptive_fusion/_bootstrap.R")

required_pkgs(c("data.table", "Matrix", "jsonlite"))
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(jsonlite)
})

paths <- module2_paths()
ensure_module2_dirs(paths)
label <- get_active_label()
cutoff <- if (label == "top1") NA_real_ else get_single_cutoff()

selected_csv <- file.path(paths$out_dir, sprintf("adaptive_neighbors_cutoff_%s.csv", label))
if (!file.exists(selected_csv)) {
  stop("Missing adaptive neighbor file. Run 02_adaptive_select_neighbors.R first: ", selected_csv)
}

cat("========== Module 2: build RNA reference latent ==========\n")
cat("  cutoff:", cutoff, "\n")
cat("  input :", selected_csv, "\n")

loaded <- load_shared_embedding(paths$shared_embedding)
emb <- as.data.table(loaded$emb)
latent_cols <- loaded$latent_cols
rna_dt <- emb[modality == "rna"]
atac_dt <- emb[modality == "atac"]

selected <- fread(selected_csv)
needed <- c("atac_cell", "rna_cell", "weight")
if (!all(needed %in% colnames(selected))) {
  stop("Selected table must contain: ", paste(needed, collapse = ", "))
}

atac_index <- setNames(seq_len(nrow(atac_dt)), atac_dt$cell)
rna_index <- setNames(seq_len(nrow(rna_dt)), rna_dt$cell)
selected[, atac_idx := atac_index[atac_cell]]
selected[, rna_idx := rna_index[rna_cell]]
selected <- selected[!is.na(atac_idx) & !is.na(rna_idx)]

weight_sum <- selected[, .(weight_sum = sum(weight)), by = atac_cell]
if (max(abs(weight_sum$weight_sum - 1), na.rm = TRUE) > 1e-6) {
  warning("Weights are not exactly normalized per ATAC cell; renormalizing.")
  selected[, weight := weight / sum(weight), by = atac_cell]
}

W <- Matrix::sparseMatrix(
  i = selected$atac_idx,
  j = selected$rna_idx,
  x = selected$weight,
  dims = c(nrow(atac_dt), nrow(rna_dt))
)
rna_mat <- as.matrix(rna_dt[, ..latent_cols])
storage.mode(rna_mat) <- "double"

rna_ref <- as.matrix(W %*% rna_mat)
colnames(rna_ref) <- paste0("rna_ref_", seq_len(ncol(rna_ref)))
rna_ref_df <- data.frame(
  atac_cell = atac_dt$cell,
  rna_ref,
  check.names = FALSE
)

out_rds <- file.path(paths$out_dir, sprintf("rna_reference_latent_cutoff_%s.rds", label))
saveRDS(rna_ref_df, out_rds, compress = FALSE)
cat("  saved:", out_rds, "\n")

manifest <- list(
  script = "03_build_rna_reference.R",
  cutoff = cutoff,
  selected_neighbors = selected_csv,
  shared_embedding = paths$shared_embedding,
  n_atac = nrow(atac_dt),
  n_rna = nrow(rna_dt),
  n_latent = length(latent_cols),
  output_rds = out_rds,
  created_at = as.character(Sys.time())
)
jsonlite::write_json(
  manifest,
  file.path(paths$out_dir, sprintf("rna_reference_manifest_cutoff_%s.json", label)),
  pretty = TRUE,
  auto_unbox = TRUE
)

cat("========== Done ==========\n")
