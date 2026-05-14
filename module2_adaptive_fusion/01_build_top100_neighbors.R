# Recompute ATAC -> RNA top-k neighbors from scSHEFT shared latent.

source("D:/RStudio/work/module2_adaptive_fusion/_bootstrap.R")

required_pkgs(c("data.table", "RANN", "jsonlite"))
suppressPackageStartupMessages({
  library(data.table)
  library(RANN)
  library(jsonlite)
})

paths <- module2_paths()
ensure_module2_dirs(paths)

k <- as.integer(Sys.getenv("MODULE2_TOPK", "1"))
if (!is.finite(k) || k < 1) stop("MODULE2_TOPK must be a positive integer")

cat("========== Module 2: build top", k, "RNA neighbors ==========\n", sep = "")
cat("  shared embedding:", paths$shared_embedding, "\n")
cat("  output dir      :", paths$out_dir, "\n")

loaded <- load_shared_embedding(paths$shared_embedding)
emb <- as.data.table(loaded$emb)
latent_cols <- loaded$latent_cols

rna_dt <- emb[modality == "rna"]
atac_dt <- emb[modality == "atac"]
if (nrow(rna_dt) == 0 || nrow(atac_dt) == 0) {
  stop("Empty RNA or ATAC rows in shared embedding")
}

rna_mat <- as.matrix(rna_dt[, ..latent_cols])
atac_mat <- as.matrix(atac_dt[, ..latent_cols])
storage.mode(rna_mat) <- "double"
storage.mode(atac_mat) <- "double"

cat("  RNA matrix :", nrow(rna_mat), "x", ncol(rna_mat), "\n")
cat("  ATAC matrix:", nrow(atac_mat), "x", ncol(atac_mat), "\n")

t0 <- Sys.time()
nn <- RANN::nn2(data = rna_mat, query = atac_mat, k = k)
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("  RANN elapsed: %.1f sec\n", elapsed))

topk_dt <- data.table(
  atac_cell = rep(atac_dt$cell, each = k),
  rank = rep(seq_len(k), times = nrow(atac_dt)),
  rna_cell = rna_dt$cell[as.vector(t(nn$nn.idx))],
  distance = as.vector(t(nn$nn.dists))
)

out_csv <- file.path(paths$out_dir, sprintf("atac_top%d_rna_neighbors.csv", k))
data.table::fwrite(topk_dt, out_csv)
cat("  saved:", out_csv, "\n")

summary <- list(
  script = "01_build_top100_neighbors.R",
  k = k,
  distance = "Euclidean distance from RANN::nn2 on raw exported scSHEFT latent",
  shared_embedding = paths$shared_embedding,
  n_rna = nrow(rna_dt),
  n_atac = nrow(atac_dt),
  n_latent = length(latent_cols),
  elapsed_seconds = elapsed,
  output_csv = out_csv,
  created_at = as.character(Sys.time())
)
jsonlite::write_json(
  summary,
  file.path(paths$out_dir, sprintf("top%d_neighbors_manifest.json", k)),
  pretty = TRUE,
  auto_unbox = TRUE
)

cat("========== Done ==========\n")

