source("D:/RStudio/work/module2_robustness/_bootstrap.R")

required_pkgs(c("SummarizedExperiment", "Matrix", "data.table", "jsonlite"))
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(Matrix)
  library(data.table)
  library(jsonlite)
})

p <- robustness_paths()
ensure_robustness_dirs(p)

rates <- parse_int_list(Sys.getenv("MODULE2_ROBUSTNESS_RATES", "0,30"), c(0L, 30L))
seeds <- parse_int_list(Sys.getenv("MODULE2_ROBUSTNESS_SEEDS", "1,2,3"), 1:3)
rates <- sort(unique(rates))
seeds <- sort(unique(seeds))
if (any(!rates %in% c(0L, 10L, 20L, 30L, 40L))) {
  stop("Rates should be selected from 0,10,20,30,40")
}

cat("========== Robustness: generate sparse ATAC SE ==========\n")
cat("  input :", p$original_atac_se, "\n")
cat("  rates :", paste(rates, collapse = ","), "\n")
cat("  seeds :", paste(seeds, collapse = ","), "\n")
cat("  out   :", p$sparse_dir, "\n")

se <- readRDS(p$original_atac_se)
mat <- assays(se)[["counts"]]
if (is.null(mat)) stop("ATAC SE must contain counts assay")
mat <- as(mat, "dgCMatrix")
nnz <- length(mat@x)
cat("  matrix:", nrow(mat), "peaks x", ncol(mat), "cells; nnz =", nnz, "\n")

outputs <- list()

make_sparse_se <- function(rate, seed) {
  label <- sprintf("dropout_%02d_seed%d", rate, seed)
  out_rds <- file.path(p$sparse_dir, paste0("atac_", label, ".rds"))
  if (file.exists(out_rds) && !tolower(Sys.getenv("MODULE2_ROBUSTNESS_OVERWRITE", "false")) %in% c("true", "1", "yes")) {
    cat("  exists, skip:", out_rds, "\n")
    return(data.table(
      dropout_pct = rate,
      seed = seed,
      label = label,
      input_rds = out_rds,
      nnz_original = nnz,
      nnz_after = NA_integer_,
      dropped_n = NA_integer_,
      skipped_existing = TRUE
    ))
  }

  se2 <- se
  mat2 <- mat
  drop_n <- floor(as.numeric(nnz) * rate / 100)
  if (drop_n > 0) {
    set.seed(seed)
    drop_idx <- sample.int(nnz, drop_n, replace = FALSE)
    mat2@x[drop_idx] <- 0
    mat2 <- Matrix::drop0(mat2)
  }
  assays(se2)[["counts"]] <- mat2
  saveRDS(se2, out_rds, compress = FALSE)
  cat("  saved:", out_rds, " nnz_after =", length(mat2@x), "\n")
  data.table(
    dropout_pct = rate,
    seed = seed,
    label = label,
    input_rds = out_rds,
    nnz_original = nnz,
    nnz_after = length(mat2@x),
    dropped_n = nnz - length(mat2@x),
    skipped_existing = FALSE
  )
}

if (0L %in% rates) {
  outputs[["dropout_00_seed0"]] <- make_sparse_se(0L, 0L)
}
for (rate in setdiff(rates, 0L)) {
  for (seed in seeds) {
    outputs[[sprintf("dropout_%02d_seed%d", rate, seed)]] <- make_sparse_se(rate, seed)
  }
}

manifest_dt <- rbindlist(outputs, use.names = TRUE, fill = TRUE)
manifest_csv <- file.path(p$out_dir, "sparse_input_manifest.csv")
fwrite(manifest_dt, manifest_csv)
write_json(
  list(
    script = "01_generate_sparse_atac.R",
    rates = rates,
    seeds = seeds,
    original_atac_se = p$original_atac_se,
    manifest_csv = manifest_csv,
    created_at = as.character(Sys.time())
  ),
  file.path(p$out_dir, "sparse_input_manifest.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)
cat("manifest:", manifest_csv, "\n")
cat("========== Done ==========\n")
