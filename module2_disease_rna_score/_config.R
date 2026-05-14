required_pkgs <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Missing packages: ", paste(missing, collapse = ", "))
  }
}

disease_rna_paths <- function() {
  data_dir <- Sys.getenv(
    "MODULE2_DATA_DIR",
    "D:/SCAVENGE-reproducibility-main/blood_data"
  )
  scsheft_input_dir <- file.path(
    data_dir,
    "multiomics_module1_pairing",
    "method_B_scSHEFT",
    "scsheft_input"
  )
  adaptive_dir <- file.path(data_dir, "module2_adaptive_fusion")
  out_dir <- Sys.getenv(
    "MODULE2_DISEASE_OUT_DIR",
    file.path(data_dir, "module2_disease_rna_score")
  )

  list(
    data_dir = data_dir,
    rna_mtx = file.path(scsheft_input_dir, "rna_counts_shared_genes.mtx"),
    shared_genes = file.path(scsheft_input_dir, "shared_genes.csv"),
    rna_cells = file.path(scsheft_input_dir, "rna_cells.csv"),
    trait_gene_csv = file.path(
      data_dir,
      "module2_gwas_gene_risk",
      "trait_gene_risk",
      "trait_gene_risk_long.csv"
    ),
    top1_neighbors = file.path(adaptive_dir, "atac_top1_rna_neighbors.csv"),
    adaptive_dir = adaptive_dir,
    out_dir = out_dir
  )
}

ensure_disease_dirs <- function(p) {
  dir.create(p$out_dir, showWarnings = FALSE, recursive = TRUE)
}

auto_select_pca_dim <- function(sdev, method = "elbow", min_dim = 2L) {
  vals <- as.numeric(sdev)^2
  vals <- vals[is.finite(vals) & vals > 0]
  n <- length(vals)
  if (n == 0) stop("Cannot auto-select PCA dimension from empty sdev")
  if (n <= min_dim) return(n)

  y <- log10(vals + .Machine$double.eps)
  x <- seq_along(y)
  x <- (x - min(x)) / max(max(x) - min(x), .Machine$double.eps)
  y <- (y - min(y)) / max(max(y) - min(y), .Machine$double.eps)

  x1 <- x[1]
  y1 <- y[1]
  x2 <- x[n]
  y2 <- y[n]
  denom <- sqrt((y2 - y1)^2 + (x2 - x1)^2)
  dist <- abs((y2 - y1) * x - (x2 - x1) * y + x2 * y1 - y2 * x1) /
    max(denom, .Machine$double.eps)
  dim <- which.max(dist)
  max(min_dim, min(dim, n))
}
