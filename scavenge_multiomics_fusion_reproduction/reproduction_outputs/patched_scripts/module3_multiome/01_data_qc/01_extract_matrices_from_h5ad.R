# 01_data_qc/01_extract_matrices_from_h5ad.R
#
# Pull ATAC peak counts + RNA gene counts out of the GSE194122 multiome .h5ad
# into R-native sparse matrices, then write them as RDS for downstream use.
#
# OpenProblems multiome convention (from the dataset documentation):
#   - adata.X holds binarized ATAC + raw RNA in a single (n_cells, n_features)
#     matrix; the modality of each feature is in /var/feature_types
#   - adata.layers['counts'] may carry the unmodified raw counts for both
# We split features by feature_types and emit separate ATAC/RNA matrices.
#
# Output (under <out_dir>/01_data_qc/):
#   BMMC_atac_counts.rds       dgCMatrix [n_cells x n_atac_peaks], binarized
#   BMMC_rna_counts.rds        dgCMatrix [n_cells x n_rna_genes],   raw counts
#   BMMC_cell_meta.rds         data.table, one row per cell
#   BMMC_atac_features.rds     data.table peak_id, chr, start, end
#   BMMC_rna_features.rds      data.table gene_id, gene_symbol
#
# Memory note: the multiome .h5ad's /X is ~70k cells x 130k features sparse;
# R will need a few GB free during the readout. We avoid materialising any
# dense slice — everything stays sparse.

source("D:/RStudio/work/scavenge_multiomics_fusion_reproduction/reproduction_outputs/patched_scripts/module3_multiome/_bootstrap.R")
required_pkgs(c("rhdf5", "Matrix", "data.table"))

paths <- m3_paths()
out_dir <- file.path(paths$out_dir, "01_data_qc")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Locate h5ad
if (paths$bmmc_source == "GSE194122") {
  h5_path <- file.path(
    paths$out_dir, "00_precheck", "_geo_cache", paths$geo_accession,
    "GSE194122_openproblems_neurips2021_multiome_BMMC_processed.h5ad"
  )
} else if (paths$bmmc_source == "figshare_s4d8") {
  h5_path <- paths$figshare_h5ad
} else {
  stop("Unknown M3_BMMC_SOURCE: ", paths$bmmc_source)
}
if (!file.exists(h5_path)) stop("h5ad not found: ", h5_path)
message("Reading: ", h5_path)

ls_full <- rhdf5::h5ls(h5_path, recursive = TRUE)

# --- 1. cell metadata -----------------------------------------------------

if (!file.exists(paths$meta_raw)) {
  stop("Missing ", paths$meta_raw,
       " — run 00_precheck/01_fetch_metadata.R first.")
}
meta <- data.table::fread(paths$meta_raw)
n_cells <- nrow(meta)
message("Cells: ", n_cells)

# --- 2. /var: feature types + names ---------------------------------------

var_top <- ls_full[ls_full$group == "/var", ]

legacy_var_cats <- list()
if ("__categories" %in% var_top$name) {
  legacy_grp <- "/var/__categories"
  legacy_items <- ls_full[ls_full$group == legacy_grp, "name"]
  for (col in legacy_items) {
    legacy_var_cats[[col]] <- as.character(
      rhdf5::h5read(h5_path, paste0(legacy_grp, "/", col))
    )
  }
}

read_var_col <- function(col_name) {
  info <- var_top[var_top$name == col_name, ]
  if (nrow(info) == 0) return(NULL)
  otype <- info$otype[1]
  path <- paste0("/var/", col_name)
  if (grepl("DATASET", otype)) {
    raw <- as.vector(rhdf5::h5read(h5_path, path))
    if (col_name %in% names(legacy_var_cats)) {
      codes <- as.integer(raw)
      return(ifelse(codes < 0L, NA_character_,
                    legacy_var_cats[[col_name]][codes + 1L]))
    }
    raw
  } else if (grepl("GROUP", otype)) {
    sub_items <- ls_full[ls_full$group == path, "name"]
    if (all(c("codes", "categories") %in% sub_items)) {
      cats  <- as.character(rhdf5::h5read(h5_path, paste0(path, "/categories")))
      codes <- as.integer(rhdf5::h5read(h5_path, paste0(path, "/codes")))
      return(ifelse(codes < 0L, NA_character_, cats[codes + 1L]))
    }
    NULL
  }
}

idx_name <- intersect(c("_index", "index"), var_top$name)[1]
feat_names <- as.character(rhdf5::h5read(h5_path, paste0("/var/", idx_name)))
ftype <- read_var_col("feature_types")
gene_id <- if ("gene_id" %in% var_top$name) read_var_col("gene_id") else NA
n_feat <- length(feat_names)
message("Features total: ", n_feat,
        " | ATAC: ", sum(ftype == "ATAC"),
        " | GEX: ",  sum(ftype == "GEX"))

atac_idx <- which(ftype == "ATAC")
rna_idx  <- which(ftype == "GEX")

# --- 3. /X: load full sparse matrix, then column-slice by modality --------

x_top <- ls_full[ls_full$group == "/X", ]
if (nrow(x_top) == 0) stop("No /X group; this script expects sparse-X anndata.")
message("/X children: ", paste(x_top$name, collapse = ", "))

x_attrs <- tryCatch(rhdf5::h5readAttributes(h5_path, "/X"), error = function(e) list())
encoding_type <- x_attrs[["encoding-type"]]
shape <- x_attrs[["shape"]]
if (is.null(shape) && "h5sparse_shape" %in% names(x_attrs)) shape <- x_attrs[["h5sparse_shape"]]
message("/X encoding-type: ",
        if (is.null(encoding_type)) "(none)" else encoding_type,
        " | shape: ",
        if (is.null(shape)) "(none)" else paste(shape, collapse = " x "))

x_data    <- rhdf5::h5read(h5_path, "/X/data")
x_indices <- rhdf5::h5read(h5_path, "/X/indices")
x_indptr  <- rhdf5::h5read(h5_path, "/X/indptr")

# Decide CSR vs CSC. anndata default is CSR (rows = cells). Detect by length:
#   CSR: indptr length = n_cells + 1
#   CSC: indptr length = n_features + 1
indptr_len <- length(x_indptr)
if (indptr_len == n_cells + 1L) {
  fmt <- "CSR"
} else if (indptr_len == n_feat + 1L) {
  fmt <- "CSC"
} else {
  stop("indptr length ", indptr_len, " does not match cells (", n_cells,
       ") or features (", n_feat, ")")
}
message("/X format: ", fmt)

# Build a dgCMatrix (cells = rows). CSR with n_cells rows naturally makes
# dgRMatrix; we transpose to dgCMatrix for downstream Matrix/SCAVENGE friendlier
# operation (cells as rows, features as cols).
if (fmt == "CSR") {
  # CSR with rows=cells, cols=features. To get dgCMatrix [cells x features]
  # we use sparseMatrix with j = indices (col indices), p = indptr (row ptr),
  # x = data, then convert via t(t(.)) — but Matrix has a direct route:
  # construct dgRMatrix then coerce.
  X <- Matrix::sparseMatrix(
    j = as.integer(x_indices) + 1L,
    p = as.integer(x_indptr),
    x = as.numeric(x_data),
    dims = c(n_cells, n_feat),
    index1 = TRUE,
    repr = "R"
  )
  X <- as(X, "CsparseMatrix")
} else {  # CSC: rows=features, cols=cells (10x convention if from raw)
  X <- Matrix::sparseMatrix(
    i = as.integer(x_indices) + 1L,
    p = as.integer(x_indptr),
    x = as.numeric(x_data),
    dims = c(n_feat, n_cells),
    index1 = TRUE,
    repr = "C"
  )
  X <- Matrix::t(X)  # -> [cells x features]
}
message("X matrix built: ", nrow(X), " x ", ncol(X), ", nnz = ", length(X@x))

# --- 4. split by modality + write -----------------------------------------

atac_mat <- X[, atac_idx, drop = FALSE]
rna_mat  <- X[, rna_idx,  drop = FALSE]

rownames(atac_mat) <- meta$barcode
rownames(rna_mat)  <- meta$barcode
colnames(atac_mat) <- feat_names[atac_idx]
colnames(rna_mat)  <- feat_names[rna_idx]

# Per OpenProblems documentation, ATAC counts in adata.X are binarized.
# Force-binarize anyway so downstream can rely on it.
if (any(atac_mat@x > 1)) {
  message("ATAC matrix had non-binary values; binarising.")
  atac_mat@x[atac_mat@x > 0] <- 1
}

saveRDS(atac_mat, file.path(out_dir, "BMMC_atac_counts.rds"))
saveRDS(rna_mat,  file.path(out_dir, "BMMC_rna_counts.rds"))
message("Wrote: ", file.path(out_dir, "BMMC_atac_counts.rds"))
message("Wrote: ", file.path(out_dir, "BMMC_rna_counts.rds"))

# Feature tables — use shared parser from _config.R (handles multiple separators)
atac_features <- parse_peak_names(feat_names[atac_idx])
rna_features  <- data.table::data.table(
  gene_id     = if (length(gene_id) == n_feat) gene_id[rna_idx] else NA,
  gene_symbol = feat_names[rna_idx]
)
saveRDS(atac_features, file.path(out_dir, "BMMC_atac_features.rds"))
saveRDS(rna_features,  file.path(out_dir, "BMMC_rna_features.rds"))

saveRDS(meta, file.path(out_dir, "BMMC_cell_meta.rds"))
message("Wrote: ", file.path(out_dir, "BMMC_cell_meta.rds"))

# --- 5. quick sanity prints ----------------------------------------------

message("\nSanity:")
message("  atac dim: ", nrow(atac_mat), " x ", ncol(atac_mat))
message("  rna  dim: ", nrow(rna_mat),  " x ", ncol(rna_mat))
message("  atac mean nnz per cell: ",
        round(length(atac_mat@x) / nrow(atac_mat), 0))
message("  rna  mean nnz per cell: ",
        round(length(rna_mat@x)  / nrow(rna_mat),  0))
