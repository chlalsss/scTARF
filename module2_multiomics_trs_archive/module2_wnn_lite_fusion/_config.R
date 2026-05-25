# Shared configuration and helpers for module 2 WNN-lite fusion scripts.

required_pkgs <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Missing packages: ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

cutoff_label <- function(x) {
  gsub("\\.", "p", sprintf("%.2f", as.numeric(x)))
}

get_single_cutoff <- function(default = 0.60) {
  as.numeric(Sys.getenv("MODULE2_CUTOFF", as.character(default)))
}

env_flag <- function(name, default = FALSE) {
  value <- Sys.getenv(name, if (default) "true" else "false")
  tolower(value) %in% c("true", "1", "yes", "on")
}

env_numeric <- function(name, default) {
  value <- as.numeric(Sys.getenv(name, as.character(default)))
  if (!is.finite(value)) {
    stop(name, " must be numeric", call. = FALSE)
  }
  value
}

env_integer <- function(name, default) {
  value <- as.integer(Sys.getenv(name, as.character(default)))
  if (!is.finite(value)) {
    stop(name, " must be an integer", call. = FALSE)
  }
  value
}

wnn_lite_paths <- function() {
  data_dir <- Sys.getenv(
    "MODULE2_DATA_DIR",
    "D:/SCAVENGE-reproducibility-main/blood_data"
  )
  scsheft_dir <- Sys.getenv(
    "MODULE2_SCSHEFT_DIR",
    file.path(
      data_dir,
      "multiomics_module1_pairing",
      "method_B_scSHEFT",
      "server_full_official_anchorfix_L1mean_20260424"
    )
  )
  scsheft_base <- dirname(scsheft_dir)
  adaptive_dir <- Sys.getenv(
    "MODULE2_ADAPTIVE_OUT_DIR",
    file.path(data_dir, "module2_adaptive_fusion")
  )
  out_dir <- Sys.getenv(
    "MODULE2_WNN_LITE_OUT_DIR",
    file.path(data_dir, "module2_wnn_lite_fusion")
  )

  list(
    data_dir = data_dir,
    scsheft_base = scsheft_base,
    scsheft_dir = scsheft_dir,
    adaptive_dir = adaptive_dir,
    out_dir = out_dir,
    fig_dir = file.path(out_dir, "figures"),
    shared_embedding = file.path(scsheft_dir, "shared_embedding.rds"),
    rna_cells = file.path(scsheft_base, "scsheft_input", "rna_cells.csv"),
    atac_cells = file.path(scsheft_base, "scsheft_input", "atac_cells.csv")
  )
}

ensure_wnn_lite_dirs <- function(paths) {
  dir.create(paths$out_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(paths$fig_dir, showWarnings = FALSE, recursive = TRUE)
}

adaptive_by_cell_path <- function(paths, cutoff) {
  file.path(
    paths$adaptive_dir,
    sprintf("adaptive_selection_by_cell_cutoff_%s.csv", cutoff_label(cutoff))
  )
}

rna_reference_path <- function(paths, cutoff) {
  file.path(
    paths$adaptive_dir,
    sprintf("rna_reference_latent_cutoff_%s.rds", cutoff_label(cutoff))
  )
}

wnn_weights_path <- function(paths, cutoff) {
  file.path(
    paths$out_dir,
    sprintf("wnn_lite_weights_cutoff_%s.csv", cutoff_label(cutoff))
  )
}

load_shared_embedding <- function(shared_path) {
  if (!file.exists(shared_path)) {
    stop("Missing shared embedding: ", shared_path, call. = FALSE)
  }
  emb <- readRDS(shared_path)
  latent_cols <- grep("^latent_", colnames(emb), value = TRUE)
  if (!all(c("cell", "modality") %in% colnames(emb)) || length(latent_cols) == 0) {
    stop("shared_embedding.rds must contain cell, modality, and latent_* columns", call. = FALSE)
  }
  list(emb = emb, latent_cols = latent_cols)
}

scale01 <- function(x, lower_q = 0, upper_q = 1, fallback = 0.5) {
  x <- as.numeric(x)
  ok <- is.finite(x)
  out <- rep(NA_real_, length(x))
  if (!any(ok)) {
    out[] <- fallback
    return(out)
  }

  lower <- as.numeric(stats::quantile(x[ok], lower_q, na.rm = TRUE, names = FALSE))
  upper <- as.numeric(stats::quantile(x[ok], upper_q, na.rm = TRUE, names = FALSE))
  if (!is.finite(lower) || !is.finite(upper) || upper <= lower) {
    out[ok] <- fallback
  } else {
    out[ok] <- (x[ok] - lower) / (upper - lower)
    out[ok] <- pmin(pmax(out[ok], 0), 1)
  }
  out[!ok] <- fallback
  out
}

l2_normalize_rows <- function(mat) {
  norms <- sqrt(rowSums(mat^2))
  norms[!is.finite(norms) | norms <= 0] <- 1
  mat / norms
}

safe_delta_scale <- function(distances) {
  d_min <- min(distances, na.rm = TRUE)
  delta <- distances - d_min
  scale <- stats::median(delta, na.rm = TRUE)
  if (!is.finite(scale) || scale <= 1e-8) {
    positives <- delta[is.finite(delta) & delta > 1e-8]
    if (length(positives) > 0) scale <- stats::median(positives)
  }
  if (!is.finite(scale) || scale <= 1e-8) scale <- 1e-8
  scale
}

distance_to_similarity <- function(distances) {
  d_min <- min(distances, na.rm = TRUE)
  scale <- safe_delta_scale(distances)
  exp(-(distances - d_min) / scale)
}

candidate_distances <- function(mat, row_norm2, i, idx) {
  sims <- as.vector(mat[idx, , drop = FALSE] %*% mat[i, ])
  d2 <- row_norm2[idx] + row_norm2[i] - 2 * sims
  sqrt(pmax(d2, 0))
}

strip_self_neighbors <- function(nn_idx, nn_dist, k) {
  n <- nrow(nn_idx)
  out_idx <- matrix(NA_integer_, nrow = n, ncol = k)
  out_dist <- matrix(NA_real_, nrow = n, ncol = k)
  for (i in seq_len(n)) {
    keep <- nn_idx[i, ] != i
    idx_i <- nn_idx[i, keep]
    dist_i <- nn_dist[i, keep]
    if (length(idx_i) < k) {
      stop("Could not find ", k, " non-self neighbors for row ", i, call. = FALSE)
    }
    out_idx[i, ] <- idx_i[seq_len(k)]
    out_dist[i, ] <- dist_i[seq_len(k)]
  }
  list(idx = out_idx, dist = out_dist)
}

load_atac_cell_types <- function(paths) {
  if (!file.exists(paths$atac_cells)) return(NULL)
  dt <- data.table::fread(paths$atac_cells)
  if (!all(c("cell", "cell_type") %in% colnames(dt))) return(NULL)
  setNames(as.character(dt$cell_type), dt$cell)
}

