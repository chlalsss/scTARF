# Shared configuration for module 2 adaptive fusion scripts.

required_pkgs <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Missing packages: ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

cutoff_label <- function(x) {
  gsub("\\.", "p", sprintf("%.2f", as.numeric(x)))
}

get_cutoffs <- function(default = c(0.50, 0.60)) {
  value <- Sys.getenv("MODULE2_CUTOFFS", "")
  if (!nzchar(value)) return(default)
  as.numeric(strsplit(value, ",", fixed = TRUE)[[1]])
}

get_single_cutoff <- function(default = 0.60) {
  as.numeric(Sys.getenv("MODULE2_CUTOFF", as.character(default)))
}

module2_paths <- function() {
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
  out_dir <- Sys.getenv(
    "MODULE2_OUT_DIR",
    file.path(data_dir, "module2_adaptive_fusion")
  )

  list(
    data_dir = data_dir,
    scsheft_base = scsheft_base,
    scsheft_dir = scsheft_dir,
    shared_embedding = file.path(scsheft_dir, "shared_embedding.rds"),
    rna_cells = file.path(scsheft_base, "scsheft_input", "rna_cells.csv"),
    atac_cells = file.path(scsheft_base, "scsheft_input", "atac_cells.csv"),
    out_dir = out_dir,
    fig_dir = file.path(out_dir, "figures")
  )
}

ensure_module2_dirs <- function(paths) {
  dir.create(paths$out_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(paths$fig_dir, showWarnings = FALSE, recursive = TRUE)
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

safe_scale_delta <- function(delta) {
  scale <- stats::median(delta, na.rm = TRUE)
  if (!is.finite(scale) || scale <= 1e-8) {
    positives <- delta[is.finite(delta) & delta > 1e-8]
    if (length(positives) > 0) scale <- stats::median(positives)
  }
  if (!is.finite(scale) || scale <= 1e-8) scale <- 1e-8
  scale
}

get_active_label <- function() {
  forced <- Sys.getenv("MODULE2_CUTOFF_LABEL", "")
  if (nzchar(forced)) return(forced)
  topk <- as.integer(Sys.getenv("MODULE2_TOPK", "1"))
  if (topk == 1) return("top1")
  cutoff_label(get_single_cutoff())
}

auto_select_pca_dim <- function(sdev, method = NULL, cumvar_threshold = NULL, min_dim = 2L) {
  vals <- as.numeric(sdev)^2
  vals <- vals[is.finite(vals) & vals > 0]
  n <- length(vals)
  if (n == 0) stop("Cannot auto-select PCA dimension from empty sdev")
  if (n <= min_dim) return(n)

  if (is.null(method)) {
    method <- tolower(Sys.getenv("MODULE2_PCA_AUTO_METHOD", "cumvar"))
  }
  if (is.null(cumvar_threshold)) {
    cumvar_threshold <- as.numeric(Sys.getenv("MODULE2_PCA_CUMVAR", "0.80"))
  }

  if (method == "cumvar") {
    cumvar <- cumsum(vals) / sum(vals)
    dim <- which(cumvar >= cumvar_threshold)[1]
    if (is.na(dim)) dim <- n
    return(max(min_dim, min(dim, n)))
  }

  # Kneedle-style scree elbow: max distance from the line joining endpoints.
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

parse_int_sequence_env <- function(name, default) {
  value <- Sys.getenv(name, "")
  if (!nzchar(value)) return(default)
  if (grepl(":", value, fixed = TRUE)) {
    bounds <- as.integer(strsplit(value, ":", fixed = TRUE)[[1]])
    if (length(bounds) != 2 || anyNA(bounds)) {
      stop(name, " must be a comma list or a:b range")
    }
    return(seq(bounds[1], bounds[2]))
  }
  parsed <- as.integer(strsplit(value, ",", fixed = TRUE)[[1]])
  if (anyNA(parsed)) stop(name, " must be a comma list or a:b range")
  parsed
}

auto_select_knn_k <- function(mat, target_iso = 0.01, k_search = 10:50) {
  if (!requireNamespace("RANN", quietly = TRUE)) {
    stop("Package RANN is required for auto_select_knn_k")
  }
  k_search <- sort(unique(as.integer(k_search)))
  k_search <- k_search[is.finite(k_search) & k_search >= 1]
  if (length(k_search) == 0) stop("k_search must contain positive integers")

  n <- nrow(mat)
  k_max <- min(max(k_search), n - 1L)
  k_search <- k_search[k_search <= k_max]
  if (length(k_search) == 0) stop("k_search exceeds available cells")

  nn <- RANN::nn2(data = mat, query = mat, k = k_max + 1L)
  isolated_at_k <- function(k) {
    from <- rep(seq_len(n), each = k)
    to <- integer(n * k)
    pos <- 1L
    for (i in seq_len(n)) {
      idx_i <- nn$nn.idx[i, ]
      idx_i <- idx_i[idx_i != i]
      to[pos:(pos + k - 1L)] <- idx_i[seq_len(k)]
      pos <- pos + k
    }
    keys <- paste(pmin(from, to), pmax(from, to), sep = "_")
    mutual_keys <- names(table(keys))[table(keys) >= 2L]
    if (length(mutual_keys) == 0) return(1)
    parts <- strsplit(mutual_keys, "_", fixed = TRUE)
    touched <- unique(as.integer(unlist(parts, use.names = FALSE)))
    1 - length(touched) / n
  }

  pct_isolated <- vapply(k_search, isolated_at_k, numeric(1))
  pass <- which(pct_isolated < target_iso)
  selected <- if (length(pass) > 0) k_search[pass[1]] else k_search[which.min(pct_isolated)]
  list(
    k = selected,
    target_iso = target_iso,
    search = data.frame(k = k_search, pct_isolated = pct_isolated)
  )
}

write_auto_value <- function(path, value) {
  writeLines(as.character(value), con = path, useBytes = TRUE)
}

read_auto_value <- function(path, what) {
  if (!file.exists(path)) stop("Missing auto-selected ", what, " file: ", path)
  value <- as.integer(readLines(path, warn = FALSE)[1])
  if (!is.finite(value) || value < 1) stop("Invalid auto-selected ", what, " in: ", path)
  value
}

