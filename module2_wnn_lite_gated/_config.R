# Shared configuration and helpers for module2_wnn_lite_gated scripts.

required_pkgs <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Missing packages: ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

cutoff_label <- function(x) {
  gsub("\\.", "p", sprintf("%.2f", as.numeric(x)))
}

env_numeric <- function(name, default) {
  value <- as.numeric(Sys.getenv(name, as.character(default)))
  if (!is.finite(value)) stop(name, " must be numeric", call. = FALSE)
  value
}

env_integer <- function(name, default) {
  value <- as.integer(Sys.getenv(name, as.character(default)))
  if (!is.finite(value)) stop(name, " must be an integer", call. = FALSE)
  value
}

env_flag <- function(name, default = FALSE) {
  value <- Sys.getenv(name, if (default) "true" else "false")
  tolower(value) %in% c("true", "1", "yes", "on")
}

get_single_cutoff <- function(default = 0.60) {
  env_numeric("MODULE2_CUTOFF", default)
}

gated_paths <- function(version = NULL) {
  data_dir <- Sys.getenv(
    "MODULE2_DATA_DIR",
    "D:/SCAVENGE-reproducibility-main/blood_data"
  )
  ungated_dir <- Sys.getenv(
    "MODULE2_WNN_LITE_OUT_DIR",
    file.path(data_dir, "module2_wnn_lite_fusion")
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
  gated_root <- Sys.getenv(
    "MODULE2_WNN_LITE_GATED_OUT_DIR",
    file.path(data_dir, "module2_wnn_lite_gated")
  )
  original_dir <- Sys.getenv(
    "MODULE2_ORIGINAL_SCAVENGE_DIR",
    file.path(data_dir, "scavenge_out_small_all_traits")
  )

  out_dir <- if (is.null(version)) gated_root else file.path(gated_root, version)
  list(
    data_dir = data_dir,
    ungated_dir = ungated_dir,
    scsheft_base = scsheft_base,
    scsheft_dir = scsheft_dir,
    adaptive_dir = adaptive_dir,
    gated_root = gated_root,
    out_dir = out_dir,
    fig_dir = file.path(out_dir, "figures"),
    shared_embedding = file.path(scsheft_dir, "shared_embedding.rds"),
    rna_cells = file.path(scsheft_base, "scsheft_input", "rna_cells.csv"),
    atac_cells = file.path(scsheft_base, "scsheft_input", "atac_cells.csv"),
    original_dir = original_dir,
    atac_se = Sys.getenv(
      "MODULE2_ATAC_SE",
      file.path(data_dir, "scATAC-Healthy-Hematopoiesis-191120.rds")
    ),
    base_weights = Sys.getenv(
      "MODULE2_WNN_LITE_BASE_WEIGHTS",
      file.path(ungated_dir, "wnn_lite_weights_cutoff_0p60.csv")
    ),
    base_same_type_by_celltype = Sys.getenv(
      "MODULE2_WNN_LITE_SAME_TYPE_BY_CELLTYPE",
      file.path(ungated_dir, "wnn_lite_same_type_by_celltype_cutoff_0p60_k30.csv")
    ),
    seed_signal_by_cell = Sys.getenv(
      "MODULE2_SEED_SIGNAL_BY_CELL",
      file.path(data_dir, "module2_gwas_gene_risk", "seed_signal_diagnostics", "seed_signal_by_cell.csv.gz")
    ),
    original_knn = Sys.getenv(
      "MODULE2_ORIGINAL_KNN",
      file.path(original_dir, "mutualknn_k30.rds")
    ),
    original_trs = Sys.getenv(
      "MODULE2_ORIGINAL_TRS",
      file.path(original_dir, "TRS_matrix_small_all_traits_fullcells.rds")
    )
  )
}

ensure_gated_dirs <- function(paths) {
  dir.create(paths$gated_root, showWarnings = FALSE, recursive = TRUE)
  dir.create(paths$out_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(paths$fig_dir, showWarnings = FALSE, recursive = TRUE)
}

rna_reference_path <- function(paths, cutoff) {
  file.path(
    paths$adaptive_dir,
    sprintf("rna_reference_latent_cutoff_%s.rds", cutoff_label(cutoff))
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

as_binary_graph <- function(adj) {
  bin <- as(adj, "dgCMatrix")
  bin@x <- rep(1, length(bin@x))
  bin
}

is_transition_celltype <- function(cell_type) {
  grepl("HSC|CMP|LMPP|CLP|GMP|Pre|Early|Unk", as.character(cell_type), ignore.case = TRUE)
}

target_definition <- function(trait) {
  if (trait == "MONO") {
    "BioClassification contains Mono"
  } else if (trait == "PLT") {
    "BioClassification contains Eryth or Baso; proxy for platelet/megakaryocytic region"
  } else if (trait == "LYM") {
    "BioClassification contains CD4, CD8, B, NK, CLP, or Pre.B"
  } else {
    "No target rule defined"
  }
}

is_target_type <- function(trait, celltype) {
  if (length(trait) > 1) {
    return(as.logical(mapply(is_target_type, trait, celltype, USE.NAMES = FALSE)))
  }
  ct <- as.character(celltype)
  if (trait == "MONO") {
    grepl("Mono", ct, ignore.case = TRUE)
  } else if (trait == "PLT") {
    grepl("Eryth|Baso", ct, ignore.case = TRUE)
  } else if (trait == "LYM") {
    grepl("CD4|CD8|\\bB\\b|NK|CLP|Pre\\.B", ct, ignore.case = TRUE)
  } else {
    rep(FALSE, length(ct))
  }
}

load_atac_metadata <- function(atac_se_path) {
  if (!file.exists(atac_se_path)) stop("Missing ATAC SE: ", atac_se_path, call. = FALSE)
  se <- readRDS(atac_se_path)
  meta <- data.table::as.data.table(as.data.frame(SummarizedExperiment::colData(se)))
  meta[, cell := colnames(se)]
  if (!"BioClassification" %in% names(meta)) {
    stop("ATAC metadata must contain BioClassification", call. = FALSE)
  }
  meta[, .(cell, BioClassification)]
}
