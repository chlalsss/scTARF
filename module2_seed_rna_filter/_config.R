# Shared configuration and helpers for module2_seed_rna_filter scripts.
#
# This module implements direction 1 from 27_module2_post-fail_plan:
# subtract RNA-incompatible cells from the ATAC z-score top 5% seed set.
# Independent verification only; no full TRS run.

required_pkgs <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Missing packages: ", paste(missing, collapse = ", "), call. = FALSE)
  }
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

parse_csv_env <- function(name, default) {
  raw <- Sys.getenv(name, default)
  parts <- strsplit(raw, ",", fixed = TRUE)[[1]]
  parts <- trimws(parts)
  parts[nzchar(parts)]
}

seed_filter_paths <- function() {
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
  out_root <- Sys.getenv(
    "MODULE2_SEED_FILTER_OUT_DIR",
    file.path(data_dir, "module2_seed_rna_filter")
  )
  original_dir <- Sys.getenv(
    "MODULE2_ORIGINAL_SCAVENGE_DIR",
    file.path(data_dir, "scavenge_out_small_all_traits")
  )

  list(
    data_dir = data_dir,
    scsheft_dir = scsheft_dir,
    scsheft_base = scsheft_base,
    out_dir = out_root,
    fig_dir = file.path(out_root, "figures"),
    rna_cells = file.path(scsheft_base, "scsheft_input", "rna_cells.csv"),
    atac_cells = file.path(scsheft_base, "scsheft_input", "atac_cells.csv"),
    topk_neighbors_csv = file.path(scsheft_dir, "atac_topk_rna_neighbors.csv"),
    topk_neighbors_rds = file.path(scsheft_dir, "atac_topk_rna_neighbors.rds"),
    seed_signal_by_cell = Sys.getenv(
      "MODULE2_SEED_SIGNAL_BY_CELL",
      file.path(data_dir, "module2_gwas_gene_risk", "seed_signal_diagnostics", "seed_signal_by_cell.csv.gz")
    ),
    atac_se = Sys.getenv(
      "MODULE2_ATAC_SE",
      file.path(data_dir, "scATAC-Healthy-Hematopoiesis-191120.rds")
    ),
    original_knn = Sys.getenv(
      "MODULE2_ORIGINAL_KNN",
      file.path(original_dir, "mutualknn_k30.rds")
    ),
    original_dir = original_dir
  )
}

ensure_seed_filter_dirs <- function(paths) {
  dir.create(paths$out_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(paths$fig_dir, showWarnings = FALSE, recursive = TRUE)
}

# Trait -> target cell type rule. Identical to module2_wnn_lite_gated to keep
# acceptance criteria consistent across modules.

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

as_binary_graph <- function(adj) {
  bin <- methods::as(adj, "dgCMatrix")
  bin@x <- rep(1, length(bin@x))
  bin
}

# Load topk RNA neighbors. Two on-disk formats are accepted:
#   .csv  : long-format with columns atac_cell, rank, rna_cell, distance
#   .rds  : named list, names = atac_cell, each element data.frame(rna_cell, distance)
#           (this is the format produced by module1_method_{A,B,C}_*.R)
# Both are normalized to a data.table with columns atac_cell, rank, rna_cell, distance.

load_topk_neighbors <- function(paths, k_max = NULL) {
  source_used <- NA_character_
  if (file.exists(paths$topk_neighbors_csv)) {
    nb <- data.table::fread(paths$topk_neighbors_csv)
    source_used <- paths$topk_neighbors_csv
  } else if (file.exists(paths$topk_neighbors_rds)) {
    raw <- readRDS(paths$topk_neighbors_rds)
    source_used <- paths$topk_neighbors_rds
    if (is.list(raw) && !is.data.frame(raw)) {
      atac_ids <- names(raw)
      if (is.null(atac_ids)) {
        stop("topk neighbors .rds is an unnamed list; cannot recover atac_cell ids", call. = FALSE)
      }
      pieces <- vector("list", length(raw))
      for (i in seq_along(raw)) {
        df <- raw[[i]]
        if (!all(c("rna_cell", "distance") %in% names(df))) {
          stop("topk neighbors .rds element ", i, " missing rna_cell/distance columns", call. = FALSE)
        }
        pieces[[i]] <- data.table::data.table(
          atac_cell = atac_ids[i],
          rank = seq_len(nrow(df)),
          rna_cell = as.character(df$rna_cell),
          distance = as.numeric(df$distance)
        )
      }
      nb <- data.table::rbindlist(pieces, use.names = TRUE)
    } else {
      nb <- data.table::as.data.table(raw)
    }
  } else {
    stop("Missing topk neighbors at: ", paths$topk_neighbors_csv, " or ", paths$topk_neighbors_rds, call. = FALSE)
  }

  required <- c("atac_cell", "rank", "rna_cell")
  missing_cols <- setdiff(required, names(nb))
  if (length(missing_cols) > 0) {
    stop("topk neighbors missing columns after normalization: ",
         paste(missing_cols, collapse = ", "),
         " (source: ", source_used, ")", call. = FALSE)
  }
  if (!is.null(k_max)) {
    nb <- nb[nb$rank <= k_max, ]
  }
  attr(nb, "source") <- source_used
  nb
}
