# Shared configuration for module2_adaptive_gamma.
#
# Direction 2 from 27_module2_post-fail_plan: per-cell heterogeneous restart
# probability gamma_i in the SCAVENGE RWR step (E step), instead of the
# original gamma = 0.05 uniform.
#
# Two ablation versions:
#   Version A (discrete) : gamma_high on progenitor/transition cells,
#                          gamma_low elsewhere
#   Version B (continuous): gamma_i is a linear function of (1 - rank(rna_reliability))
#                          so low-reliability cells get high gamma.

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

adaptive_gamma_paths <- function() {
  data_dir <- Sys.getenv(
    "MODULE2_DATA_DIR",
    "D:/SCAVENGE-reproducibility-main/blood_data"
  )
  out_root <- Sys.getenv(
    "MODULE2_ADAPTIVE_GAMMA_OUT_DIR",
    file.path(data_dir, "module2_adaptive_gamma")
  )
  original_dir <- Sys.getenv(
    "MODULE2_ORIGINAL_SCAVENGE_DIR",
    file.path(data_dir, "scavenge_out_small_all_traits")
  )
  wnn_lite_dir <- Sys.getenv(
    "MODULE2_WNN_LITE_OUT_DIR",
    file.path(data_dir, "module2_wnn_lite_fusion")
  )
  seed_filter_dir <- Sys.getenv(
    "MODULE2_SEED_FILTER_OUT_DIR",
    file.path(data_dir, "module2_seed_rna_filter")
  )

  list(
    data_dir = data_dir,
    out_dir = out_root,
    fig_dir = file.path(out_root, "figures"),
    original_dir = original_dir,
    wnn_lite_dir = wnn_lite_dir,
    seed_filter_dir = seed_filter_dir,
    atac_se = Sys.getenv(
      "MODULE2_ATAC_SE",
      file.path(data_dir, "scATAC-Healthy-Hematopoiesis-191120.rds")
    ),
    original_knn = Sys.getenv(
      "MODULE2_ORIGINAL_KNN",
      file.path(original_dir, "mutualknn_k30.rds")
    ),
    seed_signal_by_cell = Sys.getenv(
      "MODULE2_SEED_SIGNAL_BY_CELL",
      file.path(data_dir, "module2_gwas_gene_risk", "seed_signal_diagnostics", "seed_signal_by_cell.csv.gz")
    ),
    wnn_lite_weights = Sys.getenv(
      "MODULE2_WNN_LITE_BASE_WEIGHTS",
      file.path(wnn_lite_dir, "wnn_lite_weights_cutoff_0p60.csv")
    ),
    direction1_filtered_seed_template = file.path(
      seed_filter_dir,
      "filtered_seed_cells_%s_t%s.csv"
    ),
    direction1_trs_dir_template = file.path(
      seed_filter_dir,
      "trs_filtered_t%s"
    )
  )
}

ensure_adaptive_gamma_dirs <- function(paths) {
  dir.create(paths$out_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(paths$fig_dir, showWarnings = FALSE, recursive = TRUE)
}

# Cell type rules (consistent with module2_wnn_lite_gated)

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

as_binary_graph <- function(adj) {
  bin <- methods::as(adj, "dgCMatrix")
  bin@x <- rep(1, length(bin@x))
  bin
}

# Threshold label for direction 1 filename
seed_filter_label <- function(threshold) {
  gsub("\\.", "p", sprintf("%.2f", as.numeric(threshold)))
}
