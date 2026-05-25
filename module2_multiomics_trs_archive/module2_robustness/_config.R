required_pkgs <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Missing packages: ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

robustness_paths <- function() {
  data_dir <- Sys.getenv(
    "MODULE2_DATA_DIR",
    "D:/SCAVENGE-reproducibility-main/blood_data"
  )
  out_dir <- Sys.getenv(
    "MODULE2_ROBUSTNESS_OUT_DIR",
    file.path(data_dir, "module2_robustness")
  )
  adaptive_dir <- file.path(data_dir, "module2_adaptive_fusion")
  list(
    data_dir = data_dir,
    adaptive_dir = adaptive_dir,
    original_atac_se = Sys.getenv(
      "MODULE2_ATAC_SE_ORIGINAL",
      file.path(data_dir, "scATAC-Healthy-Hematopoiesis-191120.rds")
    ),
    out_dir = out_dir,
    sparse_dir = file.path(out_dir, "sparse_inputs"),
    result_dir = file.path(out_dir, "results"),
    figure_dir = file.path(out_dir, "figures")
  )
}

ensure_robustness_dirs <- function(p) {
  dir.create(p$out_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(p$sparse_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(p$result_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(p$figure_dir, showWarnings = FALSE, recursive = TRUE)
}

parse_int_list <- function(value, default) {
  if (!nzchar(value)) return(default)
  as.integer(strsplit(value, ",", fixed = TRUE)[[1]])
}

target_type <- function(trait, celltype) {
  ct <- as.character(celltype)
  if (trait == "MONO") {
    grepl("Mono", ct, ignore.case = TRUE)
  } else if (trait == "PLT") {
    grepl("Eryth|Baso|Mk|Mega|Platelet|PLT", ct, ignore.case = TRUE)
  } else if (trait == "LYM") {
    grepl("CD4|CD8|\\bB\\b|NK|CLP|Pre\\.B", ct, ignore.case = TRUE)
  } else {
    rep(FALSE, length(ct))
  }
}
