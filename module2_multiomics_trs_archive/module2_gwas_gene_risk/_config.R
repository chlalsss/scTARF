# Shared configuration and helpers for module2_gwas_gene_risk scripts.

required_pkgs <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Missing packages: ", paste(missing, collapse = ", "), call. = FALSE)
  }
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

cutoff_label <- function(x) {
  gsub("\\.", "p", sprintf("%.2f", as.numeric(x)))
}

gwas_gene_paths <- function() {
  data_dir <- Sys.getenv(
    "MODULE2_DATA_DIR",
    "D:/SCAVENGE-reproducibility-main/blood_data"
  )
  out_dir <- Sys.getenv(
    "MODULE2_GWAS_GENE_OUT_DIR",
    file.path(data_dir, "module2_gwas_gene_risk")
  )

  list(
    data_dir = data_dir,
    out_dir = out_dir,
    fig_dir = file.path(out_dir, "figures"),
    atac_se = Sys.getenv(
      "MODULE2_ATAC_SE",
      file.path(data_dir, "scATAC-Healthy-Hematopoiesis-191120.rds")
    ),
    gwas = Sys.getenv(
      "MODULE2_GWAS",
      file.path(data_dir, "CS.trans_all.PP01.strong_atac_overlap.tsv")
    ),
    original_scavenge_dir = Sys.getenv(
      "MODULE2_ORIGINAL_SCAVENGE_DIR",
      file.path(data_dir, "scavenge_out_small_all_traits")
    ),
    bg_peaks = Sys.getenv(
      "MODULE2_BG_PEAKS",
      file.path(data_dir, "scavenge_out_small_all_traits", "bg_peaks_n200.rds")
    ),
    original_knn = Sys.getenv(
      "MODULE2_ORIGINAL_KNN",
      file.path(data_dir, "scavenge_out_small_all_traits", "mutualknn_k30.rds")
    ),
    wnn_lite_knn = Sys.getenv(
      "MODULE2_WNN_LITE_KNN",
      file.path(data_dir, "module2_wnn_lite_fusion", "wnn_lite_mutualknn_cutoff_0p60_k30.rds")
    ),
    rna_genes = Sys.getenv(
      "MODULE2_RNA_GENES",
      file.path(data_dir, "multiomics_output_v2", "glue_input", "rna_genes.csv")
    ),
    atac_peaks = Sys.getenv(
      "MODULE2_ATAC_PEAKS",
      file.path(data_dir, "multiomics_output_v2", "glue_input", "atac_peaks.csv")
    ),
    gene_activity = Sys.getenv(
      "MODULE2_GENE_ACTIVITY",
      file.path(data_dir, "multiomics_output_v2", "gene_activity_matrix.rds")
    )
  )
}

ensure_gwas_gene_dirs <- function(paths) {
  dir.create(paths$out_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(paths$fig_dir, showWarnings = FALSE, recursive = TRUE)
}

parse_var_to_bed <- function(v) {
  v <- as.character(v)
  m <- regexec("^([0-9]+|X|Y|MT|M):(\\d+)_([^_]+)_([^_]+)$", v, ignore.case = TRUE)
  r <- regmatches(v, m)
  ok <- lengths(r) == 5
  chr <- start <- end <- rep(NA_character_, length(v))

  if (any(ok)) {
    raw_chr <- toupper(vapply(r[ok], `[[`, "", 2))
    raw_chr[raw_chr == "M"] <- "MT"
    chr[ok] <- paste0("chr", raw_chr)
    pos <- as.integer(vapply(r[ok], `[[`, "", 3))
    start[ok] <- as.character(pmax(pos - 1L, 0L))
    end[ok] <- as.character(pos)
  }

  data.table::data.table(
    chr = chr,
    start = as.integer(start),
    end = as.integer(end),
    ok = !is.na(chr) & !is.na(start) & !is.na(end)
  )
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
  bin <- as(adj, "dgCMatrix")
  bin@x <- rep(1, length(bin@x))
  bin
}
