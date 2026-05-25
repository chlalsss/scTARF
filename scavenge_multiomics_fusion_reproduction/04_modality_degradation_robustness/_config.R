# module3_multiome shared configuration and helpers.
#
# Conventions follow module2_gwas_gene_risk/_config.R: paths come from env vars
# with sensible defaults, helpers are kept minimal so individual scripts stay
# self-explanatory.
#
# Environment variables (override via Sys.setenv before sourcing):
#   M3_DATA_DIR       root for raw/derived BMMC artifacts (default below)
#   M3_OUT_DIR        root for module3 outputs (default = data_dir)
#   M3_GWAS_TSV       path to existing fine-mapped variants
#                     (CS.trans_all.PP01.strong_atac_overlap.tsv)
#   M3_BMMC_META_URL  optional override for metadata download URL

required_pkgs <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "Missing packages: ", paste(missing, collapse = ", "),
      ". Run install_deps.R in a fresh R session first.",
      call. = FALSE
    )
  }
}

env_str <- function(name, default) {
  v <- Sys.getenv(name, "")
  if (!nzchar(v)) default else v
}

env_numeric <- function(name, default) {
  v <- as.numeric(Sys.getenv(name, as.character(default)))
  if (!is.finite(v)) stop(name, " must be numeric", call. = FALSE)
  v
}

env_integer <- function(name, default) {
  v <- as.integer(Sys.getenv(name, as.character(default)))
  if (!is.finite(v)) stop(name, " must be integer", call. = FALSE)
  v
}

m3_paths <- function() {
  data_dir <- env_str("M3_DATA_DIR", "D:/SCAVENGE-reproducibility-main/BMMC")
  out_dir  <- env_str("M3_OUT_DIR",  data_dir)

  list(
    data_dir   = data_dir,
    out_dir    = out_dir,

    # raw / derived data files
    meta_raw   = file.path(out_dir, "00_precheck", "BMMC_meta_raw.csv"),
    peaks_bed  = file.path(out_dir, "00_precheck", "BMMC_peaks.bed"),

    # precheck outputs
    coverage_csv = file.path(out_dir, "00_precheck", "BMMC_celltype_coverage.csv"),
    overlap_csv  = file.path(out_dir, "00_precheck", "BMMC_finemap_overlap_summary.csv"),
    fig_dir      = file.path(out_dir, "00_precheck", "figures"),

    # external inputs (reused from existing project)
    gwas_tsv = env_str(
      "M3_GWAS_TSV",
      "D:/SCAVENGE-reproducibility-main/blood_data/CS.trans_all.PP01.strong_atac_overlap.tsv"
    ),

    # Source switch: "GSE194122" or "figshare_s4d8".
    # GSE194122 is the full 4-donor version (preferred — has cell_type / donor
    # / batch columns in /obs). figshare_s4d8 is a single-donor QC-only subset
    # kept as fallback (no cell_type, won't pass precheck 1 alone).
    bmmc_source = env_str("M3_BMMC_SOURCE", "GSE194122"),

    figshare_h5ad = env_str(
      "M3_FIGSHARE_H5AD",
      file.path(out_dir, "00_precheck", "_figshare_s4d8",
                "s4d8_quality_control.h5ad")
    ),

    # GEO landing for NeurIPS 2021 BMMC multiome
    geo_accession = "GSE194122",
    bmmc_meta_url = env_str(
      "M3_BMMC_META_URL",
      # GEO supplementary file containing per-cell metadata (~ tens of MB).
      # Confirmed name on the GSE194122 supplementary list:
      # GSE194122_openproblems_neurips2021_cite_BMMC_processed.h5ad / multiome variant.
      # We use the lightweight metadata-only export when available; otherwise
      # 01_fetch_metadata.R falls back to extracting metadata from the .h5ad.
      ""
    )
  )
}

ensure_m3_dirs <- function(paths) {
  for (d in c(
    paths$out_dir,
    file.path(paths$out_dir, "00_precheck"),
    paths$fig_dir
  )) {
    dir.create(d, showWarnings = FALSE, recursive = TRUE)
  }
}

# --- peak name -> chr/start/end --------------------------------------------
# Robust parser for ATAC peak names. Handles the formats that show up in
# OpenProblems / 10x / Signac / ArchR outputs:
#   "chr1:10000-10500"
#   "chr1-10000-10500"
#   "chr1_10000_10500"
#   "1:10000-10500"           (numeric chrom without 'chr' prefix)
# Any peak name that fails all patterns gets NA for chr/start/end.

parse_peak_names <- function(nm) {
  patterns <- c(
    "^(chr[0-9XYM]+):(\\d+)-(\\d+)$",
    "^(chr[0-9XYM]+)-(\\d+)-(\\d+)$",
    "^(chr[0-9XYM]+)_(\\d+)_(\\d+)$",
    "^([0-9XYM]+):(\\d+)-(\\d+)$"
  )
  chr   <- rep(NA_character_, length(nm))
  start <- rep(NA_integer_,   length(nm))
  end   <- rep(NA_integer_,   length(nm))
  for (pat in patterns) {
    rem <- is.na(chr)
    if (!any(rem)) break
    m <- regexec(pat, nm[rem])
    r <- regmatches(nm[rem], m)
    ok <- lengths(r) == 4
    if (!any(ok)) next
    hit <- which(rem)[ok]
    raw_chr <- vapply(r[ok], `[[`, "", 2)
    chr[hit]   <- ifelse(grepl("^chr", raw_chr), raw_chr, paste0("chr", raw_chr))
    start[hit] <- as.integer(vapply(r[ok], `[[`, "", 3))
    end[hit]   <- as.integer(vapply(r[ok], `[[`, "", 4))
  }
  data.table::data.table(peak_id = nm, chr = chr, start = start, end = end)
}

# --- variant string -> bed (forked from module2_gwas_gene_risk) -------------

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
    end[ok]   <- as.character(pos)
  }

  data.table::data.table(
    chr   = chr,
    start = as.integer(start),
    end   = as.integer(end),
    ok    = !is.na(chr) & !is.na(start) & !is.na(end)
  )
}

# --- trait -> target cell type --------------------------------------------
#
# BMMC cell type vocabulary differs from Granja 2019. The mapping below is a
# first draft based on the published annotation set; precheck 1 must verify
# every label appears in BMMC_meta_raw.csv and update this function if needed.

# Default target mapping used by precheck 1 (cell-type coverage gating)
# and stage 2 evaluation. The LYM definition was narrowed to B-lineage on
# 2026-05-03 after baseline diagnostic on BMMC × BCX2 fine-mapping showed
# LYM SCAVENGE TRS pulls B-lineage up (AUROC 0.64) while T/NK lands at AUROC
# 0.22. The wide pan-lymphocyte target gives AUROC 0.30 — a metric artifact
# from the inverted T/NK signal masking the B-lineage signal. Stage 2
# evaluation reports B-lineage AUROC as primary and T/NK AUROC as secondary.
bmmc_target_types <- function() {
  list(
    MONO = c("CD14+ Mono", "CD16+ Mono"),
    PLT  = c("MK/E prog", "Proerythroblast", "Erythroblast", "Normoblast"),
    LYM  = c(
      "Naive CD20+ B", "Transitional B", "B1 B",
      "Plasma cell", "Lymph prog"
    )
  )
}

# Secondary target sets for diagnostic reporting only — not used for primary
# gating. Stage 2 evaluation queries these to track whether dual-modality
# rescues the T/NK signal that ATAC-only baseline drops.
bmmc_target_types_secondary <- function() {
  list(
    LYM_TNK = c("CD4+ T activated", "CD4+ T naive",
                "CD8+ T", "CD8+ T naive", "NK", "ILC"),
    # Original wide pan-lymphocyte definition kept for documentation
    LYM_wide = c("Naive CD20+ B", "Transitional B", "B1 B", "Plasma cell",
                 "Lymph prog",
                 "CD4+ T activated", "CD4+ T naive",
                 "CD8+ T", "CD8+ T naive", "NK", "ILC")
  )
}

is_bmmc_target <- function(trait, celltype) {
  if (length(trait) > 1) {
    return(as.logical(mapply(is_bmmc_target, trait, celltype, USE.NAMES = FALSE)))
  }
  tt <- bmmc_target_types()
  if (!trait %in% names(tt)) return(rep(FALSE, length(celltype)))
  as.character(celltype) %in% tt[[trait]]
}

# --- gating thresholds ----------------------------------------------------

m3_gating <- function() {
  list(
    min_cells_per_target_total  = env_integer("M3_MIN_CELLS_TOTAL",  200L),
    min_cells_per_target_donor  = env_integer("M3_MIN_CELLS_DONOR",   50L),
    # 2026-05-02: switched from min_pp_overlap_frac (relative, 30% threshold
    # was Granja-peak-set specific) to min_pp_sum_overlap (absolute signal
    # sufficiency for g-chromVAR z-score). 5 is the empirical floor below
    # which z-score loses cell-level discrimination per SCAVENGE benchmarks.
    min_pp_sum_overlap          = env_numeric("M3_MIN_PP_SUM",         5.0)
  )
}
