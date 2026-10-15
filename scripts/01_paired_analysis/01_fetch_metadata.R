# 00_precheck/01_fetch_metadata.R
#
# Pull NeurIPS 2021 BMMC multiome cell metadata from GEO GSE194122.
# Stays metadata-only — does NOT download counts (the multiome h5ad is ~6 GB).
#
# Input:
#   GEO accession GSE194122 (supplementary files list)
#
# Output:
#   <out_dir>/00_precheck/BMMC_meta_raw.csv
#     columns: barcode, cell_type, donor, batch, site, modality
#
# Strategy: prefer the per-cell metadata table that ships as a separate
# supplementary file. Fall back to extracting /obs from the multiome .h5ad via
# rhdf5 if needed (still avoids loading X matrices).
#
# Run from any working directory; _bootstrap.R locates the project root.

source(file.path("scripts", "_bootstrap.R"))
required_pkgs(c("data.table"))

paths <- m3_paths()
ensure_m3_dirs(paths)

# --- source switch --------------------------------------------------------
#
# Two routes:
#   "figshare_s4d8" : single-donor s4d8 subset, pre-downloaded h5ad on disk
#                     (no network step). This is the default.
#   "GSE194122"     : full 4-donor h5ad.gz from NCBI GEO; will attempt
#                     download.file if not already cached.

if (paths$bmmc_source == "figshare_s4d8") {
  if (!file.exists(paths$figshare_h5ad)) {
    stop(
      "M3_BMMC_SOURCE=figshare_s4d8 but file not found at:\n  ",
      paths$figshare_h5ad,
      "\nDownload s4d8_quality_control.h5ad from ",
      "https://figshare.com/articles/dataset/22548547 and place it there, ",
      "or switch source via Sys.setenv(M3_BMMC_SOURCE = 'GSE194122')."
    )
  }
  local_path <- paths$figshare_h5ad
  message("Using figshare s4d8 source: ", local_path,
          " (", round(file.size(local_path) / 1e6, 1), " MB)")
} else if (paths$bmmc_source == "GSE194122") {

geo_acc <- paths$geo_accession
work_dir <- file.path(paths$out_dir, "00_precheck", "_geo_cache", geo_acc)
dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)

# We learned from a prior GEOquery::getGEOSuppFiles run that GSE194122
# only exposes two supplementary files; the multiome one is what we need.
# Hardcoding the URL avoids GEOquery's flaky list query and lets us show a
# progress bar via download.file().
multiome_fname <- "GSE194122_openproblems_neurips2021_multiome_BMMC_processed.h5ad.gz"
multiome_url   <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE194nnn/",
  geo_acc, "/suppl/", multiome_fname
)

local_gz   <- file.path(work_dir, multiome_fname)
local_h5ad <- sub("\\.gz$", "", local_gz)

# Anything below this is treated as a corrupt/truncated leftover from a
# failed previous run and removed before deciding whether to download.
MIN_VALID_BYTES <- 100e6  # 100 MB

is_valid_cache <- function(p) {
  file.exists(p) && file.size(p) >= MIN_VALID_BYTES
}
purge_corrupt <- function(p) {
  if (file.exists(p) && file.size(p) < MIN_VALID_BYTES) {
    message(
      "Removing corrupt/truncated cache (",
      file.size(p), " bytes): ", p
    )
    file.remove(p)
  }
}
purge_corrupt(local_h5ad)
purge_corrupt(local_gz)

if (is_valid_cache(local_h5ad)) {
  message("Decompressed h5ad already exists, skipping download: ", local_h5ad)
  local_path <- local_h5ad
} else if (is_valid_cache(local_gz)) {
  message("Compressed file already exists, skipping download: ", local_gz)
  local_path <- local_gz
} else {
  size_limit_mb <- env_numeric("M3_AUTO_DOWNLOAD_LIMIT_MB", 200)
  if (size_limit_mb < 1500) {
    stop(
      "Multiome h5ad.gz is expected to be ~1-3 GB. ",
      "Set Sys.setenv(M3_AUTO_DOWNLOAD_LIMIT_MB = 5000) and rerun."
    )
  }
  message("Downloading ", multiome_fname, " ...")
  message("URL: ", multiome_url)
  message("Target: ", local_gz)
  message("(this is a multi-GB download; progress bar follows)\n")

  # mode="wb" is required for binary on Windows.
  # method="libcurl" supports progress reporting cross-platform on R 4+.
  prev_method <- getOption("download.file.method")
  prev_timeout <- getOption("timeout")
  options(download.file.method = "libcurl", timeout = 60 * 60)
  on.exit({
    options(download.file.method = prev_method)
    options(timeout = prev_timeout)
  }, add = TRUE)

  status <- tryCatch(
    utils::download.file(multiome_url, destfile = local_gz, mode = "wb"),
    error = function(e) {
      message("libcurl failed: ", conditionMessage(e),
              " — retrying with method='wininet'")
      options(download.file.method = "wininet")
      utils::download.file(multiome_url, destfile = local_gz, mode = "wb")
    }
  )
  if (!file.exists(local_gz) || file.size(local_gz) < 1e7) {
    stop(
      "Download appears incomplete (file size ",
      if (file.exists(local_gz)) file.size(local_gz) else "0",
      " bytes). Check network and rerun."
    )
  }
  message("\nDownload finished, size: ",
          round(file.size(local_gz) / 1e6, 1), " MB")
  local_path <- local_gz
}

} else {
  stop("Unknown M3_BMMC_SOURCE: ", paths$bmmc_source,
       ". Use 'figshare_s4d8' or 'GSE194122'.")
}

# --- decompress + read /obs ------------------------------------------------

required_pkgs(c("rhdf5", "R.utils"))

if (grepl("\\.gz$", local_path)) {
  decomp_path <- sub("\\.gz$", "", local_path)
  if (!file.exists(decomp_path)) {
    message("Decompressing ", basename(local_path),
            " (this may take a few minutes) ...")
    R.utils::gunzip(local_path, destname = decomp_path,
                    remove = FALSE, overwrite = TRUE)
  } else {
    message("Decompressed file already exists: ", decomp_path)
  }
  local_path <- decomp_path
}

# anndata stores /obs as a group with categorical encoding for string cols:
#   /obs/_index            barcode dataset
#   /obs/<col>             either a numeric dataset or a categorical group
#   /obs/<col>/codes       int codes (0-indexed) when col is categorical
#   /obs/<col>/categories  unique category strings
read_anndata_obs <- function(h5_path) {
  ls_full <- rhdf5::h5ls(h5_path, recursive = TRUE)
  obs_top <- ls_full[ls_full$group == "/obs", ]
  if (nrow(obs_top) == 0) {
    stop("No /obs group found in ", h5_path)
  }

  # Detect legacy categorical layout. Older anndata (0.6.x / early 0.7.x) stores
  # categorical columns as a plain int dataset at /obs/<col> with the string
  # categories collected together at /obs/__categories/<col>. Newer versions
  # store /obs/<col> as a group with codes + categories children. Support both.
  legacy_cats <- list()
  if ("__categories" %in% obs_top$name) {
    legacy_grp_path <- "/obs/__categories"
    legacy_items <- ls_full[ls_full$group == legacy_grp_path, "name"]
    for (col in legacy_items) {
      legacy_cats[[col]] <- as.character(
        rhdf5::h5read(h5_path, paste0(legacy_grp_path, "/", col))
      )
    }
    message("Detected legacy categorical layout; columns with categories: ",
            paste(names(legacy_cats), collapse = ", "))
  }

  items <- obs_top$name
  out <- list()

  idx_name <- intersect(c("_index", "index"), items)
  if (length(idx_name) > 0) {
    out$barcode <- as.character(rhdf5::h5read(h5_path, paste0("/obs/", idx_name[1])))
  }

  decode_legacy <- function(raw, col) {
    codes <- as.integer(raw)
    ifelse(codes < 0L, NA_character_, legacy_cats[[col]][codes + 1L])
  }

  for (col in setdiff(items, c("_index", "index", "__categories"))) {
    info <- obs_top[obs_top$name == col, ]
    otype <- info$otype[1]
    path  <- paste0("/obs/", col)
    val <- tryCatch({
      if (grepl("DATASET", otype)) {
        raw <- as.vector(rhdf5::h5read(h5_path, path))
        if (col %in% names(legacy_cats)) decode_legacy(raw, col) else raw
      } else if (grepl("GROUP", otype)) {
        sub_items <- ls_full[ls_full$group == path, "name"]
        if (all(c("codes", "categories") %in% sub_items)) {
          cats  <- as.character(rhdf5::h5read(h5_path, paste0(path, "/categories")))
          codes <- as.integer(rhdf5::h5read(h5_path, paste0(path, "/codes")))
          ifelse(codes < 0L, NA_character_, cats[codes + 1L])
        } else {
          NA
        }
      } else {
        NA
      }
    }, error = function(e) {
      message("  failed to read /obs/", col, ": ", conditionMessage(e))
      NA
    })
    out[[col]] <- val
  }

  out <- out[vapply(out, function(x) length(x) > 1 || is.character(x), logical(1))]
  data.table::as.data.table(out)
}

message("Reading /obs from ", local_path, " ...")
meta <- read_anndata_obs(local_path)
message("Recovered obs: ", nrow(meta), " rows, ",
        ncol(meta), " columns -> ",
        paste(colnames(meta), collapse = ", "))

message("Raw metadata rows: ", nrow(meta))
message("Columns: ", paste(colnames(meta), collapse = ", "))

# --- normalize column names ------------------------------------------------
# NeurIPS 2021 obs columns (multiome): cell_type, batch, donor, site, technology

col_map <- list(
  barcode    = c("barcode", "cell_id", "obs_names", "index"),
  cell_type  = c("cell_type", "celltype", "CellType"),
  donor      = c("donor", "DonorID", "donor_id"),
  batch      = c("batch", "Batch"),
  site       = c("site", "Site"),
  modality   = c("technology", "Modality", "modality")
)

pick_col <- function(dt, options) {
  hit <- intersect(options, colnames(dt))
  if (length(hit) == 0) return(NULL)
  hit[1]
}

picked <- lapply(col_map, pick_col, dt = meta)
missing <- names(picked)[vapply(picked, is.null, logical(1))]
if (length(missing) > 0) {
  message(
    "Warning: missing expected columns -> ", paste(missing, collapse = ", "),
    ". Will leave them NA in BMMC_meta_raw.csv; inspect the raw file at ",
    local_path, " and adjust col_map if needed."
  )
}

out <- data.table::data.table(
  barcode   = if (!is.null(picked$barcode))   meta[[picked$barcode]]   else seq_len(nrow(meta)),
  cell_type = if (!is.null(picked$cell_type)) meta[[picked$cell_type]] else NA_character_,
  donor     = if (!is.null(picked$donor))     meta[[picked$donor]]     else NA_character_,
  batch     = if (!is.null(picked$batch))     meta[[picked$batch]]     else NA_character_,
  site      = if (!is.null(picked$site))      meta[[picked$site]]      else NA_character_,
  modality  = if (!is.null(picked$modality))  meta[[picked$modality]]  else NA_character_
)

# NOTE: previously filtered by modality keyword here, but the GSE194122
# multiome.h5ad file is already multiome-only, and its Modality column uses
# different vocabulary (e.g. GEX/ATAC) which broke the filter. Removed.

data.table::fwrite(out, paths$meta_raw)
message("Wrote ", paths$meta_raw, " (", nrow(out), " rows)")

# Quick sanity print
message("\nUnique cell_type (top 30 by count):")
print(head(out[, .N, by = cell_type][order(-N)], 30))
message("\nDonor x batch distribution:")
print(out[, .N, by = .(donor, batch)][order(donor, batch)])
