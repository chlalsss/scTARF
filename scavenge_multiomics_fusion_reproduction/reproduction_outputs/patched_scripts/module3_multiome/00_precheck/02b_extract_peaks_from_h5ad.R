# 00_precheck/02b_extract_peaks_from_h5ad.R
#
# Extract BMMC ATAC peak coordinates from the GSE194122 multiome .h5ad's /var,
# write a 3-column bed file at paths$peaks_bed for use by 03_variant_peak_overlap.R.
#
# Why a separate script: the same .h5ad already pulled in 01_fetch_metadata.R
# carries peak coordinates as /var feature names. No re-download needed.
#
# Output:
#   <out_dir>/00_precheck/BMMC_peaks.bed   (chr <tab> start <tab> end)

source("D:/RStudio/work/scavenge_multiomics_fusion_reproduction/reproduction_outputs/patched_scripts/module3_multiome/_bootstrap.R")
required_pkgs(c("rhdf5", "data.table"))

paths <- m3_paths()

# Locate the .h5ad picked by current M3_BMMC_SOURCE
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

if (!file.exists(h5_path)) {
  stop("h5ad not found at: ", h5_path,
       "\nRun 01_fetch_metadata.R first to ensure file is in place.")
}

# --- explore /var ---------------------------------------------------------

ls_full <- rhdf5::h5ls(h5_path, recursive = TRUE)
var_top <- ls_full[ls_full$group == "/var", ]
if (nrow(var_top) == 0) stop("No /var group found in ", h5_path)
message("/var children: ", paste(var_top$name, collapse = ", "))

# --- legacy categorical decode helper (same logic as 01) ------------------

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
  } else {
    NULL
  }
}

# --- feature names + modality flag ----------------------------------------

idx_name <- intersect(c("_index", "index"), var_top$name)[1]
if (is.na(idx_name)) stop("No _index/index dataset under /var")
feat_names <- as.character(rhdf5::h5read(h5_path, paste0("/var/", idx_name)))
message("Total features: ", length(feat_names))

ftype_col <- intersect(
  c("feature_types", "feature_type", "modality", "Modality"),
  var_top$name
)[1]
if (is.na(ftype_col)) {
  stop("No feature_types/modality column in /var. Available: ",
       paste(var_top$name, collapse = ", "))
}
ftype <- read_var_col(ftype_col)
if (is.null(ftype)) stop("Failed to decode /var/", ftype_col)

message("Feature type counts:")
print(table(ftype, useNA = "ifany"))

# --- pick ATAC peaks ------------------------------------------------------

atac_mask <- grepl("ATAC|peak", ftype, ignore.case = TRUE)
n_atac <- sum(atac_mask)
message("ATAC features: ", n_atac, " / ", length(ftype))
if (n_atac == 0) {
  stop("No ATAC features identified by feature_type filter. ",
       "Inspect unique values above and adjust mask.")
}
peak_names <- feat_names[atac_mask]
message("Sample peak names: ", paste(head(peak_names, 5), collapse = " | "))

# --- parse "chr1:10000-10500" / "chr1-10000-10500" / "chr1_10000_10500" ---

patterns <- c(
  "^(chr[0-9XYM]+):(\\d+)-(\\d+)$",
  "^(chr[0-9XYM]+)-(\\d+)-(\\d+)$",
  "^(chr[0-9XYM]+)_(\\d+)_(\\d+)$",
  # numeric chromosome without 'chr' prefix
  "^([0-9XYM]+):(\\d+)-(\\d+)$"
)
parsed_chr <- rep(NA_character_, n_atac)
parsed_st  <- rep(NA_integer_,   n_atac)
parsed_en  <- rep(NA_integer_,   n_atac)
for (pat in patterns) {
  rem <- is.na(parsed_chr)
  if (!any(rem)) break
  m <- regexec(pat, peak_names[rem])
  r <- regmatches(peak_names[rem], m)
  ok <- lengths(r) == 4
  if (!any(ok)) next
  hit <- which(rem)[ok]
  parsed_chr[hit] <- vapply(r[ok], `[[`, "", 2)
  parsed_st[hit]  <- as.integer(vapply(r[ok], `[[`, "", 3))
  parsed_en[hit]  <- as.integer(vapply(r[ok], `[[`, "", 4))
}

n_ok <- sum(!is.na(parsed_chr))
n_bad <- n_atac - n_ok
message("Parsed: ", n_ok, " / ", n_atac,
        " (", n_bad, " unparsed; sample: ",
        paste(head(peak_names[is.na(parsed_chr)], 5), collapse = " | "), ")")

if (n_ok == 0) {
  stop("No peak names parseable. Sample: ",
       paste(head(peak_names, 10), collapse = " | "))
}

bed <- data.table::data.table(
  chr   = ifelse(grepl("^chr", parsed_chr), parsed_chr, paste0("chr", parsed_chr)),
  start = parsed_st,
  end   = parsed_en
)
bed <- bed[!is.na(start) & !is.na(end)]

# half-open BED convention: start should be 0-based; assume 1-based input
# (anndata peak names are typically 1-based). Convert to 0-based start.
bed[, start := start - 1L]

data.table::fwrite(bed, paths$peaks_bed, sep = "\t", col.names = FALSE)
message("Wrote ", paths$peaks_bed, " (", nrow(bed), " peaks)")
message("\nFirst 5 lines:")
print(head(bed, 5))
