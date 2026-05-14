# module4 / external RNA QC: NeurIPS 2021 BMMC CITE-seq RNA
#
# Purpose:
#   Download the processed CITE-seq BMMC h5ad from GEO and check whether its
#   GEX/RNA component is suitable as an unpaired RNA reference for Granja ATAC.
#   This script intentionally does not run scDRS.

ROOT <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
THIS_DIR <- file.path(ROOT, "01_data_prep_external_rna")
RAW_DIR <- file.path(THIS_DIR, "raw")
OUT_DIR <- file.path(THIS_DIR, "output")
dir.create(RAW_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(Matrix)
})

PYTHON_EXE <- "D:/python/python.exe"
CHECK_SCRIPT <- file.path(THIS_DIR, "02_check_neurips2021_cite_rna.py")

URL <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE194nnn/GSE194122/suppl/",
  "GSE194122_openproblems_neurips2021_cite_BMMC_processed.h5ad.gz"
)
GZ_PATH <- file.path(RAW_DIR, "GSE194122_openproblems_neurips2021_cite_BMMC_processed.h5ad.gz")
H5AD_PATH <- file.path(RAW_DIR, "GSE194122_openproblems_neurips2021_cite_BMMC_processed.h5ad")

GRANJA_GA_RDS <- "D:/SCAVENGE-reproducibility-main/blood_data/multiomics_output_v2/gene_activity_matrix.rds"
GRANJA_RNA_FEATURES_RDS <- file.path(ROOT, "01_data_prep/output/granja_rna_features.rds")
GRANJA_ATAC_META_RDS <- file.path(ROOT, "01_data_prep/output/granja_atac_eval_cell_meta.rds")
GS_FILE <- file.path(ROOT, "01_data_prep/03_rna_scdrs_python/granja_3traits.gs")

GA_GENES_TXT <- file.path(OUT_DIR, "granja_ga_genes.txt")
RNA_GENES_TXT <- file.path(OUT_DIR, "granja_rna_genes.txt")
ATAC_BARCODES_TXT <- file.path(OUT_DIR, "granja_atac_barcodes.txt")
TRAIT_GENES_TXT <- file.path(OUT_DIR, "granja_3traits_genes_long.csv")

cat(">>> External RNA QC: NeurIPS 2021 BMMC CITE-seq RNA\n\n")

file_size_mb <- function(path) {
  if (!file.exists(path)) return(NA_real_)
  round(file.info(path)$size / 1024^2, 1)
}

is_valid_h5ad <- function(path) {
  if (!file.exists(path) || file.info(path)$size < 100 * 1024^2) {
    return(FALSE)
  }
  con <- file(path, open = "rb")
  on.exit(close(con), add = TRUE)
  header <- readBin(con, what = "raw", n = 8)
  identical(header, as.raw(c(0x89, 0x48, 0x44, 0x46, 0x0d, 0x0a, 0x1a, 0x0a)))
}

download_if_needed <- function(url, dest, min_mb = 500) {
  if (file.exists(dest) && file.info(dest)$size > min_mb * 1024^2) {
    cat("  found existing file: ", dest, " (", file_size_mb(dest), " MB)\n", sep = "")
    return(invisible(TRUE))
  }

  cat("  downloading from GEO:\n    ", url, "\n", sep = "")
  cat("  destination:\n    ", dest, "\n", sep = "")

  curl_bin <- Sys.which("curl.exe")
  if (nzchar(curl_bin)) {
    status <- system2(
      curl_bin,
      args = c("-L", "--fail", "--retry", "3", "-C", "-", "-o", dest, url)
    )
    if (!identical(status, 0L)) {
      stop("curl download failed with status ", status, call. = FALSE)
    }
  } else {
    download.file(url, destfile = dest, mode = "wb", method = "libcurl", quiet = FALSE)
  }

  if (!file.exists(dest) || file.info(dest)$size <= min_mb * 1024^2) {
    stop("Downloaded file is missing or unexpectedly small: ", dest, call. = FALSE)
  }
  invisible(TRUE)
}

decompress_gzip_if_needed <- function(gz_path, out_path) {
  if (is_valid_h5ad(out_path)) {
    cat("  found existing decompressed h5ad: ", out_path, " (", file_size_mb(out_path), " MB)\n", sep = "")
    return(invisible(TRUE))
  }

  cat("  decompressing gzip to h5ad ...\n")
  in_con <- gzfile(gz_path, open = "rb")
  on.exit(close(in_con), add = TRUE)
  out_con <- file(out_path, open = "wb")
  on.exit(close(out_con), add = TRUE)

  repeat {
    chunk <- readBin(in_con, what = "raw", n = 16 * 1024^2)
    if (length(chunk) == 0L) break
    writeBin(chunk, out_con)
  }

  if (!is_valid_h5ad(out_path)) {
    stop("Decompressed h5ad is missing or is not a valid HDF5/h5ad file: ", out_path, call. = FALSE)
  }
  invisible(TRUE)
}

parse_gs_to_long <- function(gs_file, out_csv) {
  lines <- readLines(gs_file, warn = FALSE)
  lines <- lines[nzchar(lines)]
  if (length(lines) < 2L) stop("Gene set file has no records: ", gs_file, call. = FALSE)

  records <- lines[-1L]
  trait <- character()
  gene <- character()
  weight <- numeric()

  for (line in records) {
    parts <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(parts) < 2L) next
    geneset <- strsplit(parts[2], ",", fixed = TRUE)[[1]]
    for (item in geneset) {
      pair <- strsplit(item, ":", fixed = TRUE)[[1]]
      if (length(pair) < 1L || !nzchar(pair[1])) next
      trait <- c(trait, parts[1])
      gene <- c(gene, pair[1])
      weight <- c(weight, if (length(pair) >= 2L) suppressWarnings(as.numeric(pair[2])) else NA_real_)
    }
  }

  out <- data.frame(trait = trait, gene = gene, weight = weight, stringsAsFactors = FALSE)
  write.csv(out, out_csv, row.names = FALSE, quote = TRUE)
  out
}

extract_feature_names <- function(x, preferred_cols = c("gene", "gene_name", "gene_symbol", "feature", "features")) {
  if (is.null(x)) {
    return(character())
  }
  if (is.character(x)) {
    return(as.character(x))
  }
  if (is.data.frame(x)) {
    for (col in preferred_cols) {
      if (col %in% colnames(x)) {
        return(as.character(x[[col]]))
      }
    }
    if (!is.null(rownames(x))) {
      return(as.character(rownames(x)))
    }
  }
  rn <- rownames(x)
  if (!is.null(rn)) {
    return(as.character(rn))
  }
  nm <- names(x)
  if (!is.null(nm)) {
    return(as.character(nm))
  }
  as.character(x)
}

cat("[1/5] Locating/downloading CITE h5ad.gz if needed ...\n")
if (is_valid_h5ad(H5AD_PATH)) {
  cat("  found existing decompressed h5ad, skip gz download:\n    ", H5AD_PATH, "\n", sep = "")
} else {
  download_if_needed(URL, GZ_PATH)
  cat("  gz size: ", file_size_mb(GZ_PATH), " MB\n", sep = "")
}
cat("\n")

cat("[2/5] Decompressing to h5ad if needed ...\n")
decompress_gzip_if_needed(GZ_PATH, H5AD_PATH)
cat("  h5ad size: ", file_size_mb(H5AD_PATH), " MB\n\n", sep = "")

cat("[3/5] Exporting Granja reference gene/barcode lists ...\n")
if (!file.exists(GRANJA_GA_RDS)) stop("Missing Granja GA RDS: ", GRANJA_GA_RDS, call. = FALSE)
if (!file.exists(GRANJA_ATAC_META_RDS)) stop("Missing Granja ATAC meta RDS: ", GRANJA_ATAC_META_RDS, call. = FALSE)
if (!file.exists(GS_FILE)) stop("Missing 3-trait gene set file: ", GS_FILE, call. = FALSE)

ga <- readRDS(GRANJA_GA_RDS)
ga_genes <- extract_feature_names(ga)
if (length(ga_genes) < 10000L) {
  stop("Too few Granja GA genes extracted: ", length(ga_genes), call. = FALSE)
}
writeLines(ga_genes, GA_GENES_TXT, useBytes = TRUE)
rm(ga)
gc(verbose = FALSE)

if (file.exists(GRANJA_RNA_FEATURES_RDS)) {
  granja_rna_genes <- readRDS(GRANJA_RNA_FEATURES_RDS)
  writeLines(extract_feature_names(granja_rna_genes), RNA_GENES_TXT, useBytes = TRUE)
} else {
  writeLines(character(), RNA_GENES_TXT, useBytes = TRUE)
}

atac_meta <- readRDS(GRANJA_ATAC_META_RDS)
atac_barcodes <- rownames(as.data.frame(atac_meta))
writeLines(atac_barcodes, ATAC_BARCODES_TXT, useBytes = TRUE)
rm(atac_meta)
gc(verbose = FALSE)

trait_genes <- parse_gs_to_long(GS_FILE, TRAIT_GENES_TXT)
cat("  Granja GA genes:      ", length(readLines(GA_GENES_TXT, warn = FALSE)), "\n", sep = "")
cat("  Granja RNA genes:     ", length(readLines(RNA_GENES_TXT, warn = FALSE)), "\n", sep = "")
cat("  Granja ATAC cells:    ", length(readLines(ATAC_BARCODES_TXT, warn = FALSE)), "\n", sep = "")
cat("  3-trait genes rows:   ", nrow(trait_genes), "\n\n", sep = "")

cat("[4/5] Running h5ad usability QC, without scDRS ...\n")
if (!file.exists(PYTHON_EXE)) stop("Missing Python executable: ", PYTHON_EXE, call. = FALSE)
if (!file.exists(CHECK_SCRIPT)) stop("Missing QC script: ", CHECK_SCRIPT, call. = FALSE)

args <- c(
  CHECK_SCRIPT,
  "--h5ad", H5AD_PATH,
  "--out-dir", OUT_DIR,
  "--granja-ga-genes", GA_GENES_TXT,
  "--granja-rna-genes", RNA_GENES_TXT,
  "--granja-atac-barcodes", ATAC_BARCODES_TXT,
  "--trait-genes", TRAIT_GENES_TXT
)
status <- system2(PYTHON_EXE, args = args)
if (!identical(status, 0L)) {
  stop("Python QC script failed with status ", status, call. = FALSE)
}

cat("\n[5/5] Main QC outputs:\n")
outputs <- c(
  "neurips2021_cite_rna_qc_summary.csv",
  "neurips2021_cite_celltype_counts.csv",
  "neurips2021_cite_candidate_target_summary.csv",
  "neurips2021_cite_marker_program_by_celltype.csv",
  "neurips2021_cite_gene_overlap.csv",
  "neurips2021_cite_trait_gene_overlap.csv",
  "neurips2021_cite_usable_report.txt"
)
for (f in outputs) {
  cat("  ", file.path(OUT_DIR, f), "\n", sep = "")
}

cat("\nExternal RNA QC OK. scDRS was not run.\n")
