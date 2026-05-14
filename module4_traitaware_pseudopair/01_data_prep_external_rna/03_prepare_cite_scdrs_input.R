# module4 / external RNA: prepare NeurIPS 2021 CITE-seq GEX for scDRS
#
# This is the pre-scDRS preprocessing step only. It does not call scDRS.
# Output is a GEX-only raw-count h5ad with unique gene symbols.

ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
THIS_DIR <- file.path(ROOT, "01_data_prep_external_rna")
SCDRS_DIR <- file.path(THIS_DIR, "03_cite_scdrs_python")
INPUT_DIR <- file.path(SCDRS_DIR, "input")
OUTPUT_DIR <- file.path(SCDRS_DIR, "output")
dir.create(INPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

PYTHON_EXE <- "D:/python/python.exe"
PY_SCRIPT <- file.path(SCDRS_DIR, "prepare_cite_scdrs_input.py")

RAW_H5AD <- file.path(
  THIS_DIR, "raw",
  "GSE194122_openproblems_neurips2021_cite_BMMC_processed.h5ad"
)
GA_GENES <- file.path(THIS_DIR, "output", "granja_ga_genes.txt")
TRAIT_GENES <- file.path(THIS_DIR, "output", "granja_3traits_genes_long.csv")
GS_IN <- file.path(ROOT, "01_data_prep/03_rna_scdrs_python/granja_3traits.gs")
GS_OUT <- file.path(SCDRS_DIR, "cite_3traits.gs")
OUT_H5AD <- file.path(INPUT_DIR, "cite_rna_raw_counts_gex_only.h5ad")

cat(">>> Pre-scDRS preprocessing: NeurIPS 2021 CITE RNA\n\n")

required <- c(RAW_H5AD, GA_GENES, TRAIT_GENES, GS_IN, PY_SCRIPT, PYTHON_EXE)
missing <- required[!file.exists(required)]
if (length(missing) > 0L) {
  stop("Missing required file(s):\n  ", paste(missing, collapse = "\n  "), call. = FALSE)
}

cat("[1/3] Copying 3-trait gene set for CITE scDRS workspace ...\n")
file.copy(GS_IN, GS_OUT, overwrite = TRUE)
cat("  ", GS_OUT, "\n\n", sep = "")

cat("[2/3] Building GEX-only raw-count h5ad ...\n")
args <- c(
  PY_SCRIPT,
  "--input-h5ad", RAW_H5AD,
  "--out-h5ad", OUT_H5AD,
  "--out-dir", SCDRS_DIR,
  "--granja-ga-genes", GA_GENES,
  "--trait-genes", TRAIT_GENES
)
status <- system2(PYTHON_EXE, args = args)
if (!identical(status, 0L)) {
  stop("Python preprocessing failed with status ", status, call. = FALSE)
}

cat("\n[3/3] Outputs ready:\n")
outputs <- c(
  OUT_H5AD,
  GS_OUT,
  file.path(SCDRS_DIR, "cite_scdrs_input_summary.csv"),
  file.path(SCDRS_DIR, "cite_scdrs_gene_overlap.csv"),
  file.path(SCDRS_DIR, "cite_scdrs_trait_gene_overlap.csv"),
  file.path(SCDRS_DIR, "cite_scdrs_celltype_counts.csv"),
  file.path(SCDRS_DIR, "cite_scdrs_ready_report.txt")
)
for (f in outputs) {
  cat("  ", f, "\n", sep = "")
}

cat("\nPre-scDRS preprocessing OK. scDRS was not run.\n")
