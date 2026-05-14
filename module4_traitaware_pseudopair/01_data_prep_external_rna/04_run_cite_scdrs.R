# module4 / external RNA: run scDRS on prepared CITE-seq GEX h5ad
#
# This runs official scDRS. It can take a long time because the CITE reference
# has 90k cells and scDRS writes full control scores for three traits.

ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
RUN_PY <- file.path(ROOT, "01_data_prep_external_rna/03_cite_scdrs_python/run_scdrs_cite.py")
PYTHON_EXE <- "D:/python/python.exe"

cat(">>> Run CITE RNA scDRS\n\n")
stopifnot(file.exists(RUN_PY), file.exists(PYTHON_EXE))

status <- system2(PYTHON_EXE, args = c(RUN_PY))
if (!identical(status, 0L)) {
  stop("CITE scDRS failed with status ", status, call. = FALSE)
}

cat("\nCITE scDRS finished.\n")
cat("Next, summarize with:\n")
cat("  source(\"01_data_prep_external_rna/05_summarize_cite_scdrs.R\")\n")
