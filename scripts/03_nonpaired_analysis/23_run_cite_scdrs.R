# module4 / external RNA: run scDRS on prepared CITE-seq GEX h5ad
#
# This runs official scDRS. It can take a long time because the CITE reference
# has 90k cells and scDRS writes full control scores for three traits.

source(file.path("scripts", "_paths.R"))
ROOT <- sctarf_nonpaired_root()
CODE_DIR <- file.path(sctarf_project_root(), "scripts", "03_nonpaired_analysis")
RUN_PY <- file.path(CODE_DIR, "22_run_scdrs_cite.py")
PYTHON_EXE <- sctarf_python()

cat(">>> Run CITE RNA scDRS\n\n")
stopifnot(file.exists(RUN_PY))
if (!file.exists(PYTHON_EXE) && !nzchar(Sys.which(PYTHON_EXE))) {
  stop("Missing Python executable: ", PYTHON_EXE, call. = FALSE)
}

status <- system2(PYTHON_EXE, args = shQuote(RUN_PY))
if (!identical(status, 0L)) {
  stop("CITE scDRS failed with status ", status, call. = FALSE)
}

cat("\nCITE scDRS finished.\n")
cat("Next, summarize with:\n")
cat("  source(\"01_data_prep_external_rna/05_summarize_cite_scdrs.R\")\n")
