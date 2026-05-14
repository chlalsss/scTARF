# module4 trait screening: run official scDRS on CITE RNA candidate traits
# -----------------------------------------------------------------------------
# Purpose:
#   Run the Python scDRS wrapper for the candidate blood traits prepared by:
#     01_data_prep_external_rna/06_prepare_trait_screening_gs.R
#
# This can take a long time because full control scores are written.
#
# Run in RStudio:
#   setwd("D:/RStudio/work/module4_traitaware_pseudopair")
#   source("01_data_prep_external_rna/07_run_cite_scdrs_trait_screening.R")
#
# Optional before source:
#   N_CTRL <- 1000L
#   FORCE_COMPUTE <- FALSE
#   FORCE_DOWNSTREAM <- FALSE
# -----------------------------------------------------------------------------

ROOT <- "D:/RStudio/work/module4_traitaware_pseudopair"
RUN_PY <- file.path(
  ROOT,
  "01_data_prep_external_rna/03_cite_scdrs_python/run_scdrs_cite_trait_screening.py"
)
GS_FILE <- file.path(
  ROOT,
  "01_data_prep_external_rna/03_cite_scdrs_python/cite_trait_screening.gs"
)
PYTHON_EXE <- "D:/python/python.exe"

if (!exists("N_CTRL")) N_CTRL <- 1000L
if (!exists("FORCE_COMPUTE")) FORCE_COMPUTE <- FALSE
if (!exists("FORCE_DOWNSTREAM")) FORCE_DOWNSTREAM <- FALSE
if (!exists("SKIP_DOWNSTREAM")) SKIP_DOWNSTREAM <- FALSE

cat(">>> Run CITE RNA scDRS for trait screening\n\n")
cat("Gene set: ", GS_FILE, "\n", sep = "")
cat("n_ctrl  : ", N_CTRL, "\n\n", sep = "")

missing <- c(RUN_PY, GS_FILE, PYTHON_EXE)[!file.exists(c(RUN_PY, GS_FILE, PYTHON_EXE))]
if (length(missing) > 0L) {
  stop("Missing required file(s):\n  ", paste(missing, collapse = "\n  "),
       call. = FALSE)
}

args <- c(RUN_PY, "--n-ctrl", as.character(N_CTRL))
if (isTRUE(FORCE_COMPUTE)) args <- c(args, "--force-compute")
if (isTRUE(FORCE_DOWNSTREAM)) args <- c(args, "--force-downstream")
if (isTRUE(SKIP_DOWNSTREAM)) args <- c(args, "--skip-downstream")

status <- system2(PYTHON_EXE, args = args)
if (!identical(status, 0L)) {
  stop("CITE scDRS trait screening failed with status ", status, call. = FALSE)
}

cat("\nCITE scDRS trait screening finished.\n")
cat("Next summarize with:\n")
cat("  source(\"01_data_prep_external_rna/08_summarize_cite_scdrs_trait_screening.R\")\n")
