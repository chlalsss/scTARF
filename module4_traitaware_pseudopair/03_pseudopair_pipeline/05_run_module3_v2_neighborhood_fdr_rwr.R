# module4: run module3 v2 RWR/TRS using neighborhood-FDR seeds
# -----------------------------------------------------------------------------
# Purpose:
#   Reuse module3's existing 05_rwr_trs.R. This script only patches the runtime
#   OUT_DIR, TRAITS, and ATAC_W_LIST so module3 reads the neighborhood-FDR seed
#   inputs prepared by:
#     03_pseudopair_pipeline/04_prepare_neighborhood_fdr_seed_inputs.R
#
# Important:
#   This does not rebuild mapping.
#   This does not rebuild KNN.
#   It only reruns RWR/TRS with the new seed vectors.
#
# Run in RStudio:
#   setwd("D:/RStudio/work/module4_traitaware_pseudopair")
#   source("03_pseudopair_pipeline/05_run_module3_v2_neighborhood_fdr_rwr.R")
# -----------------------------------------------------------------------------

ROOT <- "D:/RStudio/work/module4_traitaware_pseudopair"
MODULE3_DIR <- "D:/RStudio/work/module3_multiome/04_dual_modality"
if (!exists("OUT_DIR")) {
  OUT_DIR <- file.path(ROOT, "03_pseudopair_pipeline/output/module3_v2_neighborhood_fdr_ge0p80")
}
PATCHED_DIR <- file.path(OUT_DIR, "_patched_module3_scripts")

# Use neighborhood-specific variable names so old variables left in the RStudio
# global environment by earlier scripts cannot accidentally pull PLT back in.
if (!exists("NBD_TRAITS_TO_RUN")) NBD_TRAITS_TO_RUN <- c("MONO", "LYM")
if (!exists("NBD_ATAC_W_LIST")) NBD_ATAC_W_LIST <- c(0.7, 1.0)

MODULE3_RWR_SCRIPT <- file.path(MODULE3_DIR, "05_rwr_trs.R")

dir.create(PATCHED_DIR, recursive = TRUE, showWarnings = FALSE)

stop_missing <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0L) {
    stop("Missing required file(s):\n  ", paste(missing, collapse = "\n  "),
         call. = FALSE)
  }
}

r_literal_chr <- function(x) {
  paste(sprintf('"%s"', x), collapse = ", ")
}

r_literal_num <- function(x) {
  paste(format(x, trim = TRUE, scientific = FALSE), collapse = ", ")
}

patch_module3_rwr_script <- function(script_path, patched_path) {
  txt <- readLines(script_path, warn = FALSE)

  txt <- sub(
    "^ATAC_W_LIST\\s*<-.*$",
    paste0("ATAC_W_LIST <- c(", r_literal_num(NBD_ATAC_W_LIST), ")"),
    txt
  )
  txt <- sub(
    "^OUT_DIR\\s*<-.*$",
    paste0('OUT_DIR <- "', normalizePath(OUT_DIR, winslash = "/", mustWork = TRUE), '"'),
    txt
  )
  txt <- sub(
    "^TRAITS\\s*<-.*$",
    paste0("TRAITS <- c(", r_literal_chr(NBD_TRAITS_TO_RUN), ")"),
    txt
  )

  writeLines(txt, patched_path, useBytes = TRUE)
  patched_path
}

cat(">>> Run module3 v2 RWR/TRS with neighborhood-FDR seeds\n\n")
cat("OUT_DIR     : ", OUT_DIR, "\n", sep = "")
cat("TRAITS      : ", paste(NBD_TRAITS_TO_RUN, collapse = ", "), "\n", sep = "")
cat("ATAC_W_LIST : ", paste(NBD_ATAC_W_LIST, collapse = ", "), "\n\n", sep = "")

required_inputs <- c(
  file.path(OUT_DIR, paste0("seed_v2_", NBD_TRAITS_TO_RUN, ".rds")),
  as.vector(sapply(NBD_TRAITS_TO_RUN, function(trait) {
    file.path(OUT_DIR, sprintf("knn_v2_%s_w%.1f.rds", trait, NBD_ATAC_W_LIST))
  })),
  MODULE3_RWR_SCRIPT
)
stop_missing(required_inputs)

patched_rwr <- file.path(PATCHED_DIR, "module3_05_rwr_trs.neighborhood_fdr_patched.R")

cat("[1/2] Creating patched module3 RWR script ...\n")
patch_module3_rwr_script(MODULE3_RWR_SCRIPT, patched_rwr)
cat("  wrote: ", patched_rwr, "\n\n", sep = "")

cat("[2/2] Running module3 RWR/TRS script ...\n")
source(patched_rwr, local = new.env(parent = globalenv()))

cat("\n=== module3 v2 RWR/TRS finished with neighborhood-FDR seeds ===\n")
cat("Expected outputs in: ", OUT_DIR, "\n", sep = "")
cat("  TRS_v2_<trait>_w<weight>.tsv.gz\n")
cat("  trs_summary_v2.csv\n")
