# Reverse Step 3c: run existing module3-v2 KNN + RWR/TRS.
#
# Thin wrapper around:
#   module4_traitaware_pseudopair/03_pseudopair_pipeline/03_run_module3_v2_knn_rwr.R
#
# In this reverse run, module3's variable name "ATAC_W" means:
#   weight on the main block, which is Granja RNA PCA via Step 3b.

if (!exists("ROOT")) ROOT <- Sys.getenv("SCTARF_WORK_ROOT", normalizePath(".", winslash = "/", mustWork = FALSE))
if (!exists("MODULE4_ROOT")) MODULE4_ROOT <- Sys.getenv("SCTARF_NONPAIRED_ROOT", file.path(ROOT, "outputs", "nonpaired"))
if (!exists("REV_ROOT")) REV_ROOT <- Sys.getenv("SCTARF_REVERSE_ROOT", file.path(ROOT, "outputs", "reverse"))
NONPAIRED_SCRIPT_DIR <- if (exists("SCTARF_NONPAIRED_SCRIPT_DIR")) {
  SCTARF_NONPAIRED_SCRIPT_DIR
} else {
  file.path(ROOT, "scripts", "03_nonpaired_analysis")
}

SCRIPT_ENV <- environment()
has_local <- function(x) exists(x, envir = SCRIPT_ENV, inherits = FALSE)

if (!has_local("TRAITS_TO_RUN")) TRAITS_TO_RUN <- c("PLT", "MCH", "RBC")
if (!has_local("REVERSE_TAG")) REVERSE_TAG <- "reverse_granja_rna_main_plt_mch_rbc"
if (!has_local("MODULE3_COMPAT_OUT_DIR")) {
  MODULE3_COMPAT_OUT_DIR <- file.path(REV_ROOT, "outputs", REVERSE_TAG, "module3_v2_compat")
}
if (!has_local("ATAC_W_LIST")) ATAC_W_LIST <- c(0.7, 1.0)
if (!has_local("ALLOW_OUTSIDE_REV_ROOT")) ALLOW_OUTSIDE_REV_ROOT <- FALSE

assert_write_under_rev <- function(paths) {
  if (isTRUE(ALLOW_OUTSIDE_REV_ROOT)) return(invisible(TRUE))
  root <- normalizePath(REV_ROOT, winslash = "/", mustWork = TRUE)
  for (path in paths) {
    resolved <- normalizePath(path, winslash = "/", mustWork = FALSE)
    if (!startsWith(resolved, paste0(root, "/")) && !identical(resolved, root)) {
      stop("Refusing to write outside reverse workspace: ", path,
           "\nSet ALLOW_OUTSIDE_REV_ROOT <- TRUE only if this is intentional.",
           call. = FALSE)
    }
  }
  invisible(TRUE)
}

cat(">>> Reverse Step 3c: run module3 KNN/RWR/TRS\n\n")
cat("Compat dir : ", MODULE3_COMPAT_OUT_DIR, "\n", sep = "")
cat("Traits     : ", paste(TRAITS_TO_RUN, collapse = ", "), "\n", sep = "")
cat("Main RNA weights: ", paste(ATAC_W_LIST, collapse = ", "), "\n\n", sep = "")

assert_write_under_rev(MODULE3_COMPAT_OUT_DIR)
dir.create(MODULE3_COMPAT_OUT_DIR, recursive = TRUE, showWarnings = FALSE)

source(
  file.path(NONPAIRED_SCRIPT_DIR, "45_run_module3_v2_knn_rwr.R"),
  local = environment()
)
