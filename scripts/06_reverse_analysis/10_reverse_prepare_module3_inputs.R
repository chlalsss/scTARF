# Reverse Step 3b: prepare module3-v2 compatible inputs.
#
# This is a thin wrapper around:
#   module4_traitaware_pseudopair/03_pseudopair_pipeline/02_prepare_module3_v2_inputs.R
#
# Reuse trick:
#   - existing "atac_z_block_v2.rds" slot receives Granja RNA main PCA
#     via ATAC_Z_BLOCK_OVERRIDE.
#   - existing "rna_z_block_v2_<trait>.rds" slot receives PCA of pseudo ATAC
#     profile from reverse transfer.

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
if (!has_local("MODULE3_TRANSFER_DIR")) {
  MODULE3_TRANSFER_DIR <- file.path(REV_ROOT, "outputs", REVERSE_TAG, "transfer")
}
if (!has_local("MODULE3_COMPAT_OUT_DIR")) {
  MODULE3_COMPAT_OUT_DIR <- file.path(REV_ROOT, "outputs", REVERSE_TAG, "module3_v2_compat")
}
if (!has_local("RNA_MAIN_Z_RDS")) {
  RNA_MAIN_Z_RDS <- file.path(REV_ROOT, "outputs", REVERSE_TAG, "_inputs", "granja_rna_main_z_block_v2.rds")
}
if (!has_local("SEED_TOP_PCT")) SEED_TOP_PCT <- 0.01
if (!has_local("ATAC_Z_BLOCK_OVERRIDE")) ATAC_Z_BLOCK_OVERRIDE <- RNA_MAIN_Z_RDS
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

cat(">>> Reverse Step 3b: prepare module3-compatible inputs\n\n")
cat("Transfer dir : ", MODULE3_TRANSFER_DIR, "\n", sep = "")
cat("Compat dir   : ", MODULE3_COMPAT_OUT_DIR, "\n", sep = "")
cat("RNA main z   : ", ATAC_Z_BLOCK_OVERRIDE, "\n", sep = "")
cat("Traits       : ", paste(TRAITS_TO_RUN, collapse = ", "), "\n\n", sep = "")

assert_write_under_rev(MODULE3_COMPAT_OUT_DIR)
dir.create(MODULE3_COMPAT_OUT_DIR, recursive = TRUE, showWarnings = FALSE)

source(
  file.path(NONPAIRED_SCRIPT_DIR, "44_prepare_module3_v2_inputs.R"),
  local = environment()
)
