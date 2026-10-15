# Copy this file to paths.local.R in the repository root and edit it.
# scripts/run_workflow.R loads paths.local.R in the parent and each R child.
# Edit and source this file from the repository root before running workflows,
# or set the same environment variables manually.

PROJECT_ROOT <- normalizePath(".", winslash = "/", mustWork = FALSE)
DATA_ROOT <- file.path(PROJECT_ROOT, "local_data")
OUTPUT_ROOT <- file.path(PROJECT_ROOT, "outputs")

Sys.setenv(
  SCTARF_PROJECT_ROOT = PROJECT_ROOT,
  SCTARF_DATA_ROOT = DATA_ROOT,
  SCTARF_OUTPUT_ROOT = OUTPUT_ROOT,
  SCTARF_BLOOD_ROOT = file.path(DATA_ROOT, "Granja"),
  SCTARF_NONPAIRED_ROOT = file.path(OUTPUT_ROOT, "nonpaired"),
  SCTARF_REVERSE_ROOT = file.path(OUTPUT_ROOT, "reverse"),
  M3_DATA_DIR = file.path(DATA_ROOT, "BMMC"),
  M3_OUT_DIR = file.path(DATA_ROOT, "BMMC"),
  SCTARF_PYTHON = "python",
  SCTARF_SCDRS = "scdrs"
)

message("scTARF review paths configured for: ", PROJECT_ROOT)

# Optional paired sparsification overrides (defaults run all manuscript depths):
# P_DROP_LIST <- seq(0, 0.9, by = 0.1)
# SEED_LIST <- c(1L, 2L, 3L)
# ATAC_W_LIST <- c(0.7, 1.0)
# For a fresh efficiency benchmark (otherwise existing summaries may be reused):
# Sys.setenv(EFFICIENCY_FORCE_RERUN = "true")
