# Shared path helpers for the anonymous review repository.
#
# Reviewers can override these with environment variables instead of editing
# individual scripts.

sctarf_env <- function(name, default = "") {
  value <- Sys.getenv(name, "")
  if (nzchar(value)) value else default
}

sctarf_project_root <- function() {
  normalizePath(
    sctarf_env("SCTARF_PROJECT_ROOT", "."),
    winslash = "/",
    mustWork = FALSE
  )
}

sctarf_data_root <- function() {
  normalizePath(
    sctarf_env("SCTARF_DATA_ROOT", file.path(sctarf_project_root(), "local_data")),
    winslash = "/",
    mustWork = FALSE
  )
}

sctarf_output_root <- function() {
  normalizePath(
    sctarf_env("SCTARF_OUTPUT_ROOT", file.path(sctarf_project_root(), "outputs")),
    winslash = "/",
    mustWork = FALSE
  )
}

sctarf_bmmc_root <- function() {
  normalizePath(
    sctarf_env(
      "M3_OUT_DIR",
      sctarf_env("M3_DATA_DIR", file.path(sctarf_data_root(), "BMMC"))
    ),
    winslash = "/",
    mustWork = FALSE
  )
}

sctarf_blood_root <- function() {
  normalizePath(
    sctarf_env("SCTARF_BLOOD_ROOT", file.path(sctarf_data_root(), "Granja")),
    winslash = "/",
    mustWork = FALSE
  )
}

sctarf_nonpaired_root <- function() {
  normalizePath(
    sctarf_env("SCTARF_NONPAIRED_ROOT", file.path(sctarf_output_root(), "nonpaired")),
    winslash = "/",
    mustWork = FALSE
  )
}

sctarf_reverse_root <- function() {
  normalizePath(
    sctarf_env("SCTARF_REVERSE_ROOT", file.path(sctarf_output_root(), "reverse")),
    winslash = "/",
    mustWork = FALSE
  )
}

sctarf_executable <- function(command) {
  resolved <- if (file.exists(command)) command else unname(Sys.which(command))
  if (!nzchar(resolved)) stop("Executable not found on PATH: ", command, call. = FALSE)
  normalizePath(resolved, winslash = "/", mustWork = TRUE)
}

sctarf_python <- function() {
  sctarf_executable(sctarf_env("SCTARF_PYTHON", "python"))
}

sctarf_scdrs <- function() {
  sctarf_executable(sctarf_env("SCTARF_SCDRS", "scdrs"))
}
