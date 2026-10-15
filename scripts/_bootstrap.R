# Sourced at the top of paired-analysis scripts.
#
# Responsibility: locate the project root (so scripts work whether launched from
# RStudio, command-line, or another working directory), source _config.R, and
# create output directories. Does NOT install packages.

`%||%` <- function(a, b) if (is.null(a) || !nzchar(as.character(a))) b else a

.module3_root <- local({
  # Resolve module3 root by walking up from this file's location.
  this_file <- tryCatch(
    normalizePath(sys.frame(1)$ofile %||% ""),
    error = function(e) ""
  )
  if (!nzchar(this_file)) {
    # Fallback: assume current working directory is either module3 root or a
    # subdirectory beneath it.
    cand <- getwd()
    while (nzchar(cand) && !file.exists(file.path(cand, "_config.R"))) {
      parent <- dirname(cand)
      if (parent == cand) break
      cand <- parent
    }
    cand
  } else {
    dirname(this_file)
  }
})

.config_path <- file.path(.module3_root, "_config.R")
if (!file.exists(.config_path)) {
  stop("Cannot locate _config.R. Run paired scripts from the repository root or keep _bootstrap.R and _config.R together.",
       call. = FALSE)
}
source(.config_path, local = FALSE)

paths   <- m3_paths()
gating  <- m3_gating()
ensure_m3_dirs(paths)

invisible(list(paths = paths, gating = gating))
