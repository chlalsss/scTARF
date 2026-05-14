# Locate and load _config.R from both Rscript and RStudio source() contexts.

find_wnn_lite_code_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  script_file <- sub("--file=", "", args[grepl("--file=", args, fixed = TRUE)])

  candidates <- character()
  if (length(script_file) > 0) {
    candidates <- c(candidates, dirname(normalizePath(script_file[1], winslash = "/", mustWork = FALSE)))
  }

  frame_files <- vapply(sys.frames(), function(frame) {
    if (exists("ofile", envir = frame, inherits = FALSE)) {
      get("ofile", envir = frame, inherits = FALSE)
    } else {
      NA_character_
    }
  }, character(1))
  frame_files <- frame_files[!is.na(frame_files) & nzchar(frame_files)]
  if (length(frame_files) > 0) {
    candidates <- c(candidates, dirname(normalizePath(frame_files, winslash = "/", mustWork = FALSE)))
  }

  candidates <- c(
    candidates,
    getwd(),
    file.path(getwd(), "module2_wnn_lite_fusion"),
    "D:/RStudio/work/module2_wnn_lite_fusion"
  )

  candidates <- unique(candidates)
  hits <- candidates[file.exists(file.path(candidates, "_config.R"))]
  if (length(hits) == 0) {
    stop(
      "Cannot locate module2_wnn_lite_fusion/_config.R. ",
      "Set working directory to D:/RStudio/work or source scripts by absolute path.",
      call. = FALSE
    )
  }
  hits[1]
}

code_dir <- find_wnn_lite_code_dir()
source(file.path(code_dir, "_config.R"))

