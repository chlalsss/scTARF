# Bootstrap loader for module2_wnn_lite_gated scripts.

find_wnn_gated_code_dir <- function() {
  candidates <- unique(c(
    normalizePath(getwd(), winslash = "/", mustWork = FALSE),
    file.path(normalizePath(getwd(), winslash = "/", mustWork = FALSE), "module2_wnn_lite_gated"),
    "D:/RStudio/work/module2_wnn_lite_gated"
  ))

  for (candidate in candidates) {
    config <- file.path(candidate, "_config.R")
    if (file.exists(config)) return(candidate)
  }

  stop(
    "Cannot locate module2_wnn_lite_gated/_config.R. ",
    "Run from the project root or use the absolute script path.",
    call. = FALSE
  )
}

wnn_gated_code_dir <- find_wnn_gated_code_dir()
source(file.path(wnn_gated_code_dir, "_config.R"))
