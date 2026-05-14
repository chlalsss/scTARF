# Bootstrap loader for module2_adaptive_gamma scripts.

find_adaptive_gamma_code_dir <- function() {
  candidates <- unique(c(
    normalizePath(getwd(), winslash = "/", mustWork = FALSE),
    file.path(normalizePath(getwd(), winslash = "/", mustWork = FALSE), "module2_adaptive_gamma"),
    "D:/RStudio/work/module2_adaptive_gamma"
  ))

  for (candidate in candidates) {
    config <- file.path(candidate, "_config.R")
    if (file.exists(config)) return(candidate)
  }

  stop(
    "Cannot locate module2_adaptive_gamma/_config.R. ",
    "Run from the project root or use the absolute script path.",
    call. = FALSE
  )
}

adaptive_gamma_code_dir <- find_adaptive_gamma_code_dir()
source(file.path(adaptive_gamma_code_dir, "_config.R"))
source(file.path(adaptive_gamma_code_dir, "_rwr_hetgamma.R"))
