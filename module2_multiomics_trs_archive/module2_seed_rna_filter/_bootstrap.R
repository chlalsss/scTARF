# Bootstrap loader for module2_seed_rna_filter scripts.

find_seed_filter_code_dir <- function() {
  candidates <- unique(c(
    normalizePath(getwd(), winslash = "/", mustWork = FALSE),
    file.path(normalizePath(getwd(), winslash = "/", mustWork = FALSE), "module2_seed_rna_filter"),
    "D:/RStudio/work/module2_seed_rna_filter"
  ))

  for (candidate in candidates) {
    config <- file.path(candidate, "_config.R")
    if (file.exists(config)) return(candidate)
  }

  stop(
    "Cannot locate module2_seed_rna_filter/_config.R. ",
    "Run from the project root or use the absolute script path.",
    call. = FALSE
  )
}

seed_filter_code_dir <- find_seed_filter_code_dir()
source(file.path(seed_filter_code_dir, "_config.R"))
