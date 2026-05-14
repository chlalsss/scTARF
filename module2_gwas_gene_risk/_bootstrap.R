# Bootstrap loader for module2_gwas_gene_risk scripts.

find_gwas_gene_code_dir <- function() {
  candidates <- unique(c(
    normalizePath(getwd(), winslash = "/", mustWork = FALSE),
    file.path(normalizePath(getwd(), winslash = "/", mustWork = FALSE), "module2_gwas_gene_risk"),
    "D:/RStudio/work/module2_gwas_gene_risk"
  ))

  for (candidate in candidates) {
    config <- file.path(candidate, "_config.R")
    if (file.exists(config)) return(candidate)
  }

  stop(
    "Cannot locate module2_gwas_gene_risk/_config.R. ",
    "Run from the project root or use the absolute script path.",
    call. = FALSE
  )
}

gwas_gene_code_dir <- find_gwas_gene_code_dir()
source(file.path(gwas_gene_code_dir, "_config.R"))
