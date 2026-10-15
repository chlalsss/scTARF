# Run this script ONCE in a clean R session before running anything else in
# module3_multiome. Per project convention, install.packages() and library()
# must not share an R session — restart R after this finishes.

cran_pkgs <- c(
  "Seurat", "Signac", "pROC", "igraph", "uwot", "reticulate", "patchwork",
  "data.table",
  "ggplot2",
  "Matrix",
  "R.utils",
  "jsonlite",
  "scales",
  "tibble",
  "dplyr",
  "tidyr",
  "readr",
  "stringr",
  "hdf5r",
  "remotes",
  "irlba",     # used by SCAVENGE::do_lsi
  "RANN",      # used by SCAVENGE / general k-NN
  "matrixStats" # used by 04_dual_modality §2.3 / §2.4 row-wise stats

)

bioc_pkgs <- c(
  "BiocNeighbors", "ensembldb", "EnsDb.Hsapiens.v75",
  "BiocParallel",
  "GEOquery",
  "GenomicRanges",
  "rtracklayer",
  "rhdf5",
  "bluster",
  "SummarizedExperiment",
  "chromVAR",                       # SCAVENGE upstream: GC bias + bg_peaks + weighted deviations
  "gchromVAR",
  "BSgenome.Hsapiens.UCSC.hg19",
  "BSgenome.Hsapiens.UCSC.hg38",   # BMMC NeurIPS 2021 is hg38; replace blood_small's hg19
  # added 2026-05-03 for module 4 (dual modality) §2.2 G_trait construction
  "GenomicFeatures",
  "AnnotationDbi",
  "GenomeInfoDb",
  "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "org.Hs.eg.db"
)

install_cran <- function(pkgs) {
  todo <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(todo) == 0) {
    message("CRAN: all present (", length(pkgs), ")")
    return(invisible())
  }
  message("CRAN install: ", paste(todo, collapse = ", "))
  install.packages(todo, repos = "https://cloud.r-project.org")
}

install_bioc <- function(pkgs) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  todo <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(todo) == 0) {
    message("Bioc: all present (", length(pkgs), ")")
    return(invisible())
  }
  message("Bioc install: ", paste(todo, collapse = ", "))
  BiocManager::install(todo, ask = FALSE, update = FALSE)
}

install_cran(cran_pkgs)
install_bioc(bioc_pkgs)

install_scavenge <- function() {
  if (requireNamespace("SCAVENGE", quietly = TRUE)) {
    message("SCAVENGE: present")
    return(invisible())
  }

  repo <- Sys.getenv("SCTARF_SCAVENGE_GITHUB", "")
  if (nzchar(repo)) {
    if (!requireNamespace("remotes", quietly = TRUE)) {
      install.packages("remotes", repos = "https://cloud.r-project.org")
    }
    message("Installing SCAVENGE from GitHub repo: ", repo)
    remotes::install_github(repo, upgrade = "never")
  }

  if (!requireNamespace("SCAVENGE", quietly = TRUE)) {
    stop(
      "SCAVENGE is required by paired baseline and sparsification scripts. ",
      "Install SCAVENGE in this R library, or set SCTARF_SCAVENGE_GITHUB=",
      "<owner/repo> before running 00_install_deps.R.",
      call. = FALSE
    )
  }
}

install_scavenge()

message("\nDone. Now restart R before running any other module3 script.")
