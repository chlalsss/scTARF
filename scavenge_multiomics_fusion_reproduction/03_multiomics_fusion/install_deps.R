# Run this script ONCE in a clean R session before running anything else in
# module3_multiome. Per project convention, install.packages() and library()
# must not share an R session — restart R after this finishes.

cran_pkgs <- c(
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
  "irlba",     # used by SCAVENGE::do_lsi
  "RANN",      # used by SCAVENGE / general k-NN
  "BiocParallel",
  "matrixStats" # used by 04_dual_modality §2.3 / §2.4 row-wise stats

)

bioc_pkgs <- c(
  "GEOquery",
  "GenomicRanges",
  "rtracklayer",
  "rhdf5",
  "SummarizedExperiment",
  "chromVAR",                       # SCAVENGE upstream: GC bias + bg_peaks + weighted deviations
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

message("\nDone. Now restart R before running any other module3 script.")
