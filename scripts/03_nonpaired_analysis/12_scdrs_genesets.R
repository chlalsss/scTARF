# module4 / Phase 1.2a: prepare Granja scDRS gene sets
# ------------------------------------------------------------------
# Filter the module3 UKB MAGMA .gs file to the three traits used by
# module4: MONO, PLT, LYM. Also report overlap with Granja RNA and ATAC GA.

suppressPackageStartupMessages({
  library(Matrix)
  library(SummarizedExperiment)
})

source(file.path("scripts", "_paths.R"))
ROOT <- sctarf_nonpaired_root()
OUT_DIR <- file.path(ROOT, "01_data_prep/03_rna_scdrs_python")
DATA_OUT <- file.path(ROOT, "01_data_prep/output")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_OUT, recursive = TRUE, showWarnings = FALSE)

GS_FULL <- file.path("scripts", "01_paired_analysis", "magma_10kb_top1000_zscore.74_traits.rv1.gs")
RNA_RSE <- file.path(sctarf_blood_root(), "scRNA-Healthy-Hematopoiesis-191120.rds")
GA_RDS <- file.path(sctarf_blood_root(), "multiomics_output_v2", "gene_activity_matrix.rds")

TRAITS <- c(
  MONO = "UKB_460K.blood_MONOCYTE_COUNT",
  PLT  = "UKB_460K.blood_PLATELET_COUNT",
  LYM  = "UKB_460K.blood_LYMPHOCYTE_COUNT"
)

parse_gs <- function(path) {
  lines <- readLines(path, warn = FALSE)
  header <- strsplit(lines[1], "\t", fixed = TRUE)[[1]]
  if (!identical(header[1:2], c("TRAIT", "GENESET"))) {
    stop("Unexpected .gs header: ", paste(header, collapse = "\t"))
  }
  out <- list()
  for (line in lines[-1]) {
    parts <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(parts) < 2) next
    pairs <- strsplit(strsplit(parts[2], ",", fixed = TRUE)[[1]], ":", fixed = TRUE)
    genes <- vapply(pairs, `[`, character(1), 1)
    weights <- as.numeric(vapply(
      pairs,
      function(x) if (length(x) >= 2) x[2] else "1",
      character(1)
    ))
    out[[parts[1]]] <- list(raw = parts[2], genes = genes, weights = weights)
  }
  out
}

cat(">>> Phase 1.2a Granja scDRS gene set preparation\n\n")
stopifnot(file.exists(GS_FULL), file.exists(RNA_RSE), file.exists(GA_RDS))

cat("[1/4] Loading gene sets ...\n")
gs <- parse_gs(GS_FULL)
missing_traits <- setdiff(unname(TRAITS), names(gs))
if (length(missing_traits) > 0) {
  stop("Missing trait(s) in .gs file: ", paste(missing_traits, collapse = ", "))
}

cat("[2/4] Loading feature names for overlap QC ...\n")
rna <- readRDS(RNA_RSE)
rna_genes <- rownames(rna)
ga <- readRDS(GA_RDS)
ga_genes <- rownames(ga)

cat("[3/4] Writing filtered .gs ...\n")
out_gs <- file.path(OUT_DIR, "granja_3traits.gs")
gs_lines <- c(
  "TRAIT\tGENESET",
  vapply(names(TRAITS), function(short) {
    trait_id <- TRAITS[[short]]
    paste(trait_id, gs[[trait_id]]$raw, sep = "\t")
  }, character(1))
)
writeLines(gs_lines, out_gs, useBytes = TRUE)

cat("[4/4] Writing overlap summary ...\n")
summary <- do.call(rbind, lapply(names(TRAITS), function(short) {
  trait_id <- TRAITS[[short]]
  genes <- gs[[trait_id]]$genes
  data.frame(
    trait = short,
    trait_id = trait_id,
    n_genes = length(genes),
    n_overlap_rna = sum(genes %in% rna_genes),
    n_overlap_ga = sum(genes %in% ga_genes),
    n_overlap_both = sum(genes %in% intersect(rna_genes, ga_genes)),
    weight_min = min(gs[[trait_id]]$weights, na.rm = TRUE),
    weight_median = median(gs[[trait_id]]$weights, na.rm = TRUE),
    weight_max = max(gs[[trait_id]]$weights, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))

write.csv(summary, file.path(DATA_OUT, "granja_scdrs_geneset_summary.csv"),
          row.names = FALSE)

cat("  wrote: ", out_gs, "\n", sep = "")
cat("  wrote: ", file.path(DATA_OUT, "granja_scdrs_geneset_summary.csv"), "\n", sep = "")
print(summary, row.names = FALSE)

cat("\nPhase 1.2a OK.\n")
