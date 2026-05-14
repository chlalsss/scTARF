# module4 / Phase 1.1: Granja input inventory and QC
# ------------------------------------------------------------------
# This script standardizes lightweight metadata and writes a manifest for
# the large Granja RNA / ATAC-GeneActivity inputs. It intentionally does not
# duplicate the large matrices inside module4.

suppressPackageStartupMessages({
  library(Matrix)
  library(SummarizedExperiment)
  library(data.table)
})

ROOT <- "D:/RStudio/work/module4_traitaware_pseudopair"
BLOOD_DIR <- "D:/SCAVENGE-reproducibility-main/blood_data"
OUT_DIR <- file.path(ROOT, "01_data_prep/output")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PATHS <- list(
  rna_rse = file.path(BLOOD_DIR, "scRNA-Healthy-Hematopoiesis-191120.rds"),
  gene_activity = file.path(BLOOD_DIR, "multiomics_output_v2/gene_activity_matrix.rds"),
  atac_meta = file.path(BLOOD_DIR, "multiomics_output_v2/atac_metadata.rds"),
  eval_meta = file.path(ROOT, "04_evaluation/inputs/granja_cell_meta.rds"),
  baseline_trs = file.path(ROOT, "04_evaluation/inputs/granja_baseline_TRS_matrix.rds"),
  scavenge_meta_tsv = file.path(
    BLOOD_DIR,
    "scavenge_out_small_all_traits/SCAVENGE_small_MONO.tsv.gz"
  )
)

stop_missing <- function(paths) {
  missing <- names(paths)[!file.exists(unlist(paths))]
  if (length(missing) > 0) {
    stop("Missing required input(s): ", paste(missing, collapse = ", "))
  }
}

matrix_nnz <- function(x) {
  if (inherits(x, "sparseMatrix")) return(length(x@x))
  if (is.matrix(x)) return(sum(x != 0))
  NA_integer_
}

qc_one_matrix <- function(name, x, orientation, path) {
  data.frame(
    object = name,
    path = normalizePath(path, winslash = "/", mustWork = TRUE),
    class = paste(class(x), collapse = ";"),
    orientation = orientation,
    n_rows = nrow(x),
    n_cols = ncol(x),
    n_nonzero = matrix_nnz(x),
    first_rows = paste(head(rownames(x), 3), collapse = ";"),
    first_cols = paste(head(colnames(x), 3), collapse = ";"),
    stringsAsFactors = FALSE
  )
}

cat(">>> Phase 1.1 Granja input preparation\n\n")
stop_missing(PATHS)

cat("[1/5] Loading RNA RangedSummarizedExperiment ...\n")
rna <- readRDS(PATHS$rna_rse)
if (!inherits(rna, "SummarizedExperiment")) {
  stop("RNA object is not a SummarizedExperiment: ", paste(class(rna), collapse = ", "))
}
if (!"counts" %in% SummarizedExperiment::assayNames(rna)) {
  stop("RNA object has no 'counts' assay. Assays: ",
       paste(SummarizedExperiment::assayNames(rna), collapse = ", "))
}
rna_counts <- SummarizedExperiment::assay(rna, "counts")
rna_counts <- as(rna_counts, "CsparseMatrix")
rna_meta <- as.data.frame(SummarizedExperiment::colData(rna), stringsAsFactors = FALSE)
rna_meta$barcode <- colnames(rna)
rownames(rna_meta) <- rna_meta$barcode
if ("BioClassification" %in% colnames(rna_meta)) {
  rna_meta$cell_type <- as.character(rna_meta$BioClassification)
} else {
  warning("RNA meta has no BioClassification column; cell_type set to NA")
  rna_meta$cell_type <- NA_character_
}
rna_feature <- as.data.frame(SummarizedExperiment::rowData(rna), stringsAsFactors = FALSE)
rna_feature$gene <- rownames(rna)
cat("  RNA counts: ", paste(dim(rna_counts), collapse = " x "),
    " (gene x RNA cell)\n", sep = "")
cat("  RNA meta:   ", paste(dim(rna_meta), collapse = " x "), "\n", sep = "")

cat("\n[2/5] Loading ATAC Gene Activity matrix ...\n")
ga <- readRDS(PATHS$gene_activity)
ga <- as(ga, "CsparseMatrix")
cat("  Gene Activity: ", paste(dim(ga), collapse = " x "),
    " (gene x ATAC cell)\n", sep = "")

cat("\n[3/5] Loading ATAC metadata / evaluation metadata ...\n")
atac_meta <- readRDS(PATHS$atac_meta)
if (inherits(atac_meta, "data.table")) atac_meta <- as.data.frame(atac_meta)
eval_meta <- readRDS(PATHS$eval_meta)
if (inherits(eval_meta, "data.table")) eval_meta <- as.data.frame(eval_meta)
baseline_trs <- readRDS(PATHS$baseline_trs)
if (inherits(baseline_trs, "data.table")) baseline_trs <- as.data.frame(baseline_trs)

if (is.null(rownames(eval_meta)) || identical(rownames(eval_meta), as.character(seq_len(nrow(eval_meta))))) {
  if ("barcode" %in% colnames(eval_meta)) rownames(eval_meta) <- eval_meta$barcode
}
if (is.null(rownames(baseline_trs)) || identical(rownames(baseline_trs), as.character(seq_len(nrow(baseline_trs))))) {
  stop("baseline TRS has no meaningful rownames")
}
if (is.null(rownames(atac_meta)) || identical(rownames(atac_meta), as.character(seq_len(nrow(atac_meta))))) {
  bc_col <- intersect(c("barcode", "cell", "cell_id", "Cell", "CellID"), colnames(atac_meta))[1]
  if (!is.na(bc_col)) rownames(atac_meta) <- as.character(atac_meta[[bc_col]])
}
cat("  ATAC meta:    ", paste(dim(atac_meta), collapse = " x "), "\n", sep = "")
cat("  eval meta:    ", paste(dim(eval_meta), collapse = " x "), "\n", sep = "")
cat("  baseline TRS: ", paste(dim(baseline_trs), collapse = " x "), "\n", sep = "")

cat("\n[4/5] Alignment checks ...\n")
atac_cells <- colnames(ga)
rna_cells <- colnames(rna_counts)
ga_eval_common <- intersect(atac_cells, rownames(eval_meta))
ga_trs_common <- intersect(atac_cells, rownames(baseline_trs))
eval_trs_common <- intersect(rownames(eval_meta), rownames(baseline_trs))
rna_eval_overlap <- intersect(rna_cells, rownames(eval_meta))
rna_ga_overlap <- intersect(rna_cells, atac_cells)

alignment <- data.frame(
  check = c(
    "ATAC_GA_cells",
    "RNA_cells",
    "eval_meta_cells",
    "baseline_TRS_cells",
    "GA_intersect_eval_meta",
    "GA_intersect_baseline_TRS",
    "eval_meta_intersect_baseline_TRS",
    "RNA_intersect_eval_meta",
    "RNA_intersect_GA"
  ),
  n = c(
    length(atac_cells),
    length(rna_cells),
    nrow(eval_meta),
    nrow(baseline_trs),
    length(ga_eval_common),
    length(ga_trs_common),
    length(eval_trs_common),
    length(rna_eval_overlap),
    length(rna_ga_overlap)
  ),
  stringsAsFactors = FALSE
)
print(alignment, row.names = FALSE)

if (length(ga_eval_common) < 30000) {
  stop("Too few overlapping ATAC cells between Gene Activity and eval meta")
}
if (length(rna_ga_overlap) > 0) {
  warning("RNA and ATAC cell barcodes overlap; expected non-paired datasets.")
}

cat("\n[5/5] Writing outputs ...\n")
qc_summary <- rbind(
  qc_one_matrix("rna_counts", rna_counts, "gene x RNA_cell", PATHS$rna_rse),
  qc_one_matrix("gene_activity", ga, "gene x ATAC_cell", PATHS$gene_activity),
  qc_one_matrix("baseline_trs", baseline_trs, "ATAC_cell x trait", PATHS$baseline_trs)
)

rna_celltype_summary <- as.data.frame(sort(table(rna_meta$cell_type), decreasing = TRUE),
                                      stringsAsFactors = FALSE)
colnames(rna_celltype_summary) <- c("cell_type", "n_cells")
rna_celltype_summary$cell_type <- as.character(rna_celltype_summary$cell_type)
rna_celltype_summary$pct <- round(100 * rna_celltype_summary$n_cells / sum(rna_celltype_summary$n_cells), 2)

atac_eval_celltype_summary <- as.data.frame(sort(table(eval_meta$cell_type), decreasing = TRUE),
                                            stringsAsFactors = FALSE)
colnames(atac_eval_celltype_summary) <- c("cell_type", "n_cells")
atac_eval_celltype_summary$cell_type <- as.character(atac_eval_celltype_summary$cell_type)
atac_eval_celltype_summary$pct <- round(
  100 * atac_eval_celltype_summary$n_cells / sum(atac_eval_celltype_summary$n_cells), 2
)

manifest <- data.frame(
  key = names(PATHS),
  path = normalizePath(unlist(PATHS), winslash = "/", mustWork = TRUE),
  stringsAsFactors = FALSE
)

saveRDS(rna_meta, file.path(OUT_DIR, "granja_rna_cell_meta.rds"))
saveRDS(rna_feature, file.path(OUT_DIR, "granja_rna_features.rds"))
saveRDS(eval_meta, file.path(OUT_DIR, "granja_atac_eval_cell_meta.rds"))
saveRDS(alignment, file.path(OUT_DIR, "granja_cell_alignment.rds"))

write.csv(manifest, file.path(OUT_DIR, "granja_inputs_manifest.csv"), row.names = FALSE)
write.csv(qc_summary, file.path(OUT_DIR, "granja_input_qc_summary.csv"), row.names = FALSE)
write.csv(alignment, file.path(OUT_DIR, "granja_cell_alignment.csv"), row.names = FALSE)
write.csv(rna_celltype_summary, file.path(OUT_DIR, "granja_rna_cell_types_summary.csv"), row.names = FALSE)
write.csv(atac_eval_celltype_summary,
          file.path(OUT_DIR, "granja_atac_eval_cell_types_summary.csv"),
          row.names = FALSE)

cat("  wrote: ", file.path(OUT_DIR, "granja_inputs_manifest.csv"), "\n", sep = "")
cat("  wrote: ", file.path(OUT_DIR, "granja_input_qc_summary.csv"), "\n", sep = "")
cat("  wrote: ", file.path(OUT_DIR, "granja_cell_alignment.csv"), "\n", sep = "")
cat("  wrote: ", file.path(OUT_DIR, "granja_rna_cell_meta.rds"), "\n", sep = "")
cat("  wrote: ", file.path(OUT_DIR, "granja_rna_features.rds"), "\n", sep = "")
cat("  wrote: ", file.path(OUT_DIR, "granja_atac_eval_cell_meta.rds"), "\n", sep = "")

cat("\nPhase 1.1 OK.\n")
