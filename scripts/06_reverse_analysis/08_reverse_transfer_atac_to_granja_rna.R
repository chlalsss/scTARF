# Reverse Step 2: transfer ATAC signal to Granja RNA cells.
#
# Uses reverse RNA->ATAC mapping from Step 1.
# For each RNA cell, aggregate its top-k ATAC neighbors:
#   - pseudo ATAC profile = weighted Granja ATAC gene activity
#   - pseudo z_atac      = weighted Granja ATAC SCAVENGE TRS
#
# Output filenames intentionally keep module3-v2 compatibility with the
# existing prepare script:
#   granja_pseudo_rna_profile_<trait>.rds
#   granja_pseudo_zrna_<trait>.rds
# The object metadata states that the content is pseudo ATAC.

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
})

if (!exists("ROOT")) ROOT <- Sys.getenv("SCTARF_WORK_ROOT", normalizePath(".", winslash = "/", mustWork = FALSE))
if (!exists("MODULE4_ROOT")) MODULE4_ROOT <- Sys.getenv("SCTARF_NONPAIRED_ROOT", file.path(ROOT, "outputs", "nonpaired"))
if (!exists("REV_ROOT")) REV_ROOT <- Sys.getenv("SCTARF_REVERSE_ROOT", file.path(ROOT, "outputs", "reverse"))

SCRIPT_ENV <- environment()
has_local <- function(x) exists(x, envir = SCRIPT_ENV, inherits = FALSE)

if (!has_local("TRAITS_TO_RUN")) TRAITS_TO_RUN <- c("PLT", "MCH", "RBC")
if (!has_local("REVERSE_TAG")) REVERSE_TAG <- "reverse_granja_rna_main_plt_mch_rbc"
if (!has_local("MAP_DIR")) MAP_DIR <- file.path(REV_ROOT, "outputs", REVERSE_TAG, "mapping")
if (!has_local("TRANSFER_OUT_DIR")) {
  TRANSFER_OUT_DIR <- file.path(REV_ROOT, "outputs", REVERSE_TAG, "transfer")
}

if (!has_local("GA_PATH")) {
  source(file.path("scripts", "_paths.R"))
  GA_PATH <- file.path(sctarf_blood_root(), "multiomics_output_v2", "gene_activity_matrix.rds")
}
if (!has_local("ATAC_TRS_RDS")) {
  ATAC_TRS_RDS <- file.path(MODULE4_ROOT, "04_evaluation/inputs/granja_baseline_TRS_matrix.rds")
}
if (!has_local("FORCE_RETRANSFER")) FORCE_RETRANSFER <- FALSE
if (!has_local("ALLOW_OUTSIDE_REV_ROOT")) ALLOW_OUTSIDE_REV_ROOT <- FALSE

stop_missing <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0L) {
    stop("Missing required file(s):\n  ", paste(missing, collapse = "\n  "), call. = FALSE)
  }
}

assert_write_under_rev <- function(paths) {
  if (isTRUE(ALLOW_OUTSIDE_REV_ROOT)) return(invisible(TRUE))
  root <- normalizePath(REV_ROOT, winslash = "/", mustWork = TRUE)
  for (path in paths) {
    resolved <- normalizePath(path, winslash = "/", mustWork = FALSE)
    if (!startsWith(resolved, paste0(root, "/")) && !identical(resolved, root)) {
      stop("Refusing to write outside reverse workspace: ", path,
           "\nSet ALLOW_OUTSIDE_REV_ROOT <- TRUE only if this is intentional.",
           call. = FALSE)
    }
  }
  invisible(TRUE)
}

make_weight_matrix <- function(mapping, atac_cells) {
  n_rna <- length(mapping$rna_barcodes)
  k <- mapping$k
  neighbor_idx <- match(as.vector(mapping$atac_neighbor_barcode), atac_cells)
  if (anyNA(neighbor_idx)) {
    stop("Some ATAC neighbor barcodes are not present in the ATAC TRS/GA cell list.",
         call. = FALSE)
  }
  sparseMatrix(
    i = as.integer(neighbor_idx),
    j = rep(seq_len(n_rna), times = k),
    x = as.numeric(mapping$neighbor_weight),
    dims = c(length(atac_cells), n_rna),
    dimnames = list(atac_cells, mapping$rna_barcodes)
  )
}

vector_stats <- function(x) {
  c(
    min = min(x, na.rm = TRUE),
    q05 = as.numeric(quantile(x, 0.05, na.rm = TRUE)),
    median = median(x, na.rm = TRUE),
    mean = mean(x, na.rm = TRUE),
    q95 = as.numeric(quantile(x, 0.95, na.rm = TRUE)),
    max = max(x, na.rm = TRUE)
  )
}

cat(">>> Reverse Step 2: Granja ATAC -> Granja RNA transfer\n\n")
cat("Traits: ", paste(TRAITS_TO_RUN, collapse = ", "), "\n", sep = "")
cat("Mapping dir: ", MAP_DIR, "\n", sep = "")
cat("Output dir : ", TRANSFER_OUT_DIR, "\n\n", sep = "")

assert_write_under_rev(TRANSFER_OUT_DIR)
dir.create(TRANSFER_OUT_DIR, recursive = TRUE, showWarnings = FALSE)
required <- c(
  GA_PATH,
  ATAC_TRS_RDS,
  file.path(MAP_DIR, paste0("granja_rna_to_atac_mapping_", TRAITS_TO_RUN, ".rds"))
)
stop_missing(required)

cat("[1/3] Loading Granja ATAC Gene Activity and baseline TRS ...\n")
ga <- readRDS(GA_PATH)
ga <- as(ga, "CsparseMatrix")
atac_trs <- readRDS(ATAC_TRS_RDS)
if (inherits(atac_trs, "data.table")) atac_trs <- as.data.frame(atac_trs)
if (is.null(rownames(atac_trs))) {
  bc_col <- intersect(c("barcode", "cell", "cell_id"), colnames(atac_trs))[1]
  if (!is.na(bc_col)) {
    rownames(atac_trs) <- as.character(atac_trs[[bc_col]])
    atac_trs[[bc_col]] <- NULL
  }
}
atac_cells <- intersect(colnames(ga), rownames(atac_trs))
ga <- ga[, atac_cells, drop = FALSE]
atac_trs <- atac_trs[atac_cells, , drop = FALSE]
cat("  GA      : ", paste(dim(ga), collapse = " x "), " (gene x ATAC cell)\n", sep = "")
cat("  ATAC TRS: ", paste(dim(atac_trs), collapse = " x "), " (ATAC cell x trait)\n", sep = "")

summary_rows <- list()

cat("[2/3] Running transfer per trait ...\n")
for (trait in TRAITS_TO_RUN) {
  cat("\n=== Trait ", trait, " ===\n", sep = "")
  profile_out <- file.path(TRANSFER_OUT_DIR, paste0("granja_pseudo_rna_profile_", trait, ".rds"))
  z_out <- file.path(TRANSFER_OUT_DIR, paste0("granja_pseudo_zrna_", trait, ".rds"))
  z_csv <- file.path(TRANSFER_OUT_DIR, paste0("granja_pseudo_zrna_", trait, ".csv"))
  if (!FORCE_RETRANSFER && file.exists(profile_out) && file.exists(z_out) && file.exists(z_csv)) {
    cat("  existing transfer outputs found; skipping\n")
    obj <- readRDS(profile_out)
    z_obj <- readRDS(z_out)
    summary_rows[[trait]] <- data.frame(
      trait = trait,
      n_rna_cells = length(obj$reference_cells),
      n_g_trait_used = nrow(obj$matrix),
      k = NA_integer_,
      skipped_existing = TRUE,
      stringsAsFactors = FALSE
    )
    rm(obj, z_obj)
    next
  }

  map_file <- file.path(MAP_DIR, paste0("granja_rna_to_atac_mapping_", trait, ".rds"))
  mapping <- readRDS(map_file)
  genes <- mapping$g_trait_used
  missing_genes <- setdiff(genes, rownames(ga))
  if (length(missing_genes) > 0L) {
    stop("ATAC GA matrix is missing mapping genes for ", trait, ": ",
         paste(head(missing_genes, 10), collapse = ", "), call. = FALSE)
  }
  if (!trait %in% colnames(atac_trs)) {
    stop("ATAC baseline TRS is missing trait column: ", trait, call. = FALSE)
  }

  cat("  building ATAC-neighbor weight matrix ...\n")
  W <- make_weight_matrix(mapping, atac_cells)
  w_colsum <- Matrix::colSums(W)
  if (!all(abs(w_colsum - 1) < 1e-6)) {
    stop("Neighbor weights do not sum to 1 for all RNA cells.", call. = FALSE)
  }

  cat("  transferring pseudo ATAC profile: gene x RNA ...\n")
  pseudo_atac <- ga[genes, rownames(W), drop = FALSE] %*% W
  pseudo_atac <- as(pseudo_atac, "CsparseMatrix")
  rownames(pseudo_atac) <- genes
  colnames(pseudo_atac) <- mapping$rna_barcodes

  cat("  transferring pseudo z_atac from Granja SCAVENGE TRS ...\n")
  atac_score <- as.numeric(atac_trs[rownames(W), trait])
  if (anyNA(atac_score) || any(!is.finite(atac_score))) {
    replacement <- median(atac_score[is.finite(atac_score)], na.rm = TRUE)
    warning("Some ATAC ", trait, " TRS values are missing; replacing with median = ",
            signif(replacement, 4))
    atac_score[!is.finite(atac_score) | is.na(atac_score)] <- replacement
  }
  pseudo_z <- as.numeric(Matrix::crossprod(W, atac_score))
  names(pseudo_z) <- mapping$rna_barcodes

  max_weight <- apply(mapping$neighbor_weight, 1L, max)
  entropy <- -rowSums(mapping$neighbor_weight * log(pmax(mapping$neighbor_weight, 1e-12)))
  pseudo_z_dt <- data.table(
    barcode = mapping$rna_barcodes,
    trait = trait,
    pseudo_z_rna = pseudo_z,
    pseudo_z_atac = pseudo_z,
    top1_distance = mapping$neighbor_distance[, 1],
    top1_cos_dist = mapping$neighbor_distance[, 1]^2 / 2,
    max_neighbor_weight = max_weight,
    neighbor_weight_entropy = entropy
  )

  pseudo_profile_obj <- list(
    trait = trait,
    trait_id = mapping$trait_id,
    matrix = pseudo_atac,
    orientation = "gene x RNA_cell",
    genes = genes,
    atac_barcodes = mapping$rna_barcodes,
    rna_barcodes = mapping$rna_barcodes,
    reference_cells = mapping$rna_barcodes,
    source = "weighted average of top-k Granja ATAC gene-activity profiles",
    compatibility_note = "Saved under granja_pseudo_rna_profile_* for reuse by existing module3 v2 input preparer.",
    mapping_file = map_file,
    no_celltype_used = TRUE,
    created_at = as.character(Sys.time())
  )
  pseudo_z_obj <- list(
    trait = trait,
    trait_id = mapping$trait_id,
    scores = pseudo_z_dt,
    source = "weighted average of top-k Granja ATAC SCAVENGE TRS",
    compatibility_note = "Column pseudo_z_rna mirrors pseudo_z_atac for reuse by existing seed builder.",
    mapping_file = map_file,
    no_celltype_used = TRUE,
    created_at = as.character(Sys.time())
  )

  saveRDS(pseudo_profile_obj, profile_out)
  saveRDS(pseudo_z_obj, z_out)
  fwrite(pseudo_z_dt, z_csv)
  cat("  wrote: ", profile_out, "\n", sep = "")
  cat("  wrote: ", z_out, "\n", sep = "")
  cat("  wrote: ", z_csv, "\n", sep = "")

  z_stats <- vector_stats(pseudo_z)
  summary_rows[[trait]] <- data.frame(
    trait = trait,
    trait_id = mapping$trait_id,
    n_rna_cells = length(mapping$rna_barcodes),
    n_atac_cells = length(atac_cells),
    n_g_trait_used = length(genes),
    k = mapping$k,
    pseudo_atac_nnz = length(pseudo_atac@x),
    pseudo_atac_density = length(pseudo_atac@x) / max(prod(dim(pseudo_atac)), 1),
    pseudo_z_min = z_stats[["min"]],
    pseudo_z_q05 = z_stats[["q05"]],
    pseudo_z_median = z_stats[["median"]],
    pseudo_z_mean = z_stats[["mean"]],
    pseudo_z_q95 = z_stats[["q95"]],
    pseudo_z_max = z_stats[["max"]],
    weight_colsum_min = min(w_colsum),
    weight_colsum_max = max(w_colsum),
    skipped_existing = FALSE,
    no_celltype_used = TRUE,
    stringsAsFactors = FALSE
  )

  rm(mapping, W, pseudo_atac, pseudo_profile_obj, pseudo_z, pseudo_z_dt, pseudo_z_obj)
  invisible(gc(verbose = FALSE))
}

cat("\n[3/3] Writing transfer summary ...\n")
summary_df <- rbindlist(summary_rows, fill = TRUE)
fwrite(summary_df, file.path(TRANSFER_OUT_DIR, "reverse_transfer_summary.csv"))

cat("\n=== Reverse transfer finished ===\n")
print(summary_df)
cat("\nOutputs are in: ", TRANSFER_OUT_DIR, "\n", sep = "")
