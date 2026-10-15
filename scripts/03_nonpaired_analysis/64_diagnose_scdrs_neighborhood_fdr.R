# module4 diagnostic: scDRS neighborhood-level FDR
# -----------------------------------------------------------------------------
# Purpose:
#   Use the already-computed scDRS full_score files to test whether unsupervised
#   RNA neighborhoods are trait-associated. This is the second layer after
#   per-cell FDR.
#
# Key idea:
#   For each unsupervised CITE RNA neighborhood/cluster:
#     observed = mean(real scDRS norm_score)
#     null_b   = mean(control norm_score b)
#     p = (1 + sum(null_b >= observed)) / (1 + n_control)
#   Then BH-adjust p-values across neighborhoods for each trait.
#
# Method boundary:
#   - Cell type is not used to define neighborhoods, FDR, mapping weights, or
#     candidate seeds.
#   - Cell type is used only for post hoc diagnostics.
#   - This script does not replace the main pipeline yet.
#
# Run in RStudio:
#   Run from the repository root with optional SCTARF_* path overrides.
#   source("04_evaluation/diagnose_scdrs_neighborhood_fdr.R")
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(ggplot2)
})

if (!requireNamespace("irlba", quietly = TRUE)) {
  stop("Missing irlba. In a separate R session run: install.packages('irlba')")
}
if (!requireNamespace("igraph", quietly = TRUE)) {
  stop("Missing igraph. In a separate R session run: install.packages('igraph')")
}
if (!requireNamespace("RANN", quietly = TRUE)) {
  stop("Missing RANN. In a separate R session run: install.packages('RANN')")
}
if (!requireNamespace("pROC", quietly = TRUE)) {
  stop("Missing pROC. In a separate R session run: install.packages('pROC')")
}

if (!exists("ROOT")) ROOT <- Sys.getenv("SCTARF_NONPAIRED_ROOT", file.path("outputs", "nonpaired"))
if (!exists("MAP_DIR")) MAP_DIR <- file.path(ROOT, "02_trait_mapping/output")
if (!exists("ATAC_META_FILE")) ATAC_META_FILE <- file.path(ROOT, "04_evaluation/inputs/granja_cell_meta.rds")
if (!exists("CITE_SCORE_RDS")) {
  CITE_SCORE_RDS <- file.path(ROOT, "01_data_prep_external_rna/output/cite_rna_scdrs_scores.rds")
}
if (!exists("CITE_MTX")) CITE_MTX <- file.path(MAP_DIR, "cite_extract/cite_trait_counts.mtx")
if (!exists("CITE_GENES")) CITE_GENES <- file.path(MAP_DIR, "cite_extract/cite_trait_counts_genes.txt")
if (!exists("CITE_CELLS")) CITE_CELLS <- file.path(MAP_DIR, "cite_extract/cite_trait_counts_cells.txt")
if (!exists("FULL_SCORE_DIR")) {
  FULL_SCORE_DIR <- file.path(ROOT, "01_data_prep_external_rna/03_cite_scdrs_python/output")
}
if (!exists("OUT_DIR")) OUT_DIR <- file.path(ROOT, "04_evaluation/module4_v2_diagnostics/neighborhood_fdr")
FIG_DIR <- file.path(OUT_DIR, "figures")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

if (!exists("TRAITS")) TRAITS <- c("MONO", "PLT", "LYM")
if (!exists("TRAIT_IDS")) {
  TRAIT_IDS <- c(
    MONO = "UKB_460K.blood_MONOCYTE_COUNT",
    PLT  = "UKB_460K.blood_PLATELET_COUNT",
    LYM  = "UKB_460K.blood_LYMPHOCYTE_COUNT"
  )
}

if (!exists("RNA_REFERENCE_LABEL")) RNA_REFERENCE_LABEL <- "CITE RNA"
if (!exists("RNA_PCA_DIMS")) RNA_PCA_DIMS <- 30L
if (!exists("RNA_KNN_K")) RNA_KNN_K <- 30L
if (!exists("LEIDEN_RESOLUTION")) LEIDEN_RESOLUTION <- 0.005
if (!exists("GROUP_FDR_THRESHOLDS")) GROUP_FDR_THRESHOLDS <- c(0.05, 0.10, 0.20)
if (!exists("TOP_N_ATAC")) TOP_N_ATAC <- 350L

if (!exists("ATAC_TARGETS")) ATAC_TARGETS <- list(
  MONO = c("11_CD14.Mono.1", "12_CD14.Mono.2"),
  PLT  = c("02_Early.Eryth", "03_Late.Eryth", "04_Early.Baso"),
  LYM  = c(
    "16_Pre.B", "17_B", "18_Plasma",
    "19_CD8.N", "20_CD4.N1", "21_CD4.N2", "22_CD4.M",
    "23_CD8.EM", "24_CD8.CM", "25_NK"
  )
)

if (!exists("CITE_TARGETS")) CITE_TARGETS <- list(
  MONO = c("CD14+ Mono", "CD16+ Mono"),
  PLT = c("MK/E prog", "Proerythroblast", "Erythroblast", "Normoblast", "Reticulocyte"),
  LYM = c(
    "CD4+ T activated", "CD4+ T naive", "CD4+ T activated integrinB7+",
    "CD4+ T CD314+ CD45RA+", "CD8+ T naive",
    "CD8+ T naive CD127+ CD26- CD101-", "CD8+ T CD49f+",
    "CD8+ T CD57+ CD45RA+", "CD8+ T CD57+ CD45RO+",
    "CD8+ T CD69+ CD45RA+", "CD8+ T CD69+ CD45RO+",
    "CD8+ T TIGIT+ CD45RA+", "CD8+ T TIGIT+ CD45RO+",
    "T reg", "T prog cycling", "MAIT", "gdT CD158b+", "gdT TCRVD2+",
    "dnT", "NK", "NK CD158e1+", "ILC", "ILC1", "Lymph prog",
    "Naive CD20+ B IGKC+", "Naive CD20+ B IGKC-",
    "B1 B IGKC+", "B1 B IGKC-", "Transitional B",
    "Plasma cell IGKC+", "Plasma cell IGKC-",
    "Plasmablast IGKC+", "Plasmablast IGKC-"
  )
)

if (!exists("CITE_MYELOID")) {
  CITE_MYELOID <- c("CD14+ Mono", "CD16+ Mono", "G/M prog", "cDC1", "cDC2", "pDC")
}

stop_missing <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0L) {
    stop("Missing required file(s):\n  ", paste(missing, collapse = "\n  "),
         call. = FALSE)
  }
}

compute_auc <- function(score, label) {
  keep <- is.finite(score) & !is.na(label)
  score <- score[keep]
  label <- label[keep]
  if (sum(label == 1L) < 10L || sum(label == 0L) < 10L) return(NA_real_)
  if (length(unique(score)) < 2L) return(NA_real_)
  as.numeric(pROC::auc(pROC::roc(
    label, score, quiet = TRUE, direction = "<", levels = c(0L, 1L)
  )))
}

safe_div <- function(num, den) {
  if (!is.finite(den) || den <= 0) return(NA_real_)
  num / den
}

top_table_label <- function(x, n = 10L) {
  tab <- sort(table(x), decreasing = TRUE)
  if (length(tab) == 0L) return("")
  paste(paste0(names(head(tab, n)), ":", as.integer(head(tab, n))), collapse = "; ")
}

top_weight_label <- function(dt, type_col = "cite_cell_type",
                             weight_col = "neighbor_weight", n = 10L) {
  if (nrow(dt) == 0L) return("")
  s <- dt[, .(total_weight = sum(get(weight_col), na.rm = TRUE)),
          by = type_col][order(-total_weight)]
  paste(paste0(head(s[[type_col]], n), ":", round(head(s$total_weight, n), 1)),
        collapse = "; ")
}

load_atac_meta <- function(path) {
  meta <- readRDS(path)
  if (inherits(meta, "data.table")) meta <- as.data.frame(meta)
  if (is.null(rownames(meta)) ||
      identical(rownames(meta), as.character(seq_len(nrow(meta))))) {
    if ("barcode" %in% colnames(meta)) rownames(meta) <- meta$barcode
  }
  if (!"cell_type" %in% colnames(meta)) {
    stop("ATAC meta missing cell_type column", call. = FALSE)
  }
  meta
}

log_normalize_counts <- function(counts_g_by_c, lib_size, scale_factor = 1e4) {
  lib <- lib_size[colnames(counts_g_by_c)]
  if (anyNA(lib) || any(!is.finite(lib))) {
    stop("Missing CITE library sizes for extracted RNA matrix.", call. = FALSE)
  }
  lib[lib < 1] <- 1
  norm <- Matrix::t(Matrix::t(counts_g_by_c) / lib) * scale_factor
  log1p(norm)
}

build_rna_pca <- function(cite_counts, cite_meta, cite_cells, dims = 30L) {
  lib_size <- cite_meta$scdrs_input_total_counts[match(cite_cells, cite_meta$cell)]
  names(lib_size) <- cite_cells
  rna_log <- log_normalize_counts(cite_counts, lib_size)

  # Avoid materializing the full cell x gene dense matrix. Centering/scaling is
  # still handled by irlba, but the sparse matrix can stay sparse.
  x <- Matrix::t(rna_log)
  gene_mean <- Matrix::colMeans(x)
  gene_m2 <- Matrix::colMeans(x * x)
  gene_sd <- sqrt(pmax(gene_m2 - gene_mean * gene_mean, 0))
  gene_sd[!is.finite(gene_sd) | gene_sd < 1e-8] <- 1

  dims <- min(dims, nrow(x) - 1L, ncol(x) - 1L)
  pca <- irlba::prcomp_irlba(x, n = dims, center = gene_mean, scale. = gene_sd)
  z <- pca$x
  rownames(z) <- cite_cells
  colnames(z) <- paste0("RNA_PC", seq_len(ncol(z)))
  z
}

build_knn_graph <- function(z, k = 30L) {
  k <- min(k + 1L, nrow(z))
  nn <- RANN::nn2(z, z, k = k)
  idx <- nn$nn.idx[, -1, drop = FALSE]
  d <- nn$nn.dists[, -1, drop = FALSE]
  n <- nrow(z)
  edges <- data.table(
    from = rep(seq_len(n), times = ncol(idx)),
    to = as.vector(idx),
    dist = as.vector(d)
  )
  edges <- edges[from != to]
  edges[, a := pmin(from, to)]
  edges[, b := pmax(from, to)]
  edges <- unique(edges[, .(from = a, to = b)])
  igraph::graph_from_data_frame(edges, directed = FALSE, vertices = data.frame(name = seq_len(n)))
}

cluster_graph <- function(g, resolution = 1.0) {
  if ("cluster_leiden" %in% getNamespaceExports("igraph")) {
    cl <- igraph::cluster_leiden(g, resolution = resolution)
  } else {
    cl <- igraph::cluster_louvain(g)
  }
  as.integer(igraph::membership(cl))
}

full_score_file <- function(trait) {
  file.path(FULL_SCORE_DIR, paste0(TRAIT_IDS[[trait]], ".full_score.gz"))
}

read_full_score_needed <- function(path) {
  dt <- fread(path)
  if (!"cell" %in% colnames(dt) && "barcode" %in% colnames(dt)) {
    setnames(dt, "barcode", "cell")
  }
  if (!"pval" %in% colnames(dt) && "mc_pval" %in% colnames(dt)) {
    dt[, pval := mc_pval]
  }
  needed_base <- c("cell", "norm_score", "pval")
  ctrl_cols <- grep("^ctrl_norm_score_", colnames(dt), value = TRUE)
  if (!all(needed_base %in% colnames(dt)) || length(ctrl_cols) == 0L) {
    stop("Unexpected full_score format: ", path, call. = FALSE)
  }
  dt[, c(needed_base, ctrl_cols), with = FALSE]
}

test_groups_with_controls <- function(full_dt, group_id) {
  ctrl_cols <- grep("^ctrl_norm_score_", colnames(full_dt), value = TRUE)
  score <- full_dt$norm_score
  group_levels <- sort(unique(group_id))
  rows <- vector("list", length(group_levels))

  for (i in seq_along(group_levels)) {
    g <- group_levels[i]
    ix <- which(group_id == g)
    observed <- mean(score[ix], na.rm = TRUE)
    ctrl_mean <- colMeans(as.matrix(full_dt[ix, ..ctrl_cols]), na.rm = TRUE)
    p_mc <- (1 + sum(ctrl_mean >= observed, na.rm = TRUE)) /
      (1 + sum(is.finite(ctrl_mean)))
    rows[[i]] <- data.frame(
      group_id = as.integer(g),
      n_scored_cells = length(ix),
      observed_mean_score = observed,
      null_mean = mean(ctrl_mean, na.rm = TRUE),
      null_sd = stats::sd(ctrl_mean, na.rm = TRUE),
      p_mc = p_mc,
      stringsAsFactors = FALSE
    )
  }

  out <- rbindlist(rows)
  out[, qval := p.adjust(p_mc, method = "BH")]
  out[, z_vs_null := (observed_mean_score - null_mean) / pmax(null_sd, 1e-8)]
  out
}

cat(">>> Diagnose scDRS neighborhood-level FDR\n\n")
cat("Neighborhood source: ", RNA_REFERENCE_LABEL,
    " PCA/kNN/unsupervised graph cluster\n", sep = "")
cat("RNA_PCA_DIMS = ", RNA_PCA_DIMS,
    "; RNA_KNN_K = ", RNA_KNN_K,
    "; LEIDEN_RESOLUTION = ", LEIDEN_RESOLUTION, "\n\n", sep = "")

stop_missing(c(ATAC_META_FILE, CITE_SCORE_RDS, CITE_MTX, CITE_GENES, CITE_CELLS))

cat("[1/5] Loading metadata and CITE extracted counts ...\n")
atac_meta <- load_atac_meta(ATAC_META_FILE)
cite_meta <- readRDS(CITE_SCORE_RDS)
cite_meta <- as.data.frame(cite_meta)
if (!all(c("cell", "cell_type", "scdrs_input_total_counts") %in% colnames(cite_meta))) {
  stop("CITE score RDS must contain cell, cell_type, and scdrs_input_total_counts.",
       call. = FALSE)
}
rownames(cite_meta) <- cite_meta$cell

cite_counts <- Matrix::readMM(CITE_MTX)
cite_counts <- as(cite_counts, "CsparseMatrix")
cite_genes <- readLines(CITE_GENES, warn = FALSE)
cite_cells <- readLines(CITE_CELLS, warn = FALSE)
rownames(cite_counts) <- cite_genes
colnames(cite_counts) <- cite_cells
cite_meta <- cite_meta[match(cite_cells, cite_meta$cell), , drop = FALSE]
if (anyNA(cite_meta$cell)) stop("CITE metadata does not match extracted CITE cells.", call. = FALSE)

cat("[2/5] Building unsupervised CITE RNA neighborhoods ...\n")
rna_z <- build_rna_pca(cite_counts, cite_meta, cite_cells, dims = RNA_PCA_DIMS)
g <- build_knn_graph(rna_z, k = RNA_KNN_K)
group_id <- cluster_graph(g, resolution = LEIDEN_RESOLUTION)
if (length(group_id) != length(cite_cells)) {
  stop("Cluster membership length does not match CITE cells.", call. = FALSE)
}
names(group_id) <- cite_cells

group_base <- data.table(
  cell = cite_cells,
  group_id = group_id,
  cell_type_posthoc = cite_meta$cell_type
)
group_size <- group_base[, .(n_cells = .N), by = group_id][order(group_id)]
fwrite(group_base, file.path(OUT_DIR, "rna_unsupervised_group_membership.csv.gz"))
fwrite(group_size, file.path(OUT_DIR, "rna_unsupervised_group_sizes.csv"))
cat("  groups: ", length(unique(group_id)),
    "; median group size: ", median(group_size$n_cells), "\n", sep = "")

summary_rows <- list()
group_rows <- list()
cell_rows <- list()
source_rows <- list()
atac_by_celltype_rows <- list()

cat("[3/5] Testing group-level scDRS enrichment per trait ...\n")
for (trait in TRAITS) {
  cat("=== Trait ", trait, " ===\n", sep = "")
  map_file <- file.path(MAP_DIR, paste0("granja_cite_mapping_", trait, ".rds"))
  fs_file <- full_score_file(trait)
  stop_missing(c(map_file, fs_file))

  full_dt <- read_full_score_needed(fs_file)
  full_ix <- match(cite_cells, full_dt$cell)
  has_full_score <- !is.na(full_ix)
  if (sum(has_full_score) < 1000L) {
    stop("Too few CITE cells have full_score rows for ", trait, ": ",
         sum(has_full_score), call. = FALSE)
  }
  if (any(!has_full_score)) {
    cat("  note: ", sum(!has_full_score),
        " extracted CITE cells are absent from full_score and will be ignored in group tests.\n",
        sep = "")
  }
  full_dt <- full_dt[full_ix[has_full_score]]
  group_id_for_score <- group_id[has_full_score]

  group_test <- test_groups_with_controls(full_dt, group_id_for_score)
  group_test[, trait := trait]
  group_test <- merge(group_test, group_size, by = "group_id", all.x = TRUE)

  group_posthoc <- group_base[, .(
    top_celltypes_posthoc = top_table_label(cell_type_posthoc, 8L),
    target_fraction_posthoc = mean(cell_type_posthoc %in% CITE_TARGETS[[trait]]),
    myeloid_fraction_posthoc = mean(cell_type_posthoc %in% CITE_MYELOID)
  ), by = group_id]
  group_test <- merge(group_test, group_posthoc, by = "group_id", all.x = TRUE)
  setcolorder(group_test, c("trait", "group_id"))
  group_rows[[trait]] <- group_test

  mapping <- readRDS(map_file)
  atac_bc <- mapping$atac_barcodes
  nn_bc <- mapping$rna_neighbor_barcode
  w <- mapping$neighbor_weight
  k <- ncol(nn_bc)
  nn_group <- matrix(group_id[as.vector(nn_bc)], nrow = nrow(nn_bc), ncol = k)
  nn_type <- matrix(cite_meta$cell_type[match(as.vector(nn_bc), cite_meta$cell)],
                    nrow = nrow(nn_bc), ncol = k)
  nn_target <- nn_type %in% CITE_TARGETS[[trait]]
  nn_myeloid <- nn_type %in% CITE_MYELOID

  source_target_weight <- rowSums(w * nn_target, na.rm = TRUE)
  source_myeloid_weight <- rowSums(w * nn_myeloid, na.rm = TRUE)
  atac_cell_type <- as.character(atac_meta[atac_bc, "cell_type"])
  atac_target <- atac_cell_type %in% ATAC_TARGETS[[trait]]
  atac_target_rate <- mean(atac_target)

  edge_dt_base <- data.table(
    trait = trait,
    atac_barcode = rep(atac_bc, times = k),
    rank = rep(seq_len(k), each = length(atac_bc)),
    cite_barcode = as.vector(nn_bc),
    neighbor_weight = as.vector(w),
    cite_group_id = as.vector(nn_group),
    cite_cell_type = as.vector(nn_type)
  )

  for (thr in GROUP_FDR_THRESHOLDS) {
    sig_groups <- group_test$group_id[group_test$qval <= thr & group_test$observed_mean_score > 0]
    cell_group_sig <- group_id %in% sig_groups
    group_sig_weight <- rowSums(w * (nn_group %in% sig_groups), na.rm = TRUE)
    group_z <- group_test$z_vs_null[match(as.vector(nn_group), group_test$group_id)]
    group_z <- matrix(group_z, nrow = nrow(nn_group), ncol = k)
    group_z[!is.finite(group_z) | group_z < 0] <- 0
    group_z_evidence <- rowSums(w * group_z, na.rm = TRUE)

    score_auc <- compute_auc(group_sig_weight, as.integer(atac_target))
    z_auc <- compute_auc(group_z_evidence, as.integer(atac_target))

    candidate_top <- rep(FALSE, length(atac_bc))
    ok <- which(is.finite(group_sig_weight) & group_sig_weight > 0)
    if (length(ok) > 0L) {
      ord <- ok[order(group_sig_weight[ok], group_z_evidence[ok], decreasing = TRUE)]
      candidate_top[head(ord, min(TOP_N_ATAC, length(ord)))] <- TRUE
    }

    candidate_any <- is.finite(group_sig_weight) & group_sig_weight > 0
    candidate_ge50 <- is.finite(group_sig_weight) & group_sig_weight >= 0.5
    candidate_ge80 <- is.finite(group_sig_weight) & group_sig_weight >= 0.8
    candidates <- list(
      group_weight_gt0 = candidate_any,
      group_weight_ge0p50 = candidate_ge50,
      group_weight_ge0p80 = candidate_ge80,
      group_top350_positive = candidate_top
    )

    cell_dt <- data.table(
      trait = trait,
      group_fdr_threshold = thr,
      barcode = atac_bc,
      atac_cell_type_posthoc = atac_cell_type,
      is_atac_target_posthoc = atac_target,
      group_sig_source_weight = group_sig_weight,
      group_z_evidence = group_z_evidence,
      source_target_weight_posthoc = source_target_weight,
      source_myeloid_weight_posthoc = source_myeloid_weight,
      no_celltype_used_to_select = TRUE
    )

    for (strategy in names(candidates)) {
      candidate <- candidates[[strategy]]
      cell_dt[[paste0("candidate_", strategy)]] <- candidate

      n_candidate <- sum(candidate)
      n_candidate_target <- sum(candidate & atac_target)
      precision <- safe_div(n_candidate_target, n_candidate)
      recall <- safe_div(n_candidate_target, sum(atac_target))
      enrichment <- precision / max(atac_target_rate, 1e-12)

      edge_dt <- copy(edge_dt_base)
      edge_dt[, is_candidate_atac := candidate[match(atac_barcode, atac_bc)]]
      top_sources <- top_weight_label(edge_dt[is_candidate_atac == TRUE])

      summary_rows[[paste(trait, thr, strategy, sep = "_")]] <- data.frame(
        trait = trait,
        group_fdr_threshold = thr,
        strategy = strategy,
        n_rna_groups = nrow(group_test),
        n_sig_groups = length(sig_groups),
        n_sig_group_cells = sum(cell_group_sig),
        sig_group_cell_rate = mean(cell_group_sig),
        sig_group_target_fraction_posthoc =
          mean(group_base$cell_type_posthoc[cell_group_sig] %in% CITE_TARGETS[[trait]]),
        sig_group_myeloid_fraction_posthoc =
          mean(group_base$cell_type_posthoc[cell_group_sig] %in% CITE_MYELOID),
        sig_groups_top_celltypes_posthoc =
          top_table_label(group_base$cell_type_posthoc[cell_group_sig], 10L),
        n_atac_cells = length(atac_bc),
        atac_target_rate = atac_target_rate,
        group_sig_source_weight_auc_posthoc = score_auc,
        group_z_evidence_auc_posthoc = z_auc,
        mean_group_sig_source_weight_all = mean(group_sig_weight, na.rm = TRUE),
        mean_group_sig_source_weight_atac_target = mean(group_sig_weight[atac_target], na.rm = TRUE),
        mean_group_sig_source_weight_atac_rest = mean(group_sig_weight[!atac_target], na.rm = TRUE),
        n_candidate_atac = n_candidate,
        candidate_precision_posthoc = precision,
        candidate_recall_posthoc = recall,
        candidate_enrichment_posthoc = enrichment,
        candidate_mean_source_target_weight_posthoc =
          mean(source_target_weight[candidate], na.rm = TRUE),
        candidate_mean_source_myeloid_weight_posthoc =
          mean(source_myeloid_weight[candidate], na.rm = TRUE),
        candidate_top_atac_celltypes_posthoc =
          top_table_label(atac_cell_type[candidate], 10L),
        candidate_top_rna_sources_by_weight_posthoc = top_sources,
        no_celltype_used_to_select = TRUE,
        stringsAsFactors = FALSE
      )

      by_ct <- data.table(
        trait = trait,
        group_fdr_threshold = thr,
        strategy = strategy,
        atac_cell_type = atac_cell_type,
        is_target_celltype_posthoc = atac_target,
        is_candidate = candidate,
        group_sig_source_weight = group_sig_weight,
        group_z_evidence = group_z_evidence,
        source_target_weight_posthoc = source_target_weight,
        source_myeloid_weight_posthoc = source_myeloid_weight
      )[, .(
        n_cells = .N,
        n_candidate = sum(is_candidate),
        candidate_rate = mean(is_candidate),
        mean_group_sig_source_weight = mean(group_sig_source_weight, na.rm = TRUE),
        mean_group_z_evidence = mean(group_z_evidence, na.rm = TRUE),
        mean_source_target_weight_posthoc = mean(source_target_weight_posthoc, na.rm = TRUE),
        mean_source_myeloid_weight_posthoc = mean(source_myeloid_weight_posthoc, na.rm = TRUE),
        is_target_celltype_posthoc = unique(is_target_celltype_posthoc)[1]
      ), by = .(trait, group_fdr_threshold, strategy, atac_cell_type)][
        order(trait, group_fdr_threshold, strategy, -candidate_rate, -n_candidate)
      ]
      atac_by_celltype_rows[[paste(trait, thr, strategy, sep = "_")]] <- by_ct

      src <- edge_dt[is_candidate_atac == TRUE, .(
        total_weight = sum(neighbor_weight, na.rm = TRUE),
        n_edges = .N
      ), by = cite_cell_type][order(-total_weight)]
      src[, `:=`(
        trait = trait,
        group_fdr_threshold = thr,
        strategy = strategy
      )]
      setcolorder(src, c("trait", "group_fdr_threshold", "strategy",
                         "cite_cell_type", "total_weight", "n_edges"))
      source_rows[[paste(trait, thr, strategy, sep = "_")]] <- src
    }

    cell_rows[[paste(trait, thr, sep = "_")]] <- cell_dt
  }

  cat("  groups tested: ", nrow(group_test),
      "; min group qval: ", signif(min(group_test$qval, na.rm = TRUE), 4),
      "; q<=0.10 groups: ", sum(group_test$qval <= 0.10 & group_test$observed_mean_score > 0),
      "\n", sep = "")
}

cat("[4/5] Writing diagnostics ...\n")
summary_dt <- rbindlist(summary_rows, use.names = TRUE, fill = TRUE)
group_dt <- rbindlist(group_rows, use.names = TRUE, fill = TRUE)
cell_dt <- rbindlist(cell_rows, use.names = TRUE, fill = TRUE)
atac_by_ct_dt <- rbindlist(atac_by_celltype_rows, use.names = TRUE, fill = TRUE)
source_dt <- rbindlist(source_rows, use.names = TRUE, fill = TRUE)

fwrite(summary_dt, file.path(OUT_DIR, "neighborhood_fdr_seed_summary.csv"))
fwrite(group_dt, file.path(OUT_DIR, "neighborhood_fdr_group_tests.csv"))
fwrite(cell_dt, file.path(OUT_DIR, "neighborhood_fdr_cell_level.csv.gz"))
fwrite(atac_by_ct_dt, file.path(OUT_DIR, "neighborhood_fdr_by_atac_celltype_posthoc.csv"))
fwrite(source_dt, file.path(OUT_DIR, "neighborhood_fdr_seed_source_by_cite_celltype_posthoc.csv"))

cat("[5/5] Writing figures ...\n")
plot_dt <- summary_dt[strategy %in% c("group_weight_gt0", "group_weight_ge0p50", "group_top350_positive")]

p1_dt <- plot_dt[is.finite(candidate_enrichment_posthoc)]
p1_line_dt <- p1_dt[, if (.N >= 2L) .SD, by = .(trait, strategy)]
p1 <- ggplot(p1_dt, aes(x = factor(group_fdr_threshold),
                        y = candidate_enrichment_posthoc,
                        color = strategy, group = strategy)) +
  geom_hline(yintercept = 1, linetype = 2, color = "gray40") +
  geom_point(size = 2) +
  facet_wrap(~ trait, scales = "free_y") +
  labs(title = "Post hoc enrichment of neighborhood-FDR ATAC seed candidates",
       x = "RNA neighborhood FDR threshold",
       y = "seed enrichment vs ATAC target rate",
       color = "strategy") +
  theme_bw()
if (nrow(p1_line_dt) > 0L) {
  p1 <- p1 + geom_line(data = p1_line_dt, linewidth = 0.7)
}
ggsave(file.path(FIG_DIR, "neighborhood_fdr_seed_enrichment.png"), p1,
       width = 8, height = 4.5, dpi = 150)

p2_dt <- unique(summary_dt[, .(
  trait,
  group_fdr_threshold,
  group_sig_source_weight_auc_posthoc,
  group_z_evidence_auc_posthoc
)])
p2_dt <- melt(
  p2_dt,
  id.vars = c("trait", "group_fdr_threshold"),
  variable.name = "evidence",
  value.name = "auroc"
)
p2_dt <- p2_dt[is.finite(auroc)]
p2_line_dt <- p2_dt[, if (.N >= 2L) .SD, by = .(trait, evidence)]
p2 <- ggplot(p2_dt, aes(x = factor(group_fdr_threshold), y = auroc,
                        color = evidence, group = evidence)) +
  geom_hline(yintercept = 0.5, linetype = 2, color = "gray40") +
  geom_point(size = 2) +
  facet_wrap(~ trait) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "ATAC AUROC of neighborhood-level RNA evidence",
       x = "RNA neighborhood FDR threshold", y = "AUROC", color = "evidence") +
  theme_bw()
if (nrow(p2_line_dt) > 0L) {
  p2 <- p2 + geom_line(data = p2_line_dt, linewidth = 0.7)
}
ggsave(file.path(FIG_DIR, "neighborhood_fdr_evidence_auroc.png"), p2,
       width = 8, height = 4.5, dpi = 150)

p3_dt <- plot_dt[is.finite(candidate_mean_source_myeloid_weight_posthoc)]
p3_line_dt <- p3_dt[, if (.N >= 2L) .SD, by = .(trait, strategy)]
p3 <- ggplot(p3_dt, aes(x = factor(group_fdr_threshold),
                        y = candidate_mean_source_myeloid_weight_posthoc,
                        color = strategy, group = strategy)) +
  geom_point(size = 2) +
  facet_wrap(~ trait, scales = "free_y") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Post hoc myeloid RNA-source weight in neighborhood-FDR seeds",
       x = "RNA neighborhood FDR threshold",
       y = "mean myeloid source weight",
       color = "strategy") +
  theme_bw()
if (nrow(p3_line_dt) > 0L) {
  p3 <- p3 + geom_line(data = p3_line_dt, linewidth = 0.7)
}
ggsave(file.path(FIG_DIR, "neighborhood_fdr_seed_myeloid_source_weight.png"), p3,
       width = 8, height = 4.5, dpi = 150)

cat("\n=== Neighborhood FDR diagnostics ready ===\n")
print(summary_dt[strategy == "group_top350_positive" & group_fdr_threshold == 0.10,
                 .(trait, group_fdr_threshold, strategy, n_sig_groups,
                   n_sig_group_cells, sig_group_target_fraction_posthoc,
                   sig_group_myeloid_fraction_posthoc,
                   group_sig_source_weight_auc_posthoc,
                   n_candidate_atac,
                   candidate_precision_posthoc,
                   candidate_enrichment_posthoc,
                   candidate_mean_source_myeloid_weight_posthoc)])
cat("\nOutputs are in: ", OUT_DIR, "\n", sep = "")
