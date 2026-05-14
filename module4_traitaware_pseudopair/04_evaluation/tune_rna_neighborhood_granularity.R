# module4 diagnostic: tune unsupervised CITE RNA neighborhood granularity
# -----------------------------------------------------------------------------
# Purpose:
#   Before running neighborhood-level scDRS FDR, choose a reasonable RNA group
#   granularity. The previous run used a clustering setting that was too fine:
#   almost every cell became its own group.
#
# This script does NOT read scDRS full_score files and does NOT use cell type to
# build groups. Cell type is only summarized post hoc so we can inspect whether
# the unsupervised groups are biologically interpretable.
#
# Run in RStudio:
#   setwd("D:/RStudio/work/module4_traitaware_pseudopair")
#   source("04_evaluation/tune_rna_neighborhood_granularity.R")
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(ggplot2)
})

if (!requireNamespace("irlba", quietly = TRUE)) {
  stop("Missing irlba. In a separate R session run: install.packages('irlba')")
}
if (!requireNamespace("RANN", quietly = TRUE)) {
  stop("Missing RANN. In a separate R session run: install.packages('RANN')")
}
if (!requireNamespace("igraph", quietly = TRUE)) {
  stop("Missing igraph. In a separate R session run: install.packages('igraph')")
}

ROOT <- "D:/RStudio/work/module4_traitaware_pseudopair"
CITE_SCORE_RDS <- file.path(ROOT, "01_data_prep_external_rna/output/cite_rna_scdrs_scores.rds")
CITE_MTX <- file.path(ROOT, "02_trait_mapping/output/cite_extract/cite_trait_counts.mtx")
CITE_GENES <- file.path(ROOT, "02_trait_mapping/output/cite_extract/cite_trait_counts_genes.txt")
CITE_CELLS <- file.path(ROOT, "02_trait_mapping/output/cite_extract/cite_trait_counts_cells.txt")
OUT_DIR <- file.path(ROOT, "04_evaluation/module4_v2_diagnostics/neighborhood_granularity")
FIG_DIR <- file.path(OUT_DIR, "figures")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

if (!exists("RNA_PCA_DIMS")) RNA_PCA_DIMS <- 30L
if (!exists("KNN_K_LIST")) KNN_K_LIST <- c(30L, 50L)
if (!exists("RESOLUTION_LIST")) {
  RESOLUTION_LIST <- c(1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 0.1, 0.2, 0.5, 1.0)
}

stop_missing <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0L) {
    stop("Missing required file(s):\n  ", paste(missing, collapse = "\n  "),
         call. = FALSE)
  }
}

top_table_label <- function(x, n = 8L) {
  tab <- sort(table(x), decreasing = TRUE)
  if (length(tab) == 0L) return("")
  paste(paste0(names(head(tab, n)), ":", as.integer(head(tab, n))), collapse = "; ")
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
  n <- nrow(z)
  edges <- data.table(
    from = rep(seq_len(n), times = ncol(idx)),
    to = as.vector(idx)
  )
  edges <- edges[from != to]
  edges[, a := pmin(from, to)]
  edges[, b := pmax(from, to)]
  edges <- unique(edges[, .(from = a, to = b)])
  igraph::graph_from_data_frame(edges, directed = FALSE, vertices = data.frame(name = seq_len(n)))
}

cluster_graph <- function(g, resolution) {
  if ("cluster_leiden" %in% getNamespaceExports("igraph")) {
    cl <- igraph::cluster_leiden(g, resolution = resolution)
  } else {
    # cluster_louvain has no resolution parameter in igraph. If this branch is
    # used, all rows for different resolutions will be identical.
    cl <- igraph::cluster_louvain(g)
  }
  as.integer(igraph::membership(cl))
}

summarize_sizes <- function(group_id) {
  sizes <- as.integer(table(group_id))
  data.frame(
    n_groups = length(sizes),
    min_group_size = min(sizes),
    q01_group_size = as.numeric(quantile(sizes, 0.01)),
    q05_group_size = as.numeric(quantile(sizes, 0.05)),
    median_group_size = median(sizes),
    mean_group_size = mean(sizes),
    q95_group_size = as.numeric(quantile(sizes, 0.95)),
    q99_group_size = as.numeric(quantile(sizes, 0.99)),
    max_group_size = max(sizes),
    pct_singleton_groups = mean(sizes == 1) * 100,
    pct_groups_lt10 = mean(sizes < 10) * 100,
    pct_groups_lt30 = mean(sizes < 30) * 100,
    stringsAsFactors = FALSE
  )
}

cat(">>> Tune CITE RNA neighborhood granularity\n\n")
cat("RNA_PCA_DIMS = ", RNA_PCA_DIMS, "\n", sep = "")
cat("KNN_K_LIST = ", paste(KNN_K_LIST, collapse = ", "), "\n", sep = "")
cat("RESOLUTION_LIST = ", paste(RESOLUTION_LIST, collapse = ", "), "\n\n", sep = "")

stop_missing(c(CITE_SCORE_RDS, CITE_MTX, CITE_GENES, CITE_CELLS))

cat("[1/4] Loading CITE extracted counts and metadata ...\n")
cite_meta <- readRDS(CITE_SCORE_RDS)
cite_meta <- as.data.frame(cite_meta)
if (!all(c("cell", "cell_type", "scdrs_input_total_counts") %in% colnames(cite_meta))) {
  stop("CITE score RDS must contain cell, cell_type, and scdrs_input_total_counts.",
       call. = FALSE)
}
cite_counts <- Matrix::readMM(CITE_MTX)
cite_counts <- as(cite_counts, "CsparseMatrix")
cite_genes <- readLines(CITE_GENES, warn = FALSE)
cite_cells <- readLines(CITE_CELLS, warn = FALSE)
rownames(cite_counts) <- cite_genes
colnames(cite_counts) <- cite_cells
cite_meta <- cite_meta[match(cite_cells, cite_meta$cell), , drop = FALSE]
if (anyNA(cite_meta$cell)) stop("CITE metadata does not match extracted CITE cells.", call. = FALSE)

cat("[2/4] Building sparse PCA on CITE RNA trait-gene matrix ...\n")
rna_z <- build_rna_pca(cite_counts, cite_meta, cite_cells, dims = RNA_PCA_DIMS)
saveRDS(
  list(
    z = rna_z,
    cells = cite_cells,
    genes = cite_genes,
    source = "CITE extracted trait-gene log-normalized PCA",
    no_celltype_used = TRUE,
    created_at = as.character(Sys.time())
  ),
  file.path(OUT_DIR, "cite_rna_pca_for_granularity.rds")
)

summary_rows <- list()
membership_rows <- list()
by_group_rows <- list()

cat("[3/4] Sweeping graph k and clustering resolution ...\n")
for (k in KNN_K_LIST) {
  cat("  k = ", k, ": building graph ...\n", sep = "")
  g <- build_knn_graph(rna_z, k = k)
  for (resolution in RESOLUTION_LIST) {
    cat("    resolution = ", resolution, " ...\n", sep = "")
    group_id <- cluster_graph(g, resolution)
    sizes <- summarize_sizes(group_id)
    config_id <- paste0("k", k, "_res", format(resolution, scientific = FALSE))
    config_id <- gsub("[^A-Za-z0-9]+", "p", config_id)

    summary_rows[[paste(k, resolution, sep = "_")]] <- cbind(
      data.frame(
        config_id = config_id,
        knn_k = k,
        resolution = resolution,
        no_celltype_used_for_grouping = TRUE,
        stringsAsFactors = FALSE
      ),
      sizes
    )

    membership_rows[[paste(k, resolution, sep = "_")]] <- data.table(
      config_id = config_id,
      knn_k = k,
      resolution = resolution,
      cell = cite_cells,
      group_id = group_id
    )

    by_group <- data.table(
      config_id = config_id,
      knn_k = k,
      resolution = resolution,
      group_id = group_id,
      cell_type_posthoc = cite_meta$cell_type
    )[, .(
      n_cells = .N,
      top_celltypes_posthoc = top_table_label(cell_type_posthoc, 6L)
    ), by = .(config_id, knn_k, resolution, group_id)][
      order(config_id, -n_cells)
    ]
    by_group_rows[[paste(k, resolution, sep = "_")]] <- by_group
  }
}

cat("[4/4] Writing outputs ...\n")
summary_dt <- rbindlist(summary_rows, use.names = TRUE, fill = TRUE)
membership_dt <- rbindlist(membership_rows, use.names = TRUE, fill = TRUE)
by_group_dt <- rbindlist(by_group_rows, use.names = TRUE, fill = TRUE)

fwrite(summary_dt, file.path(OUT_DIR, "rna_group_granularity_summary.csv"))
fwrite(membership_dt, file.path(OUT_DIR, "rna_group_membership_by_resolution.csv.gz"))
fwrite(by_group_dt, file.path(OUT_DIR, "rna_group_celltype_posthoc_by_resolution.csv"))

p1 <- ggplot(summary_dt, aes(x = resolution, y = n_groups, color = factor(knn_k))) +
  geom_point(size = 2) +
  geom_line(linewidth = 0.7) +
  scale_x_log10() +
  scale_y_log10() +
  labs(title = "Number of unsupervised RNA groups",
       x = "Leiden resolution", y = "number of groups", color = "kNN k") +
  theme_bw()
ggsave(file.path(FIG_DIR, "n_groups_by_resolution.png"), p1,
       width = 6.5, height = 4, dpi = 150)

p2 <- ggplot(summary_dt, aes(x = resolution, y = median_group_size, color = factor(knn_k))) +
  geom_point(size = 2) +
  geom_line(linewidth = 0.7) +
  scale_x_log10() +
  scale_y_log10() +
  labs(title = "Median RNA group size",
       x = "Leiden resolution", y = "median group size", color = "kNN k") +
  theme_bw()
ggsave(file.path(FIG_DIR, "median_group_size_by_resolution.png"), p2,
       width = 6.5, height = 4, dpi = 150)

p3 <- ggplot(summary_dt, aes(x = resolution, y = pct_groups_lt30, color = factor(knn_k))) +
  geom_point(size = 2) +
  geom_line(linewidth = 0.7) +
  scale_x_log10() +
  coord_cartesian(ylim = c(0, 100)) +
  labs(title = "Small unsupervised RNA groups",
       x = "Leiden resolution", y = "% groups with <30 cells", color = "kNN k") +
  theme_bw()
ggsave(file.path(FIG_DIR, "small_groups_by_resolution.png"), p3,
       width = 6.5, height = 4, dpi = 150)

cat("\n=== RNA neighborhood granularity sweep finished ===\n")
print(summary_dt[order(knn_k, resolution)])
cat("\nSuggested target range: n_groups about 50-500, median_group_size about 100-2000,\n")
cat("and pct_groups_lt30 preferably low. Use this summary to choose k/resolution\n")
cat("before rerunning neighborhood-level scDRS FDR.\n")
cat("\nOutputs are in: ", OUT_DIR, "\n", sep = "")
