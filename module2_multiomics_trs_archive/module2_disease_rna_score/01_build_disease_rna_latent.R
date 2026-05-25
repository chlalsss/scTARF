source("D:/RStudio/work/module2_disease_rna_score/_bootstrap.R")

required_pkgs(c("data.table", "Matrix", "irlba", "jsonlite"))
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(irlba)
  library(jsonlite)
})

p <- disease_rna_paths()
ensure_disease_dirs(p)

target_traits <- trimws(strsplit(Sys.getenv("MODULE2_TRAITS", "MONO,PLT,LYM"),
                                 ",", fixed = TRUE)[[1]])
n_pcs_raw <- Sys.getenv("MODULE2_DISEASE_PCS", "auto")
n_pcs_auto <- tolower(n_pcs_raw) == "auto"
n_pcs_max <- as.integer(Sys.getenv("MODULE2_DISEASE_PCS_MAX", "80"))
n_pcs <- if (n_pcs_auto) n_pcs_max else as.integer(n_pcs_raw)
if (!is.finite(n_pcs) || n_pcs < 1) stop("MODULE2_DISEASE_PCS must be 'auto' or a positive integer")

cat("traits:", paste(target_traits, collapse = ","),
    " | PCs:", if (n_pcs_auto) paste0("auto(max=", n_pcs, ")") else n_pcs, "\n")

gene_risk <- fread(p$trait_gene_csv)
if (!all(c("trait", "gene_id") %in% colnames(gene_risk))) {
  stop("trait_gene_csv must contain trait and gene_id columns")
}

rna_mat <- readMM(p$rna_mtx)
rna_genes <- fread(p$shared_genes)$gene
rna_cells <- fread(p$rna_cells)
if (!"cell" %in% colnames(rna_cells)) {
  stop("rna_cells.csv must contain a cell column")
}
if (length(rna_genes) != nrow(rna_mat)) {
  stop("shared_genes row count does not match RNA matrix rows")
}
if (nrow(rna_cells) != ncol(rna_mat)) {
  stop("rna_cells row count does not match RNA matrix columns")
}
rownames(rna_mat) <- rna_genes
colnames(rna_mat) <- rna_cells$cell

build_trait_latent <- function(tr) {
  G_trait <- unique(gene_risk[trait == tr, gene_id])
  if (length(G_trait) == 0) {
    stop("No risk genes found for trait: ", tr)
  }
  G_use <- intersect(G_trait, rna_genes)
  cat("[", tr, "] |G_trait| =", length(G_trait),
      "|G_use| =", length(G_use),
      "(missing", length(G_trait) - length(G_use), "genes)\n")
  if (length(G_use) < 2) {
    stop("Need at least 2 shared risk genes for PCA in ", tr, "; got ", length(G_use))
  }

  rna_sub <- rna_mat[G_use, , drop = FALSE]
  lib_size <- Matrix::colSums(rna_sub)
  lib_size[lib_size == 0] <- 1
  rna_norm <- t(t(rna_sub) / lib_size) * 10000
  rna_log <- log1p(rna_norm)

  rna_log_t <- as.matrix(t(rna_log))
  gene_means <- colMeans(rna_log_t)
  gene_sds <- apply(rna_log_t, 2, sd)
  gene_sds[gene_sds < 1e-8] <- 1
  rna_scaled <- sweep(rna_log_t, 2, gene_means)
  rna_scaled <- sweep(rna_scaled, 2, gene_sds, "/")

  pca_n <- min(n_pcs, length(G_use) - 1L, nrow(rna_scaled) - 1L)
  if (!is.finite(pca_n) || pca_n < 1) {
    stop("Invalid MODULE2_DISEASE_PCS after dimension checks for ", tr, ": ", pca_n)
  }
  pca <- irlba::prcomp_irlba(
    rna_scaled,
    n = pca_n,
    center = FALSE,
    scale. = FALSE
  )
  n_pcs_use <- if (n_pcs_auto) auto_select_pca_dim(pca$sdev) else pca_n
  cat("[", tr, "] selected PCs =", n_pcs_use, "\n")

  disease_latent <- data.frame(
    cell = rownames(rna_log_t),
    pca$x[, seq_len(n_pcs_use), drop = FALSE],
    check.names = FALSE
  )
  colnames(disease_latent)[-1] <- paste0("disease_", seq_len(n_pcs_use))

  out_rds <- file.path(p$out_dir, sprintf("rna_disease_latent_%s.rds", tr))
  saveRDS(disease_latent, out_rds, compress = FALSE)
  cat("[", tr, "] saved:", out_rds, "\n")

  list(
    trait = tr,
    G_trait_size = length(G_trait),
    G_use_size = length(G_use),
    n_pcs = n_pcs_use,
    n_pcs_mode = if (n_pcs_auto) "auto_elbow" else "manual",
    n_pcs_requested = n_pcs_raw,
    variance_first10 = head(pca$sdev^2 / sum(pca$sdev^2), 10),
    output = out_rds
  )
}

trait_results <- lapply(target_traits, build_trait_latent)
names(trait_results) <- target_traits

manifest <- list(
  script = "01_build_disease_rna_latent.R",
  mode = "per_trait",
  target_traits = target_traits,
  n_pcs_mode = if (n_pcs_auto) "auto_elbow" else "manual",
  n_pcs_requested = n_pcs_raw,
  traits = trait_results,
  created_at = as.character(Sys.time())
)
write_json(
  manifest,
  file.path(p$out_dir, "rna_disease_latent_manifest.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

summary_dt <- rbindlist(lapply(trait_results, function(x) {
  data.table(
    trait = x$trait,
    G_trait_size = x$G_trait_size,
    G_use_size = x$G_use_size,
    n_pcs = x$n_pcs,
    n_pcs_mode = x$n_pcs_mode,
    n_pcs_requested = x$n_pcs_requested,
    output = x$output
  )
}))
fwrite(summary_dt, file.path(p$out_dir, "rna_disease_latent_summary.csv"))
cat("saved manifest:", file.path(p$out_dir, "rna_disease_latent_manifest.json"), "\n")
