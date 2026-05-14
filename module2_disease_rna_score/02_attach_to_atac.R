source("D:/RStudio/work/module2_disease_rna_score/_bootstrap.R")

required_pkgs(c("data.table", "jsonlite"))
suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

p <- disease_rna_paths()
ensure_disease_dirs(p)

target_traits <- trimws(strsplit(Sys.getenv("MODULE2_TRAITS", "MONO,PLT,LYM"),
                                 ",", fixed = TRUE)[[1]])

top1 <- fread(p$top1_neighbors)
if (!all(c("atac_cell", "rank", "rna_cell") %in% colnames(top1))) {
  stop("top1 neighbor table must contain atac_cell, rank, and rna_cell columns")
}
stopifnot(max(top1$rank) == 1L)

attach_trait <- function(tr) {
  latent_rds <- file.path(p$out_dir, sprintf("rna_disease_latent_%s.rds", tr))
  if (!file.exists(latent_rds)) {
    stop("Missing disease latent for ", tr, ": ", latent_rds)
  }
  disease_latent <- readRDS(latent_rds)
  lookup <- as.matrix(disease_latent[, -1])
  rownames(lookup) <- disease_latent$cell

  missing_rna <- setdiff(top1$rna_cell, rownames(lookup))
  if (length(missing_rna) > 0) {
    stop("Missing disease latent rows for top1 RNA cells in ", tr, ": ", length(missing_rna))
  }

  ref_mat <- lookup[top1$rna_cell, , drop = FALSE]
  rna_ref_df <- data.frame(
    atac_cell = top1$atac_cell,
    ref_mat,
    check.names = FALSE
  )
  colnames(rna_ref_df)[-1] <- paste0("rna_ref_", seq_len(ncol(ref_mat)))

  out_label <- sprintf("top1disease_%s", tr)
  out_rds <- file.path(
    p$adaptive_dir,
    sprintf("rna_reference_latent_cutoff_%s.rds", out_label)
  )
  saveRDS(rna_ref_df, out_rds, compress = FALSE)
  cat("[", tr, "] saved:", out_rds, "\n")

  list(
    trait = tr,
    output_label = out_label,
    latent_rds = latent_rds,
    output_rds = out_rds,
    n_atac = nrow(rna_ref_df),
    ref_dim = ncol(ref_mat)
  )
}

trait_results <- lapply(target_traits, attach_trait)
names(trait_results) <- target_traits

write_json(
  list(
    script = "02_attach_to_atac.R",
    mode = "per_trait",
    traits = trait_results,
    created_at = as.character(Sys.time())
  ),
  file.path(p$out_dir, "attach_manifest.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

summary_dt <- rbindlist(lapply(trait_results, function(x) {
  data.table(
    trait = x$trait,
    output_label = x$output_label,
    ref_dim = x$ref_dim,
    n_atac = x$n_atac,
    output_rds = x$output_rds
  )
}))
fwrite(summary_dt, file.path(p$out_dir, "attach_summary.csv"))
cat("saved manifest:", file.path(p$out_dir, "attach_manifest.json"), "\n")
