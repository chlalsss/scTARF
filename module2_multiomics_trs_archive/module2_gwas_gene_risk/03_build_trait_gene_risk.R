# Aggregate weighted GWAS variant -> peak -> gene links into trait gene-risk scores.

source("D:/RStudio/work/module2_gwas_gene_risk/_bootstrap.R")

required_pkgs(c("data.table", "jsonlite"))
suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

paths <- gwas_gene_paths()
ensure_gwas_gene_dirs(paths)

map_dir <- file.path(paths$out_dir, "gwas_peak_gene_map")
links_path <- Sys.getenv(
  "MODULE2_GWAS_PEAK_GENE_LINKS",
  file.path(map_dir, "gwas_variant_peak_gene_links.csv.gz")
)
out_dir <- file.path(paths$out_dir, "trait_gene_risk")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

timestamped_path <- function(path) {
  ext <- tools::file_ext(path)
  stem <- tools::file_path_sans_ext(basename(path))
  suffix <- format(Sys.time(), "%Y%m%d_%H%M%S")
  file.path(
    dirname(path),
    paste0(stem, "_", suffix, if (nzchar(ext)) paste0(".", ext) else "")
  )
}

safe_fwrite <- function(x, path) {
  tryCatch(
    {
      data.table::fwrite(x, path)
      path
    },
    error = function(e) {
      alt <- timestamped_path(path)
      warning(
        "Could not write ", path, ": ", conditionMessage(e),
        ". Writing ", alt, " instead.",
        call. = FALSE
      )
      data.table::fwrite(x, alt)
      alt
    }
  )
}

cat("========== Module 2 GWAS gene-risk: aggregate trait gene risk ==========\n")
cat("  links     :", links_path, "\n")
cat("  output dir:", out_dir, "\n")

if (!file.exists(links_path)) {
  stop("Missing GWAS peak-gene links. Run 02_build_gwas_peak_gene_map.R first: ", links_path, call. = FALSE)
}

links <- data.table::fread(links_path)
required_cols <- c("trait", "var", "PP", "weighted_PP", "peak_id", "gene_id", "link_type", "link_weight")
if (!all(required_cols %in% colnames(links))) {
  stop("Link table must contain: ", paste(required_cols, collapse = ", "), call. = FALSE)
}

links[, `:=`(
  PP = as.numeric(PP),
  weighted_PP = as.numeric(weighted_PP),
  link_weight = as.numeric(link_weight)
)]

gene_risk <- links[, .(
  gene_risk_score = sum(weighted_PP, na.rm = TRUE),
  max_weighted_PP = max(weighted_PP, na.rm = TRUE),
  max_PP = max(PP, na.rm = TRUE),
  n_links = .N,
  n_variants = uniqueN(var),
  n_peaks = uniqueN(peak_id),
  mean_link_weight = mean(link_weight, na.rm = TRUE),
  min_distance_to_tss = min(distance_to_tss, na.rm = TRUE),
  has_body_promoter_link = any(link_type == "body_promoter")
), by = .(trait, gene_id)]

gene_risk[!is.finite(min_distance_to_tss), min_distance_to_tss := NA_real_]
gene_risk[, gene_risk_rank := frank(gene_risk_score, ties.method = "average"), by = trait]
gene_risk[, gene_risk_rank01 := (gene_risk_rank - 0.5) / .N, by = trait]
gene_risk[, gene_risk_score01 := {
  max_score <- max(gene_risk_score, na.rm = TRUE)
  if (is.finite(max_score) && max_score > 0) gene_risk_score / max_score else 0
}, by = trait]
setorder(gene_risk, trait, -gene_risk_score, gene_id)

gene_risk_out <- file.path(out_dir, "trait_gene_risk_long.csv")
gene_risk_out <- safe_fwrite(gene_risk, gene_risk_out)

trait_summary <- gene_risk[, .(
  n_genes = .N,
  n_genes_body_promoter = sum(has_body_promoter_link),
  n_genes_nearest_only = sum(!has_body_promoter_link),
  total_gene_risk_score = sum(gene_risk_score, na.rm = TRUE),
  median_gene_risk_score = as.numeric(median(gene_risk_score, na.rm = TRUE)),
  max_gene_risk_score = max(gene_risk_score, na.rm = TRUE),
  median_n_variants_per_gene = as.numeric(median(n_variants)),
  median_n_peaks_per_gene = as.numeric(median(n_peaks))
), by = trait][order(trait)]
trait_summary_out <- file.path(out_dir, "trait_gene_risk_summary.csv")
trait_summary_out <- safe_fwrite(trait_summary, trait_summary_out)

top_n <- env_integer("MODULE2_TOP_RISK_GENES_N", 50)
top_genes <- gene_risk[, head(.SD, top_n), by = trait]
top_genes_out <- file.path(out_dir, sprintf("trait_gene_risk_top%d.csv", top_n))
top_genes_out <- safe_fwrite(top_genes, top_genes_out)

link_type_gene_summary <- links[, .(
  n_links = .N,
  n_variants = uniqueN(var),
  n_peaks = uniqueN(peak_id),
  n_genes = uniqueN(gene_id),
  sum_weighted_PP = sum(weighted_PP, na.rm = TRUE)
), by = .(trait, link_type)][order(trait, link_type)]
link_type_gene_summary_out <- file.path(out_dir, "trait_gene_risk_link_type_summary.csv")
link_type_gene_summary_out <- safe_fwrite(link_type_gene_summary, link_type_gene_summary_out)

manifest <- list(
  script = "03_build_trait_gene_risk.R",
  links = links_path,
  gene_risk_long = gene_risk_out,
  trait_summary = trait_summary_out,
  top_genes = top_genes_out,
  link_type_summary = link_type_gene_summary_out,
  top_n = top_n,
  created_at = as.character(Sys.time())
)
jsonlite::write_json(
  manifest,
  file.path(out_dir, "trait_gene_risk_manifest.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

cat("\n========== Done ==========\n")
cat("  gene risk :", gene_risk_out, "\n")
cat("  summary   :", trait_summary_out, "\n")
cat("  top genes :", top_genes_out, "\n")
