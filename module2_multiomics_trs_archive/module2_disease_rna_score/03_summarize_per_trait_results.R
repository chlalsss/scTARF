source("D:/RStudio/work/module2_disease_rna_score/_bootstrap.R")

required_pkgs(c("data.table", "Matrix", "SummarizedExperiment", "jsonlite"))
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(SummarizedExperiment)
  library(jsonlite)
})

p <- disease_rna_paths()
ensure_disease_dirs(p)

traits <- trimws(strsplit(Sys.getenv("MODULE2_TRAITS", "MONO,PLT,LYM"),
                          ",", fixed = TRUE)[[1]])
npcs_raw <- Sys.getenv("MODULE2_FUSED_NPCS", "50")
k_raw <- Sys.getenv("MODULE2_KNN_K", "30")
resolve_npcs <- function(label) {
  if (tolower(npcs_raw) != "auto") return(as.integer(npcs_raw))
  txt <- file.path(p$adaptive_dir, sprintf("fused_embedding_auto_npcs_cutoff_%s.txt", label))
  as.integer(readLines(txt, warn = FALSE)[1])
}
resolve_k <- function(label, npcs) {
  if (tolower(k_raw) != "auto") return(as.integer(k_raw))
  txt <- file.path(p$adaptive_dir, sprintf("fused_knn_auto_k_cutoff_%s_pca%d.txt", label, npcs))
  as.integer(readLines(txt, warn = FALSE)[1])
}
se_path <- file.path(p$data_dir, "scATAC-Healthy-Hematopoiesis-191120.rds")
atac_only_trs_path <- file.path(
  p$adaptive_dir,
  "trs_atac_only_k30",
  "TRS_matrix_fusedKNN_selected_traits_fullcells.rds"
)

SE <- readRDS(se_path)
meta <- as.data.table(as.data.frame(colData(SE)))
meta[, cell := colnames(SE)]
celltype_map <- setNames(as.character(meta$BioClassification), meta$cell)

is_target_type <- function(trait, celltype) {
  ct <- as.character(celltype)
  if (trait == "MONO") {
    grepl("Mono", ct, ignore.case = TRUE)
  } else if (trait == "PLT") {
    grepl("Eryth|Baso|Mk|Mega|Platelet|PLT", ct, ignore.case = TRUE)
  } else if (trait == "LYM") {
    grepl("CD4|CD8|\\bB\\b|NK|CLP|Pre\\.B", ct, ignore.case = TRUE)
  } else {
    rep(FALSE, length(ct))
  }
}

top5_fraction <- function(trs_path, trait, method) {
  trs <- as.data.table(readRDS(trs_path))
  if (!all(c("cell", trait) %in% colnames(trs))) {
    stop("Missing columns in TRS file: ", trs_path)
  }
  dt <- merge(
    trs[, .(cell, score = get(trait))],
    meta[, .(cell, BioClassification)],
    by = "cell",
    all.x = TRUE
  )
  dt[, target_type := is_target_type(trait, BioClassification)]
  n_top <- max(1L, ceiling(nrow(dt) * 0.05))
  top_idx <- order(dt$score, decreasing = TRUE)[seq_len(n_top)]
  data.table(
    trait = trait,
    method = method,
    n_top = n_top,
    target_fraction_top = mean(dt$target_type[top_idx], na.rm = TRUE),
    target_fraction_all = mean(dt$target_type, na.rm = TRUE),
    enrichment_ratio = mean(dt$target_type[top_idx], na.rm = TRUE) /
      mean(dt$target_type, na.rm = TRUE),
    source = trs_path
  )
}

change2_reported <- data.table(
  trait = c("MONO", "PLT", "LYM"),
  method = "top1_change2_reported",
  n_top = NA_integer_,
  target_fraction_top = c(0.514269406392694, 0.802511415525114, 0.875570776255708),
  target_fraction_all = c(0.175951823734231, 0.0740339060448656, 0.424624693190251),
  source = "reported baseline from completed change 2 before weighted/no-op experiments"
)
change2_reported[, enrichment_ratio := target_fraction_top / target_fraction_all]

top_rows <- list()
for (tr in traits) {
  label <- sprintf("top1disease_%s", tr)
  npcs <- resolve_npcs(label)
  k <- resolve_k(label, npcs)
  top_rows[[paste(tr, "atac", sep = "_")]] <- top5_fraction(
    atac_only_trs_path,
    tr,
    "atac_only_scSHEFT"
  )
  per_trait_trs <- file.path(
    p$adaptive_dir,
    sprintf("trs_fused_knn_cutoff_top1disease_%s_pca%d_k%d", tr, npcs, k),
    "TRS_matrix_fusedKNN_selected_traits_fullcells.rds"
  )
  top_rows[[paste(tr, "pertrait", sep = "_")]] <- top5_fraction(
    per_trait_trs,
    tr,
    "top1disease_per_trait"
  )
}
top_dt <- rbindlist(top_rows, use.names = TRUE)
top_dt <- rbind(
  top_dt,
  change2_reported[trait %in% traits],
  fill = TRUE,
  use.names = TRUE
)
top_dt[, method_order := match(method, c(
  "atac_only_scSHEFT",
  "top1_change2_reported",
  "top1disease_per_trait"
))]
setorder(top_dt, trait, method_order)
top_dt[, method_order := NULL]
top_out <- file.path(p$out_dir, "per_trait_trs_top5_comparison.csv")
fwrite(top_dt, top_out)

same_type_stats <- function(adj_path, method, trait = NA_character_) {
  adj <- readRDS(adj_path)
  cells <- rownames(adj)
  types <- celltype_map[cells]
  types[is.na(types)] <- "Unknown"
  adj_t <- as(adj, "TsparseMatrix")
  i <- adj_t@i + 1L
  j <- adj_t@j + 1L
  keep <- i != j
  i <- i[keep]
  j <- j[keep]
  same <- types[i] == types[j]
  total <- tabulate(i, nbins = length(types))
  same_n <- tabulate(i[same], nbins = length(types))
  frac <- ifelse(total > 0, same_n / total, NA_real_)
  data.table(
    trait = trait,
    method = method,
    n_cells = length(types),
    n_cells_with_edges = sum(total > 0),
    mean_same_type_frac = mean(frac, na.rm = TRUE),
    median_same_type_frac = median(frac, na.rm = TRUE),
    n_isolated = sum(total == 0),
    source = adj_path
  )
}

same_rows <- list(
  atac_only = same_type_stats(
    file.path(p$adaptive_dir, "atac_only_mutualknn_k30.rds"),
    "atac_only_scSHEFT"
  ),
  top1 = same_type_stats(
    {
      top1_npcs <- resolve_npcs("top1")
      top1_k <- resolve_k("top1", top1_npcs)
      file.path(p$adaptive_dir, sprintf("fused_mutualknn_cutoff_top1_pca%d_k%d.rds", top1_npcs, top1_k))
    },
    "top1_full_latent"
  )
)
for (tr in traits) {
  label <- sprintf("top1disease_%s", tr)
  npcs <- resolve_npcs(label)
  k <- resolve_k(label, npcs)
  same_rows[[tr]] <- same_type_stats(
    file.path(
      p$adaptive_dir,
      sprintf("fused_mutualknn_cutoff_top1disease_%s_pca%d_k%d.rds", tr, npcs, k)
    ),
    "top1disease_per_trait",
    tr
  )
}
same_dt <- rbindlist(same_rows, use.names = TRUE, fill = TRUE)
same_out <- file.path(p$out_dir, "per_trait_knn_same_type_summary.csv")
fwrite(same_dt, same_out)

manifest <- list(
  script = "03_summarize_per_trait_results.R",
  traits = traits,
  top5_csv = top_out,
  same_type_csv = same_out,
  created_at = as.character(Sys.time())
)
write_json(
  manifest,
  file.path(p$out_dir, "per_trait_summary_manifest.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

cat("top5:", top_out, "\n")
print(top_dt[, .(trait, method, target_fraction_top, target_fraction_all, enrichment_ratio)])
cat("same-type:", same_out, "\n")
print(same_dt[, .(trait, method, mean_same_type_frac, median_same_type_frac, n_isolated)])
