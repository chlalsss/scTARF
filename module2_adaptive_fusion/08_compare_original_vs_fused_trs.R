# Compare original SCAVENGE TRS with fused-KNN TRS.
#
# Main focus: trait specificity by BioClassification, not just correlation.

source("D:/RStudio/work/module2_adaptive_fusion/_bootstrap.R")

required_pkgs(c("data.table", "dplyr", "ggplot2", "SummarizedExperiment", "jsonlite"))
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(SummarizedExperiment)
  library(jsonlite)
})

paths <- module2_paths()
ensure_module2_dirs(paths)

label <- get_active_label()
cutoff <- if (label == "top1" || startsWith(label, "top1disease")) NA_real_ else get_single_cutoff()
npcs_raw <- Sys.getenv("MODULE2_FUSED_NPCS", "auto")
npcs <- if (tolower(npcs_raw) == "auto") {
  read_auto_value(
    file.path(paths$out_dir, sprintf("fused_embedding_auto_npcs_cutoff_%s.txt", label)),
    "fused PCA dimension"
  )
} else {
  as.integer(npcs_raw)
}
k_raw <- Sys.getenv("MODULE2_KNN_K", "auto")
k <- if (tolower(k_raw) == "auto") {
  read_auto_value(
    file.path(paths$out_dir, sprintf("fused_knn_auto_k_cutoff_%s_pca%d.txt", label, npcs)),
    "KNN k"
  )
} else {
  as.integer(k_raw)
}
traits <- strsplit(Sys.getenv("MODULE2_TRAITS", "MONO,PLT,LYM"), ",", fixed = TRUE)[[1]]
traits <- trimws(traits)

data_dir <- paths$data_dir
se_path <- Sys.getenv(
  "MODULE2_ATAC_SE",
  file.path(data_dir, "scATAC-Healthy-Hematopoiesis-191120.rds")
)
original_trs_path <- Sys.getenv(
  "MODULE2_ORIGINAL_TRS",
  file.path(data_dir, "scavenge_out_small_all_traits", "TRS_matrix_small_all_traits_fullcells.rds")
)
default_run_label <- sprintf("trs_fused_knn_cutoff_%s_pca%d_k%d", label, npcs, k)
run_label <- Sys.getenv("MODULE2_RUN_LABEL", default_run_label)
fused_trs_dir <- file.path(paths$out_dir, run_label)
fused_trs_path <- file.path(fused_trs_dir, "TRS_matrix_fusedKNN_selected_traits_fullcells.rds")
out_dir <- file.path(fused_trs_dir, "comparison_original_vs_fused")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("========== Module 2: compare original vs fused TRS ==========\n")
cat("  traits       :", paste(traits, collapse = ", "), "\n")
cat("  original TRS :", original_trs_path, "\n")
cat("  fused TRS    :", fused_trs_path, "\n")
cat("  output dir   :", out_dir, "\n")

stopifnot(file.exists(se_path), file.exists(original_trs_path), file.exists(fused_trs_path))

SE_atac <- readRDS(se_path)
meta <- as.data.table(as.data.frame(colData(SE_atac)))
meta[, cell := colnames(SE_atac)]
if (!"BioClassification" %in% colnames(meta)) {
  stop("ATAC metadata must contain BioClassification")
}

original <- as.data.table(readRDS(original_trs_path))
fused <- as.data.table(readRDS(fused_trs_path))
if (!"cell" %in% colnames(original) || !"cell" %in% colnames(fused)) {
  stop("TRS matrices must contain cell column")
}

available_traits <- intersect(traits, intersect(colnames(original), colnames(fused)))
if (length(available_traits) == 0) {
  stop("No requested traits are present in both original and fused TRS matrices")
}
if (length(setdiff(traits, available_traits)) > 0) {
  warning("Missing traits skipped: ", paste(setdiff(traits, available_traits), collapse = ", "))
}
traits <- available_traits

merged <- merge(
  original[, c("cell", traits), with = FALSE],
  fused[, c("cell", traits), with = FALSE],
  by = "cell",
  suffixes = c("_original", "_fused")
)
merged <- merge(
  merged,
  meta[, .(cell, BioClassification, UMAP1, UMAP2)],
  by = "cell",
  all.x = TRUE
)

align_summary <- data.table(
  n_original = nrow(original),
  n_fused = nrow(fused),
  n_metadata = nrow(meta),
  n_merged = nrow(merged),
  n_missing_metadata = sum(is.na(merged$BioClassification))
)
print(align_summary)
fwrite(align_summary, file.path(out_dir, "trs_compare_alignment_summary.csv"))

long_list <- list()
summary_list <- list()
for (tr in traits) {
  ocol <- paste0(tr, "_original")
  fcol <- paste0(tr, "_fused")
  dt <- merged[, .(
    cell,
    BioClassification,
    UMAP1,
    UMAP2,
    trait = tr,
    original = get(ocol),
    fused = get(fcol)
  )]
  dt[, delta := fused - original]
  long_list[[tr]] <- dt

  summary_list[[tr]] <- data.table(
    trait = tr,
    n_cells = nrow(dt),
    pearson = suppressWarnings(cor(dt$original, dt$fused, method = "pearson")),
    spearman = suppressWarnings(cor(dt$original, dt$fused, method = "spearman")),
    mean_original = mean(dt$original),
    mean_fused = mean(dt$fused),
    mean_delta = mean(dt$delta),
    median_delta = median(dt$delta),
    sd_delta = sd(dt$delta)
  )
}
long_dt <- rbindlist(long_list)
summary_dt <- rbindlist(summary_list)
fwrite(summary_dt, file.path(out_dir, "TRS_original_vs_fused_summary.csv"))
fwrite(long_dt, file.path(out_dir, "TRS_original_vs_fused_long.csv"))

# Target cell-type groups for trait specificity.
is_target_type <- function(trait, celltype) {
  if (length(trait) > 1) {
    return(as.logical(mapply(is_target_type, trait, celltype, USE.NAMES = FALSE)))
  }
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

specificity_list <- list()
for (tr in traits) {
  dt <- long_dt[trait == tr]
  dt[, target_type := is_target_type(tr, BioClassification)]
  celltype_stats <- dt[, .(
    n_cells = .N,
    target_type = target_type[1],
    median_original = median(original),
    median_fused = median(fused),
    mean_original = mean(original),
    mean_fused = mean(fused),
    median_delta = median(delta),
    mean_delta = mean(delta)
  ), by = .(trait, BioClassification)]
  specificity_list[[tr]] <- celltype_stats
}
specificity_dt <- rbindlist(specificity_list)
fwrite(specificity_dt, file.path(out_dir, "TRS_celltype_specificity_summary.csv"))

top_enrichment <- list()
for (tr in traits) {
  dt <- copy(long_dt[trait == tr])
  dt[, target_type := is_target_type(tr, BioClassification)]
  for (method in c("original", "fused")) {
    score <- dt[[method]]
    for (top_frac in c(0.05, 0.10)) {
      n_top <- max(1L, ceiling(nrow(dt) * top_frac))
      top_idx <- order(score, decreasing = TRUE)[seq_len(n_top)]
      top_enrichment[[paste(tr, method, top_frac, sep = "_")]] <- data.table(
        trait = tr,
        method = method,
        top_frac = top_frac,
        n_top = n_top,
        target_fraction_top = mean(dt$target_type[top_idx], na.rm = TRUE),
        target_fraction_all = mean(dt$target_type, na.rm = TRUE),
        enrichment_ratio = mean(dt$target_type[top_idx], na.rm = TRUE) /
          mean(dt$target_type, na.rm = TRUE)
      )
    }
  }
}
top_dt <- rbindlist(top_enrichment)
fwrite(top_dt, file.path(out_dir, "TRS_celltype_enrichment_top5_top10.csv"))

# Plot 1: trait specificity boxplot by cell type, original vs fused.
plot_dt <- melt(
  long_dt,
  id.vars = c("cell", "BioClassification", "trait"),
  measure.vars = c("original", "fused"),
  variable.name = "method",
  value.name = "TRS"
)
plot_dt[, target_type := is_target_type(trait, BioClassification)]

cell_order <- plot_dt[, .(median_fused = median(TRS[method == "fused"])), by = BioClassification][
  order(median_fused)
]$BioClassification
plot_dt[, BioClassification := factor(BioClassification, levels = cell_order)]

p_specificity <- ggplot(plot_dt, aes(x = BioClassification, y = TRS, fill = method)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.2) +
  facet_wrap(~ trait, scales = "free_y", ncol = 1) +
  coord_flip() +
  scale_fill_manual(values = c(original = "#8C8C8C", fused = "#4C78A8")) +
  labs(
    title = "TRS cell-type specificity: original vs fused-KNN",
    x = "ATAC cell type",
    y = "TRS"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "top")
ggsave(
  file.path(out_dir, "Fig_TRS_celltype_specificity_MONO_PLT_LYM.png"),
  p_specificity,
  width = 10,
  height = 12,
  dpi = 300
)

# Plot 2: correlation scatter.
scatter_dt <- rbindlist(lapply(traits, function(tr) {
  dt <- long_dt[trait == tr]
  if (nrow(dt) > 6000) {
    set.seed(9527)
    dt <- dt[sample.int(.N, 6000)]
  }
  dt
}))
p_corr <- ggplot(scatter_dt, aes(x = original, y = fused)) +
  geom_point(alpha = 0.25, size = 0.35) +
  facet_wrap(~ trait, scales = "free") +
  labs(
    title = "Original vs fused-KNN TRS",
    x = "Original TRS",
    y = "Fused-KNN TRS"
  ) +
  theme_bw(base_size = 12)
ggsave(
  file.path(out_dir, "Fig_TRS_original_vs_fused_correlation.png"),
  p_corr,
  width = 9,
  height = 4.8,
  dpi = 300
)

# Plot 3: top enrichment ratio.
p_top <- ggplot(top_dt, aes(x = method, y = enrichment_ratio, fill = method)) +
  geom_col(width = 0.65) +
  facet_grid(trait ~ top_frac, labeller = label_both) +
  scale_fill_manual(values = c(original = "#8C8C8C", fused = "#4C78A8")) +
  labs(
    title = "Target cell-type enrichment in top TRS cells",
    x = NULL,
    y = "Target fraction enrichment ratio"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")
ggsave(
  file.path(out_dir, "Fig_TRS_target_enrichment_top5_top10.png"),
  p_top,
  width = 8,
  height = 6,
  dpi = 300
)

manifest <- list(
  script = "08_compare_original_vs_fused_trs.R",
  traits = traits,
  original_trs = original_trs_path,
  fused_trs = fused_trs_path,
  output_dir = out_dir,
  summary_csv = file.path(out_dir, "TRS_original_vs_fused_summary.csv"),
  specificity_csv = file.path(out_dir, "TRS_celltype_specificity_summary.csv"),
  enrichment_csv = file.path(out_dir, "TRS_celltype_enrichment_top5_top10.csv"),
  created_at = as.character(Sys.time())
)
jsonlite::write_json(
  manifest,
  file.path(out_dir, "compare_original_vs_fused_manifest.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

cat("\n========== Done ==========\n")
cat("  summary    :", file.path(out_dir, "TRS_original_vs_fused_summary.csv"), "\n")
cat("  specificity:", file.path(out_dir, "TRS_celltype_specificity_summary.csv"), "\n")
cat("  enrichment :", file.path(out_dir, "TRS_celltype_enrichment_top5_top10.csv"), "\n")
cat("  figures    :", out_dir, "\n")

