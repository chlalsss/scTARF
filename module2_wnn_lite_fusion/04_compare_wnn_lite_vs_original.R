# Compare original SCAVENGE TRS with WNN-lite TRS.
#
# Main focus:
#   1. TRS correlation and global shifts.
#   2. Top TRS cell-type composition.
#   3. Target-region enrichment for MONO / LYM and PLT.
#
# Note for PLT:
#   This blood ATAC annotation has no explicit Mk/Mega/Platelet class.
#   Previous project notes describe PLT signal as concentrated in the
#   red/megakaryocytic region, so this script uses Eryth/Baso labels as the
#   available annotated proxy and also reports raw top cell-type composition.

source("D:/RStudio/work/module2_wnn_lite_fusion/_bootstrap.R")

required_pkgs(c("data.table", "dplyr", "ggplot2", "SummarizedExperiment", "jsonlite"))
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(SummarizedExperiment)
  library(jsonlite)
})

paths <- wnn_lite_paths()
ensure_wnn_lite_dirs(paths)

cutoff <- get_single_cutoff()
label <- cutoff_label(cutoff)
k <- env_integer("WNN_LITE_KNN_K", env_integer("MODULE2_KNN_K", 30))
traits <- strsplit(Sys.getenv("MODULE2_TRAITS", "MONO,PLT,LYM"), ",", fixed = TRUE)[[1]]
traits <- trimws(traits)
traits <- traits[nzchar(traits)]

data_dir <- paths$data_dir
se_path <- Sys.getenv(
  "MODULE2_ATAC_SE",
  file.path(data_dir, "scATAC-Healthy-Hematopoiesis-191120.rds")
)
original_trs_path <- Sys.getenv(
  "MODULE2_ORIGINAL_TRS",
  file.path(data_dir, "scavenge_out_small_all_traits", "TRS_matrix_small_all_traits_fullcells.rds")
)
run_label <- Sys.getenv(
  "WNN_LITE_RUN_LABEL",
  Sys.getenv("MODULE2_RUN_LABEL", sprintf("trs_wnn_lite_cutoff_%s_k%d", label, k))
)
wnn_trs_dir <- file.path(paths$out_dir, run_label)
wnn_trs_path <- file.path(wnn_trs_dir, "TRS_matrix_wnnLite_selected_traits_fullcells.rds")
out_dir <- file.path(wnn_trs_dir, "comparison_original_vs_wnn_lite")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("========== Module 2 WNN-lite: compare original vs WNN-lite TRS ==========\n")
cat("  traits       :", paste(traits, collapse = ", "), "\n")
cat("  original TRS :", original_trs_path, "\n")
cat("  WNN-lite TRS :", wnn_trs_path, "\n")
cat("  ATAC SE      :", se_path, "\n")
cat("  output dir   :", out_dir, "\n")

stopifnot(file.exists(se_path), file.exists(original_trs_path), file.exists(wnn_trs_path))

SE_atac <- readRDS(se_path)
meta <- as.data.table(as.data.frame(colData(SE_atac)))
meta[, cell := colnames(SE_atac)]
if (!"BioClassification" %in% colnames(meta)) {
  stop("ATAC metadata must contain BioClassification", call. = FALSE)
}

original <- as.data.table(readRDS(original_trs_path))
wnn <- as.data.table(readRDS(wnn_trs_path))
if (!"cell" %in% colnames(original) || !"cell" %in% colnames(wnn)) {
  stop("TRS matrices must contain cell column", call. = FALSE)
}

available_traits <- intersect(traits, intersect(colnames(original), colnames(wnn)))
if (length(available_traits) == 0) {
  stop("No requested traits are present in both original and WNN-lite TRS matrices", call. = FALSE)
}
if (length(setdiff(traits, available_traits)) > 0) {
  warning("Missing traits skipped: ", paste(setdiff(traits, available_traits), collapse = ", "))
}
traits <- available_traits

merged <- merge(
  original[, c("cell", traits), with = FALSE],
  wnn[, c("cell", traits), with = FALSE],
  by = "cell",
  suffixes = c("_original", "_wnn_lite")
)
merged <- merge(
  merged,
  meta[, .(cell, BioClassification, UMAP1, UMAP2)],
  by = "cell",
  all.x = TRUE
)

align_summary <- data.table(
  n_original = nrow(original),
  n_wnn_lite = nrow(wnn),
  n_metadata = nrow(meta),
  n_merged = nrow(merged),
  n_missing_metadata = sum(is.na(merged$BioClassification))
)
print(align_summary)
fwrite(align_summary, file.path(out_dir, "trs_compare_alignment_summary.csv"))

is_target_type <- function(trait, celltype) {
  if (length(trait) > 1) {
    return(as.logical(mapply(is_target_type, trait, celltype, USE.NAMES = FALSE)))
  }
  ct <- as.character(celltype)
  if (trait == "MONO") {
    grepl("Mono", ct, ignore.case = TRUE)
  } else if (trait == "PLT") {
    grepl("Eryth|Baso", ct, ignore.case = TRUE)
  } else if (trait == "LYM") {
    grepl("CD4|CD8|\\bB\\b|NK|CLP|Pre\\.B", ct, ignore.case = TRUE)
  } else {
    rep(FALSE, length(ct))
  }
}

long_list <- list()
summary_list <- list()
for (tr in traits) {
  ocol <- paste0(tr, "_original")
  wcol <- paste0(tr, "_wnn_lite")
  dt <- merged[, .(
    cell,
    BioClassification,
    UMAP1,
    UMAP2,
    trait = tr,
    original = get(ocol),
    wnn_lite = get(wcol)
  )]
  dt[, delta := wnn_lite - original]
  long_list[[tr]] <- dt

  summary_list[[tr]] <- data.table(
    trait = tr,
    n_cells = nrow(dt),
    pearson = suppressWarnings(cor(dt$original, dt$wnn_lite, method = "pearson")),
    spearman = suppressWarnings(cor(dt$original, dt$wnn_lite, method = "spearman")),
    mean_original = mean(dt$original),
    mean_wnn_lite = mean(dt$wnn_lite),
    mean_delta = mean(dt$delta),
    median_delta = median(dt$delta),
    sd_delta = sd(dt$delta)
  )
}
long_dt <- rbindlist(long_list)
summary_dt <- rbindlist(summary_list)
fwrite(summary_dt, file.path(out_dir, "TRS_original_vs_wnn_lite_summary.csv"))
fwrite(long_dt, file.path(out_dir, "TRS_original_vs_wnn_lite_long.csv"))

specificity_list <- list()
for (tr in traits) {
  dt <- long_dt[trait == tr]
  dt[, target_type := is_target_type(tr, BioClassification)]
  specificity_list[[tr]] <- dt[, .(
    n_cells = .N,
    target_type = target_type[1],
    median_original = median(original),
    median_wnn_lite = median(wnn_lite),
    mean_original = mean(original),
    mean_wnn_lite = mean(wnn_lite),
    median_delta = median(delta),
    mean_delta = mean(delta)
  ), by = .(trait, BioClassification)]
}
specificity_dt <- rbindlist(specificity_list)
fwrite(specificity_dt, file.path(out_dir, "TRS_celltype_specificity_summary.csv"))

top_enrichment <- list()
top_composition <- list()
for (tr in traits) {
  dt <- copy(long_dt[trait == tr])
  dt[, target_type := is_target_type(tr, BioClassification)]
  for (method in c("original", "wnn_lite")) {
    score <- dt[[method]]
    for (top_frac in c(0.05, 0.10)) {
      n_top <- max(1L, ceiling(nrow(dt) * top_frac))
      top_idx <- order(score, decreasing = TRUE)[seq_len(n_top)]
      top_dt <- dt[top_idx]
      all_target <- mean(dt$target_type, na.rm = TRUE)
      top_target <- mean(top_dt$target_type, na.rm = TRUE)

      top_enrichment[[paste(tr, method, top_frac, sep = "_")]] <- data.table(
        trait = tr,
        method = method,
        top_frac = top_frac,
        n_top = n_top,
        target_definition = if (tr == "PLT") "Eryth/Baso proxy" else "annotated target classes",
        target_fraction_top = top_target,
        target_fraction_all = all_target,
        enrichment_ratio = if (all_target > 0) top_target / all_target else NA_real_
      )

      comp <- top_dt[, .(
        n = .N,
        fraction = .N / n_top,
        median_score = median(get(method))
      ), by = BioClassification][order(-n)]
      comp[, `:=`(trait = tr, method = method, top_frac = top_frac, rank = seq_len(.N))]
      top_composition[[paste(tr, method, top_frac, sep = "_")]] <- comp
    }
  }
}
top_dt <- rbindlist(top_enrichment)
composition_dt <- rbindlist(top_composition)
setcolorder(composition_dt, c("trait", "method", "top_frac", "rank", "BioClassification"))
fwrite(top_dt, file.path(out_dir, "TRS_target_enrichment_top5_top10.csv"))
fwrite(composition_dt, file.path(out_dir, "TRS_top_celltype_composition_top5_top10.csv"))

plot_dt <- melt(
  long_dt,
  id.vars = c("cell", "BioClassification", "trait"),
  measure.vars = c("original", "wnn_lite"),
  variable.name = "method",
  value.name = "TRS"
)

cell_order <- plot_dt[, .(median_wnn_lite = median(TRS[method == "wnn_lite"])), by = BioClassification][
  order(median_wnn_lite)
]$BioClassification
plot_dt[, BioClassification := factor(BioClassification, levels = cell_order)]

p_specificity <- ggplot(plot_dt, aes(x = BioClassification, y = TRS, fill = method)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.2) +
  facet_wrap(~ trait, scales = "free_y", ncol = 1) +
  coord_flip() +
  scale_fill_manual(values = c(original = "#8C8C8C", wnn_lite = "#4C78A8")) +
  labs(
    title = "TRS cell-type specificity: original vs WNN-lite",
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

scatter_dt <- rbindlist(lapply(traits, function(tr) {
  dt <- long_dt[trait == tr]
  if (nrow(dt) > 6000) {
    set.seed(9527)
    dt <- dt[sample.int(.N, 6000)]
  }
  dt
}))
p_corr <- ggplot(scatter_dt, aes(x = original, y = wnn_lite)) +
  geom_point(alpha = 0.25, size = 0.35) +
  facet_wrap(~ trait, scales = "free") +
  labs(
    title = "Original vs WNN-lite TRS",
    x = "Original TRS",
    y = "WNN-lite TRS"
  ) +
  theme_bw(base_size = 12)
ggsave(
  file.path(out_dir, "Fig_TRS_original_vs_wnn_lite_correlation.png"),
  p_corr,
  width = 9,
  height = 4.8,
  dpi = 300
)

p_top <- ggplot(top_dt, aes(x = method, y = enrichment_ratio, fill = method)) +
  geom_col(width = 0.65) +
  facet_grid(trait ~ top_frac, labeller = label_both) +
  scale_fill_manual(values = c(original = "#8C8C8C", wnn_lite = "#4C78A8")) +
  labs(
    title = "Target-region enrichment in top TRS cells",
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
  script = "04_compare_wnn_lite_vs_original.R",
  traits = traits,
  cutoff = cutoff,
  k = k,
  original_trs = original_trs_path,
  wnn_lite_trs = wnn_trs_path,
  output_dir = out_dir,
  summary_csv = file.path(out_dir, "TRS_original_vs_wnn_lite_summary.csv"),
  specificity_csv = file.path(out_dir, "TRS_celltype_specificity_summary.csv"),
  enrichment_csv = file.path(out_dir, "TRS_target_enrichment_top5_top10.csv"),
  top_composition_csv = file.path(out_dir, "TRS_top_celltype_composition_top5_top10.csv"),
  created_at = as.character(Sys.time())
)
jsonlite::write_json(
  manifest,
  file.path(out_dir, "compare_original_vs_wnn_lite_manifest.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

cat("\n========== Done ==========\n")
cat("  summary        :", file.path(out_dir, "TRS_original_vs_wnn_lite_summary.csv"), "\n")
cat("  specificity    :", file.path(out_dir, "TRS_celltype_specificity_summary.csv"), "\n")
cat("  enrichment     :", file.path(out_dir, "TRS_target_enrichment_top5_top10.csv"), "\n")
cat("  top composition:", file.path(out_dir, "TRS_top_celltype_composition_top5_top10.csv"), "\n")
cat("  figures        :", out_dir, "\n")

