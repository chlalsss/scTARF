# Prepare one case-study package from existing TRS outputs.
# Defaults to MONO with v3_top5_significant. Change CASE_TRAIT / CASE_ATAC_W if needed.

source("D:/RStudio/work/multiomics_fusion_final_delivery_plan_20260528/scripts/_final_delivery_common.R")
ensure_final_dirs()
need_packages(c("data.table", "ggplot2"))

library(data.table)
library(ggplot2)

CASE_TRAIT <- if (exists("CASE_TRAIT")) CASE_TRAIT else "MONO"
CASE_ATAC_W <- if (exists("CASE_ATAC_W")) CASE_ATAC_W else 1.0
CASE_TOP_FRAC <- if (exists("CASE_TOP_FRAC")) CASE_TOP_FRAC else 0.05

out_dir <- file.path(FINAL_ROOT, "06_case_study", CASE_TRAIT)
fig_dir <- file.path(out_dir, "figures")
dir_ensure(fig_dir)

meta_file <- file.path(MODULE4_ROOT, "04_evaluation/inputs/granja_cell_meta.rds")
baseline_file <- file.path(MODULE4_ROOT, "04_evaluation/inputs/granja_baseline_TRS_matrix.rds")

trs_candidates <- c(
  file.path(MODULE4_ROOT, sprintf("03_pseudopair_pipeline/output/v3_top5_significant/module3_mono_neighborhood_fdr_ge0p80/TRS_v2_%s_w%.1f.tsv.gz", CASE_TRAIT, CASE_ATAC_W)),
  file.path(MODULE4_ROOT, sprintf("03_pseudopair_pipeline/output/v3_sparse_other_traits/module3_neighborhood_fdr_ge0p80/TRS_v2_%s_w%.1f.tsv.gz", CASE_TRAIT, CASE_ATAC_W)),
  file.path(MODULE4_ROOT, sprintf("03_pseudopair_pipeline/output/v3_dense_top1/module3_neighborhood_fdr_ge0p80/TRS_v2_%s_w%.1f.tsv.gz", CASE_TRAIT, CASE_ATAC_W))
)
fusion_trs_file <- first_existing(trs_candidates)

if (!file.exists(meta_file)) stop("Missing metadata: ", meta_file, call. = FALSE)
if (!file.exists(baseline_file)) stop("Missing baseline TRS matrix: ", baseline_file, call. = FALSE)
if (is.na(fusion_trs_file)) {
  stop("No fusion TRS file found for CASE_TRAIT=", CASE_TRAIT, " and CASE_ATAC_W=", CASE_ATAC_W,
       "\nChecked:\n  ", paste(trs_candidates, collapse = "\n  "), call. = FALSE)
}

meta <- readRDS(meta_file)
if (inherits(meta, "data.table")) meta <- as.data.frame(meta)
if (is.null(rownames(meta)) && "barcode" %in% names(meta)) rownames(meta) <- meta$barcode
baseline <- readRDS(baseline_file)
if (inherits(baseline, "data.table")) baseline <- as.data.frame(baseline)
if ("barcode" %in% names(baseline)) {
  rownames(baseline) <- baseline$barcode
  baseline$barcode <- NULL
}
if (!CASE_TRAIT %in% names(baseline)) stop("Baseline lacks trait column: ", CASE_TRAIT, call. = FALSE)

fusion <- fread(fusion_trs_file)
bc_col <- intersect(c("barcode", "cell", "cell_id", "V1"), names(fusion))[1]
trs_col <- intersect(c("TRS", "trs", "score"), names(fusion))[1]
if (is.na(bc_col) || is.na(trs_col)) {
  bc_col <- names(fusion)[1]
  trs_col <- names(fusion)[2]
}
fusion_score <- fusion[[trs_col]]
names(fusion_score) <- fusion[[bc_col]]

base_score <- baseline[[CASE_TRAIT]]
names(base_score) <- rownames(baseline)

common <- Reduce(intersect, list(rownames(meta), names(base_score), names(fusion_score)))
case <- data.frame(
  barcode = common,
  cell_type = as.character(meta[common, "cell_type"]),
  baseline_trs = as.numeric(base_score[common]),
  fusion_trs = as.numeric(fusion_score[common]),
  stringsAsFactors = FALSE
)
case$baseline_rank <- rank(-case$baseline_trs, ties.method = "average")
case$fusion_rank <- rank(-case$fusion_trs, ties.method = "average")
case$rank_shift_to_higher_risk <- case$baseline_rank - case$fusion_rank
case$fusion_top <- case$fusion_rank <= ceiling(nrow(case) * CASE_TOP_FRAC)
case$baseline_top <- case$baseline_rank <= ceiling(nrow(case) * CASE_TOP_FRAC)

write_csv(case, file.path(out_dir, "case_cell_level_scores.csv"))

by_ct <- aggregate(
  cbind(n_cells = rep(1, nrow(case)), fusion_top = case$fusion_top, baseline_top = case$baseline_top),
  by = list(cell_type = case$cell_type),
  FUN = sum
)
by_ct$fusion_top_rate <- by_ct$fusion_top / by_ct$n_cells
by_ct$baseline_top_rate <- by_ct$baseline_top / by_ct$n_cells
by_ct <- by_ct[order(-by_ct$fusion_top_rate, -by_ct$fusion_top), ]
write_csv(by_ct, file.path(out_dir, "case_high_risk_celltype_table.csv"))

p1 <- ggplot(head(by_ct, 15), aes(x = reorder(cell_type, fusion_top_rate), y = fusion_top_rate)) +
  geom_col(fill = "#C81D25") +
  coord_flip() +
  labs(x = NULL, y = sprintf("Fraction in top %.0f%% fusion TRS", CASE_TOP_FRAC * 100),
       title = paste(CASE_TRAIT, "case study: high-risk cell-type enrichment")) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())
ggsave(file.path(fig_dir, "fig_case_high_risk_celltype_bar.png"), p1, width = 8, height = 5, dpi = 300)

p2 <- ggplot(case, aes(x = baseline_trs, y = fusion_trs, color = fusion_top)) +
  geom_point(alpha = 0.35, size = 0.7) +
  scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "#C81D25")) +
  labs(x = "Original SCAVENGE TRS", y = "Fusion TRS", color = "Fusion top risk",
       title = paste(CASE_TRAIT, "case study: baseline vs fusion TRS")) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())
ggsave(file.path(fig_dir, "fig_case_scavenge_vs_fusion_rank_shift.png"), p2, width = 6, height = 5, dpi = 300)

umap_cols <- intersect(c("UMAP_1", "UMAP_2", "umap_1", "umap_2", "UMAP1", "UMAP2"), names(meta))
if (length(umap_cols) >= 2L) {
  case$umap_1 <- as.numeric(meta[case$barcode, umap_cols[1]])
  case$umap_2 <- as.numeric(meta[case$barcode, umap_cols[2]])
  p3 <- ggplot(case, aes(x = umap_1, y = umap_2, color = fusion_trs)) +
    geom_point(size = 0.45, alpha = 0.8) +
    scale_color_viridis_c(option = "magma") +
    labs(x = "UMAP 1", y = "UMAP 2", color = "Fusion TRS",
         title = paste(CASE_TRAIT, "case study: fusion TRS on UMAP")) +
    theme_void(base_size = 11) +
    theme(legend.position = "right")
  ggsave(file.path(fig_dir, "fig_case_trs_umap.png"), p3, width = 6, height = 5, dpi = 300)
}

text <- c(
  paste0("# Case Study: ", CASE_TRAIT),
  "",
  paste0("- Fusion TRS file: `", fusion_trs_file, "`"),
  paste0("- ATAC weight: ", CASE_ATAC_W),
  paste0("- High-risk threshold: top ", CASE_TOP_FRAC * 100, "% cells by fusion TRS"),
  "",
  "## Key Tables",
  "",
  "- `case_cell_level_scores.csv`: cell-level baseline and fusion TRS.",
  "- `case_high_risk_celltype_table.csv`: high-risk-cell enrichment by cell type.",
  "",
  "## Interpretation Placeholder",
  "",
  "Describe which cell types are enriched among fusion high-risk cells, whether they match the expected trait biology, and how fusion changes the ranking relative to original SCAVENGE."
)
writeLines(text, file.path(out_dir, "case_study_text_template.md"), useBytes = TRUE)

cat("Wrote case-study assets under: ", out_dir, "\n", sep = "")
