#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(Matrix)
  library(patchwork)
  library(BuenColors)
})

set.seed(42)

ROOT <- "D:/RStudio/work/ppt_dataset_figures"
BMMC_DIR <- "D:/SCAVENGE-reproducibility-main/BMMC"
FIG_DIR <- file.path(ROOT, "figures")
TABLE_DIR <- file.path(ROOT, "tables")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TABLE_DIR, showWarnings = FALSE, recursive = TRUE)

stop_missing <- function(path) {
  if (!file.exists(path)) stop("Missing input: ", path, call. = FALSE)
  invisible(path)
}

message("Loading BMMC UMAP coordinates ...")
umap_file <- file.path(BMMC_DIR, "04_dual_modality", "figures", "bmmc_atac_umap_coords.rds")
lsi_file <- file.path(BMMC_DIR, "03_scavenge_baseline", "atac_lsi.rds")

if (file.exists(umap_file)) {
  umap_dt <- as.data.table(readRDS(umap_file))
} else {
  if (!requireNamespace("uwot", quietly = TRUE)) {
    stop("Package 'uwot' is required to compute UMAP when cached coordinates are absent.")
  }
  stop_missing(lsi_file)
  message("Cached UMAP not found. Computing UMAP from ATAC LSI ...")
  lsi <- readRDS(lsi_file)
  if (is.data.frame(lsi)) {
    if ("barcode" %in% names(lsi)) {
      rownames(lsi) <- lsi$barcode
      lsi$barcode <- NULL
    }
    lsi <- as.matrix(lsi[, vapply(lsi, is.numeric, logical(1)), drop = FALSE])
  }
  umap_coords <- uwot::umap(
    lsi,
    n_neighbors = 30,
    min_dist = 0.3,
    n_components = 2,
    n_threads = 4
  )
  umap_dt <- data.table(
    barcode = rownames(lsi),
    UMAP1 = umap_coords[, 1],
    UMAP2 = umap_coords[, 2]
  )
  saveRDS(umap_dt, umap_file)
}

needed_umap <- c("barcode", "UMAP1", "UMAP2")
if (!all(needed_umap %in% names(umap_dt))) {
  stop("UMAP object must contain: ", paste(needed_umap, collapse = ", "), call. = FALSE)
}

message("Loading BMMC metadata ...")
meta_file <- file.path(BMMC_DIR, "01_data_qc", "BMMC_cell_meta.rds")
meta <- as.data.table(readRDS(stop_missing(meta_file)))
needed_meta <- c("barcode", "cell_type", "donor")
if (!all(needed_meta %in% names(meta))) {
  stop("Metadata must contain: ", paste(needed_meta, collapse = ", "), call. = FALSE)
}

message("Loading or computing paired modality QC metrics ...")
qc_file <- file.path(TABLE_DIR, "bmmc_paired_multiome_qc.csv.gz")
if (file.exists(qc_file)) {
  qc <- fread(qc_file)
} else {
  atac_file <- file.path(BMMC_DIR, "01_data_qc", "BMMC_atac_counts.rds")
  rna_file <- file.path(BMMC_DIR, "01_data_qc", "BMMC_rna_counts.rds")
  atac_counts <- readRDS(stop_missing(atac_file))
  rna_counts <- readRDS(stop_missing(rna_file))

  atac_barcodes <- rownames(atac_counts)
  rna_barcodes <- rownames(rna_counts)
  if (is.null(atac_barcodes) || is.null(rna_barcodes)) {
    stop("ATAC/RNA count matrices must have rownames matching cell barcodes.", call. = FALSE)
  }

  qc <- data.table(
    barcode = atac_barcodes,
    atac_peaks_detected = as.numeric(Matrix::rowSums(atac_counts != 0)),
    rna_umi_count = as.numeric(Matrix::rowSums(rna_counts[atac_barcodes, , drop = FALSE])),
    rna_genes_detected = as.numeric(Matrix::rowSums(rna_counts[atac_barcodes, , drop = FALSE] != 0))
  )
  fwrite(qc, qc_file)
}

meta_cols <- intersect(c("barcode", "cell_type", "donor", "site"), names(meta))
plot_dt <- merge(umap_dt, meta[, ..meta_cols], by = "barcode")
plot_dt <- merge(plot_dt, qc, by = "barcode")
plot_dt[, donor := factor(donor)]
plot_dt[, log10_rna_umi := log10(rna_umi_count + 1)]

cell_order <- plot_dt[, .N, by = cell_type][order(-N)]$cell_type
plot_dt[, cell_type := factor(cell_type, levels = cell_order)]

n_cells <- nrow(plot_dt)
n_donors <- uniqueN(plot_dt$donor)
n_sites <- if ("site" %in% names(plot_dt)) uniqueN(plot_dt$site) else NA_integer_
n_cell_types <- uniqueN(plot_dt$cell_type)
summary_dt <- data.table(
  metric = c("cells", "donors", "sites", "cell_types", "median_atac_peaks", "median_rna_umi", "median_rna_genes"),
  value = c(
    n_cells,
    n_donors,
    n_sites,
    n_cell_types,
    median(plot_dt$atac_peaks_detected),
    median(plot_dt$rna_umi_count),
    median(plot_dt$rna_genes_detected)
  )
)
fwrite(summary_dt, file.path(TABLE_DIR, "bmmc_paired_multiome_summary.csv"))

plot_dt <- plot_dt[sample(.N)]

base_umap_theme <- theme_void(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5, color = "grey15"),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "grey30", margin = margin(b = 8)),
    plot.caption = element_text(face = "bold", size = 15, hjust = 0.5, color = "grey15", margin = margin(t = 10)),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    legend.key.height = unit(0.42, "cm"),
    legend.key.width = unit(0.42, "cm"),
    plot.margin = margin(8, 10, 12, 10)
  )

cell_palette <- setNames(hcl.colors(n_cell_types, "Dark 3"), levels(plot_dt$cell_type))

p_celltype <- ggplot(plot_dt, aes(UMAP1, UMAP2, color = cell_type)) +
  geom_point(size = 1, na.rm = TRUE) +
  scale_color_manual(values = cell_palette) +
  guides(color = guide_legend(override.aes = list(size = 3), ncol = 1)) +
  pretty_plot() +
  labs(
    title = "Cell type",
    subtitle = sprintf("%s cells, %s annotated cell types", format(n_cells, big.mark = ","), n_cell_types),
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.key.height = unit(0.34, "cm"),
    legend.key.width = unit(0.34, "cm"),
    plot.title = element_text(face = "bold", size = 17, hjust = 0.5, color = "grey10"),
    plot.subtitle = element_text(size = 13, hjust = 0.5, color = "grey30", margin = margin(b = 8))
  )

p_rna <- ggplot(plot_dt, aes(UMAP1, UMAP2, color = log10_rna_umi)) +
  geom_point(size = 0.32, alpha = 0.76, stroke = 0) +
  scale_color_viridis_c(
    option = "C",
    name = "log10\nRNA UMI"
  ) +
  labs(
    title = "RNA depth",
    subtitle = "RNA UMI counts projected on the same ATAC-derived UMAP"
  ) +
  base_umap_theme

celltype_png <- file.path(FIG_DIR, "fig_bmmc_celltype_umap.png")
rna_png <- file.path(FIG_DIR, "fig_bmmc_rna_depth_umap.png")

ggsave(celltype_png, p_celltype, width = 8.4, height = 6.0, dpi = 320, bg = "white")
ggsave(rna_png, p_rna, width = 8.4, height = 6.0, dpi = 320, bg = "white")

message("Wrote: ", celltype_png)
message("Wrote: ", rna_png)
message("Wrote: ", qc_file)
