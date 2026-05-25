#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

ROOT <- "D:/RStudio/work/ppt_dataset_figures"
FIG_DIR <- file.path(ROOT, "figures")
TABLE_DIR <- file.path(ROOT, "tables")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TABLE_DIR, showWarnings = FALSE, recursive = TRUE)

stop_missing <- function(path) {
  if (!file.exists(path)) stop("Missing input: ", path, call. = FALSE)
  invisible(path)
}

GWAS_TSV <- "D:/SCAVENGE-reproducibility-main/blood_data/CS.trans_all.PP01.strong_atac_overlap.tsv"
BMMC_OVERLAP_CSV <- "D:/SCAVENGE-reproducibility-main/BMMC/00_precheck/BMMC_finemap_overlap_summary.csv"

traits <- c("MONO", "PLT", "LYM")
trait_labels <- c(
  MONO = "MONO\nMonocyte count",
  PLT = "PLT\nPlatelet count",
  LYM = "LYM\nLymphocyte count"
)

message("Loading GWAS fine-mapping input ...")
gwas <- fread(stop_missing(GWAS_TSV))
if ("var" %in% names(gwas) && !"variant" %in% names(gwas)) {
  setnames(gwas, "var", "variant")
}
needed <- c("variant", "trait", "PP")
if (!all(needed %in% names(gwas))) {
  stop("GWAS TSV must contain: ", paste(needed, collapse = ", "), call. = FALSE)
}

granja_counts <- gwas[trait %in% traits, .(
  n_variants = uniqueN(variant),
  pp_sum = sum(as.numeric(PP), na.rm = TRUE)
), by = trait]
granja_counts[, dataset := "Granja scATAC peak set"]

message("Loading BMMC peak-overlap summary ...")
bmmc_summary <- fread(stop_missing(BMMC_OVERLAP_CSV))
needed_bmmc <- c("trait", "n_overlap", "pp_sum_overlap")
if (!all(needed_bmmc %in% names(bmmc_summary))) {
  stop("BMMC overlap summary must contain: ", paste(needed_bmmc, collapse = ", "), call. = FALSE)
}
bmmc_counts <- bmmc_summary[trait %in% traits, .(
  trait,
  n_variants = as.integer(n_overlap),
  pp_sum = as.numeric(pp_sum_overlap)
)]
bmmc_counts[, dataset := "BMMC multiome ATAC peak set"]

plot_dt <- rbindlist(list(granja_counts, bmmc_counts), use.names = TRUE)
plot_dt[, trait_label := trait_labels[trait]]
plot_dt[, trait_label_table := gsub("\n", " - ", trait_label, fixed = TRUE)]
plot_dt[, plot_trait := factor(trait, levels = traits, labels = trait_labels[traits])]
plot_dt[, dataset := factor(dataset, levels = c("Granja scATAC peak set", "BMMC multiome ATAC peak set"))]
plot_dt[, label := trimws(format(n_variants, big.mark = ","))]
setorder(plot_dt, plot_trait, dataset)

table_out <- file.path(TABLE_DIR, "gwas_variant_overlap_counts.csv")
fwrite(
  plot_dt[, .(trait, trait_label = trait_label_table, dataset = as.character(dataset), n_variants, pp_sum)],
  table_out
)

palette <- c(
  "Granja scATAC peak set" = "#2A6F97",
  "BMMC multiome ATAC peak set" = "#E76F51"
)

p <- ggplot(plot_dt, aes(plot_trait, n_variants, fill = dataset)) +
  geom_col(position = position_dodge(width = 0.68), width = 0.58) +
  geom_text(
    aes(label = label),
    position = position_dodge(width = 0.68),
    vjust = -0.35,
    size = 3.6,
    color = "grey15"
  ) +
  scale_fill_manual(values = palette, name = NULL) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.13)),
    labels = function(x) format(x, big.mark = ",", scientific = FALSE)
  ) +
  labs(
    title = "GWAS Fine-Mapped Variants Mapped to ATAC Peak Sets",
    subtitle = "Variants are filtered at PP > 0.01 and counted after overlap with the peak set used by each experiment",
    x = NULL,
    y = "Overlapping variants"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 10.5, hjust = 0.5, color = "grey25"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 10.5, color = "grey15"),
    axis.text.y = element_text(color = "grey20"),
    legend.position = "top",
    legend.text = element_text(size = 10.5),
    plot.margin = margin(10, 14, 8, 10)
  )

png_out <- file.path(FIG_DIR, "fig_gwas_variant_overlap_counts.png")
pdf_out <- file.path(FIG_DIR, "fig_gwas_variant_overlap_counts.pdf")
ggsave(png_out, p, width = 9.2, height = 5.2, dpi = 320, bg = "white")
ggsave(pdf_out, p, width = 9.2, height = 5.2, bg = "white")

message("Wrote: ", png_out)
message("Wrote: ", pdf_out)
message("Wrote: ", table_out)
