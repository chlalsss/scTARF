# Plot diagnostics for adaptive RNA neighbor selection.

source("D:/RStudio/work/module2_adaptive_fusion/_bootstrap.R")

required_pkgs(c("data.table", "ggplot2"))
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

paths <- module2_paths()
ensure_module2_dirs(paths)
k <- as.integer(Sys.getenv("MODULE2_TOPK", "100"))
cutoffs <- get_cutoffs()
main_cutoff <- get_single_cutoff()
main_label <- cutoff_label(main_cutoff)

topk_csv <- file.path(paths$out_dir, sprintf("atac_top%d_rna_neighbors.csv", k))
summary_csv <- file.path(paths$out_dir, "adaptive_selection_summary.csv")
if (!file.exists(topk_csv)) stop("Missing top-k neighbor file: ", topk_csv)
if (!file.exists(summary_csv)) stop("Missing adaptive selection summary: ", summary_csv)

cat("========== Module 2: diagnostics plots ==========\n")

dt <- fread(topk_csv)
setorder(dt, atac_cell, rank)
dt[, d_min := min(distance), by = atac_cell]
dt[, delta := distance - d_min]
scale_dt <- dt[, .(scale = safe_scale_delta(delta)), by = atac_cell]
dt <- merge(dt, scale_dt, by = "atac_cell", all.x = TRUE, sort = FALSE)
dt[, similarity := exp(-delta / scale)]

plot_sample_n <- as.integer(Sys.getenv("MODULE2_PLOT_SAMPLE_N", "500000"))
plot_dt <- if (nrow(dt) > plot_sample_n) {
  set.seed(9527)
  dt[sample.int(.N, plot_sample_n)]
} else {
  dt
}

p1 <- ggplot(plot_dt, aes(x = similarity)) +
  geom_histogram(bins = 80, fill = "#4C78A8", color = "white", linewidth = 0.15) +
  geom_vline(xintercept = main_cutoff, color = "#D62728", linewidth = 0.8) +
  labs(
    title = sprintf("RNA neighbor similarity distribution (cutoff %.2f)", main_cutoff),
    x = "Relative similarity",
    y = "Neighbor count"
  ) +
  theme_bw(base_size = 12)
ggsave(
  file.path(paths$fig_dir, sprintf("Fig_similarity_distribution_cutoff_%s.png", main_label)),
  p1,
  width = 7,
  height = 4.8,
  dpi = 300
)

rank_dt <- dt[, .(
  median_similarity = median(similarity),
  q25_similarity = as.numeric(quantile(similarity, 0.25)),
  q75_similarity = as.numeric(quantile(similarity, 0.75))
), by = rank]
setorder(rank_dt, rank)

p2 <- ggplot(rank_dt, aes(x = rank, y = median_similarity)) +
  geom_ribbon(aes(ymin = q25_similarity, ymax = q75_similarity),
              fill = "#A6CEE3", alpha = 0.5) +
  geom_line(color = "#1F77B4", linewidth = 0.8) +
  geom_hline(yintercept = main_cutoff, color = "#D62728", linetype = "dashed") +
  labs(
    title = "Rank-similarity curve",
    x = "RNA neighbor rank",
    y = "Relative similarity"
  ) +
  theme_bw(base_size = 12)
ggsave(
  file.path(paths$fig_dir, sprintf("Fig_rank_similarity_cutoff_%s.png", main_label)),
  p2,
  width = 7,
  height = 4.8,
  dpi = 300
)

summary_dt <- fread(summary_csv)
by_cell_list <- list()
for (cutoff in cutoffs) {
  label <- cutoff_label(cutoff)
  f <- file.path(paths$out_dir, sprintf("adaptive_selection_by_cell_cutoff_%s.csv", label))
  if (file.exists(f)) {
    x <- fread(f, select = c("atac_cell", "n_selected"))
    x[, cutoff := cutoff]
    by_cell_list[[label]] <- x
  }
}
by_cell <- rbindlist(by_cell_list, use.names = TRUE)

p3 <- ggplot(by_cell, aes(x = n_selected)) +
  geom_histogram(binwidth = 2, fill = "#59A14F", color = "white", linewidth = 0.1) +
  facet_wrap(~ cutoff, scales = "free_y") +
  labs(
    title = "Selected RNA neighbors per ATAC cell",
    x = "Selected RNA neighbor count",
    y = "ATAC cell count"
  ) +
  theme_bw(base_size = 12)
ggsave(
  file.path(paths$fig_dir, "Fig_selected_neighbor_count_by_cutoff.png"),
  p3,
  width = 9,
  height = 6,
  dpi = 300
)

p4 <- ggplot(summary_dt, aes(x = cutoff, y = median_n_selected)) +
  geom_line(color = "#222222", linewidth = 0.8) +
  geom_point(size = 2.4, color = "#222222") +
  geom_ribbon(aes(ymin = q10_n_selected, ymax = q90_n_selected),
              fill = "#BAB0AC", alpha = 0.35) +
  labs(
    title = "Selected neighbor count across cutoffs",
    x = "Similarity cutoff",
    y = "Selected RNA neighbors per ATAC cell"
  ) +
  theme_bw(base_size = 12)
ggsave(
  file.path(paths$fig_dir, "Fig_selected_neighbor_count_summary.png"),
  p4,
  width = 7,
  height = 4.8,
  dpi = 300
)

cat("  figures:", paths$fig_dir, "\n")
cat("========== Done ==========\n")
