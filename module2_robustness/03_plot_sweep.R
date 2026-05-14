source("D:/RStudio/work/module2_robustness/_bootstrap.R")

required_pkgs(c("data.table", "ggplot2", "jsonlite"))
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(jsonlite)
})

p <- robustness_paths()
ensure_robustness_dirs(p)

result_csv <- file.path(p$out_dir, "robustness_top5_results.csv")
if (!file.exists(result_csv)) stop("Missing robustness results: ", result_csv)
dt <- fread(result_csv)

summary_dt <- dt[, .(
  mean_top5 = mean(target_fraction_top, na.rm = TRUE),
  sd_top5 = if (.N > 1) sd(target_fraction_top, na.rm = TRUE) else 0,
  n = .N
), by = .(dropout_pct, method, trait)]
summary_csv <- file.path(p$out_dir, "robustness_top5_summary.csv")
fwrite(summary_dt, summary_csv)

for (tr in unique(summary_dt$trait)) {
  plot_dt <- summary_dt[trait == tr]
  fig <- ggplot(plot_dt, aes(x = dropout_pct, y = mean_top5, color = method, group = method)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = mean_top5 - sd_top5, ymax = mean_top5 + sd_top5), width = 1.5) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_x_continuous(breaks = sort(unique(plot_dt$dropout_pct))) +
    labs(
      title = sprintf("%s top5 target enrichment under ATAC dropout", tr),
      x = "Dropped nonzero ATAC entries (%)",
      y = "Top5 target fraction"
    ) +
    theme_bw(base_size = 12)
  ggsave(
    file.path(p$figure_dir, sprintf("Fig_robustness_top5_%s.png", tr)),
    fig,
    width = 6.5,
    height = 4.5,
    dpi = 300
  )
}

write_json(
  list(
    script = "03_plot_sweep.R",
    result_csv = result_csv,
    summary_csv = summary_csv,
    figure_dir = p$figure_dir,
    created_at = as.character(Sys.time())
  ),
  file.path(p$out_dir, "robustness_plot_manifest.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)
print(summary_dt)
cat("summary:", summary_csv, "\n")
