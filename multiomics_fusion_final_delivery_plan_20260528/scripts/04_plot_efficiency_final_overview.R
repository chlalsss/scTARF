# Plot final efficiency overview from the formal MONO benchmark table.

suppressPackageStartupMessages(library(ggplot2))

root <- "D:/RStudio/work/multiomics_fusion_final_delivery_plan_20260528"
src <- file.path(root, "final_outputs/04_robustness_efficiency/efficiency/efficiency_runtime_memory_table.csv")
x <- read.csv(src, stringsAsFactors = FALSE)
x <- x[x$status == "OK", , drop = FALSE]

x$method <- ifelse(x$pairing_setting == "paired", "Paired-Fusion", "Non-paired v3.3")
x$cell_scale <- factor(
  x$n_cells,
  levels = c(1000, 10000, 50000, 100000),
  labels = c("1K", "10K", "50K", "100K")
)

mk <- function(metric, value, label) {
  data.frame(
    method = x$method,
    cell_scale = x$cell_scale,
    metric = metric,
    value = value,
    label = label,
    stringsAsFactors = FALSE
  )
}

plot_dt <- rbind(
  mk("Runtime (s)", x$runtime_seconds, sprintf("%.1f", x$runtime_seconds)),
  mk("Memory (GB)", x$peak_memory_gb, sprintf("%.2f", x$peak_memory_gb)),
  mk("AUROC", x$auroc, sprintf("%.3f", x$auroc)),
  mk("ACC", x$acc_top_prevalence, sprintf("%.3f", x$acc_top_prevalence))
)
plot_dt$metric <- factor(plot_dt$metric, levels = c("Runtime (s)", "Memory (GB)", "AUROC", "ACC"))

p <- ggplot(plot_dt, aes(x = cell_scale, y = value, color = method, group = method)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  geom_text(aes(label = label), vjust = -0.75, size = 3, show.legend = FALSE) +
  facet_wrap(~ metric, scales = "free_y", ncol = 2) +
  scale_color_manual(values = c("Paired-Fusion" = "#2166AC", "Non-paired v3.3" = "#B2182B")) +
  labs(
    title = "Final efficiency benchmark (MONO)",
    subtitle = "Paired and non-paired final methods across 1K, 10K, 50K and 100K cells",
    x = "Number of cells",
    y = NULL,
    color = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

ascii_out <- file.path(root, "final_outputs/04_robustness_efficiency/efficiency/fig_efficiency_final_overview.png")
ggsave(ascii_out, p, width = 9, height = 6, dpi = 300)

final_name <- intToUtf8(c(0x6700, 0x7ec8, 0x8f93, 0x51fa))
eff_name <- paste0("04_", intToUtf8(c(0x6548, 0x7387)))
final_out <- file.path(root, "final_outputs", final_name, eff_name, "fig_efficiency_final_overview.png")
file.copy(ascii_out, final_out, overwrite = TRUE)

cat("Wrote overview figure:\n")
cat("  ", ascii_out, "\n", sep = "")
cat("  ", final_out, "\n", sep = "")
