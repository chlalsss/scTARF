#!/usr/bin/env Rscript
# Fig3 Panel B: Granja vs CITE RNA suitability score per trait, colored by gate.
# Reads the screening tool's comparison CSV so the chart stays in sync with the run.
# install.packages() intentionally lives elsewhere; this script only loads libs.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

out_dir <- "D:/RStudio/work/module4_traitaware_pseudopair/01_data_prep_rna_reference_screening/output"
src <- file.path(out_dir, "rna_suitability_trait_reference_comparison.csv")
png_out <- file.path(out_dir, "rna_suitability_scores_barplot.png")

cmp <- fread(src)

# story ordering: granja-selected PASS -> cite PASS -> cite RISK -> none/FAIL
trait_order <- c("MONO", "LYM", "MCH", "RDW", "RBC", "EOS", "WBC", "PLT")

gate_short <- function(g) {
  ifelse(grepl("^PASS", g), "PASS",
  ifelse(grepl("^RISK", g), "RISK",
  ifelse(grepl("^FAIL", g), "FAIL", "NONE")))
}

# long format: one row per (trait, reference)
mk <- function(ref) {
  data.table(
    trait     = cmp$trait,
    reference = ref,
    score     = cmp[[paste0(ref, "_rna_suitability_score")]],
    gate      = gate_short(cmp[[paste0(ref, "_final_gate")]]),
    recommend = cmp$recommended_reference == ref
  )
}
dt <- rbind(mk("granja_rna"), mk("cite_rna"))
dt[, trait := factor(trait, levels = trait_order)]
dt[, reference := factor(reference, levels = c("granja_rna", "cite_rna"),
                         labels = c("Granja RNA", "CITE RNA"))]
dt[, gate := factor(gate, levels = c("PASS", "RISK", "FAIL", "NONE"))]

gate_cols <- c(PASS = "#59A14F", RISK = "#E8A33D", FAIL = "#C44E4E", NONE = "#FFFFFF")

# bars with a score; reference with no scDRS score -> NA -> annotated text instead
bars <- dt[!is.na(score)]
noscore <- dt[is.na(score)]

p <- ggplot(bars, aes(x = trait, y = score, fill = gate, group = reference)) +
  geom_col(aes(linewidth = recommend, color = recommend),
           position = position_dodge2(width = 0.8, padding = 0.12), width = 0.78) +
  geom_text(aes(label = sprintf("%.1f", score), fontface = ifelse(recommend, "bold", "plain")),
            position = position_dodge2(width = 0.8, padding = 0.12),
            vjust = -0.4, size = 3.0, color = "#111111") +
  geom_text(data = bars[recommend == TRUE],
            aes(label = "★", y = score + 5.5),
            position = position_dodge2(width = 0.8, padding = 0.12),
            size = 4.2, color = "#111111") +
  scale_fill_manual(values = gate_cols,
                    breaks = c("PASS", "RISK", "FAIL"),
                    labels = c("PASS_PRIOR", "RISK_PRIOR", "FAIL"),
                    name = NULL) +
  scale_color_manual(values = c(`TRUE` = "#111111", `FALSE` = NA), guide = "none") +
  scale_linewidth_manual(values = c(`TRUE` = 1.0, `FALSE` = 0), guide = "none") +
  scale_y_continuous(limits = c(0, 105), expand = c(0, 0)) +
  labs(x = NULL, y = "RNA suitability score (0-100)",
       subtitle = "left bar = Granja RNA, right bar = CITE RNA; black outline + ★ = recommended") +
  theme_classic(base_size = 13) +
  theme(legend.position = "top",
        legend.text = element_text(size = 10),
        plot.subtitle = element_text(size = 9.5, color = "#555555"),
        axis.text.x = element_text(size = 12),
        panel.grid.major.y = element_line(color = "#E8E8E8", linewidth = 0.4))

# "no scDRS score" annotations under the missing (Granja) bars
if (nrow(noscore) > 0) {
  noscore[, xpos := as.numeric(trait) - 0.2]  # Granja is the left (first) bar
  p <- p + geom_text(data = noscore, inherit.aes = FALSE,
                     aes(x = xpos, y = 6, label = "no\nscDRS\nscore"),
                     size = 2.4, fontface = "italic", color = "#9A9A9A", lineheight = 0.9)
}

ggsave(png_out, p, width = 8.2, height = 4.8, dpi = 200, bg = "white")
cat("saved", png_out, "\n")
