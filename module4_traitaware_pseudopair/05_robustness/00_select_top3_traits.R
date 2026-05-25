# module4 robustness: select the top 3 traits from completed AUROC results
# -----------------------------------------------------------------------------
# Selection rule:
#   1. For each trait, keep baseline_atac and the best v2 result.
#   2. Rank by delta = best_v2_auroc - baseline_auroc, descending.
#   3. Keep the top 3 traits for robustness testing.
#
# With the current results this should select:
#   MONO, RDW, RBC
#
# Run in RStudio:
#   setwd("D:/RStudio/work/module4_traitaware_pseudopair")
#   source("05_robustness/00_select_top3_traits.R")
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

ROOT <- "D:/RStudio/work/module4_traitaware_pseudopair"
OUT_DIR <- file.path(ROOT, "05_robustness/output/top3_embedding_stress")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

AUROC_FILES <- c(
  file.path(
    ROOT,
    "04_evaluation/module4_v2_neighborhood_fdr_ge0p80/module4_v2_neighborhood_fdr_ge0p80_auroc.csv"
  ),
  file.path(
    ROOT,
    "04_evaluation/trait_screening_ge0p80/trait_screening_ge0p80_auroc.csv"
  ),
  file.path(
    ROOT,
    "04_evaluation/trait_screening_remaining_ge0p80/trait_screening_remaining_ge0p80_auroc.csv"
  )
)

stop_missing <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0L) {
    stop("Missing required file(s):\n  ", paste(missing, collapse = "\n  "),
         call. = FALSE)
  }
}

cat(">>> Select top 3 traits for robustness\n\n")
stop_missing(AUROC_FILES)

dt <- rbindlist(lapply(AUROC_FILES, fread), use.names = TRUE, fill = TRUE)
dt <- dt[is.finite(auroc)]

base <- dt[method == "baseline_atac", .(
  baseline_auroc = auroc[1],
  baseline_ci_lo = auroc_ci_lo[1],
  baseline_ci_hi = auroc_ci_hi[1],
  n_positive = n_positive[1],
  n_total = n_total[1]
), by = target]

v2_best <- dt[method == "v2", .SD[which.max(auroc)], by = target]
setnames(v2_best, c("auroc", "auroc_ci_lo", "auroc_ci_hi", "atac_w"),
         c("best_v2_auroc", "best_v2_ci_lo", "best_v2_ci_hi", "best_v2_atac_w"))
v2_best <- v2_best[, .(
  target,
  best_v2_atac_w,
  best_v2_auroc,
  best_v2_ci_lo,
  best_v2_ci_hi
)]

ranked <- merge(base, v2_best, by = "target", all = FALSE)
ranked[, delta_vs_baseline := best_v2_auroc - baseline_auroc]
setorder(ranked, -delta_vs_baseline, -best_v2_auroc)
ranked[, rank_by_delta := seq_len(.N)]

selected <- ranked[rank_by_delta <= 3L]

fwrite(ranked, file.path(OUT_DIR, "all_completed_traits_ranked.csv"))
fwrite(selected, file.path(OUT_DIR, "selected_top3_traits.csv"))

cat("All completed traits ranked by best-v2 delta:\n")
print(ranked[, .(
  rank_by_delta,
  target,
  baseline_auroc = round(baseline_auroc, 4),
  best_v2_auroc = round(best_v2_auroc, 4),
  best_v2_atac_w,
  delta_vs_baseline = round(delta_vs_baseline, 4)
)])

cat("\nSelected top 3 traits:\n")
print(selected[, .(target, best_v2_atac_w, delta_vs_baseline)])

cat("\nOutputs are in: ", OUT_DIR, "\n", sep = "")
