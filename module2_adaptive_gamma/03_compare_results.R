# 03: Compare 6-condition results, with two acceptance verdicts:
#   (a) Direction 2 alone:        C2/C3 vs C0
#   (b) Direction 1+2 stacking:   C4/C5 vs C1
#
# Also cross-checks C0 against direction-1 04 output to validate that our
# custom hetgamma RWR matches SCAVENGE's randomWalk_sparse when gamma is
# constant.
#
# Outputs:
#   acceptance_summary.csv        # per-trait verdict table
#   sanity_C0_vs_direction1.csv   # absolute TRS delta vs direction-1 04
#   compare_results_manifest.json
#   figures/Fig_condition_top5.png

source("D:/RStudio/work/module2_adaptive_gamma/_bootstrap.R")

required_pkgs(c("data.table", "ggplot2", "jsonlite"))
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(jsonlite)
})

paths <- adaptive_gamma_paths()
ensure_adaptive_gamma_dirs(paths)

traits <- parse_csv_env("MODULE2_TRAITS", "MONO,PLT,LYM")
direction1_threshold <- env_numeric("MODULE2_SEED_FILTER_PRIMARY", 0.50)
acceptance_pp <- env_numeric("MODULE2_ADAPTIVE_GAMMA_ACCEPT_PP", 1.0)

cat("========== Module 2 adaptive gamma: 03 compare results ==========\n")
cat("  acceptance threshold:", acceptance_pp, "pp\n")

top5_path <- file.path(paths$out_dir, "condition_top5_summary.csv")
if (!file.exists(top5_path)) {
  stop("Missing 02 output: ", top5_path, "\nRun 02_run_scavenge_adaptive_gamma.R first.",
       call. = FALSE)
}
top5_dt <- fread(top5_path)

# Acceptance verdicts ---------------------------------------------------

cat("\n[1/3] Acceptance check\n")

verdict_rows <- list()
for (tr in traits) {
  trait_name <- tr
  tr_dt <- top5_dt[trait == trait_name]
  get_top <- function(cid) {
    v <- tr_dt[condition == cid, top_target_fraction]
    if (length(v) == 0) NA_real_ else v[1]
  }
  C0 <- get_top("C0")
  C1 <- get_top("C1")
  C2 <- get_top("C2")
  C3 <- get_top("C3")
  C4 <- get_top("C4")
  C5 <- get_top("C5")

  # Direction 2 alone vs C0
  dir2_A_vs_baseline <- if (is.finite(C2) && is.finite(C0)) (C2 - C0) * 100 else NA_real_
  dir2_B_vs_baseline <- if (is.finite(C3) && is.finite(C0)) (C3 - C0) * 100 else NA_real_
  # Direction 1+2 vs direction 1 alone
  stack_A_vs_dir1 <- if (is.finite(C4) && is.finite(C1)) (C4 - C1) * 100 else NA_real_
  stack_B_vs_dir1 <- if (is.finite(C5) && is.finite(C1)) (C5 - C1) * 100 else NA_real_

  verdict_rows[[tr]] <- data.table(
    trait = tr,
    C0_orig_const = C0,
    C1_filt_const = C1,
    C2_orig_A = C2,
    C3_orig_B = C3,
    C4_filt_A = C4,
    C5_filt_B = C5,
    dir2_A_vs_baseline_pp = dir2_A_vs_baseline,
    dir2_B_vs_baseline_pp = dir2_B_vs_baseline,
    stack_A_vs_dir1_pp = stack_A_vs_dir1,
    stack_B_vs_dir1_pp = stack_B_vs_dir1,
    target_definition = target_definition(tr)
  )
}
verdict_dt <- rbindlist(verdict_rows, use.names = TRUE, fill = TRUE)

# Per-trait flags
verdict_dt[, accept_dir2_A_alone := dir2_A_vs_baseline_pp >= acceptance_pp]
verdict_dt[, accept_dir2_B_alone := dir2_B_vs_baseline_pp >= acceptance_pp]
verdict_dt[, accept_stack_A := stack_A_vs_dir1_pp >= acceptance_pp]
verdict_dt[, accept_stack_B := stack_B_vs_dir1_pp >= acceptance_pp]

verdict_path <- file.path(paths$out_dir, "acceptance_summary.csv")
fwrite(verdict_dt, verdict_path)
cat("\n  Acceptance verdict table:\n")
print(verdict_dt[, .(trait, C0_orig_const, C2_orig_A, dir2_A_vs_baseline_pp,
                     accept_dir2_A_alone, C3_orig_B, dir2_B_vs_baseline_pp,
                     accept_dir2_B_alone)])
cat("\n")
print(verdict_dt[, .(trait, C1_filt_const, C4_filt_A, stack_A_vs_dir1_pp,
                     accept_stack_A, C5_filt_B, stack_B_vs_dir1_pp,
                     accept_stack_B)])

# Aggregate: how many traits pass each criterion
agg <- data.table(
  criterion = c(
    "direction 2 (Version A) alone vs C0",
    "direction 2 (Version B) alone vs C0",
    "direction 1+2 (Version A) stacking vs C1",
    "direction 1+2 (Version B) stacking vs C1"
  ),
  n_traits_pass = c(
    sum(verdict_dt$accept_dir2_A_alone, na.rm = TRUE),
    sum(verdict_dt$accept_dir2_B_alone, na.rm = TRUE),
    sum(verdict_dt$accept_stack_A, na.rm = TRUE),
    sum(verdict_dt$accept_stack_B, na.rm = TRUE)
  ),
  n_traits_total = nrow(verdict_dt)
)
cat("\n  Aggregate pass count:\n")
print(agg)
agg_path <- file.path(paths$out_dir, "acceptance_aggregate.csv")
fwrite(agg, agg_path)

# Sanity: C0 vs direction-1 04 baseline ---------------------------------

cat("\n[2/3] Sanity: C0 (our hetgamma RWR with const gamma) vs direction-1 04 baseline\n")
sanity_rows <- list()
dir1_threshold_label <- seed_filter_label(direction1_threshold)
dir1_top5_path <- file.path(
  sprintf(paths$direction1_trs_dir_template, dir1_threshold_label),
  "trs_top5_enrichment_summary.csv"
)
if (file.exists(dir1_top5_path)) {
  dir1_top5 <- fread(dir1_top5_path)
  for (tr in traits) {
    trait_name <- tr
    C0_val <- top5_dt[condition == "C0" & trait == trait_name, top_target_fraction]
    dir1_baseline <- dir1_top5[trait == trait_name, top_original_target_fraction]
    if (length(C0_val) > 0 && length(dir1_baseline) > 0) {
      sanity_rows[[tr]] <- data.table(
        trait = tr,
        C0_top5 = C0_val[1],
        direction1_04_top_original = dir1_baseline[1],
        abs_diff_pp = abs(C0_val[1] - dir1_baseline[1]) * 100
      )
    }
  }
} else {
  warning("Missing direction-1 04 summary: ", dir1_top5_path)
}
if (length(sanity_rows) > 0) {
  sanity_dt <- rbindlist(sanity_rows, use.names = TRUE, fill = TRUE)
  sanity_path <- file.path(paths$out_dir, "sanity_C0_vs_direction1.csv")
  fwrite(sanity_dt, sanity_path)
  cat("\n  C0 hetgamma-RWR top5 vs direction-1 04 baseline:\n")
  print(sanity_dt)
  max_diff <- max(sanity_dt$abs_diff_pp, na.rm = TRUE)
  if (is.finite(max_diff) && max_diff > 0.5) {
    warning(sprintf(
      "C0 differs from direction-1 04 baseline by %.2f pp on at least one trait. ",
      max_diff
    ),
    "Hetgamma RWR may not be matching SCAVENGE's randomWalk_sparse exactly.")
  } else {
    cat("\n  Sanity OK: max |C0 - dir1_baseline| =", round(max_diff, 3), "pp\n")
  }
} else {
  sanity_dt <- data.table()
  sanity_path <- NA_character_
}

# Plot ------------------------------------------------------------------

cat("\n[3/3] Plot\n")
plot_dt <- top5_dt[, .(condition, trait, top_target_fraction)]
plot_dt[, condition := factor(condition, levels = c("C0", "C1", "C2", "C3", "C4", "C5"))]

p <- ggplot(plot_dt, aes(x = condition, y = top_target_fraction * 100, fill = condition)) +
  geom_col() +
  facet_wrap(~ trait, ncol = 1, scales = "free_y") +
  geom_text(aes(label = sprintf("%.1f", top_target_fraction * 100)),
            vjust = -0.3, size = 3) +
  labs(
    title = "Top 5% TRS target enrichment across 6 conditions",
    subtitle = "C0=orig+const  C1=filt+const  C2=orig+A  C3=orig+B  C4=filt+A  C5=filt+B",
    x = NULL,
    y = "Top 5% target cell fraction (%)"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

fig_path <- file.path(paths$fig_dir, "Fig_condition_top5.png")
ggsave(fig_path, p, width = 8, height = 8, dpi = 300)
cat("  wrote:", fig_path, "\n")

# Manifest --------------------------------------------------------------

manifest <- list(
  script = "03_compare_results.R",
  acceptance_pp = acceptance_pp,
  outputs = list(
    acceptance_summary = verdict_path,
    acceptance_aggregate = agg_path,
    sanity = sanity_path,
    figure = fig_path
  ),
  notes = c(
    "Direction 2 alone passes if C2 or C3 beats C0 by >= 1pp on >= 2 traits.",
    "Direction 1+2 stacking adds value if C4 or C5 beats C1 by >= 1pp on >= 2 traits.",
    "Sanity check requires |C0 - direction1_04 baseline| < 0.5pp per trait."
  ),
  created_at = as.character(Sys.time())
)
manifest_out <- file.path(paths$out_dir, "compare_results_manifest.json")
write_json(manifest, manifest_out, pretty = TRUE, auto_unbox = TRUE)

cat("\n========== Done ==========\n")
cat("  manifest:", manifest_out, "\n")
