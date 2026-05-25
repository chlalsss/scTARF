suppressPackageStartupMessages({
  library(data.table)
})

bmmc_dir <- "D:/SCAVENGE-reproducibility-main/BMMC"
meta_file <- file.path(bmmc_dir, "01_data_qc", "BMMC_cell_meta.rds")
dual_dir <- file.path(bmmc_dir, "04_dual_modality")
robust_dir <- file.path(bmmc_dir, "05_robustness")
out_dir <- "D:/RStudio/work/paper_metric_supplement"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

meta <- readRDS(meta_file)
if (inherits(meta, "data.table")) meta <- as.data.frame(meta)
if (is.null(rownames(meta)) || identical(rownames(meta), as.character(seq_len(nrow(meta))))) {
  if (!"barcode" %in% names(meta)) stop("meta has no rownames or barcode column")
  rownames(meta) <- meta$barcode
}
if (!"cell_type" %in% names(meta)) stop("meta missing cell_type")

targets <- list(
  MONO = c("CD14+ Mono", "CD16+ Mono"),
  PLT = c("MK/E prog"),
  PLT_broad = c("MK/E prog", "Proerythroblast", "Erythroblast", "Normoblast"),
  LYM = c("B1 B", "Naive CD20+ B", "Transitional B", "Plasma cell"),
  LYM_TNK = c("NK", "ILC", "CD4+ T activated", "CD4+ T naive", "CD8+ T", "CD8+ T naive"),
  LYM_wide = c(
    "NK", "ILC", "CD4+ T activated", "CD4+ T naive", "CD8+ T", "CD8+ T naive",
    "Lymph prog", "B1 B", "Naive CD20+ B", "Transitional B", "Plasma cell"
  )
)
pred_trait <- c(
  MONO = "MONO", PLT = "PLT", PLT_broad = "PLT",
  LYM = "LYM", LYM_TNK = "LYM", LYM_wide = "LYM"
)

metric_from_vector <- function(score, target_label, section, method, trait, target,
                               p_drop = NA_real_, seed = NA_integer_, atac_w = NA_real_) {
  common <- intersect(names(score), rownames(meta))
  score <- score[common]
  lbl <- as.character(meta[common, "cell_type"]) %in% targets[[target_label]]
  ord <- order(score, decreasing = TRUE, na.last = NA)
  n_total <- length(ord)
  n_top5 <- floor(n_total * 0.05)
  if (n_top5 < 1) stop("no top5 cells for ", method, " ", trait, " ", target)
  top_idx <- ord[seq_len(n_top5)]
  target_frac_top5 <- mean(lbl[top_idx])
  target_frac_bg <- mean(lbl)
  data.table(
    section = section,
    method = method,
    trait = trait,
    target = target,
    p_drop = p_drop,
    seed = seed,
    atac_w = atac_w,
    n_total = n_total,
    n_top5 = n_top5,
    n_top5_in_target = sum(lbl[top_idx]),
    n_target_total = sum(lbl),
    target_frac_top5 = target_frac_top5,
    target_frac_bg = target_frac_bg,
    enrichment_ratio = target_frac_top5 / target_frac_bg
  )
}

read_trs_file <- function(path) {
  dt <- fread(path)
  bc_col <- intersect(c("barcode", "cell", "rn", "V1"), names(dt))[1]
  trs_col <- intersect(c("TRS", "trs", "score"), names(dt))[1]
  if (is.na(bc_col) || is.na(trs_col)) {
    bc_col <- names(dt)[1]
    trs_col <- names(dt)[2]
  }
  score <- dt[[trs_col]]
  names(score) <- dt[[bc_col]]
  score
}

# 4.4: paper Table 4.5 method set.
baseline_file <- file.path(bmmc_dir, "03_scavenge_baseline", "TRS_matrix_baseline.rds")
baseline_mat <- readRDS(baseline_file)
if (inherits(baseline_mat, "data.table")) baseline_mat <- as.data.frame(baseline_mat)
if ("cell" %in% names(baseline_mat)) {
  rownames(baseline_mat) <- baseline_mat$cell
  baseline_mat$cell <- NULL
}

main_targets <- c("MONO", "PLT", "LYM_wide")
supplement_targets <- c("MONO", "PLT", "PLT_broad", "LYM_wide")
res44 <- rbindlist(lapply(main_targets, function(tgt) {
  tr <- unname(pred_trait[tgt])
  score <- baseline_mat[[tr]]
  names(score) <- rownames(baseline_mat)
  metric_from_vector(score, tgt, "4.4", "ATAC-only baseline", tr, tgt)
}))

for (w in c(1.0, 0.7)) {
  for (tgt in main_targets) {
    tr <- unname(pred_trait[tgt])
    f <- file.path(dual_dir, sprintf("TRS_v2_%s_w%.1f.tsv.gz", tr, w))
    res44 <- rbind(
      res44,
      metric_from_vector(
        read_trs_file(f), tgt, "4.4", sprintf("v2 atac_w=%.1f", w), tr, tgt, atac_w = w
      ),
      fill = TRUE
    )
  }
}
setorder(res44, target, method)
fwrite(res44, file.path(out_dir, "bmmc_4_4_top5_enrichment_metrics.csv"))

res44_with_broad_plt <- rbindlist(lapply(supplement_targets, function(tgt) {
  tr <- unname(pred_trait[tgt])
  score <- baseline_mat[[tr]]
  names(score) <- rownames(baseline_mat)
  metric_from_vector(score, tgt, "4.4", "ATAC-only baseline", tr, tgt)
}))
for (w in c(1.0, 0.7)) {
  for (tgt in supplement_targets) {
    tr <- unname(pred_trait[tgt])
    f <- file.path(dual_dir, sprintf("TRS_v2_%s_w%.1f.tsv.gz", tr, w))
    res44_with_broad_plt <- rbind(
      res44_with_broad_plt,
      metric_from_vector(
        read_trs_file(f), tgt, "4.4", sprintf("v2 atac_w=%.1f", w), tr, tgt, atac_w = w
      ),
      fill = TRUE
    )
  }
}
setorder(res44_with_broad_plt, target, method)
fwrite(res44_with_broad_plt, file.path(out_dir, "bmmc_4_4_top5_enrichment_metrics_with_broad_plt.csv"))

# 4.5: robustness TRS files, all seeds/dropout for the three main paper targets.
robust_files <- list.files(
  robust_dir,
  pattern = "^TRS_p[0-9.]+_seed[0-9]+_(MONO|PLT|LYM)_w[0-9.]+\\.tsv\\.gz$",
  full.names = TRUE
)
parse_robust <- function(f) {
  m <- regmatches(
    basename(f),
    regexec("^TRS_p([0-9.]+)_seed([0-9]+)_([A-Z]+)_w([0-9.]+)\\.tsv\\.gz$", basename(f))
  )[[1]]
  list(p_drop = as.numeric(m[2]), seed = as.integer(m[3]), trait = m[4], atac_w = as.numeric(m[5]), file = f)
}
robust_meta <- lapply(robust_files, parse_robust)

res45 <- rbindlist(lapply(robust_meta, function(x) {
  target_set <- names(pred_trait)[pred_trait == x$trait]
  target_set <- intersect(target_set, supplement_targets)
  rbindlist(lapply(target_set, function(tgt) {
    metric_from_vector(
      read_trs_file(x$file), tgt, "4.5", sprintf("atac_w=%.1f", x$atac_w), x$trait, tgt,
      p_drop = x$p_drop, seed = x$seed, atac_w = x$atac_w
    )
  }))
}), fill = TRUE)
setorder(res45, target, p_drop, atac_w, seed)
fwrite(res45, file.path(out_dir, "bmmc_4_5_top5_enrichment_metrics_by_seed.csv"))

res45_summary <- res45[, .(
  n = .N,
  mean_top5_target_frac = mean(target_frac_top5),
  sd_top5_target_frac = sd(target_frac_top5),
  lo95_top5_target_frac = mean(target_frac_top5) - 1.96 * sd(target_frac_top5) / sqrt(.N),
  hi95_top5_target_frac = mean(target_frac_top5) + 1.96 * sd(target_frac_top5) / sqrt(.N),
  background_target_frac = mean(target_frac_bg),
  mean_enrichment_ratio = mean(enrichment_ratio),
  sd_enrichment_ratio = sd(enrichment_ratio),
  lo95_enrichment_ratio = mean(enrichment_ratio) - 1.96 * sd(enrichment_ratio) / sqrt(.N),
  hi95_enrichment_ratio = mean(enrichment_ratio) + 1.96 * sd(enrichment_ratio) / sqrt(.N)
), by = .(section, trait, target, p_drop, atac_w)]
setorder(res45_summary, target, p_drop, atac_w)
fwrite(res45_summary, file.path(out_dir, "bmmc_4_5_top5_enrichment_metrics_summary.csv"))

res45_summary_main <- res45_summary[target %in% main_targets]
fwrite(res45_summary_main, file.path(out_dir, "bmmc_4_5_top5_enrichment_metrics_summary_main_targets.csv"))

cat("Saved:\n")
cat(" - ", file.path(out_dir, "bmmc_4_4_top5_enrichment_metrics.csv"), "\n", sep = "")
cat(" - ", file.path(out_dir, "bmmc_4_4_top5_enrichment_metrics_with_broad_plt.csv"), "\n", sep = "")
cat(" - ", file.path(out_dir, "bmmc_4_5_top5_enrichment_metrics_by_seed.csv"), "\n", sep = "")
cat(" - ", file.path(out_dir, "bmmc_4_5_top5_enrichment_metrics_summary.csv"), "\n", sep = "")
cat(" - ", file.path(out_dir, "bmmc_4_5_top5_enrichment_metrics_summary_main_targets.csv"), "\n", sep = "")
cat("\n4.4 main metrics:\n")
print(res44[, .(
  method, trait, target,
  top5 = sprintf("%.1f%%", 100 * target_frac_top5),
  background = sprintf("%.1f%%", 100 * target_frac_bg),
  enrichment = sprintf("%.2fx", enrichment_ratio)
)])

cat("\n4.5 p=0 and p=0.9 summary:\n")
print(res45_summary[p_drop %in% c(0, 0.9), .(
  trait, target, p_drop, atac_w,
  top5_mean = sprintf("%.1f%%", 100 * mean_top5_target_frac),
  background = sprintf("%.1f%%", 100 * background_target_frac),
  enrichment_mean = sprintf("%.2fx", mean_enrichment_ratio)
)])
