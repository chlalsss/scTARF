source("D:/RStudio/work/module2_robustness/_bootstrap.R")

required_pkgs(c("data.table", "SummarizedExperiment", "jsonlite"))
suppressPackageStartupMessages({
  library(data.table)
  library(SummarizedExperiment)
  library(jsonlite)
})

p <- robustness_paths()
ensure_robustness_dirs(p)

manifest_csv <- file.path(p$out_dir, "sparse_input_manifest.csv")
if (!file.exists(manifest_csv)) {
  stop("Missing sparse input manifest. Run 01_generate_sparse_atac.R first: ", manifest_csv)
}
inputs <- fread(manifest_csv)
traits <- trimws(strsplit(Sys.getenv("MODULE2_TRAITS", "MONO,PLT,LYM"), ",", fixed = TRUE)[[1]])

baseline_knn_path <- Sys.getenv(
  "MODULE2_ROBUSTNESS_BASELINE_KNN",
  file.path(p$adaptive_dir, "atac_only_mutualknn_k23.rds")
)
fused_knn_path <- Sys.getenv(
  "MODULE2_ROBUSTNESS_FUSED_KNN",
  file.path(p$adaptive_dir, "fused_mutualknn_cutoff_top1_pca23_k16.rds")
)
stopifnot(file.exists(baseline_knn_path), file.exists(fused_knn_path))

baseline_k <- as.integer(Sys.getenv("MODULE2_ROBUSTNESS_BASELINE_K", "23"))
fused_k <- as.integer(Sys.getenv("MODULE2_ROBUSTNESS_FUSED_K", "16"))
fused_npcs <- as.integer(Sys.getenv("MODULE2_ROBUSTNESS_FUSED_NPCS", "23"))

top5_fraction <- function(trs_path, se_path, method, row) {
  se <- readRDS(se_path)
  meta <- data.table(cell = colnames(se), BioClassification = as.character(colData(se)$BioClassification))
  trs <- as.data.table(readRDS(trs_path))
  rbindlist(lapply(traits, function(tr) {
    if (!tr %in% colnames(trs)) stop("Missing trait ", tr, " in ", trs_path)
    dt <- merge(trs[, .(cell, score = get(tr))], meta, by = "cell", all.x = TRUE)
    dt[, target := target_type(tr, BioClassification)]
    n_top <- max(1L, ceiling(nrow(dt) * 0.05))
    idx <- order(dt$score, decreasing = TRUE)[seq_len(n_top)]
    data.table(
      dropout_pct = row$dropout_pct,
      seed = row$seed,
      label = row$label,
      method = method,
      trait = tr,
      n_top = n_top,
      target_fraction_top = mean(dt$target[idx], na.rm = TRUE),
      target_fraction_all = mean(dt$target, na.rm = TRUE),
      enrichment_ratio = mean(dt$target[idx], na.rm = TRUE) / mean(dt$target, na.rm = TRUE),
      trs_rds = trs_path
    )
  }))
}

run_07 <- function(row, method, knn_path, k, npcs) {
  run_label <- file.path("robustness_runs", row$label, method)
  Sys.setenv(
    MODULE2_OUT_DIR = p$result_dir,
    MODULE2_ATAC_SE = row$input_rds,
    MODULE2_KNN_PATH = knn_path,
    MODULE2_RUN_LABEL = run_label,
    MODULE2_TOPK = "1",
    MODULE2_FUSED_NPCS = as.character(npcs),
    MODULE2_KNN_K = as.character(k),
    MODULE2_TRAITS = paste(traits, collapse = ",")
  )
  Sys.unsetenv("MODULE2_CUTOFF_LABEL")
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd("D:/RStudio/work/module2_adaptive_fusion")
  source("D:/RStudio/work/module2_adaptive_fusion/_config.R", local = .GlobalEnv)
  run_env <- new.env(parent = .GlobalEnv)
  source("D:/RStudio/work/module2_adaptive_fusion/07_run_scavenge_fused_knn.R", local = run_env)
  file.path(p$result_dir, run_label, "TRS_matrix_fusedKNN_selected_traits_fullcells.rds")
}

all_results <- list()
for (ii in seq_len(nrow(inputs))) {
  row <- inputs[ii]
  cat("\n========== Robustness input:", row$label, "==========\n")
  base_trs <- run_07(row, "baseline_atac_only_auto", baseline_knn_path, baseline_k, fused_npcs)
  fused_trs <- run_07(row, "fused_top1_auto", fused_knn_path, fused_k, fused_npcs)
  all_results[[paste0(row$label, "_baseline")]] <- top5_fraction(
    base_trs, row$input_rds, "baseline_atac_only_auto", row
  )
  all_results[[paste0(row$label, "_fused")]] <- top5_fraction(
    fused_trs, row$input_rds, "fused_top1_auto", row
  )
}

result_dt <- rbindlist(all_results, use.names = TRUE, fill = TRUE)
result_csv <- file.path(p$out_dir, "robustness_top5_results.csv")
fwrite(result_dt, result_csv)
write_json(
  list(
    script = "02_run_mini_sweep.R",
    sparse_manifest = manifest_csv,
    baseline_knn = baseline_knn_path,
    fused_knn = fused_knn_path,
    baseline_k = baseline_k,
    fused_k = fused_k,
    fused_npcs = fused_npcs,
    result_csv = result_csv,
    created_at = as.character(Sys.time())
  ),
  file.path(p$out_dir, "robustness_run_manifest.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)
print(result_dt)
cat("results:", result_csv, "\n")
