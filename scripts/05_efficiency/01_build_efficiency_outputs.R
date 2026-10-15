# Build data-scale efficiency outputs for final delivery.
#
# The formal benchmark is intentionally scoped to one trait by default
# (MONO), because it rebuilds the fused KNN and RWR/TRS for each scale.
# It reports paired and non-paired settings separately, using:
#   - paired: BMMC paired multiome fusion
#   - non-paired: module4 v3.3 sparse mapping + top5 significant completion

source(file.path("scripts", "05_efficiency", "00_final_delivery_common.R"))
source(file.path("scripts", "_paths.R"))
ensure_final_dirs()
need_packages(c("ggplot2", "Matrix", "BiocNeighbors", "data.table", "pROC"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(Matrix)
  library(BiocNeighbors)
  library(data.table)
  library(pROC)
})

out_dir <- file.path(FINAL_ROOT, "04_robustness_efficiency/efficiency")
dir_ensure(out_dir)

trait_to_run <- Sys.getenv("EFFICIENCY_TRAIT", "MONO")
cell_sizes_env <- Sys.getenv("EFFICIENCY_CELL_SIZES", "")
cell_sizes <- if (nzchar(cell_sizes_env)) {
  as.integer(strsplit(cell_sizes_env, "[,; ]+")[[1]])
} else {
  c(1000L, 10000L, 50000L, 100000L)
}
cell_sizes <- cell_sizes[is.finite(cell_sizes) & cell_sizes > 0L]
if (length(cell_sizes) == 0L) stop("No valid EFFICIENCY_CELL_SIZES supplied.", call. = FALSE)

force_rerun <- tolower(Sys.getenv("EFFICIENCY_FORCE_RERUN", "false")) %in% c("1", "true", "yes")
atac_w_default <- as.numeric(Sys.getenv("EFFICIENCY_ATAC_W", "0.7"))
if (!is.finite(atac_w_default) || atac_w_default < 0 || atac_w_default > 1) {
  stop("EFFICIENCY_ATAC_W must be between 0 and 1.", call. = FALSE)
}

formal_file <- file.path(out_dir, "efficiency_runtime_memory_table.csv")
step_file <- file.path(out_dir, "efficiency_step_times_formal_benchmark.csv")

make_elapsed_demo_outputs <- function() {
  elapsed_files <- list.files(
    ACCEPTANCE_ROOT,
    pattern = "^elapsed_time_summary\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  )

  elapsed_rows <- lapply(elapsed_files, function(f) {
    x <- read_csv_if_exists(f)
    if (is.null(x) || nrow(x) == 0L) return(NULL)
    x$run_id <- basename(dirname(f))
    x$source_file <- f
    x
  })
  elapsed <- bind_rows_fill(elapsed_rows)

  if (nrow(elapsed) == 0L) return(invisible(NULL))

  elapsed$elapsed_seconds <- suppressWarnings(as.numeric(elapsed$elapsed_seconds))
  write_csv(elapsed, file.path(out_dir, "efficiency_elapsed_step_times_from_demo_runs.csv"))

  elapsed_plot <- elapsed[is.finite(elapsed$elapsed_seconds), , drop = FALSE]
  total <- aggregate(elapsed_seconds ~ run_id + source_file, elapsed_plot, sum, na.rm = TRUE)
  names(total)[names(total) == "elapsed_seconds"] <- "total_elapsed_seconds"
  total$benchmark_type <- "acceptance_demo"
  total$n_cells <- NA_integer_
  total$peak_memory_gb <- NA_real_
  total$auroc <- NA_real_
  write_csv(total, file.path(out_dir, "efficiency_total_runtime_from_demo_runs.csv"))

  p1 <- ggplot(elapsed_plot, aes(x = reorder(step_id, elapsed_seconds), y = elapsed_seconds, fill = status)) +
    geom_col() +
    coord_flip() +
    facet_wrap(~ run_id, scales = "free_y") +
    labs(x = NULL, y = "Elapsed seconds", fill = NULL, title = "Step runtime in existing demo runs") +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank())
  ggsave(file.path(out_dir, "fig_efficiency_step_time_demo.png"), p1, width = 12, height = 7, dpi = 300)

  p2 <- ggplot(total, aes(x = run_id, y = total_elapsed_seconds)) +
    geom_col(fill = "#2B6CB0") +
    coord_flip() +
    labs(x = NULL, y = "Total elapsed seconds", title = "Total runtime in existing demo runs") +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank())
  ggsave(file.path(out_dir, "fig_efficiency_total_time_demo.png"), p2, width = 9, height = 5, dpi = 300)
}

gc_memory_gb <- function() {
  z <- gc()
  if (!all(c("used") %in% colnames(z))) return(NA_real_)
  n_cells <- if ("Ncells" %in% rownames(z)) z["Ncells", "used"] else 0
  v_cells <- if ("Vcells" %in% rownames(z)) z["Vcells", "used"] else 0
  (as.numeric(n_cells) * 56 + as.numeric(v_cells) * 8) / 1024^3
}

object_memory_gb <- function(...) {
  objs <- list(...)
  sum(vapply(objs, object.size, numeric(1))) / 1024^3
}

normalise_seed <- function(seed, cells) {
  if (is.null(names(seed))) stop("Seed vector has no names.", call. = FALSE)
  s <- as.numeric(seed[cells])
  s[!is.finite(s) | is.na(s)] <- 0
  s
}

sample_indices <- function(labels, seed_vec, n_out, rng_seed) {
  set.seed(rng_seed)
  n_available <- length(labels)
  replace <- n_out > n_available
  idx <- sample.int(n_available, n_out, replace = replace)

  if (sum(labels[idx] == 1L, na.rm = TRUE) < 10L && sum(labels == 1L, na.rm = TRUE) >= 10L) {
    pos <- which(labels == 1L)
    idx[seq_len(min(10L, length(idx)))] <- sample(pos, min(10L, length(pos)), replace = TRUE)
  }
  if (sum(seed_vec[idx] > 0, na.rm = TRUE) == 0L && sum(seed_vec > 0, na.rm = TRUE) > 0L) {
    seed_pos <- which(seed_vec > 0)
    idx[1] <- sample(seed_pos, 1L)
  }
  idx
}

build_knn <- function(emb, k = 30L) {
  n <- nrow(emb)
  k <- min(k, max(1L, n - 1L))
  nn <- BiocNeighbors::findKNN(emb, k = k, BNPARAM = AnnoyParam())
  median_d <- apply(nn$distance, 1, median)
  median_d[median_d == 0 | !is.finite(median_d)] <- 1
  edge_w <- exp(-nn$distance / median_d)

  i_idx <- rep(seq_len(n), each = k)
  j_idx <- as.integer(t(nn$index))
  w_idx <- as.numeric(t(edge_w))
  A <- sparseMatrix(i = i_idx, j = j_idx, x = w_idx, dims = c(n, n))
  At <- t(A)
  mutual_mask <- (A != 0) * (At != 0)
  A_mut <- ((A + At) / 2) * mutual_mask
  A_mut <- drop0(A_mut)
  rs <- rowSums(A_mut)
  rs[rs == 0] <- 1
  A_norm <- Diagonal(n, 1 / rs) %*% A_mut
  rownames(A_norm) <- rownames(emb)
  colnames(A_norm) <- rownames(emb)
  list(adj = A_norm, iso_frac = mean(diff(A_mut@p) == 0))
}

rwr <- function(A_norm, seed_vec, gamma = 0.05, tol = 1e-6, max_iter = 200L) {
  if (sum(seed_vec > 0, na.rm = TRUE) == 0L) {
    stop("No positive seed cells in sampled benchmark subset.", call. = FALSE)
  }
  s <- seed_vec / sum(seed_vec)
  x <- s
  for (iter in seq_len(max_iter)) {
    x_new <- gamma * s + (1 - gamma) * as.numeric(A_norm %*% x)
    delta <- sum(abs(x_new - x))
    x <- x_new
    if (delta < tol) break
  }
  x
}

post_process_trs <- function(trs, cap_q = 0.95) {
  cap <- as.numeric(stats::quantile(trs, cap_q, na.rm = TRUE))
  trs_c <- pmin(trs, cap)
  rng <- range(trs_c, na.rm = TRUE)
  if (!all(is.finite(rng)) || rng[2] - rng[1] <= 0) return(rep(0, length(trs)))
  (trs_c - rng[1]) / (rng[2] - rng[1])
}

metric_from_scores <- function(scores, labels) {
  keep <- is.finite(scores) & !is.na(labels)
  scores <- scores[keep]
  labels <- labels[keep]
  n_pos <- sum(labels == 1L)
  n_neg <- sum(labels == 0L)
  if (length(scores) < 20L || n_pos < 10L || n_neg < 10L || stats::var(scores) < 1e-12) {
    return(list(auroc = NA_real_, acc = NA_real_, n_positive = n_pos, n_total = length(scores)))
  }

  roc_obj <- suppressMessages(pROC::roc(labels, scores, quiet = TRUE, direction = "<", levels = c(0L, 1L)))
  auroc <- as.numeric(pROC::auc(roc_obj))

  pred <- integer(length(scores))
  top_n <- max(1L, n_pos)
  pred[order(scores, decreasing = TRUE)[seq_len(top_n)]] <- 1L
  acc <- mean(pred == labels)

  list(auroc = auroc, acc = acc, n_positive = n_pos, n_total = length(scores))
}

load_setting <- function(setting) {
  if (identical(setting, "paired")) {
    root <- sctarf_bmmc_root()
    compat <- file.path(root, "04_dual_modality")
    meta <- readRDS(file.path(root, "01_data_qc/BMMC_cell_meta.rds"))
    if (inherits(meta, "data.table")) meta <- as.data.frame(meta)
    if (is.null(rownames(meta)) || identical(rownames(meta), as.character(seq_len(nrow(meta))))) {
      rownames(meta) <- meta$barcode
    }
    target_types <- c("CD14+ Mono", "CD16+ Mono")
    list(
      dataset = "BMMC paired multiome",
      pairing_setting = "paired",
      method_version = "Paired-Fusion",
      method_label = "Paired-Fusion: paired RNA-ATAC fusion",
      source_dir = compat,
      atac_z = readRDS(file.path(compat, "atac_z_block_v2.rds")),
      rna_z = readRDS(file.path(compat, paste0("rna_z_block_v2_", trait_to_run, ".rds"))),
      seed = readRDS(file.path(compat, paste0("seed_v2_", trait_to_run, ".rds"))),
      meta = meta,
      target_types = target_types
    )
  } else {
    root <- sctarf_nonpaired_root()
    compat <- file.path(
      root,
      "03_pseudopair_pipeline/output/v3_sparse_mapping_top5significant/module3_mono_neighborhood_fdr_ge0p80"
    )
    meta <- readRDS(file.path(root, "01_data_prep/output/granja_atac_eval_cell_meta.rds"))
    if (inherits(meta, "data.table")) meta <- as.data.frame(meta)
    if (is.null(rownames(meta)) || identical(rownames(meta), as.character(seq_len(nrow(meta))))) {
      rownames(meta) <- meta$barcode
    }
    target_types <- c("11_CD14.Mono.1", "12_CD14.Mono.2")
    list(
      dataset = "Granja non-paired ATAC + RNA reference",
      pairing_setting = "non_paired",
      method_version = "v3.3",
      method_label = "v3.3: sparse mapping + top5 significant fusion",
      source_dir = compat,
      atac_z = readRDS(file.path(compat, "atac_z_block_v2.rds")),
      rna_z = readRDS(file.path(compat, paste0("rna_z_block_v2_", trait_to_run, ".rds"))),
      seed = readRDS(file.path(compat, paste0("seed_v2_", trait_to_run, ".rds"))),
      meta = meta,
      target_types = target_types
    )
  }
}

benchmark_one <- function(cfg, n_cells, rng_seed) {
  step_log <- data.frame()
  peak_memory_gb <- gc_memory_gb()
  tic <- Sys.time()
  step_start <- Sys.time()

  common <- Reduce(intersect, list(rownames(cfg$atac_z), rownames(cfg$rna_z), names(cfg$seed), rownames(cfg$meta)))
  labels_all <- as.integer(as.character(cfg$meta[common, "cell_type"]) %in% cfg$target_types)
  seed_all <- normalise_seed(cfg$seed, common)

  idx <- sample_indices(labels_all, seed_all, n_cells, rng_seed)
  source_n <- length(common)
  sampled_source <- common[idx]
  bench_cells <- paste0(sampled_source, "__bench", seq_along(sampled_source))

  atac <- cfg$atac_z[sampled_source, , drop = FALSE]
  rna <- cfg$rna_z[sampled_source, , drop = FALSE]
  rownames(atac) <- bench_cells
  rownames(rna) <- bench_cells
  labels <- labels_all[idx]
  seed_vec <- seed_all[idx]
  names(seed_vec) <- bench_cells

  peak_memory_gb <- max(peak_memory_gb, gc_memory_gb(), object_memory_gb(atac, rna, seed_vec), na.rm = TRUE)
  step_log <- rbind(step_log, data.frame(
    pairing_setting = cfg$pairing_setting,
    trait = trait_to_run,
    n_cells = n_cells,
    step = "load_and_sample",
    elapsed_seconds = as.numeric(difftime(Sys.time(), step_start, units = "secs")),
    stringsAsFactors = FALSE
  ))

  step_start <- Sys.time()
  fusion <- cbind(sqrt(atac_w_default) * atac, sqrt(1 - atac_w_default) * rna)
  peak_memory_gb <- max(peak_memory_gb, gc_memory_gb(), object_memory_gb(atac, rna, fusion, seed_vec), na.rm = TRUE)
  step_log <- rbind(step_log, data.frame(
    pairing_setting = cfg$pairing_setting,
    trait = trait_to_run,
    n_cells = n_cells,
    step = "build_fusion_embedding",
    elapsed_seconds = as.numeric(difftime(Sys.time(), step_start, units = "secs")),
    stringsAsFactors = FALSE
  ))

  step_start <- Sys.time()
  knn_res <- build_knn(fusion)
  peak_memory_gb <- max(peak_memory_gb, gc_memory_gb(), object_memory_gb(fusion, knn_res$adj, seed_vec), na.rm = TRUE)
  step_log <- rbind(step_log, data.frame(
    pairing_setting = cfg$pairing_setting,
    trait = trait_to_run,
    n_cells = n_cells,
    step = "build_mutual_knn",
    elapsed_seconds = as.numeric(difftime(Sys.time(), step_start, units = "secs")),
    stringsAsFactors = FALSE
  ))

  step_start <- Sys.time()
  trs_raw <- rwr(knn_res$adj, seed_vec)
  trs <- post_process_trs(trs_raw)
  peak_memory_gb <- max(peak_memory_gb, gc_memory_gb(), object_memory_gb(knn_res$adj, trs_raw, trs), na.rm = TRUE)
  step_log <- rbind(step_log, data.frame(
    pairing_setting = cfg$pairing_setting,
    trait = trait_to_run,
    n_cells = n_cells,
    step = "rwr_trs",
    elapsed_seconds = as.numeric(difftime(Sys.time(), step_start, units = "secs")),
    stringsAsFactors = FALSE
  ))

  step_start <- Sys.time()
  met <- metric_from_scores(trs, labels)
  step_log <- rbind(step_log, data.frame(
    pairing_setting = cfg$pairing_setting,
    trait = trait_to_run,
    n_cells = n_cells,
    step = "evaluate_auroc_acc",
    elapsed_seconds = as.numeric(difftime(Sys.time(), step_start, units = "secs")),
    stringsAsFactors = FALSE
  ))

  total_elapsed <- as.numeric(difftime(Sys.time(), tic, units = "secs"))

  row <- data.frame(
    dataset = cfg$dataset,
    pairing_setting = cfg$pairing_setting,
    method_version = cfg$method_version,
    method_label = cfg$method_label,
    trait = trait_to_run,
    n_cells = n_cells,
    source_cells_available = source_n,
    sampling = ifelse(n_cells > source_n, "with_replacement", "without_replacement"),
    atac_w = atac_w_default,
    top_k = ifelse(cfg$pairing_setting == "non_paired", 5L, NA_integer_),
    runtime_seconds = total_elapsed,
    runtime_minutes = total_elapsed / 60,
    peak_memory_gb = peak_memory_gb,
    auroc = met$auroc,
    acc_top_prevalence = met$acc,
    n_positive = met$n_positive,
    n_total = met$n_total,
    iso_frac = knn_res$iso_frac,
    status = "OK",
    note = ifelse(
      n_cells > source_n,
      "Scale tested by bootstrap resampling existing cells because source data has fewer cells.",
      "Scale tested by subsampling existing cells."
    ),
    source_dir = cfg$source_dir,
    stringsAsFactors = FALSE
  )

  rm(atac, rna, fusion, knn_res, trs_raw, trs)
  gc()
  list(row = row, steps = step_log)
}

run_formal_benchmark <- function() {
  settings <- c("paired", "non_paired")
  all_rows <- list()
  all_steps <- list()

  for (setting in settings) {
    cat("Loading benchmark setting: ", setting, "\n", sep = "")
    cfg <- load_setting(setting)
    for (n_cells in cell_sizes) {
      cat("  Benchmark ", setting, " / ", trait_to_run, " / n=", n_cells, "\n", sep = "")
      res <- tryCatch(
        benchmark_one(cfg, n_cells, rng_seed = 20260528L + n_cells + nchar(setting)),
        error = function(e) {
          list(
            row = data.frame(
              dataset = cfg$dataset,
              pairing_setting = cfg$pairing_setting,
              method_version = cfg$method_version,
              method_label = cfg$method_label,
              trait = trait_to_run,
              n_cells = n_cells,
              source_cells_available = nrow(cfg$atac_z),
              sampling = ifelse(n_cells > nrow(cfg$atac_z), "with_replacement", "without_replacement"),
              atac_w = atac_w_default,
              top_k = ifelse(cfg$pairing_setting == "non_paired", 5L, NA_integer_),
              runtime_seconds = NA_real_,
              runtime_minutes = NA_real_,
              peak_memory_gb = NA_real_,
              auroc = NA_real_,
              acc_top_prevalence = NA_real_,
              n_positive = NA_integer_,
              n_total = NA_integer_,
              iso_frac = NA_real_,
              status = "FAILED",
              note = conditionMessage(e),
              source_dir = cfg$source_dir,
              stringsAsFactors = FALSE
            ),
            steps = data.frame()
          )
        }
      )
      all_rows[[length(all_rows) + 1L]] <- res$row
      if (nrow(res$steps) > 0L) all_steps[[length(all_steps) + 1L]] <- res$steps
    }
    rm(cfg)
    gc()
  }

  table <- bind_rows_fill(all_rows)
  steps <- bind_rows_fill(all_steps)
  write_csv(table, formal_file)
  write_csv(steps, step_file)
  write_csv(table[table$pairing_setting == "paired", , drop = FALSE],
            file.path(out_dir, "efficiency_runtime_memory_table_paired.csv"))
  write_csv(table[table$pairing_setting == "non_paired", , drop = FALSE],
            file.path(out_dir, "efficiency_runtime_memory_table_non_paired_v3_3.csv"))
  table
}

make_plots <- function(formal) {
  plot_dt <- formal[formal$status == "OK" & is.finite(formal$runtime_seconds), , drop = FALSE]
  if (nrow(plot_dt) == 0L) return(invisible(NULL))

  plot_dt$cell_scale <- factor(
    plot_dt$n_cells,
    levels = sort(unique(plot_dt$n_cells)),
    labels = paste0(sort(unique(plot_dt$n_cells)) / 1000, "K")
  )
  plot_dt$setting_label <- ifelse(
    plot_dt$pairing_setting == "paired",
    "Paired",
    "Non-paired v3.3"
  )

  p_runtime <- ggplot(plot_dt, aes(x = n_cells, y = runtime_seconds, color = setting_label)) +
    geom_point(size = 2.4) +
    geom_line(linewidth = 0.7) +
    scale_x_continuous(trans = "log10", breaks = sort(unique(plot_dt$n_cells))) +
    labs(x = "Number of cells", y = "Runtime seconds", color = NULL,
         title = "Runtime scaling of final fusion methods") +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank())
  ggsave(file.path(out_dir, "fig_efficiency_runtime.png"), p_runtime, width = 7, height = 5, dpi = 300)
  ggsave(file.path(out_dir, "fig_efficiency_runtime_scaling.png"), p_runtime, width = 7, height = 5, dpi = 300)

  p_memory <- ggplot(plot_dt, aes(x = n_cells, y = peak_memory_gb, color = setting_label)) +
    geom_point(size = 2.4) +
    geom_line(linewidth = 0.7) +
    scale_x_continuous(trans = "log10", breaks = sort(unique(plot_dt$n_cells))) +
    labs(x = "Number of cells", y = "Estimated peak memory GB", color = NULL,
         title = "Memory scaling of final fusion methods") +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank())
  ggsave(file.path(out_dir, "fig_efficiency_memory.png"), p_memory, width = 7, height = 5, dpi = 300)
  ggsave(file.path(out_dir, "fig_efficiency_memory_scaling.png"), p_memory, width = 7, height = 5, dpi = 300)

  p_auroc <- ggplot(plot_dt, aes(x = n_cells, y = auroc, color = setting_label)) +
    geom_point(size = 2.4) +
    geom_line(linewidth = 0.7) +
    scale_x_continuous(trans = "log10", breaks = sort(unique(plot_dt$n_cells))) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "Number of cells", y = "AUROC", color = NULL,
         title = "AUROC stability across data scales") +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank())
  ggsave(file.path(out_dir, "fig_efficiency_auroc_stability.png"), p_auroc, width = 7, height = 5, dpi = 300)
}

write_template <- function() {
  template <- expand.grid(
    pairing_setting = c("paired", "non_paired"),
    trait = trait_to_run,
    n_cells = cell_sizes,
    stringsAsFactors = FALSE
  )
  template$dataset <- ifelse(
    template$pairing_setting == "paired",
    "BMMC paired multiome",
    "Granja non-paired ATAC + RNA reference"
  )
  template$method_version <- ifelse(template$pairing_setting == "paired", "Paired-Fusion", "v3.3")
  template$top_k <- ifelse(template$pairing_setting == "paired", NA_integer_, 5L)
  template$atac_w <- atac_w_default
  template$runtime_seconds <- NA_real_
  template$peak_memory_gb <- NA_real_
  template$auroc <- NA_real_
  template$acc_top_prevalence <- NA_real_
  template$note <- "Template retained for manual reruns or external machine benchmarks."
  write_csv(template, file.path(out_dir, "manual_efficiency_runs_template.csv"))
}

write_notes <- function(formal) {
  ok <- formal[formal$status == "OK", , drop = FALSE]
  lines <- c(
    "# Efficiency Notes",
    "",
    paste0("Formal benchmark trait: ", trait_to_run, "."),
    paste0("Cell scales: ", paste(cell_sizes, collapse = ", "), "."),
    "",
    "The formal table reports paired and non-paired settings separately.",
    "Paired rows use BMMC paired multiome fusion. Non-paired rows use final method v3.3: sparse mapping plus top5-significant RNA completion.",
    "",
    "Runtime includes subsampling/bootstrap sampling, fused embedding construction, mutual KNN, RWR/TRS, and AUROC/ACC evaluation.",
    "Peak memory is estimated from R heap/object usage during the benchmark, not from an OS-level profiler.",
    "ACC is computed by calling the top N cells positive, where N equals the number of target-positive cells in that benchmark subset.",
    "Rows with sampling = with_replacement use bootstrap resampling because the local source dataset has fewer cells than the requested scale.",
    ""
  )

  if (nrow(ok) > 0L) {
    max_runtime <- aggregate(runtime_seconds ~ pairing_setting, ok, max, na.rm = TRUE)
    max_memory <- aggregate(peak_memory_gb ~ pairing_setting, ok, max, na.rm = TRUE)
    min_auroc <- aggregate(auroc ~ pairing_setting, ok, min, na.rm = TRUE)
    summary <- merge(merge(max_runtime, max_memory, by = "pairing_setting"), min_auroc, by = "pairing_setting")
    names(summary) <- c("pairing_setting", "max_runtime_seconds", "max_peak_memory_gb", "min_auroc")
    lines <- c(lines, "Current benchmark summary:", "", capture.output(print(summary, row.names = FALSE)), "")
  }

  lines <- c(
    lines,
    "Interpretation checklist:",
    "- If 50K/100K rows finish with finite runtime and memory, the method is practically runnable at those scales on this machine.",
    "- Stable AUROC/ACC across scales supports that performance is not an artifact of a small sample.",
    "- Any with_replacement 100K result should be described as a scale stress test, not as a new independent 100K-cell biological dataset."
  )
  writeLines(lines, file.path(out_dir, "efficiency_notes.md"), useBytes = TRUE)
}

copy_to_readable_efficiency_dir <- function() {
  readable_eff_dirs <- Sys.glob(file.path(FINAL_ROOT, "*", "04_*"))
  readable_eff_dirs <- readable_eff_dirs[dir.exists(readable_eff_dirs)]
  if (length(readable_eff_dirs) == 0L) return(invisible(NULL))

  files <- c(
    "efficiency_runtime_memory_table.csv",
    "efficiency_runtime_memory_table_paired.csv",
    "efficiency_runtime_memory_table_non_paired_v3_3.csv",
    "efficiency_step_times_formal_benchmark.csv",
    "fig_efficiency_runtime.png",
    "fig_efficiency_memory.png",
    "fig_efficiency_auroc_stability.png",
    "efficiency_notes.md"
  )
  for (d in readable_eff_dirs) {
    for (f in files) {
      from <- file.path(out_dir, f)
      if (file.exists(from)) file.copy(from, file.path(d, f), overwrite = TRUE)
    }
  }
  invisible(NULL)
}

make_elapsed_demo_outputs()
write_template()

if (!file.exists(formal_file) || force_rerun) {
  formal <- run_formal_benchmark()
} else {
  formal <- read_csv_if_exists(formal_file)
  if (is.null(formal) || nrow(formal) == 0L) formal <- run_formal_benchmark()
}

make_plots(formal)
write_notes(formal)
copy_to_readable_efficiency_dir()

append_manifest(data.frame(
  final_file = c(
    file.path(out_dir, "efficiency_runtime_memory_table.csv"),
    file.path(out_dir, "efficiency_runtime_memory_table_paired.csv"),
    file.path(out_dir, "efficiency_runtime_memory_table_non_paired_v3_3.csv"),
    file.path(out_dir, "efficiency_step_times_formal_benchmark.csv"),
    file.path(out_dir, "fig_efficiency_runtime.png"),
    file.path(out_dir, "fig_efficiency_memory.png"),
    file.path(out_dir, "fig_efficiency_auroc_stability.png"),
    file.path(out_dir, "efficiency_notes.md")
  ),
  source_file = c(
    "BMMC paired module3 outputs and Granja module4 v3.3 outputs",
    file.path(sctarf_bmmc_root(), "04_dual_modality"),
    file.path(sctarf_nonpaired_root(), "03_pseudopair_pipeline/output/v3_sparse_mapping_top5significant"),
    "generated by formal benchmark in 04_build_efficiency_outputs.R",
    "generated by formal benchmark in 04_build_efficiency_outputs.R",
    "generated by formal benchmark in 04_build_efficiency_outputs.R",
    "generated by formal benchmark in 04_build_efficiency_outputs.R",
    "generated by formal benchmark in 04_build_efficiency_outputs.R"
  ),
  script = "scripts/04_build_efficiency_outputs.R",
  dataset = "paired + non_paired",
  trait = trait_to_run,
  method = "Paired-Fusion and non-paired v3.3",
  status = "OK",
  note = "Formal data-scale efficiency benchmark.",
  stringsAsFactors = FALSE
))

cat("Wrote efficiency outputs under: ", out_dir, "\n", sep = "")
