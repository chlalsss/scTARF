# 02: Run SCAVENGE TRS with heterogeneous restart probability gamma_i.
#
# Six conditions (matrix of seed source × gamma policy):
#
#   id  | seed source                  | gamma policy
#   ----|------------------------------|--------------------------------
#   C0  | original (z top 5%)          | constant 0.05  (sanity baseline)
#   C1  | filtered (direction 1, X=0.5)| constant 0.05  (sanity dir1 only)
#   C2  | original                     | Version A (discrete progenitor)
#   C3  | original                     | Version B (continuous reliability)
#   C4  | filtered (direction 1, X=0.5)| Version A
#   C5  | filtered (direction 1, X=0.5)| Version B
#
# C0 / C1 cross-check our custom hetgamma RWR against SCAVENGE's randomWalk_sparse
# (results should match the direction-1 04 output to ~1e-6 in TRS).
# C2-C3 isolate direction 2 alone.
# C4-C5 evaluate direction 1 + direction 2 stacked.
#
# Outputs (in <out_dir>/runs/<condition>/):
#   SCAVENGE_<trait>.tsv.gz            # per-cell TRS
#   trait_summary.csv                  # per-trait run stats (n_seed, iter, etc.)
#
# Top-level outputs (in <out_dir>/):
#   condition_top5_summary.csv         # one row per (condition, trait)
#   condition_run_manifest.json
#
# Parameters:
#   MODULE2_TRAITS                  : default "MONO,PLT,LYM"
#   MODULE2_RW_GAMMA_BASE           : constant gamma for C0/C1 (default 0.05)
#   MODULE2_TOP_PERCENT             : top-x% for enrichment (default 0.05)
#   MODULE2_CAP_QUANTILE            : default 0.95
#   MODULE2_SCALE_PERCENT           : default 0.01
#   MODULE2_SEED_FILTER_PRIMARY     : direction 1 threshold to use (default 0.50)
#   MODULE2_RWR_MAX_ITER            : default 200
#   MODULE2_RWR_EPS                 : default 1e-5
#   MODULE2_ADAPTIVE_GAMMA_CONDITIONS: comma list (default "C0,C1,C2,C3,C4,C5")

source("D:/RStudio/work/module2_adaptive_gamma/_bootstrap.R")

required_pkgs(c("data.table", "Matrix", "jsonlite"))
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(jsonlite)
})

paths <- adaptive_gamma_paths()
ensure_adaptive_gamma_dirs(paths)

traits <- parse_csv_env("MODULE2_TRAITS", "MONO,PLT,LYM")
gamma_base <- env_numeric("MODULE2_RW_GAMMA_BASE", 0.05)
top_pct <- env_numeric("MODULE2_TOP_PERCENT", 0.05)
cap_q <- env_numeric("MODULE2_CAP_QUANTILE", 0.95)
scale_pct <- env_numeric("MODULE2_SCALE_PERCENT", 0.01)
direction1_threshold <- env_numeric("MODULE2_SEED_FILTER_PRIMARY", 0.50)
max_iter <- env_integer("MODULE2_RWR_MAX_ITER", 200L)
eps <- env_numeric("MODULE2_RWR_EPS", 1e-5)

requested_conditions <- parse_csv_env("MODULE2_ADAPTIVE_GAMMA_CONDITIONS",
                                      "C0,C1,C2,C3,C4,C5")
all_conditions <- c("C0", "C1", "C2", "C3", "C4", "C5")
unknown <- setdiff(requested_conditions, all_conditions)
if (length(unknown) > 0) {
  stop("Unknown conditions: ", paste(unknown, collapse = ", "), call. = FALSE)
}
conditions <- intersect(all_conditions, requested_conditions)

cat("========== Module 2 adaptive gamma: 02 run SCAVENGE ==========\n")
cat("  traits     :", paste(traits, collapse = ", "), "\n")
cat("  conditions :", paste(conditions, collapse = ", "), "\n")
cat("  gamma_base :", gamma_base, "\n")
cat("  cap_q      :", cap_q, "\n")
cat("  scale_pct  :", scale_pct, "\n")
cat("  top_pct    :", top_pct, "\n")
cat("  dir1 thr   :", direction1_threshold, "\n")
cat("  RWR eps    :", eps, "  max_iter:", max_iter, "\n")

# Load shared inputs ----------------------------------------------------

cat("\n[1/4] Load shared inputs\n")
stopifnot(file.exists(paths$original_knn))
stopifnot(file.exists(paths$seed_signal_by_cell))
original_knn <- readRDS(paths$original_knn)
adj_full <- as_binary_graph(original_knn)
cells_full <- rownames(adj_full)

# Remove isolated cells before RWR (matches direction 1 04 behavior)
deg <- Matrix::rowSums(adj_full)
graph_cells <- cells_full[deg > 0]
isolated_cells <- cells_full[deg == 0]
if (length(isolated_cells) > 0) {
  cat("  isolated cells removed before RWR:", length(isolated_cells), "\n")
}
adj <- adj_full[graph_cells, graph_cells]

seed_signal_dt <- fread(paths$seed_signal_by_cell)
seed_signal_dt[, target_type := as.logical(target_type)]
seed_signal_dt[, seed := as.logical(seed)]

# Gamma vectors (only loaded for relevant conditions)
load_gamma_versionA <- function() {
  fp <- file.path(paths$out_dir, "gamma_versionA.csv")
  if (!file.exists(fp)) stop("Missing gamma_versionA.csv; run 01 first.", call. = FALSE)
  dt <- fread(fp)
  setNames(as.numeric(dt$gamma), dt$cell)
}
load_gamma_versionB <- function() {
  fp <- file.path(paths$out_dir, "gamma_versionB.csv")
  if (!file.exists(fp)) stop("Missing gamma_versionB.csv; run 01 first.", call. = FALSE)
  dt <- fread(fp)
  setNames(as.numeric(dt$gamma), dt$cell)
}

# Seed loaders ----------------------------------------------------------

load_seed_original <- function(trait_name) {
  signal_tr <- seed_signal_dt[trait == trait_name]
  if (nrow(signal_tr) == 0) stop("No seed_signal rows for trait ", trait_name, call. = FALSE)
  signal_tr[seed == TRUE, cell]
}

load_seed_filtered <- function(trait, threshold) {
  fp <- sprintf(paths$direction1_filtered_seed_template,
                trait, seed_filter_label(threshold))
  if (!file.exists(fp)) {
    stop("Missing direction-1 filtered seed file: ", fp,
         "\nRun module2_seed_rna_filter/02_filter_seed_by_rna.R first.",
         call. = FALSE)
  }
  as.character(fread(fp)$cell)
}

# Per-condition spec ----------------------------------------------------

condition_spec <- function(cid) {
  switch(cid,
    "C0" = list(label = "C0_orig_seed_const_gamma",
                seed_kind = "original", gamma_kind = "const"),
    "C1" = list(label = "C1_filt_seed_const_gamma",
                seed_kind = "filtered", gamma_kind = "const"),
    "C2" = list(label = "C2_orig_seed_versionA",
                seed_kind = "original", gamma_kind = "A"),
    "C3" = list(label = "C3_orig_seed_versionB",
                seed_kind = "original", gamma_kind = "B"),
    "C4" = list(label = "C4_filt_seed_versionA",
                seed_kind = "filtered", gamma_kind = "A"),
    "C5" = list(label = "C5_filt_seed_versionB",
                seed_kind = "filtered", gamma_kind = "B"),
    stop("Unknown condition: ", cid, call. = FALSE)
  )
}

# Load gamma vectors lazily
gamma_A <- NULL
gamma_B <- NULL
need_A <- any(c("C2", "C4") %in% conditions)
need_B <- any(c("C3", "C5") %in% conditions)
if (need_A) {
  gamma_A_full <- load_gamma_versionA()
  gamma_A <- gamma_A_full[graph_cells]
  if (any(is.na(gamma_A))) {
    stop("gamma_A has NA after subsetting to graph_cells", call. = FALSE)
  }
}
if (need_B) {
  gamma_B_full <- load_gamma_versionB()
  gamma_B <- gamma_B_full[graph_cells]
  if (any(is.na(gamma_B))) {
    stop("gamma_B has NA after subsetting to graph_cells", call. = FALSE)
  }
}

# Run -------------------------------------------------------------------

cat("\n[2/4] Run conditions\n")
all_top5 <- list()
all_trait_summary <- list()
runs_root <- file.path(paths$out_dir, "runs")
dir.create(runs_root, showWarnings = FALSE, recursive = TRUE)

for (cid in conditions) {
  spec <- condition_spec(cid)
  cat("\n  ===== Condition", cid, "(", spec$label, ") =====\n")
  cond_dir <- file.path(runs_root, spec$label)
  dir.create(cond_dir, showWarnings = FALSE, recursive = TRUE)

  trait_summary_list <- list()
  for (tr in traits) {
    cat("    [", cid, "/", tr, "]\n", sep = "")
    trait_name <- tr
    signal_tr <- seed_signal_dt[trait == trait_name]
    if (nrow(signal_tr) == 0) next
    z_score <- setNames(as.numeric(signal_tr$z_score), signal_tr$cell)
    target_named <- setNames(as.logical(signal_tr$target_type), signal_tr$cell)

    if (spec$seed_kind == "original") {
      seed_total <- load_seed_original(tr)
    } else {
      seed_total <- load_seed_filtered(tr, direction1_threshold)
    }
    seed_in_graph <- intersect(seed_total, graph_cells)
    if (length(seed_in_graph) == 0) {
      warning("No seed in graph for ", cid, "/", tr); next
    }

    if (spec$gamma_kind == "const") {
      gamma_in <- gamma_base
    } else if (spec$gamma_kind == "A") {
      gamma_in <- gamma_A
    } else {
      gamma_in <- gamma_B
    }

    np_score <- random_walk_hetgamma(
      intM = adj,
      queryCells = seed_in_graph,
      gamma = gamma_in,
      eps = eps,
      max_iter = max_iter,
      verbose = FALSE
    )

    # Drop nodes with zero stationary probability (consistent with direction-1 04)
    omit <- np_score == 0
    np_use <- np_score[!omit]
    keep_cells <- names(np_use)

    scale_factor <- compute_scale_factor(z_score, scale_pct = scale_pct)
    TRS_unscaled <- cap_and_scale(np_use, q_ceiling = cap_q)
    TRS <- TRS_unscaled * scale_factor

    out_df <- data.table(
      cell = keep_cells,
      z_score = as.numeric(z_score[keep_cells]),
      seed = keep_cells %in% seed_in_graph,
      target_type = as.logical(target_named[keep_cells]),
      np_score = as.numeric(np_use),
      TRS = as.numeric(TRS)
    )
    out_file <- file.path(cond_dir, sprintf("SCAVENGE_%s.tsv.gz", tr))
    fwrite(out_df, out_file, sep = "\t")

    # top 5%
    setorder(out_df, -TRS)
    n_top <- max(1L, ceiling(nrow(out_df) * top_pct))
    top_dt <- out_df[seq_len(n_top)]
    top_target_frac <- mean(top_dt$target_type, na.rm = TRUE)

    iter_used <- attr(np_score, "n_iter")
    final_delta <- attr(np_score, "final_delta")

    all_top5[[paste(cid, tr, sep = "_")]] <- data.table(
      condition = cid,
      condition_label = spec$label,
      seed_kind = spec$seed_kind,
      gamma_kind = spec$gamma_kind,
      trait = tr,
      n_seed_total = length(seed_total),
      n_seed_in_graph = length(seed_in_graph),
      n_cells_in_TRS = nrow(out_df),
      n_top = n_top,
      top_target_fraction = top_target_frac,
      mean_TRS = mean(TRS),
      median_TRS = median(TRS),
      rwr_iters = iter_used,
      rwr_final_delta = final_delta,
      output_file = out_file,
      target_definition = target_definition(tr)
    )
    trait_summary_list[[tr]] <- all_top5[[paste(cid, tr, sep = "_")]]

    cat(sprintf("      seed=%d  iter=%d  delta=%.1e  top5=%.2f%%\n",
                length(seed_in_graph), iter_used, final_delta,
                top_target_frac * 100))
  }
  cond_summary_dt <- rbindlist(trait_summary_list, use.names = TRUE, fill = TRUE)
  cond_summary_path <- file.path(cond_dir, "trait_summary.csv")
  fwrite(cond_summary_dt, cond_summary_path)
  all_trait_summary[[cid]] <- cond_summary_dt
}

# Top-level summary -----------------------------------------------------

cat("\n[3/4] Write top-level summary\n")
top5_dt <- rbindlist(all_top5, use.names = TRUE, fill = TRUE)
setorder(top5_dt, trait, condition)
top5_path <- file.path(paths$out_dir, "condition_top5_summary.csv")
fwrite(top5_dt, top5_path)
cat("\n  condition_top5_summary.csv:\n")
print(top5_dt[, .(condition, trait, n_seed_in_graph,
                  top_target_fraction, rwr_iters, rwr_final_delta)])

# Manifest --------------------------------------------------------------

cat("\n[4/4] Write manifest\n")
manifest <- list(
  script = "02_run_scavenge_adaptive_gamma.R",
  conditions_run = conditions,
  parameters = list(
    traits = traits,
    gamma_base = gamma_base,
    cap_quantile = cap_q,
    scale_percent = scale_pct,
    top_percent = top_pct,
    direction1_threshold = direction1_threshold,
    rwr_max_iter = max_iter,
    rwr_eps = eps
  ),
  inputs = list(
    original_knn = paths$original_knn,
    seed_signal = paths$seed_signal_by_cell,
    direction1_filtered_seed_dir = paths$seed_filter_dir,
    gamma_versionA = file.path(paths$out_dir, "gamma_versionA.csv"),
    gamma_versionB = file.path(paths$out_dir, "gamma_versionB.csv")
  ),
  outputs = list(
    runs_dir = runs_root,
    top5_summary = top5_path
  ),
  notes = c(
    "C0/C1 are sanity baselines using constant gamma; their TRS should match the direction-1 04 output to within numerical precision.",
    "C2/C3 isolate direction 2 alone (no seed filter).",
    "C4/C5 evaluate direction 1 + direction 2 stacked.",
    "Acceptance: C2 or C3 vs original_seed const baseline (>= 1pp on >=2 traits) for direction 2 alone; C4 or C5 vs C1 for stacking value."
  ),
  created_at = as.character(Sys.time())
)
manifest_out <- file.path(paths$out_dir, "condition_run_manifest.json")
write_json(manifest, manifest_out, pretty = TRUE, auto_unbox = TRUE)

cat("\n========== Done ==========\n")
cat("  manifest:", manifest_out, "\n")
