# module4 v3 sparse: RNA neighborhood-FDR seed candidates from sparse mappings
# -----------------------------------------------------------------------------
# Purpose:
#   Given ATAC -> CITE RNA mappings and RNA neighborhood-level FDR results,
#   compute ATAC seed candidates by summing each ATAC cell's top5 RNA-neighbor
#   weight that lands in significant RNA neighborhoods.
#
# Default seed rule used downstream:
#   candidate_group_weight_ge0p80:
#     group_sig_source_weight >= 0.80
#
# Method boundary:
#   - Cell type is not used for mapping, RNA neighborhood FDR, or seed selection.
#   - Cell type is used only for post hoc diagnostics and AUROC target proxies.
#
# Run through RUN_V3_SPARSE_OTHER_TRAITS.R, or directly after setting MAP_DIR,
# TRAITS_TO_SCREEN, RNA_GROUP_TEST_FILES, and OUT_DIR.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

if (!exists("ROOT")) ROOT <- Sys.getenv("SCTARF_NONPAIRED_ROOT", file.path("outputs", "nonpaired"))

if (!exists("MAP_DIR")) {
  MAP_DIR <- file.path(ROOT, "02_trait_mapping/output/v3_sparse_other_traits")
}
if (!exists("RNA_GROUP_FILE")) {
  RNA_GROUP_FILE <- file.path(
    ROOT,
    "04_evaluation/module4_v2_diagnostics/neighborhood_fdr/rna_unsupervised_group_membership.csv.gz"
  )
}
if (!exists("RNA_GROUP_TEST_FILES")) {
  RNA_GROUP_TEST_FILES <- c(
    file.path(ROOT, "04_evaluation/module4_v2_diagnostics/neighborhood_fdr/neighborhood_fdr_group_tests.csv"),
    file.path(ROOT, "01_data_prep_external_rna/output/trait_screening/cite_trait_screening_neighborhood_fdr_group_tests.csv")
  )
}
if (!exists("ATAC_META_FILE")) {
  ATAC_META_FILE <- file.path(ROOT, "04_evaluation/inputs/granja_cell_meta.rds")
}
if (!exists("OUT_DIR")) {
  OUT_DIR <- file.path(ROOT, "02_trait_mapping/output/v3_sparse_other_traits/neighborhood_seed_diagnostics")
}
if (!exists("TRAITS_TO_SCREEN")) {
  TRAITS_TO_SCREEN <- c("PLT", "LYM", "MCH", "RDW", "RBC", "EOS", "WBC")
}
if (!exists("GROUP_TEST_TRAIT_MAP")) {
  GROUP_TEST_TRAIT_MAP <- c(EOS = "EOSINOPHIL", WBC = "WHITE")
}
if (!exists("GROUP_FDR_THRESHOLD")) GROUP_FDR_THRESHOLD <- 0.10
if (!exists("TOP_N_ATAC")) TOP_N_ATAC <- 350L

if (!exists("POSTHOC_TARGETS_LIST")) {
  red_proxy <- c("02_Early.Eryth", "03_Late.Eryth", "04_Early.Baso")
  lymph_proxy <- c(
    "16_Pre.B", "17_B", "18_Plasma",
    "19_CD8.N", "20_CD4.N1", "21_CD4.N2", "22_CD4.M",
    "23_CD8.EM", "24_CD8.CM", "25_NK"
  )
  granulocyte_proxy <- c("04_Early.Baso", "07_GMP", "08_GMP.Neut")
  leukocyte_proxy <- c(
    "04_Early.Baso",
    "06_CLP.1", "07_GMP", "08_GMP.Neut", "09_pDC", "10_cDC",
    "11_CD14.Mono.1", "12_CD14.Mono.2",
    "15_CLP.2", "16_Pre.B", "17_B", "18_Plasma",
    "19_CD8.N", "20_CD4.N1", "21_CD4.N2", "22_CD4.M",
    "23_CD8.EM", "24_CD8.CM", "25_NK"
  )
  POSTHOC_TARGETS_LIST <- list(
    PLT = red_proxy,
    LYM = lymph_proxy,
    MCH = red_proxy,
    RDW = red_proxy,
    RBC = red_proxy,
    EOS = granulocyte_proxy,
    WBC = leukocyte_proxy
  )
}
if (!exists("POSTHOC_TARGET_NAMES")) {
  POSTHOC_TARGET_NAMES <- c(
    PLT = "red_lineage_proxy",
    LYM = "lymph_proxy",
    MCH = "red_lineage_proxy",
    RDW = "red_lineage_proxy",
    RBC = "red_lineage_proxy",
    EOS = "granulocyte_like_proxy",
    WBC = "leukocyte_non_erythroid_proxy"
  )
}

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

stop_missing <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0L) {
    stop("Missing required file(s):\n  ", paste(missing, collapse = "\n  "),
         call. = FALSE)
  }
}

top_table_label <- function(x, n = 10L) {
  x <- x[!is.na(x)]
  tab <- sort(table(x), decreasing = TRUE)
  if (length(tab) == 0L) return("")
  paste(paste0(names(head(tab, n)), ":", as.integer(head(tab, n))), collapse = "; ")
}

lineage_atac <- function(cell_type) {
  out <- rep("other", length(cell_type))
  out[cell_type %in% c("02_Early.Eryth", "03_Late.Eryth", "04_Early.Baso")] <- "red_lineage"
  out[cell_type %in% c("11_CD14.Mono.1", "12_CD14.Mono.2")] <- "mono"
  out[cell_type %in% c("09_pDC", "10_cDC", "08_GMP.Neut")] <- "myeloid_other"
  out[cell_type %in% c("01_HSC", "05_CMP.LMPP", "06_CLP.1", "07_GMP", "15_CLP.2")] <- "progenitor"
  out[cell_type %in% c(
    "16_Pre.B", "17_B", "18_Plasma",
    "19_CD8.N", "20_CD4.N1", "21_CD4.N2", "22_CD4.M",
    "23_CD8.EM", "24_CD8.CM", "25_NK"
  )] <- "lymph"
  out[grepl("Unk", cell_type)] <- "unknown"
  out
}

safe_auc <- function(score, label) {
  if (!requireNamespace("pROC", quietly = TRUE)) return(NA_real_)
  keep <- is.finite(score) & !is.na(label)
  score <- score[keep]
  label <- as.integer(label[keep])
  if (sum(label == 1L) < 10L || sum(label == 0L) < 10L) return(NA_real_)
  if (length(unique(score)) < 2L) return(NA_real_)
  as.numeric(pROC::auc(pROC::roc(label, score, quiet = TRUE,
                                 direction = "<", levels = c(0L, 1L))))
}

safe_div <- function(num, den) {
  if (!is.finite(den) || den <= 0) return(NA_real_)
  num / den
}

load_atac_meta <- function(path) {
  meta <- readRDS(path)
  if (inherits(meta, "data.table")) meta <- as.data.frame(meta)
  if (is.null(rownames(meta)) ||
      identical(rownames(meta), as.character(seq_len(nrow(meta))))) {
    if ("barcode" %in% colnames(meta)) rownames(meta) <- meta$barcode
  }
  if (!"cell_type" %in% colnames(meta)) {
    stop("ATAC meta missing cell_type column.", call. = FALSE)
  }
  meta
}

read_group_tests <- function(paths) {
  stop_missing(paths)
  rows <- lapply(paths, function(path) {
    dt <- fread(path)
    needed <- c("trait", "group_id", "qval", "observed_mean_score")
    missing_cols <- setdiff(needed, names(dt))
    if (length(missing_cols) > 0L) {
      stop("RNA group test file missing columns in ", path, ": ",
           paste(missing_cols, collapse = ", "), call. = FALSE)
    }
    dt[, source_file := basename(path)]
    dt
  })
  rbindlist(rows, use.names = TRUE, fill = TRUE)
}

cat(">>> Sparse mapping RNA-neighborhood seed candidates\n\n")
cat("MAP_DIR             : ", MAP_DIR, "\n", sep = "")
cat("RNA_GROUP_FILE      : ", RNA_GROUP_FILE, "\n", sep = "")
cat("OUT_DIR             : ", OUT_DIR, "\n", sep = "")
cat("TRAITS_TO_SCREEN    : ", paste(TRAITS_TO_SCREEN, collapse = ", "), "\n", sep = "")
cat("GROUP_FDR_THRESHOLD : ", GROUP_FDR_THRESHOLD, "\n\n", sep = "")

stop_missing(c(RNA_GROUP_FILE, ATAC_META_FILE))
stop_missing(file.path(MAP_DIR, paste0("granja_cite_mapping_", TRAITS_TO_SCREEN, ".rds")))

rna_group <- fread(RNA_GROUP_FILE)
if (!all(c("cell", "group_id") %in% names(rna_group))) {
  stop("RNA group membership file must contain cell and group_id.", call. = FALSE)
}
group_lookup <- rna_group$group_id
names(group_lookup) <- rna_group$cell

group_tests <- read_group_tests(RNA_GROUP_TEST_FILES)
atac_meta <- load_atac_meta(ATAC_META_FILE)

summary_rows <- list()
cell_rows <- list()
by_ct_rows <- list()
source_rows <- list()

for (trait_name in TRAITS_TO_SCREEN) {
  cat("=== ", trait_name, " ===\n", sep = "")
  map_file <- file.path(MAP_DIR, paste0("granja_cite_mapping_", trait_name, ".rds"))
  mapping <- readRDS(map_file)

  group_trait_name <- if (trait_name %in% names(GROUP_TEST_TRAIT_MAP)) {
    unname(GROUP_TEST_TRAIT_MAP[[trait_name]])
  } else {
    trait_name
  }

  sig_groups <- group_tests[
    trait == group_trait_name &
      qval <= GROUP_FDR_THRESHOLD &
      observed_mean_score > 0,
    unique(group_id)
  ]
  if (length(sig_groups) == 0L) {
    warning("No significant RNA groups for trait ", trait_name)
  }

  nn_group <- matrix(
    group_lookup[as.vector(mapping$rna_neighbor_barcode)],
    nrow = nrow(mapping$rna_neighbor_barcode),
    ncol = ncol(mapping$rna_neighbor_barcode)
  )
  if (anyNA(nn_group)) {
    stop("Some mapped RNA neighbor barcodes are absent from RNA_GROUP_FILE for ",
         trait_name, call. = FALSE)
  }

  is_sig_neighbor <- nn_group %in% sig_groups
  dim(is_sig_neighbor) <- dim(nn_group)
  sig_weight <- rowSums(mapping$neighbor_weight * is_sig_neighbor, na.rm = TRUE)
  names(sig_weight) <- mapping$atac_barcodes

  atac_ct <- as.character(atac_meta[mapping$atac_barcodes, "cell_type"])
  if (anyNA(atac_ct)) {
    stop("Some ATAC barcodes are missing from metadata for ", trait_name,
         call. = FALSE)
  }
  atac_lineage <- lineage_atac(atac_ct)

  target_cells <- POSTHOC_TARGETS_LIST[[trait_name]]
  target_name <- if (trait_name %in% names(POSTHOC_TARGET_NAMES)) {
    unname(POSTHOC_TARGET_NAMES[[trait_name]])
  } else {
    "posthoc_target"
  }
  is_target <- if (is.null(target_cells)) {
    rep(NA, length(atac_ct))
  } else {
    atac_ct %in% target_cells
  }
  target_rate <- mean(is_target == TRUE, na.rm = TRUE)

  candidate_gt0 <- sig_weight > 0
  candidate_ge0p50 <- sig_weight >= 0.50
  candidate_ge0p80 <- sig_weight >= 0.80
  candidate_topn <- rep(FALSE, length(sig_weight))
  positive <- which(sig_weight > 0 & is.finite(sig_weight))
  if (length(positive) > 0L) {
    ord <- positive[order(sig_weight[positive], decreasing = TRUE)]
    candidate_topn[head(ord, min(TOP_N_ATAC, length(ord)))] <- TRUE
  }

  strategies <- list(
    group_weight_gt0 = candidate_gt0,
    group_weight_ge0p50 = candidate_ge0p50,
    group_weight_ge0p80 = candidate_ge0p80,
    group_top350_positive = candidate_topn
  )

  score_auc <- safe_auc(sig_weight, is_target)

  for (strategy in names(strategies)) {
    candidate <- strategies[[strategy]]
    n_candidate <- sum(candidate)
    n_target_candidate <- sum(candidate & is_target, na.rm = TRUE)
    precision <- safe_div(n_target_candidate, n_candidate)
    recall <- safe_div(n_target_candidate, sum(is_target == TRUE, na.rm = TRUE))
    enrichment <- precision / max(target_rate, 1e-12)

    summary_rows[[paste(trait_name, strategy, sep = "_")]] <- data.table(
      trait = trait_name,
      group_test_trait = group_trait_name,
      group_fdr_threshold = GROUP_FDR_THRESHOLD,
      strategy = strategy,
      n_sig_rna_groups = length(sig_groups),
      n_atac_cells = length(sig_weight),
      posthoc_target_name = target_name,
      n_target_atac_cells_posthoc = sum(is_target == TRUE, na.rm = TRUE),
      target_atac_rate_posthoc = target_rate,
      sig_group_weight_auc_posthoc = score_auc,
      n_candidate_atac = n_candidate,
      candidate_target_precision_posthoc = precision,
      candidate_target_recall_posthoc = recall,
      candidate_target_enrichment_posthoc = enrichment,
      candidate_weight_min = if (n_candidate > 0L) min(sig_weight[candidate]) else NA_real_,
      candidate_weight_median = if (n_candidate > 0L) median(sig_weight[candidate]) else NA_real_,
      candidate_weight_max = if (n_candidate > 0L) max(sig_weight[candidate]) else NA_real_,
      candidate_top_atac_celltypes_posthoc = top_table_label(atac_ct[candidate], 12L),
      candidate_top_atac_lineages_posthoc = top_table_label(atac_lineage[candidate], 8L),
      no_celltype_used_to_select = TRUE
    )

    by_ct <- data.table(
      trait = trait_name,
      strategy = strategy,
      atac_cell_type = atac_ct,
      atac_lineage_posthoc = atac_lineage,
      is_target_posthoc = is_target,
      is_candidate = candidate,
      sig_group_source_weight = sig_weight
    )[, .(
      n_cells = .N,
      n_candidate = sum(is_candidate),
      candidate_rate = mean(is_candidate),
      mean_sig_group_source_weight = mean(sig_group_source_weight, na.rm = TRUE),
      is_target_posthoc = unique(is_target_posthoc)[1]
    ), by = .(trait, strategy, atac_cell_type, atac_lineage_posthoc)][
      order(trait, strategy, -candidate_rate, -n_candidate)
    ]
    by_ct_rows[[paste(trait_name, strategy, sep = "_")]] <- by_ct
  }

  cell_rows[[trait_name]] <- data.table(
    trait = trait_name,
    group_test_trait = group_trait_name,
    atac_barcode = mapping$atac_barcodes,
    atac_cell_type_posthoc = atac_ct,
    atac_lineage_posthoc = atac_lineage,
    posthoc_target_name = target_name,
    is_target_posthoc = is_target,
    sig_group_source_weight = sig_weight,
    candidate_group_weight_gt0 = candidate_gt0,
    candidate_group_weight_ge0p50 = candidate_ge0p50,
    candidate_group_weight_ge0p80 = candidate_ge0p80,
    candidate_group_top350_positive = candidate_topn,
    no_celltype_used_to_select = TRUE
  )

  edge_dt <- data.table(
    atac_barcode = rep(mapping$atac_barcodes, times = ncol(mapping$rna_neighbor_barcode)),
    rna_barcode = as.vector(mapping$rna_neighbor_barcode),
    neighbor_weight = as.vector(mapping$neighbor_weight),
    rna_group_id = as.vector(nn_group),
    is_sig_rna_group = as.vector(is_sig_neighbor)
  )
  source_rows[[trait_name]] <- edge_dt[is_sig_rna_group == TRUE, .(
    total_weight = sum(neighbor_weight, na.rm = TRUE),
    n_edges = .N
  ), by = rna_group_id][order(-total_weight)][, trait := trait_name]

  cat("  sig RNA groups=", length(sig_groups),
      "; ge0p80 seeds=", sum(candidate_ge0p80),
      "; top350 positive seeds=", sum(candidate_topn),
      "\n", sep = "")
}

summary_dt <- rbindlist(summary_rows, use.names = TRUE, fill = TRUE)
cell_dt <- rbindlist(cell_rows, use.names = TRUE, fill = TRUE)
by_ct_dt <- rbindlist(by_ct_rows, use.names = TRUE, fill = TRUE)
source_dt <- rbindlist(source_rows, use.names = TRUE, fill = TRUE)

fwrite(summary_dt, file.path(OUT_DIR, "sparse_neighborhood_seed_summary.csv"))
fwrite(cell_dt, file.path(OUT_DIR, "sparse_neighborhood_seed_cell_level.csv.gz"))
fwrite(by_ct_dt, file.path(OUT_DIR, "sparse_neighborhood_seed_by_celltype_posthoc.csv"))
fwrite(source_dt, file.path(OUT_DIR, "sparse_neighborhood_seed_sig_rna_group_source_weights.csv"))

cat("\n=== Sparse neighborhood seed candidates finished ===\n")
print(summary_dt[
  strategy %in% c("group_weight_ge0p80", "group_top350_positive"),
  .(
    trait, strategy, n_sig_rna_groups, sig_group_weight_auc_posthoc,
    posthoc_target_name, n_candidate_atac,
    candidate_target_precision_posthoc,
    candidate_target_enrichment_posthoc,
    candidate_top_atac_celltypes_posthoc
  )
])
cat("\nOutputs are in: ", OUT_DIR, "\n", sep = "")
