# module4 trait screening: ATAC candidate diagnostics after mapping
# -----------------------------------------------------------------------------
# Purpose:
#   For MCH/RDW/RBC, combine:
#     1. trait-aware ATAC -> CITE RNA mapping
#     2. RNA neighborhood-level scDRS FDR results
#
#   Then ask: which ATAC cells receive mapping weight from significant RNA
#   neighborhoods, and are those ATAC cells enriched for the expected red-lineage
#   target post hoc?
#
# Method boundary:
#   - celltype is NOT used in mapping.
#   - celltype is NOT used in RNA neighborhood FDR.
#   - celltype is NOT used to define candidate ATAC cells.
#   - celltype is used only for post hoc diagnostics.
#
# Run in RStudio:
#   setwd("D:/RStudio/work/module4_traitaware_pseudopair")
#   source("02_trait_mapping/04_screen_trait_mapping_atac_candidates.R")
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

ROOT <- "D:/RStudio/work/module4_traitaware_pseudopair"
if (!exists("MAP_DIR")) MAP_DIR <- file.path(ROOT, "02_trait_mapping/output/trait_screening")
if (!exists("RNA_GROUP_FILE")) {
  RNA_GROUP_FILE <- file.path(
    ROOT,
    "04_evaluation/module4_v2_diagnostics/neighborhood_fdr/rna_unsupervised_group_membership.csv.gz"
  )
}
if (!exists("RNA_GROUP_TEST_FILE")) {
  RNA_GROUP_TEST_FILE <- file.path(
    ROOT,
    "01_data_prep_external_rna/output/trait_screening/cite_trait_screening_neighborhood_fdr_group_tests.csv"
  )
}
if (!exists("ATAC_META_FILE")) ATAC_META_FILE <- file.path(ROOT, "04_evaluation/inputs/granja_cell_meta.rds")
if (!exists("OUT_DIR")) OUT_DIR <- file.path(ROOT, "02_trait_mapping/output/trait_screening/diagnostics")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

if (!exists("TRAITS_TO_SCREEN")) TRAITS_TO_SCREEN <- c("MCH", "RDW", "RBC")
if (!exists("GROUP_FDR_THRESHOLD")) GROUP_FDR_THRESHOLD <- 0.10
if (!exists("TOP_N_ATAC")) TOP_N_ATAC <- 350L

ATAC_RED_TARGETS <- c("02_Early.Eryth", "03_Late.Eryth", "04_Early.Baso")
if (!exists("POSTHOC_TARGETS")) POSTHOC_TARGETS <- ATAC_RED_TARGETS
if (!exists("POSTHOC_TARGET_NAME")) POSTHOC_TARGET_NAME <- "red_lineage_proxy"
if (!exists("GROUP_TEST_TRAIT_MAP")) GROUP_TEST_TRAIT_MAP <- setNames(TRAITS_TO_SCREEN, TRAITS_TO_SCREEN)

lineage_atac <- function(cell_type) {
  out <- rep("other", length(cell_type))
  out[cell_type %in% ATAC_RED_TARGETS] <- "red_lineage_posthoc"
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

top_table_label <- function(x, n = 10L) {
  x <- x[!is.na(x)]
  tab <- sort(table(x), decreasing = TRUE)
  if (length(tab) == 0L) return("")
  paste(paste0(names(head(tab, n)), ":", as.integer(head(tab, n))), collapse = "; ")
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

stop_missing <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0L) {
    stop("Missing required file(s):\n  ", paste(missing, collapse = "\n  "),
         call. = FALSE)
  }
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

cat(">>> Screen ATAC candidates from trait-screening mappings\n\n")
cat("Traits: ", paste(TRAITS_TO_SCREEN, collapse = ", "), "\n", sep = "")
cat("GROUP_FDR_THRESHOLD = ", GROUP_FDR_THRESHOLD, "\n", sep = "")
cat("TOP_N_ATAC = ", TOP_N_ATAC, "\n\n", sep = "")
cat("Post hoc target proxy: ", POSTHOC_TARGET_NAME, "\n\n", sep = "")

stop_missing(c(RNA_GROUP_FILE, RNA_GROUP_TEST_FILE, ATAC_META_FILE))
stop_missing(file.path(MAP_DIR, paste0("granja_cite_mapping_", TRAITS_TO_SCREEN, ".rds")))

rna_group <- fread(RNA_GROUP_FILE)
if (!all(c("cell", "group_id") %in% names(rna_group))) {
  stop("RNA group membership file must contain cell and group_id.", call. = FALSE)
}

group_tests <- fread(RNA_GROUP_TEST_FILE)
if (!all(c("trait", "group_id", "qval", "observed_mean_score") %in% names(group_tests))) {
  stop("RNA group test file missing required columns.", call. = FALSE)
}

atac_meta <- load_atac_meta(ATAC_META_FILE)

summary_rows <- list()
cell_rows <- list()
by_ct_rows <- list()
source_rows <- list()

for (trait_name in TRAITS_TO_SCREEN) {
  cat("=== ", trait_name, " ===\n", sep = "")
  map_file <- file.path(MAP_DIR, paste0("granja_cite_mapping_", trait_name, ".rds"))
  m <- readRDS(map_file)
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

  group_lookup <- rna_group$group_id
  names(group_lookup) <- rna_group$cell

  nn_group <- matrix(
    group_lookup[as.vector(m$rna_neighbor_barcode)],
    nrow = nrow(m$rna_neighbor_barcode),
    ncol = ncol(m$rna_neighbor_barcode)
  )
  is_sig_neighbor <- nn_group %in% sig_groups
  dim(is_sig_neighbor) <- dim(nn_group)

  sig_weight <- rowSums(m$neighbor_weight * is_sig_neighbor, na.rm = TRUE)
  names(sig_weight) <- m$atac_barcodes

  atac_ct <- atac_meta[m$atac_barcodes, "cell_type"]
  if (anyNA(atac_ct)) stop("Some ATAC barcodes are missing from metadata.", call. = FALSE)
  atac_lineage <- lineage_atac(atac_ct)
  is_target <- atac_ct %in% POSTHOC_TARGETS
  target_rate <- mean(is_target)

  candidate_gt0 <- sig_weight > 0
  candidate_ge0p50 <- sig_weight >= 0.50
  candidate_ge0p80 <- sig_weight >= 0.80
  top_n <- min(TOP_N_ATAC, length(sig_weight))
  ord <- order(sig_weight, decreasing = TRUE)
  candidate_topn <- rep(FALSE, length(sig_weight))
  candidate_topn[ord[seq_len(top_n)]] <- sig_weight[ord[seq_len(top_n)]] > 0

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
    n_target_candidate <- sum(candidate & is_target)
    precision <- if (n_candidate > 0L) n_target_candidate / n_candidate else NA_real_
    recall <- if (sum(is_target) > 0L) n_target_candidate / sum(is_target) else NA_real_
    enrichment <- precision / max(target_rate, 1e-12)

    summary_rows[[paste(trait_name, strategy, sep = "_")]] <- data.table(
      trait = trait_name,
      group_test_trait = group_trait_name,
      group_fdr_threshold = GROUP_FDR_THRESHOLD,
      strategy = strategy,
      n_sig_rna_groups = length(sig_groups),
      n_atac_cells = length(sig_weight),
      posthoc_target_name = POSTHOC_TARGET_NAME,
      n_target_atac_cells_posthoc = sum(is_target),
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
    atac_barcode = m$atac_barcodes,
    atac_cell_type_posthoc = atac_ct,
    atac_lineage_posthoc = atac_lineage,
    posthoc_target_name = POSTHOC_TARGET_NAME,
    is_target_posthoc = is_target,
    sig_group_source_weight = sig_weight,
    candidate_group_weight_gt0 = candidate_gt0,
    candidate_group_weight_ge0p50 = candidate_ge0p50,
    candidate_group_weight_ge0p80 = candidate_ge0p80,
    candidate_group_top350_positive = candidate_topn,
    no_celltype_used_to_select = TRUE
  )

  edge_dt <- data.table(
    atac_barcode = rep(m$atac_barcodes, times = ncol(m$rna_neighbor_barcode)),
    rna_barcode = as.vector(m$rna_neighbor_barcode),
    neighbor_weight = as.vector(m$neighbor_weight),
    rna_group_id = as.vector(nn_group),
    is_sig_rna_group = as.vector(is_sig_neighbor)
  )
  source_rows[[trait_name]] <- edge_dt[is_sig_rna_group == TRUE, .(
    total_weight = sum(neighbor_weight, na.rm = TRUE),
    n_edges = .N
  ), by = rna_group_id][order(-total_weight)][, trait := trait_name]

  cat("  sig RNA groups=", length(sig_groups),
      "; score AUROC=", signif(score_auc, 4),
      "; top350 target precision=",
      signif(summary_rows[[paste(trait_name, "group_top350_positive", sep = "_")]]$candidate_target_precision_posthoc, 4),
      "; ge0p80 target precision=",
      signif(summary_rows[[paste(trait_name, "group_weight_ge0p80", sep = "_")]]$candidate_target_precision_posthoc, 4),
      "\n", sep = "")
}

summary_dt <- rbindlist(summary_rows, use.names = TRUE, fill = TRUE)
cell_dt <- rbindlist(cell_rows, use.names = TRUE, fill = TRUE)
by_ct_dt <- rbindlist(by_ct_rows, use.names = TRUE, fill = TRUE)
source_dt <- rbindlist(source_rows, use.names = TRUE, fill = TRUE)

fwrite(summary_dt, file.path(OUT_DIR, "trait_screening_atac_candidate_summary.csv"))
fwrite(cell_dt, file.path(OUT_DIR, "trait_screening_atac_candidate_cell_level.csv.gz"))
fwrite(by_ct_dt, file.path(OUT_DIR, "trait_screening_atac_candidate_by_celltype_posthoc.csv"))
fwrite(source_dt, file.path(OUT_DIR, "trait_screening_sig_rna_group_source_weights.csv"))

cat("\n=== ATAC candidate screening finished ===\n")
print(summary_dt[
  strategy %in% c("group_weight_ge0p80", "group_top350_positive"),
  .(
    trait, strategy, n_sig_rna_groups, sig_group_weight_auc_posthoc,
    posthoc_target_name, n_candidate_atac, candidate_target_precision_posthoc,
    candidate_target_enrichment_posthoc, candidate_top_atac_celltypes_posthoc
  )
])
cat("\nOutputs are in: ", OUT_DIR, "\n", sep = "")
