# Run the non-paired dense top1 workflow for traits missing top1 AUROC.
#
# This driver is intentionally scoped to non-paired Granja ATAC + CITE RNA
# reference data. It writes into a new module4 tag, so existing v3_dense_top1
# MONO/RBC/RDW outputs are not overwritten.

ROOT <- "D:/RStudio/work/module4_traitaware_pseudopair"

if (!exists("TOP1_TRAITS_TO_RUN")) {
  TOP1_TRAITS_TO_RUN <- c("EOS", "LYM", "MCH", "PLT", "WBC")
}
if (!exists("RUN_STEPS")) {
  RUN_STEPS <- c(
    "prepare_unified_inputs",
    "dense_top1_mapping",
    "transfer",
    "seed_candidates",
    "module3_blocks_and_seed",
    "knn_rwr"
  )
}
if (!exists("DENSE_TOP1_TAG")) DENSE_TOP1_TAG <- "v3_dense_top1_missing_traits"
if (!exists("V3_ATAC_W_LIST")) V3_ATAC_W_LIST <- c(0.7, 1.0)
if (!exists("GROUP_FDR_THRESHOLD")) GROUP_FDR_THRESHOLD <- 0.10

TRAIT_IDS_ALL <- c(
  PLT = "UKB_460K.blood_PLATELET_COUNT",
  LYM = "UKB_460K.blood_LYMPHOCYTE_COUNT",
  MCH = "UKB_460K.blood_MEAN_CORPUSCULAR_HEMOGLOBIN",
  RDW = "UKB_460K.blood_RBC_DISTRIB_WIDTH",
  RBC = "UKB_460K.blood_RED_COUNT",
  EOS = "UKB_460K.blood_EOSINOPHIL_COUNT",
  WBC = "UKB_460K.blood_WHITE_COUNT"
)

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
  EOS = granulocyte_proxy,
  WBC = leukocyte_proxy
)
POSTHOC_TARGET_NAMES <- c(
  PLT = "red_lineage_proxy",
  LYM = "lymph_proxy",
  MCH = "red_lineage_proxy",
  EOS = "granulocyte_like_proxy",
  WBC = "leukocyte_non_erythroid_proxy"
)

step_enabled <- function(x) x %in% RUN_STEPS

stop_missing <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0L) {
    stop("Missing required file(s):\n  ", paste(missing, collapse = "\n  "),
         call. = FALSE)
  }
}

validate_tag <- function(x) {
  if (!grepl("^[A-Za-z0-9_.-]+$", x)) {
    stop("DENSE_TOP1_TAG may only contain letters, numbers, underscore, dot, or hyphen.",
         call. = FALSE)
  }
}

run_script <- function(script, vars = list()) {
  e <- new.env(parent = globalenv())
  for (nm in names(vars)) assign(nm, vars[[nm]], envir = e)
  cat("\n>>> source(", script, ")\n", sep = "")
  source(file.path(ROOT, script), local = e)
  invisible(e)
}

read_gs_table <- function(path) {
  stop_missing(path)
  x <- utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("TRAIT", "GENESET") %in% names(x))) {
    stop("Gene-set file must contain TRAIT and GENESET: ", path, call. = FALSE)
  }
  x[, c("TRAIT", "GENESET")]
}

write_unified_gene_set <- function(out_file, trait_ids) {
  gs_files <- c(
    file.path(ROOT, "01_data_prep/03_rna_scdrs_python/granja_3traits.gs"),
    file.path(ROOT, "01_data_prep_external_rna/03_cite_scdrs_python/cite_trait_screening.gs")
  )
  gs <- do.call(rbind, lapply(gs_files, read_gs_table))
  keep <- gs[gs$TRAIT %in% unname(trait_ids), , drop = FALSE]
  missing_ids <- setdiff(unname(trait_ids), keep$TRAIT)
  if (length(missing_ids) > 0L) {
    stop("Missing trait IDs in gene-set files: ",
         paste(missing_ids, collapse = ", "), call. = FALSE)
  }
  keep <- keep[match(unname(trait_ids), keep$TRAIT), , drop = FALSE]
  utils::write.table(keep, out_file, sep = "\t", row.names = FALSE, quote = FALSE)
  out_file
}

merge_score_rds <- function(out_rds, out_csv, traits) {
  base_rds <- file.path(ROOT, "01_data_prep_external_rna/output/cite_rna_scdrs_scores.rds")
  screening_rds <- file.path(
    ROOT,
    "01_data_prep_external_rna/output/trait_screening/cite_rna_scdrs_scores_trait_screening.rds"
  )
  remaining_rds <- file.path(
    ROOT,
    "01_data_prep_external_rna/output/trait_screening/cite_rna_scdrs_scores_trait_screening_remaining.rds"
  )
  stop_missing(base_rds)
  out <- as.data.frame(readRDS(base_rds), stringsAsFactors = FALSE)
  for (path in c(screening_rds, remaining_rds)) {
    if (!file.exists(path)) next
    x <- as.data.frame(readRDS(path), stringsAsFactors = FALSE)
    score_cols <- grep("_norm_score$", names(x), value = TRUE)
    score_cols <- setdiff(score_cols, names(out))
    if (length(score_cols) > 0L) {
      out <- merge(out, x[, c("cell", score_cols), drop = FALSE],
                   by = "cell", all.x = TRUE, sort = FALSE)
    }
  }
  missing_cols <- setdiff(paste0(traits, "_norm_score"), names(out))
  if (length(missing_cols) > 0L) {
    stop("Unified CITE score RDS missing: ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  saveRDS(out, out_rds)
  con <- gzfile(out_csv, open = "wt")
  on.exit(close(con), add = TRUE)
  utils::write.csv(out, con, row.names = FALSE, quote = TRUE)
  out_rds
}

validate_tag(DENSE_TOP1_TAG)
if ("MONO" %in% TOP1_TRAITS_TO_RUN) {
  stop("This driver is for traits still missing top1. MONO already has a top1 run.",
       call. = FALSE)
}
unknown <- setdiff(TOP1_TRAITS_TO_RUN, names(TRAIT_IDS_ALL))
if (length(unknown) > 0L) {
  stop("Unknown TOP1_TRAITS_TO_RUN: ", paste(unknown, collapse = ", "),
       call. = FALSE)
}

TRAIT_IDS <- TRAIT_IDS_ALL[TOP1_TRAITS_TO_RUN]
BASE_OUT_DIR <- file.path(ROOT, "02_trait_mapping/output", DENSE_TOP1_TAG)
UNIFIED_INPUT_DIR <- file.path(BASE_OUT_DIR, "_inputs")
MAP_DIR <- BASE_OUT_DIR
TRANSFER_DIR <- file.path(ROOT, "03_pseudopair_pipeline/output", DENSE_TOP1_TAG, "transfer")
SEED_DIAG_DIR <- file.path(BASE_OUT_DIR, "neighborhood_seed_diagnostics")
COMPAT_DIR <- file.path(ROOT, "03_pseudopair_pipeline/output", DENSE_TOP1_TAG, "module3_neighborhood_fdr_ge0p80")

UNIFIED_GS_FILE <- file.path(UNIFIED_INPUT_DIR, paste0(DENSE_TOP1_TAG, ".gs"))
UNIFIED_SCORE_RDS <- file.path(UNIFIED_INPUT_DIR, paste0("cite_rna_scdrs_scores_", DENSE_TOP1_TAG, ".rds"))
UNIFIED_SCORE_CSV <- file.path(UNIFIED_INPUT_DIR, paste0("cite_rna_scdrs_scores_", DENSE_TOP1_TAG, ".csv.gz"))
SEED_CELL_FILE <- file.path(SEED_DIAG_DIR, "sparse_neighborhood_seed_cell_level.csv.gz")

dir.create(UNIFIED_INPUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat(">>> non-paired dense top1 missing-traits driver\n\n")
cat("Traits: ", paste(TOP1_TRAITS_TO_RUN, collapse = ", "), "\n", sep = "")
cat("Tag   : ", DENSE_TOP1_TAG, "\n", sep = "")
cat("Steps : ", paste(RUN_STEPS, collapse = ", "), "\n\n", sep = "")

if (step_enabled("prepare_unified_inputs")) {
  write_unified_gene_set(UNIFIED_GS_FILE, TRAIT_IDS)
  merge_score_rds(UNIFIED_SCORE_RDS, UNIFIED_SCORE_CSV, TOP1_TRAITS_TO_RUN)
}

if (step_enabled("dense_top1_mapping")) {
  run_script("02_trait_mapping/01_align_granja_cite_by_trait.R", list(
    GS_FILE = UNIFIED_GS_FILE,
    ALIGN_OUT_DIR = MAP_DIR,
    TRAITS_ALL = TRAIT_IDS,
    TRAITS_TO_RUN = TOP1_TRAITS_TO_RUN,
    K_NEIGHBOR = 1L,
    BATCH_SIZE = 512L,
    USE_GWAS_WEIGHT = TRUE,
    FORCE_REEXTRACT_CITE = FALSE
  ))
}

if (step_enabled("transfer")) {
  run_script("03_pseudopair_pipeline/01_transfer_dual.R", list(
    TRANSFER_MAP_DIR = MAP_DIR,
    TRANSFER_OUT_DIR = TRANSFER_DIR,
    CITE_MTX = file.path(MAP_DIR, "cite_extract/cite_trait_counts.mtx"),
    CITE_GENES = file.path(MAP_DIR, "cite_extract/cite_trait_counts_genes.txt"),
    CITE_CELLS = file.path(MAP_DIR, "cite_extract/cite_trait_counts_cells.txt"),
    CITE_SCORE_RDS = UNIFIED_SCORE_RDS,
    TRAITS_ALL = TOP1_TRAITS_TO_RUN,
    TRAITS_TO_RUN = TOP1_TRAITS_TO_RUN
  ))
}

if (step_enabled("seed_candidates")) {
  run_script("02_trait_mapping/09_screen_sparse_mapping_neighborhood_seeds.R", list(
    MAP_DIR = MAP_DIR,
    OUT_DIR = SEED_DIAG_DIR,
    TRAITS_TO_SCREEN = TOP1_TRAITS_TO_RUN,
    RNA_GROUP_FILE = file.path(
      ROOT,
      "04_evaluation/module4_v2_diagnostics/neighborhood_fdr/rna_unsupervised_group_membership.csv.gz"
    ),
    RNA_GROUP_TEST_FILES = c(
      file.path(ROOT, "04_evaluation/module4_v2_diagnostics/neighborhood_fdr/neighborhood_fdr_group_tests.csv"),
      file.path(ROOT, "01_data_prep_external_rna/output/trait_screening/cite_trait_screening_neighborhood_fdr_group_tests.csv")
    ),
    GROUP_TEST_TRAIT_MAP = c(EOS = "EOSINOPHIL", WBC = "WHITE"),
    GROUP_FDR_THRESHOLD = GROUP_FDR_THRESHOLD,
    POSTHOC_TARGETS_LIST = POSTHOC_TARGETS_LIST,
    POSTHOC_TARGET_NAMES = POSTHOC_TARGET_NAMES
  ))
}

if (step_enabled("module3_blocks_and_seed")) {
  run_script("03_pseudopair_pipeline/07_prepare_trait_screening_module3_inputs.R", list(
    MODULE3_TRANSFER_DIR = TRANSFER_DIR,
    MODULE3_COMPAT_OUT_DIR = COMPAT_DIR,
    TRAITS_TO_RUN = TOP1_TRAITS_TO_RUN,
    SEED_CELL_FILE = SEED_CELL_FILE,
    SEED_SUMMARY_FILE = file.path(COMPAT_DIR, "seed_summary_neighborhood_fdr_ge0p80.csv")
  ))
}

if (step_enabled("knn_rwr")) {
  run_script("03_pseudopair_pipeline/03_run_module3_v2_knn_rwr.R", list(
    MODULE3_COMPAT_OUT_DIR = COMPAT_DIR,
    TRAITS_TO_RUN = TOP1_TRAITS_TO_RUN,
    ATAC_W_LIST = V3_ATAC_W_LIST
  ))
}

cat("\n=== dense top1 missing-traits selected steps finished ===\n")
cat("Mapping dir: ", MAP_DIR, "\n", sep = "")
cat("Module3 dir: ", COMPAT_DIR, "\n", sep = "")
