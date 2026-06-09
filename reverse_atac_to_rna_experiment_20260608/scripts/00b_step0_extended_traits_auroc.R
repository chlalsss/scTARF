# Step 0b: extended RNA-only scDRS baseline AUROC table.
#
# Reads existing official scDRS outputs only. Extra CITE trait-screening traits
# are included when available; Granja RNA is left NA for traits without existing
# official scDRS score files.

suppressPackageStartupMessages({
  library(data.table)
  library(pROC)
})

ROOT <- "D:/RStudio/work"
REV_ROOT <- file.path(ROOT, "reverse_atac_to_rna_experiment_20260608")
OUT_DIR <- file.path(REV_ROOT, "outputs")
MANIFEST_DIR <- file.path(REV_ROOT, "inputs_manifest")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(MANIFEST_DIR, recursive = TRUE, showWarnings = FALSE)

CITE_BASE_RDS <- file.path(
  ROOT,
  "module4_traitaware_pseudopair/01_data_prep_external_rna/output/cite_rna_scdrs_scores.rds"
)
CITE_SCREEN_RDS <- file.path(
  ROOT,
  "module4_traitaware_pseudopair/01_data_prep_external_rna/output/trait_screening/cite_rna_scdrs_scores_trait_screening.rds"
)
CITE_REMAINING_RDS <- file.path(
  ROOT,
  "module4_traitaware_pseudopair/01_data_prep_external_rna/output/trait_screening/cite_rna_scdrs_scores_trait_screening_remaining.rds"
)
GRANJA_BASE_RDS <- file.path(
  ROOT,
  "module4_traitaware_pseudopair/01_data_prep/output/granja_rna_scdrs_scores.rds"
)
GRANJA_COMBINED_RDS <- file.path(
  ROOT,
  "module4_traitaware_pseudopair/01_data_prep/output/granja_rna_scdrs_scores_with_missing_traits.rds"
)

TRAIT_IDS <- c(
  LYM = "UKB_460K.blood_LYMPHOCYTE_COUNT",
  MONO = "UKB_460K.blood_MONOCYTE_COUNT",
  PLT = "UKB_460K.blood_PLATELET_COUNT",
  EOS = "UKB_460K.blood_EOSINOPHIL_COUNT",
  MCH = "UKB_460K.blood_MEAN_CORPUSCULAR_HEMOGLOBIN",
  RBC = "UKB_460K.blood_RED_COUNT",
  RDW = "UKB_460K.blood_RBC_DISTRIB_WIDTH",
  WBC = "UKB_460K.blood_WHITE_COUNT"
)

erythroid_mk_cite <- c("MK/E prog", "Proerythroblast", "Erythroblast",
                       "Normoblast", "Reticulocyte")
mono_cite <- c("CD14+ Mono", "CD16+ Mono")
myeloid_other_cite <- c("G/M prog", "cDC1", "cDC2", "pDC")
lymph_cite <- c(
  "CD4+ T activated", "CD4+ T naive", "CD4+ T activated integrinB7+",
  "CD4+ T CD314+ CD45RA+", "CD8+ T naive",
  "CD8+ T naive CD127+ CD26- CD101-", "CD8+ T CD49f+",
  "CD8+ T CD57+ CD45RA+", "CD8+ T CD57+ CD45RO+",
  "CD8+ T CD69+ CD45RA+", "CD8+ T CD69+ CD45RO+",
  "CD8+ T TIGIT+ CD45RA+", "CD8+ T TIGIT+ CD45RO+",
  "T reg", "T prog cycling", "MAIT", "gdT CD158b+", "gdT TCRVD2+",
  "dnT", "NK", "NK CD158e1+", "ILC", "ILC1", "Lymph prog",
  "Naive CD20+ B IGKC+", "Naive CD20+ B IGKC-",
  "B1 B IGKC+", "B1 B IGKC-", "Transitional B",
  "Plasma cell IGKC+", "Plasma cell IGKC-",
  "Plasmablast IGKC+", "Plasmablast IGKC-"
)

CITE_TARGETS <- list(
  LYM = lymph_cite,
  MONO = mono_cite,
  PLT = erythroid_mk_cite,
  EOS = c(mono_cite, myeloid_other_cite, lymph_cite),
  MCH = erythroid_mk_cite,
  RBC = erythroid_mk_cite,
  RDW = erythroid_mk_cite,
  WBC = c(mono_cite, myeloid_other_cite, lymph_cite)
)
CITE_TARGET_NAMES <- c(
  LYM = "lymphoid_broad",
  MONO = "monocyte",
  PLT = "erythroid_mk",
  EOS = "leukocyte_non_erythroid_proxy",
  MCH = "erythroid_mk",
  RBC = "erythroid_mk",
  RDW = "erythroid_mk",
  WBC = "leukocyte_non_erythroid_proxy"
)

GRANJA_TARGETS <- list(
  LYM = c(
    "16_Pre.B", "17_B", "18_Plasma",
    "19_CD8.N", "23_CD8.EM", "24_CD8.CM",
    "20_CD4.N1", "21_CD4.N2", "22_CD4.M",
    "25_NK"
  ),
  MONO = c("11_CD14.Mono.1", "12_CD14.Mono.2"),
  PLT = c("02_Early.Eryth", "03_Late.Eryth", "04_Early.Baso"),
  EOS = c(
    "01_HSC", "05_CMP.LMPP", "06_CLP.1", "07_GMP", "08_GMP.Neut",
    "09_pDC", "10_cDC", "11_CD14.Mono.1", "12_CD14.Mono.2",
    "13_CD16.Mono", "15_CLP.2", "16_Pre.B", "17_B", "18_Plasma",
    "19_CD8.N", "20_CD4.N1", "21_CD4.N2", "22_CD4.M",
    "23_CD8.EM", "24_CD8.CM", "25_NK"
  ),
  MCH = c("02_Early.Eryth", "03_Late.Eryth", "04_Early.Baso"),
  RBC = c("02_Early.Eryth", "03_Late.Eryth", "04_Early.Baso"),
  RDW = c("02_Early.Eryth", "03_Late.Eryth", "04_Early.Baso"),
  WBC = c(
    "01_HSC", "05_CMP.LMPP", "06_CLP.1", "07_GMP", "08_GMP.Neut",
    "09_pDC", "10_cDC", "11_CD14.Mono.1", "12_CD14.Mono.2",
    "13_CD16.Mono", "15_CLP.2", "16_Pre.B", "17_B", "18_Plasma",
    "19_CD8.N", "20_CD4.N1", "21_CD4.N2", "22_CD4.M",
    "23_CD8.EM", "24_CD8.CM", "25_NK"
  )
)
GRANJA_TARGET_NAMES <- c(
  LYM = "lymphoid_broad",
  MONO = "monocyte",
  PLT = "erythroid_baso",
  EOS = "leukocyte_non_erythroid_proxy",
  MCH = "erythroid_baso",
  RBC = "erythroid_baso",
  RDW = "erythroid_baso",
  WBC = "leukocyte_non_erythroid_proxy"
)

compute_auc <- function(score, label) {
  keep <- is.finite(score) & !is.na(label)
  score <- score[keep]
  label <- label[keep]
  if (length(score) < 20 || sum(label == 1L) < 10 || sum(label == 0L) < 10) {
    return(NA_real_)
  }
  as.numeric(pROC::auc(pROC::roc(
    label, score,
    quiet = TRUE,
    direction = "<",
    levels = c(0L, 1L)
  )))
}

merge_score_rds <- function(paths) {
  out <- as.data.table(readRDS(paths[1]))
  for (path in paths[-1]) {
    if (!file.exists(path)) next
    x <- as.data.table(readRDS(path))
    score_cols <- grep("_norm_score$", names(x), value = TRUE)
    score_cols <- setdiff(score_cols, names(out))
    if (length(score_cols) == 0) next
    out <- merge(
      out,
      x[, c("cell", score_cols), with = FALSE],
      by = "cell",
      all.x = TRUE,
      sort = FALSE
    )
  }
  out
}

score_one_dataset <- function(dataset, dt, targets, target_names, trait_order) {
  rows <- lapply(trait_order, function(trait) {
    score_col <- paste0(trait, "_norm_score")
    if (!(score_col %in% names(dt)) || !(trait %in% names(targets))) {
      return(data.table(
        dataset = dataset,
        trait = trait,
        trait_id = TRAIT_IDS[[trait]],
        target_definition = NA_character_,
        baseline_scdrs_auroc = NA_real_,
        n_scored_cells = 0L,
        n_positive_target = NA_integer_,
        note = "No existing official scDRS score for this dataset/trait"
      ))
    }

    score <- as.numeric(dt[[score_col]])
    label <- as.integer(dt$cell_type %in% targets[[trait]])
    keep <- is.finite(score) & !is.na(label)
    data.table(
      dataset = dataset,
      trait = trait,
      trait_id = TRAIT_IDS[[trait]],
      target_definition = target_names[[trait]],
      baseline_scdrs_auroc = compute_auc(score, label),
      n_scored_cells = sum(keep),
      n_positive_target = sum(label[keep] == 1L),
      note = fifelse(
        trait %in% c("EOS", "WBC"),
        "CITE has no eosinophil-specific cell type; target uses leukocyte_non_erythroid proxy.",
        "Original module4 target definition."
      )
    )
  })
  rbindlist(rows, use.names = TRUE, fill = TRUE)
}

cat(">>> Step 0b extended RNA-only scDRS baseline AUROC\n\n")
required <- c(CITE_BASE_RDS, CITE_SCREEN_RDS, CITE_REMAINING_RDS, GRANJA_BASE_RDS)
missing <- required[!file.exists(required)]
if (length(missing) > 0) {
  stop("Missing required input(s):\n  ", paste(missing, collapse = "\n  "))
}

cat("[1/3] Reading CITE base + trait-screening scDRS RDS files ...\n")
cite_dt <- merge_score_rds(c(CITE_BASE_RDS, CITE_SCREEN_RDS, CITE_REMAINING_RDS))

cat("[2/3] Reading Granja scDRS RDS file ...\n")
granja_source_rds <- if (file.exists(GRANJA_COMBINED_RDS)) GRANJA_COMBINED_RDS else GRANJA_BASE_RDS
cat("  ", granja_source_rds, "\n", sep = "")
granja_dt <- as.data.table(readRDS(granja_source_rds))

cat("[3/3] Computing extended AUROC table ...\n")
trait_order <- names(TRAIT_IDS)
long <- rbindlist(list(
  score_one_dataset("CITE_RNA", cite_dt, CITE_TARGETS, CITE_TARGET_NAMES, trait_order),
  score_one_dataset("GRANJA_RNA", granja_dt, GRANJA_TARGETS, GRANJA_TARGET_NAMES, trait_order)
), use.names = TRUE, fill = TRUE)

long[, baseline_scdrs_auroc := round(baseline_scdrs_auroc, 6)]

wide <- dcast(
  long,
  trait + trait_id ~ dataset,
  value.var = "baseline_scdrs_auroc"
)
wide[, trait_order_index := match(trait, trait_order)]
setorder(wide, trait_order_index)
wide[, weaker_dataset := fifelse(
  is.na(CITE_RNA) & is.na(GRANJA_RNA), "NO_SCORE",
  fifelse(is.na(CITE_RNA), "GRANJA_RNA_ONLY",
          fifelse(is.na(GRANJA_RNA), "CITE_RNA_ONLY",
                  fifelse(CITE_RNA <= GRANJA_RNA, "CITE_RNA", "GRANJA_RNA")))
)]
wide[, weaker_auroc := pmin(CITE_RNA, GRANJA_RNA, na.rm = TRUE)]
wide[is.infinite(weaker_auroc), weaker_auroc := NA_real_]
wide[, auroc_gap_cite_minus_granja := round(CITE_RNA - GRANJA_RNA, 6)]

target_lookup <- unique(long[dataset == "CITE_RNA", .(trait, target_definition, note)])
wide <- merge(wide, target_lookup, by = "trait", all.x = TRUE, sort = FALSE)
setcolorder(wide, c(
  "trait", "trait_id", "target_definition",
  "CITE_RNA", "GRANJA_RNA", "weaker_dataset",
  "weaker_auroc", "auroc_gap_cite_minus_granja", "note",
  "trait_order_index"
))
setorder(wide, trait_order_index)
wide[, trait_order_index := NULL]

manifest <- data.table(
  role = c(
    "CITE base official scDRS scores",
    "CITE trait-screening official scDRS scores",
    "CITE remaining trait-screening official scDRS scores",
    "Granja base official scDRS scores",
    "Granja combined official scDRS scores after missing-trait rerun"
  ),
  path = c(CITE_BASE_RDS, CITE_SCREEN_RDS, CITE_REMAINING_RDS, GRANJA_BASE_RDS, GRANJA_COMBINED_RDS),
  exists = file.exists(c(CITE_BASE_RDS, CITE_SCREEN_RDS, CITE_REMAINING_RDS, GRANJA_BASE_RDS, GRANJA_COMBINED_RDS))
)

fwrite(long, file.path(OUT_DIR, "step0b_extended_rna_only_scdrs_baseline_auroc_long.csv"))
fwrite(wide, file.path(OUT_DIR, "step0b_extended_rna_only_scdrs_baseline_auroc_wide.csv"))
fwrite(manifest, file.path(MANIFEST_DIR, "step0b_extended_inputs_manifest.csv"))

cat("\nExtended per-trait AUROC:\n")
print(wide)
cat("\nStep 0b OK.\n")
