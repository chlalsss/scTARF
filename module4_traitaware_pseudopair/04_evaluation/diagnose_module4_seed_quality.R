# module4 / diagnostic: seed and pseudo-zRNA quality after mapping + transfer
# -----------------------------------------------------------------------------
# 目的:
#   现在 module4 v2 AUROC 不理想, 最需要先判断问题是不是出在 seed。
#
#   本脚本检查:
#     1. pseudo_z_rna 本身能不能区分目标 ATAC cell type
#     2. top 1% seed 是否富集到目标 ATAC cell type
#     3. seed 主要落在哪些 cell type
#     4. 如果 TRS 已存在, 顺便比较 seed score 和最终 TRS 的 AUROC
#
# 重要:
#   cell_type 只在这里用于事后评估, 不参与 mapping / transfer / KNN / RWR。
#
# 运行:
#   setwd("D:/RStudio/work/module4_traitaware_pseudopair")
#   source("04_evaluation/diagnose_module4_seed_quality.R")
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

if (!requireNamespace("pROC", quietly = TRUE)) {
  stop("Missing pROC. In a separate R session run: install.packages('pROC')")
}

ROOT <- "D:/RStudio/work/module4_traitaware_pseudopair"
META_FILE <- file.path(ROOT, "04_evaluation/inputs/granja_cell_meta.rds")
TRANSFER_DIR <- file.path(ROOT, "03_pseudopair_pipeline/output/transfer")
MODULE3_COMPAT_DIR <- file.path(ROOT, "03_pseudopair_pipeline/output/module3_v2_compat")
OUT_DIR <- file.path(ROOT, "04_evaluation/module4_v2_diagnostics")
FIG_DIR <- file.path(OUT_DIR, "figures")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

TRAITS <- c("MONO", "PLT", "LYM")

TARGETS <- list(
  MONO = c("11_CD14.Mono.1", "12_CD14.Mono.2"),
  PLT  = c("02_Early.Eryth", "03_Late.Eryth", "04_Early.Baso"),
  LYM  = c(
    "16_Pre.B", "17_B", "18_Plasma",
    "19_CD8.N", "23_CD8.EM", "24_CD8.CM",
    "20_CD4.N1", "21_CD4.N2", "22_CD4.M",
    "25_NK"
  )
)

MYELOID_DISTRACTOR <- c(
  "08_GMP.Neut", "09_pDC", "10_cDC",
  "11_CD14.Mono.1", "12_CD14.Mono.2"
)

compute_auc <- function(score, label) {
  keep <- is.finite(score) & !is.na(label)
  score <- score[keep]
  label <- label[keep]
  if (sum(label == 1L) < 10L || sum(label == 0L) < 10L) return(NA_real_)
  as.numeric(pROC::auc(pROC::roc(
    label, score, quiet = TRUE, direction = "<", levels = c(0L, 1L)
  )))
}

load_meta <- function(path) {
  meta <- readRDS(path)
  if (inherits(meta, "data.table")) meta <- as.data.frame(meta)
  if (is.null(rownames(meta)) ||
      identical(rownames(meta), as.character(seq_len(nrow(meta))))) {
    if ("barcode" %in% colnames(meta)) rownames(meta) <- meta$barcode
  }
  if (!"cell_type" %in% colnames(meta)) {
    stop("meta missing cell_type column", call. = FALSE)
  }
  meta
}

load_trs_optional <- function(trait, atac_w) {
  f <- file.path(MODULE3_COMPAT_DIR, sprintf("TRS_v2_%s_w%.1f.tsv.gz", trait, atac_w))
  if (!file.exists(f)) return(NULL)
  dt <- fread(f)
  bc_col <- intersect(c("barcode", "cell", "cell_id", "V1"), colnames(dt))[1]
  trs_col <- intersect(c("TRS", "trs", "score"), colnames(dt))[1]
  if (is.na(bc_col) || is.na(trs_col)) {
    bc_col <- colnames(dt)[1]
    trs_col <- colnames(dt)[2]
  }
  out <- dt[[trs_col]]
  names(out) <- dt[[bc_col]]
  out
}

top_table_label <- function(x, n = 8L) {
  tab <- sort(table(x), decreasing = TRUE)
  if (length(tab) == 0L) return("")
  paste(paste0(names(head(tab, n)), ":", as.integer(head(tab, n))), collapse = "; ")
}

cat(">>> Diagnose module4 seed / pseudo-zRNA quality\n\n")

meta <- load_meta(META_FILE)
summary_rows <- list()
celltype_rows <- list()
cell_rows <- list()
plot_rows <- list()

for (trait in TRAITS) {
  cat("=== Trait ", trait, " ===\n", sep = "")
  seed_file <- file.path(MODULE3_COMPAT_DIR, paste0("seed_v2_", trait, ".rds"))
  z_file <- file.path(TRANSFER_DIR, paste0("granja_pseudo_zrna_", trait, ".rds"))
  if (!file.exists(seed_file)) stop("Missing seed file: ", seed_file, call. = FALSE)
  if (!file.exists(z_file)) stop("Missing pseudo zRNA file: ", z_file, call. = FALSE)

  seed <- readRDS(seed_file)
  z_obj <- readRDS(z_file)
  z_dt <- as.data.frame(z_obj$scores)
  if (!all(c("barcode", "pseudo_z_rna") %in% colnames(z_dt))) {
    stop("pseudo zRNA object has unexpected format: ", z_file, call. = FALSE)
  }

  common <- intersect(names(seed), rownames(meta))
  common <- intersect(common, z_dt$barcode)
  if (length(common) == 0L) stop("No common cells for trait ", trait, call. = FALSE)

  z <- z_dt$pseudo_z_rna[match(common, z_dt$barcode)]
  names(z) <- common
  seed_vec <- as.numeric(seed[common] > 0)
  cell_type <- as.character(meta[common, "cell_type"])
  target <- as.integer(cell_type %in% TARGETS[[trait]])
  myeloid <- as.integer(cell_type %in% MYELOID_DISTRACTOR)

  n_total <- length(common)
  n_target <- sum(target)
  n_seed <- sum(seed_vec)
  n_seed_target <- sum(seed_vec == 1L & target == 1L)
  target_rate <- n_target / n_total
  seed_precision <- n_seed_target / max(n_seed, 1)
  seed_recall <- n_seed_target / max(n_target, 1)
  seed_enrichment <- seed_precision / max(target_rate, 1e-12)

  pseudo_auc <- compute_auc(z, target)
  pseudo_delta <- mean(z[target == 1L], na.rm = TRUE) -
    mean(z[target == 0L], na.rm = TRUE)
  myeloid_auc <- compute_auc(z, myeloid)

  trs_07 <- load_trs_optional(trait, 0.7)
  trs_10 <- load_trs_optional(trait, 1.0)
  trs_07_auc <- if (is.null(trs_07)) NA_real_ else {
    cc <- intersect(common, names(trs_07))
    compute_auc(trs_07[cc], as.integer(meta[cc, "cell_type"] %in% TARGETS[[trait]]))
  }
  trs_10_auc <- if (is.null(trs_10)) NA_real_ else {
    cc <- intersect(common, names(trs_10))
    compute_auc(trs_10[cc], as.integer(meta[cc, "cell_type"] %in% TARGETS[[trait]]))
  }

  seed_celltypes <- top_table_label(cell_type[seed_vec == 1L], n = 10L)

  summary_rows[[trait]] <- data.frame(
    trait = trait,
    n_total = n_total,
    n_target = n_target,
    target_rate = target_rate,
    n_seed = n_seed,
    n_seed_target = n_seed_target,
    seed_target_precision = seed_precision,
    seed_target_recall = seed_recall,
    seed_target_enrichment = seed_enrichment,
    pseudo_z_target_mean = mean(z[target == 1L], na.rm = TRUE),
    pseudo_z_rest_mean = mean(z[target == 0L], na.rm = TRUE),
    pseudo_z_delta_target_minus_rest = pseudo_delta,
    pseudo_z_target_auroc = pseudo_auc,
    pseudo_z_myeloid_auroc = myeloid_auc,
    trs_w0p7_auroc = trs_07_auc,
    trs_w1p0_auroc = trs_10_auc,
    top_seed_celltypes = seed_celltypes,
    no_celltype_used_in_pipeline = TRUE,
    stringsAsFactors = FALSE
  )

  cell_dt <- data.table(
    trait = trait,
    barcode = common,
    cell_type = cell_type,
    is_target = target == 1L,
    is_myeloid_distractor = myeloid == 1L,
    pseudo_z_rna = z,
    is_seed = seed_vec == 1L
  )
  cell_rows[[trait]] <- cell_dt

  by_ct <- cell_dt[, .(
    n_cells = .N,
    n_seed = sum(is_seed),
    seed_rate = mean(is_seed),
    mean_pseudo_z = mean(pseudo_z_rna, na.rm = TRUE),
    median_pseudo_z = median(pseudo_z_rna, na.rm = TRUE),
    q95_pseudo_z = as.numeric(quantile(pseudo_z_rna, 0.95, na.rm = TRUE)),
    is_trait_target_celltype = unique(is_target)[1],
    is_myeloid_distractor_celltype = unique(is_myeloid_distractor)[1]
  ), by = .(trait, cell_type)][order(trait, -seed_rate, -n_seed)]
  celltype_rows[[trait]] <- by_ct

  plot_rows[[trait]] <- cell_dt[, .(
    trait,
    cell_type_group = fifelse(is_target, "target",
                              fifelse(is_myeloid_distractor, "myeloid_distractor", "other")),
    pseudo_z_rna,
    is_seed
  )]

  cat(sprintf(
    "  pseudo_z AUROC=%.4f, seed precision=%.3f, recall=%.3f, enrichment=%.2f\n",
    pseudo_auc, seed_precision, seed_recall, seed_enrichment
  ))
  cat("  top seed cell types: ", seed_celltypes, "\n\n", sep = "")
}

summary_dt <- rbindlist(summary_rows, use.names = TRUE, fill = TRUE)
cell_dt <- rbindlist(cell_rows, use.names = TRUE, fill = TRUE)
celltype_dt <- rbindlist(celltype_rows, use.names = TRUE, fill = TRUE)
plot_dt <- rbindlist(plot_rows, use.names = TRUE, fill = TRUE)

fwrite(summary_dt, file.path(OUT_DIR, "seed_diagnostics_summary.csv"))
fwrite(cell_dt, file.path(OUT_DIR, "seed_cell_level_scores.csv.gz"))
fwrite(celltype_dt, file.path(OUT_DIR, "seed_by_celltype.csv"))

p1 <- ggplot(summary_dt, aes(x = trait, y = pseudo_z_target_auroc, fill = trait)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0.5, linetype = 2, color = "gray40") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Pseudo zRNA target AUROC", x = NULL, y = "AUROC") +
  theme_bw() +
  theme(legend.position = "none")
ggsave(file.path(FIG_DIR, "pseudo_z_target_auroc.png"), p1,
       width = 5.5, height = 4, dpi = 150)

p2 <- ggplot(summary_dt, aes(x = trait, y = seed_target_precision, fill = trait)) +
  geom_col(width = 0.7) +
  geom_point(aes(y = target_rate), color = "black", size = 2) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Seed target precision vs background target rate",
       subtitle = "black point = target rate among all cells",
       x = NULL, y = "fraction") +
  theme_bw() +
  theme(legend.position = "none")
ggsave(file.path(FIG_DIR, "seed_precision_vs_background.png"), p2,
       width = 6, height = 4, dpi = 150)

p3 <- ggplot(plot_dt, aes(x = cell_type_group, y = pseudo_z_rna, fill = cell_type_group)) +
  geom_boxplot(outlier.size = 0.2) +
  facet_wrap(~ trait, scales = "free_y") +
  labs(title = "Pseudo zRNA distribution by diagnostic cell-type group",
       x = NULL, y = "pseudo zRNA") +
  theme_bw() +
  theme(legend.position = "none")
ggsave(file.path(FIG_DIR, "pseudo_z_by_diagnostic_group.png"), p3,
       width = 8, height = 4.5, dpi = 150)

cat("\n=== Seed diagnostics finished ===\n")
print(summary_dt)
cat("\nOutputs are in: ", OUT_DIR, "\n", sep = "")
