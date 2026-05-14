# module4 / Phase 1.4: check trait-aware alignment with cell-type labels
# -----------------------------------------------------------------------------
# 目的:
#   这一步只做事后检查. 它会读取 01 脚本产生的 mapping, 用 CITE RNA 的
#   cell_type 给每个 ATAC 细胞做加权投票标签, 再和 Granja ATAC 自己的
#   cell_type 做粗粒度比较.
#
# 重要:
#   cell_type 没有参与 01 脚本的对齐过程. 这里只用于检查对齐是否合理.
#
# 运行方式:
#   source("02_trait_mapping/02_check_alignment_by_celltype.R")
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

ROOT <- "D:/RStudio/work/module4_traitaware_pseudopair"
MAP_DIR <- file.path(ROOT, "02_trait_mapping/output")
OUT_DIR <- file.path(MAP_DIR, "diagnostics")
FIG_DIR <- file.path(OUT_DIR, "figures")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

ATAC_META_PATH <- file.path(ROOT, "04_evaluation/inputs/granja_cell_meta.rds")
CITE_SCORE_RDS <- file.path(ROOT, "01_data_prep_external_rna/output/cite_rna_scdrs_scores.rds")

if (!exists("TRAITS_TO_CHECK")) TRAITS_TO_CHECK <- c("MONO", "PLT", "LYM")

ATAC_TARGETS <- list(
  MONO = c("11_CD14.Mono.1", "12_CD14.Mono.2"),
  PLT = c("02_Early.Eryth", "03_Late.Eryth", "04_Early.Baso"),
  LYM = c(
    "16_Pre.B", "17_B", "18_Plasma",
    "19_CD8.N", "20_CD4.N1", "21_CD4.N2", "22_CD4.M",
    "23_CD8.EM", "24_CD8.CM", "25_NK"
  )
)

CITE_TARGETS <- list(
  MONO = c("CD14+ Mono", "CD16+ Mono"),
  PLT = c("MK/E prog", "Proerythroblast", "Erythroblast", "Normoblast", "Reticulocyte"),
  LYM = c(
    "CD4+ T activated", "CD4+ T naive", "CD4+ T activated integrinB7+",
    "CD4+ T CD314+ CD45RA+", "CD8+ T naive",
    "CD8+ T naive CD127+ CD26- CD101-", "CD8+ T CD49f+",
    "CD8+ T CD57+ CD45RA+", "CD8+ T CD57+ CD45RO+",
    "CD8+ T CD69+ CD45RA+", "CD8+ T CD69+ CD45RO+",
    "CD8+ T TIGIT+ CD45RA+", "CD8+ T TIGIT+ CD45RO+",
    "T reg", "T prog cycling", "MAIT", "gdT CD158b+", "gdT TCRVD2+",
    "dnT", "NK", "NK CD158e1+", "ILC", "ILC1", "Lymph prog",
    "Naive CD20+ B IGKC+", "Naive CD20+ B IGKC-", "B1 B IGKC+",
    "B1 B IGKC-", "Transitional B", "Plasma cell IGKC+",
    "Plasma cell IGKC-", "Plasmablast IGKC+", "Plasmablast IGKC-"
  )
)

coarse_atac <- function(x) {
  out <- rep("OTHER", length(x))
  out[x %in% ATAC_TARGETS$MONO] <- "MONO"
  out[x %in% ATAC_TARGETS$PLT] <- "ERY_MK_PLT"
  out[x %in% ATAC_TARGETS$LYM] <- "LYM"
  out[x %in% c("09_pDC", "10_cDC", "08_GMP.Neut")] <- "MYELOID_OTHER"
  out[x %in% c("01_HSC", "05_CMP.LMPP", "06_CLP.1", "07_GMP", "15_CLP.2")] <- "PROGENITOR"
  out[grepl("Unk", x)] <- "UNKNOWN"
  out
}

coarse_cite <- function(x) {
  out <- rep("OTHER", length(x))
  out[x %in% CITE_TARGETS$MONO] <- "MONO"
  out[x %in% CITE_TARGETS$PLT] <- "ERY_MK_PLT"
  out[x %in% CITE_TARGETS$LYM] <- "LYM"
  out[x %in% c("G/M prog", "cDC1", "cDC2", "pDC")] <- "MYELOID_OTHER"
  out[x %in% c("HSC")] <- "PROGENITOR"
  out
}

weighted_vote <- function(types, weights) {
  n <- nrow(types)
  vapply(seq_len(n), function(i) {
    s <- tapply(weights[i, ], types[i, ], sum)
    names(s)[which.max(s)]
  }, character(1))
}

top_pred_label <- function(x) {
  tab <- sort(table(x), decreasing = TRUE)
  if (length(tab) == 0L) return(NA_character_)
  paste0(names(tab)[1], " (", as.integer(tab[1]), ")")
}

cat(">>> Phase 1.4 alignment diagnostics by cell type\n\n")

atac_meta <- readRDS(ATAC_META_PATH)
if (inherits(atac_meta, "data.table")) atac_meta <- as.data.frame(atac_meta)
if (is.null(rownames(atac_meta)) || identical(rownames(atac_meta), as.character(seq_len(nrow(atac_meta))))) {
  rownames(atac_meta) <- atac_meta$barcode
}
cite_meta <- readRDS(CITE_SCORE_RDS)
cite_meta <- as.data.frame(cite_meta)
rownames(cite_meta) <- cite_meta$cell

summary_rows <- list()
cell_rows <- list()
by_atac_rows <- list()
conf_rows <- list()
distance_rows <- list()

for (trait in TRAITS_TO_CHECK) {
  cat("=== Checking ", trait, " ===\n", sep = "")
  map_file <- file.path(MAP_DIR, paste0("granja_cite_mapping_", trait, ".rds"))
  if (!file.exists(map_file)) stop("Missing mapping file: ", map_file, call. = FALSE)
  m <- readRDS(map_file)

  atac_bc <- m$atac_barcodes
  nn_bc <- m$rna_neighbor_barcode
  w <- m$neighbor_weight
  d <- m$neighbor_distance

  true_ct <- atac_meta[atac_bc, "cell_type"]
  if (any(is.na(true_ct))) stop("Some ATAC barcodes are missing in ATAC meta")

  nn_type <- matrix(
    cite_meta[as.vector(nn_bc), "cell_type"],
    nrow = nrow(nn_bc),
    ncol = ncol(nn_bc)
  )
  if (any(is.na(nn_type))) stop("Some neighbor CITE barcodes are missing in CITE meta")

  pred_ct <- weighted_vote(nn_type, w)
  true_coarse <- coarse_atac(true_ct)
  pred_coarse <- coarse_cite(pred_ct)

  target_true <- true_ct %in% ATAC_TARGETS[[trait]]
  target_pred <- pred_ct %in% CITE_TARGETS[[trait]]
  target_recall <- if (sum(target_true) > 0) mean(target_pred[target_true]) else NA_real_
  target_specificity <- if (sum(!target_true) > 0) mean(!target_pred[!target_true]) else NA_real_
  target_balanced <- mean(c(target_recall, target_specificity), na.rm = TRUE)
  coarse_agreement <- mean(true_coarse == pred_coarse)

  top1 <- d[, 1]
  iso_thr <- as.numeric(quantile(top1, 0.95, na.rm = TRUE)) * 1.05
  pct_iso <- mean(top1 > iso_thr) * 100
  row_ok <- all(abs(rowSums(w) - 1) < 1e-6)

  summary_rows[[trait]] <- data.frame(
    trait = trait,
    n_atac_cells = length(atac_bc),
    n_target_atac_cells = sum(target_true),
    n_g_trait_used = length(m$g_trait_used),
    k = m$k,
    no_celltype_used_in_alignment = TRUE,
    weight_rowsum_ok = row_ok,
    pct_isolated = pct_iso,
    cos_dist_median = median(d^2 / 2),
    cos_dist_q05 = as.numeric(quantile(d^2 / 2, 0.05)),
    cos_dist_q95 = as.numeric(quantile(d^2 / 2, 0.95)),
    coarse_lineage_agreement = coarse_agreement,
    trait_target_recall = target_recall,
    trait_target_specificity = target_specificity,
    trait_target_balanced_accuracy = target_balanced,
    pred_target_rate_overall = mean(target_pred),
    stringsAsFactors = FALSE
  )

  cell_dt <- data.table(
    trait = trait,
    atac_barcode = atac_bc,
    atac_cell_type = true_ct,
    atac_coarse = true_coarse,
    predicted_cite_cell_type = pred_ct,
    predicted_cite_coarse = pred_coarse,
    is_trait_target_atac = target_true,
    is_trait_target_predicted = target_pred,
    top1_distance = top1,
    top1_cos_dist = top1^2 / 2
  )
  cell_rows[[trait]] <- cell_dt

  by_ct <- cell_dt[, .(
    n_cells = .N,
    coarse_agreement = mean(atac_coarse == predicted_cite_coarse),
    trait_target_recall = if (sum(is_trait_target_atac) > 0) {
      mean(is_trait_target_predicted[is_trait_target_atac])
    } else {
      NA_real_
    },
    pred_target_rate = mean(is_trait_target_predicted),
    most_common_predicted_cite_type = top_pred_label(predicted_cite_cell_type)
  ), by = .(trait, atac_cell_type, atac_coarse)]
  by_atac_rows[[trait]] <- by_ct

  conf <- cell_dt[, .N, by = .(trait, atac_coarse, predicted_cite_coarse)]
  conf[, frac := N / sum(N), by = .(trait, atac_coarse)]
  conf_rows[[trait]] <- conf

  distance_rows[[trait]] <- data.table(
    trait = trait,
    cos_dist = as.vector(d^2 / 2)
  )

  cat(sprintf(
    "  coarse agreement = %.3f; target recall = %.3f; target specificity = %.3f\n",
    coarse_agreement, target_recall, target_specificity
  ))
}

summary_dt <- rbindlist(summary_rows, use.names = TRUE, fill = TRUE)
cell_dt <- rbindlist(cell_rows, use.names = TRUE, fill = TRUE)
by_atac_dt <- rbindlist(by_atac_rows, use.names = TRUE, fill = TRUE)
conf_dt <- rbindlist(conf_rows, use.names = TRUE, fill = TRUE)
dist_dt <- rbindlist(distance_rows, use.names = TRUE, fill = TRUE)

fwrite(summary_dt, file.path(OUT_DIR, "alignment_diagnostics_summary.csv"))
fwrite(cell_dt, file.path(OUT_DIR, "alignment_cell_level_labels.csv.gz"))
fwrite(by_atac_dt, file.path(OUT_DIR, "alignment_by_atac_celltype.csv"))
fwrite(conf_dt, file.path(OUT_DIR, "alignment_coarse_confusion.csv"))

p1 <- ggplot(summary_dt, aes(x = trait, y = trait_target_balanced_accuracy, fill = trait)) +
  geom_col(width = 0.7) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = NULL, y = "target balanced accuracy",
       title = "Trait-target consistency after alignment") +
  theme_bw() +
  theme(legend.position = "none")
ggsave(file.path(FIG_DIR, "alignment_target_balanced_accuracy.png"),
       p1, width = 5.5, height = 4, dpi = 150)

p2 <- ggplot(conf_dt, aes(x = predicted_cite_coarse, y = atac_coarse, fill = frac)) +
  geom_tile(color = "white") +
  facet_wrap(~ trait, ncol = 1) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(x = "predicted CITE coarse label", y = "ATAC coarse label",
       fill = "row fraction", title = "Coarse lineage agreement after alignment") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
ggsave(file.path(FIG_DIR, "alignment_coarse_confusion_heatmap.png"),
       p2, width = 7.5, height = 8.5, dpi = 150)

p3 <- ggplot(dist_dt, aes(x = cos_dist, color = trait, fill = trait)) +
  geom_density(alpha = 0.25) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(x = "cosine distance", y = "density",
       title = "Top-k neighbor distance distribution") +
  theme_bw()
ggsave(file.path(FIG_DIR, "alignment_distance_distribution.png"),
       p3, width = 6.5, height = 4, dpi = 150)

cat("\n=== Alignment diagnostics finished ===\n")
print(summary_dt)
cat("\nOutputs are in: ", OUT_DIR, "\n", sep = "")
