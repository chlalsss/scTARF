# module4 / diagnostic: where transferred RNA signal comes from
# -----------------------------------------------------------------------------
# 目的:
#   上一步已经看到 pseudo-zRNA / seed 质量不稳定, 特别是 PLT 明显被带偏.
#   本脚本继续往前追一层:
#     1. CITE RNA 自己的 scDRS 分数在源数据里是否指向正确 cell type
#     2. 每个 ATAC cell 在 mapping 后主要借到了哪些 CITE RNA neighbor
#     3. top seed 的 RNA 来源是目标谱系, myeloid, 还是其它细胞
#
# 重要:
#   cell_type 只用于事后诊断. mapping / transfer / KNN / RWR 中没有使用 cell_type.
#
# 运行:
#   setwd("D:/RStudio/work/module4_traitaware_pseudopair")
#   source("04_evaluation/diagnose_module4_rna_sources.R")
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

if (!requireNamespace("pROC", quietly = TRUE)) {
  stop("Missing pROC. In a separate R session run: install.packages('pROC')")
}

ROOT <- "D:/RStudio/work/module4_traitaware_pseudopair"
MAP_DIR <- file.path(ROOT, "02_trait_mapping/output")
TRANSFER_DIR <- file.path(ROOT, "03_pseudopair_pipeline/output/transfer")
MODULE3_COMPAT_DIR <- file.path(ROOT, "03_pseudopair_pipeline/output/module3_v2_compat")
ATAC_META_FILE <- file.path(ROOT, "04_evaluation/inputs/granja_cell_meta.rds")
CITE_SCORE_RDS <- file.path(ROOT, "01_data_prep_external_rna/output/cite_rna_scdrs_scores.rds")
OUT_DIR <- file.path(ROOT, "04_evaluation/module4_v2_diagnostics/rna_sources")
FIG_DIR <- file.path(OUT_DIR, "figures")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

TRAITS <- c("MONO", "PLT", "LYM")

ATAC_TARGETS <- list(
  MONO = c("11_CD14.Mono.1", "12_CD14.Mono.2"),
  PLT  = c("02_Early.Eryth", "03_Late.Eryth", "04_Early.Baso"),
  LYM  = c(
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
    "Naive CD20+ B IGKC+", "Naive CD20+ B IGKC-",
    "B1 B IGKC+", "B1 B IGKC-", "Transitional B",
    "Plasma cell IGKC+", "Plasma cell IGKC-",
    "Plasmablast IGKC+", "Plasmablast IGKC-"
  )
)

CITE_MYELOID <- c("CD14+ Mono", "CD16+ Mono", "G/M prog", "cDC1", "cDC2", "pDC")

coarse_cite <- function(x) {
  out <- rep("OTHER", length(x))
  out[x %in% CITE_TARGETS$MONO] <- "MONO"
  out[x %in% CITE_TARGETS$PLT] <- "ERY_MK_PLT"
  out[x %in% CITE_TARGETS$LYM] <- "LYM"
  out[x %in% c("G/M prog", "cDC1", "cDC2", "pDC")] <- "MYELOID_OTHER"
  out[x %in% c("HSC")] <- "PROGENITOR"
  out
}

compute_auc <- function(score, label) {
  keep <- is.finite(score) & !is.na(label)
  score <- score[keep]
  label <- label[keep]
  if (sum(label == 1L) < 10L || sum(label == 0L) < 10L) return(NA_real_)
  as.numeric(pROC::auc(pROC::roc(
    label, score, quiet = TRUE, direction = "<", levels = c(0L, 1L)
  )))
}

top_table_label <- function(x, n = 10L) {
  tab <- sort(table(x), decreasing = TRUE)
  if (length(tab) == 0L) return("")
  paste(paste0(names(head(tab, n)), ":", as.integer(head(tab, n))), collapse = "; ")
}

weighted_top_label <- function(types, weights) {
  vapply(seq_len(nrow(types)), function(i) {
    s <- tapply(weights[i, ], types[i, ], sum)
    names(s)[which.max(s)]
  }, character(1))
}

load_atac_meta <- function(path) {
  meta <- readRDS(path)
  if (inherits(meta, "data.table")) meta <- as.data.frame(meta)
  if (is.null(rownames(meta)) ||
      identical(rownames(meta), as.character(seq_len(nrow(meta))))) {
    if ("barcode" %in% colnames(meta)) rownames(meta) <- meta$barcode
  }
  if (!"cell_type" %in% colnames(meta)) {
    stop("ATAC meta missing cell_type column", call. = FALSE)
  }
  meta
}

source_score_summary <- function(cite_meta, trait) {
  score_col <- paste0(trait, "_norm_score")
  if (!score_col %in% colnames(cite_meta)) {
    stop("CITE score RDS missing column: ", score_col, call. = FALSE)
  }

  score <- as.numeric(cite_meta[[score_col]])
  label <- as.integer(cite_meta$cell_type %in% CITE_TARGETS[[trait]])
  q <- as.numeric(stats::quantile(score, 0.99, na.rm = TRUE))
  top <- is.finite(score) & score >= q
  target_rate <- mean(label == 1L, na.rm = TRUE)
  top_precision <- mean(label[top] == 1L, na.rm = TRUE)

  data.frame(
    trait = trait,
    cite_n_cells = nrow(cite_meta),
    cite_target_rate = target_rate,
    cite_score_target_auroc = compute_auc(score, label),
    cite_top1pct_threshold = q,
    cite_top1pct_precision = top_precision,
    cite_top1pct_enrichment = top_precision / max(target_rate, 1e-12),
    cite_top1pct_celltypes = top_table_label(cite_meta$cell_type[top], 10L),
    stringsAsFactors = FALSE
  )
}

cat(">>> Diagnose module4 RNA source composition\n\n")

atac_meta <- load_atac_meta(ATAC_META_FILE)
cite_meta <- readRDS(CITE_SCORE_RDS)
cite_meta <- as.data.frame(cite_meta)
if (!all(c("cell", "cell_type") %in% colnames(cite_meta))) {
  stop("CITE score RDS must contain cell and cell_type columns.", call. = FALSE)
}
rownames(cite_meta) <- cite_meta$cell

summary_rows <- list()
cell_rows <- list()
celltype_rows <- list()
edge_rows <- list()
plot_rows <- list()

for (trait in TRAITS) {
  cat("=== Trait ", trait, " ===\n", sep = "")

  map_file <- file.path(MAP_DIR, paste0("granja_cite_mapping_", trait, ".rds"))
  seed_file <- file.path(MODULE3_COMPAT_DIR, paste0("seed_v2_", trait, ".rds"))
  z_file <- file.path(TRANSFER_DIR, paste0("granja_pseudo_zrna_", trait, ".rds"))
  if (!file.exists(map_file)) stop("Missing mapping file: ", map_file, call. = FALSE)
  if (!file.exists(seed_file)) stop("Missing seed file: ", seed_file, call. = FALSE)
  if (!file.exists(z_file)) stop("Missing pseudo zRNA file: ", z_file, call. = FALSE)

  mapping <- readRDS(map_file)
  seed <- readRDS(seed_file)
  z_obj <- readRDS(z_file)
  z_dt <- as.data.frame(z_obj$scores)

  atac_bc <- mapping$atac_barcodes
  nn_bc <- mapping$rna_neighbor_barcode
  w <- mapping$neighbor_weight
  d <- mapping$neighbor_distance
  k <- ncol(nn_bc)

  if (!all(atac_bc %in% rownames(atac_meta))) {
    stop("Some ATAC barcodes are missing from ATAC meta for ", trait, call. = FALSE)
  }
  if (!all(as.vector(nn_bc) %in% rownames(cite_meta))) {
    stop("Some mapped CITE neighbors are missing from CITE meta for ", trait, call. = FALSE)
  }

  score_col <- paste0(trait, "_norm_score")
  cite_target_types <- CITE_TARGETS[[trait]]
  atac_target_types <- ATAC_TARGETS[[trait]]

  nn_type <- matrix(
    cite_meta[as.vector(nn_bc), "cell_type"],
    nrow = nrow(nn_bc),
    ncol = ncol(nn_bc)
  )
  nn_score <- matrix(
    as.numeric(cite_meta[as.vector(nn_bc), score_col]),
    nrow = nrow(nn_bc),
    ncol = ncol(nn_bc)
  )

  source_target_mass <- rowSums(w * (nn_type %in% cite_target_types), na.rm = TRUE)
  source_myeloid_mass <- rowSums(w * (nn_type %in% CITE_MYELOID), na.rm = TRUE)
  source_score_weighted <- rowSums(w * nn_score, na.rm = TRUE)
  pred_source_type <- weighted_top_label(nn_type, w)
  pred_source_coarse <- coarse_cite(pred_source_type)

  pseudo_z <- z_dt$pseudo_z_rna[match(atac_bc, z_dt$barcode)]
  seed_vec <- as.logical(seed[atac_bc] > 0)
  seed_vec[is.na(seed_vec)] <- FALSE
  atac_type <- as.character(atac_meta[atac_bc, "cell_type"])
  atac_target <- atac_type %in% atac_target_types

  cell_dt <- data.table(
    trait = trait,
    barcode = atac_bc,
    atac_cell_type = atac_type,
    is_atac_trait_target = atac_target,
    is_seed = seed_vec,
    pseudo_z_rna = pseudo_z,
    weighted_source_score = source_score_weighted,
    source_target_weight = source_target_mass,
    source_myeloid_weight = source_myeloid_mass,
    predicted_source_cell_type = pred_source_type,
    predicted_source_coarse = pred_source_coarse,
    top1_source_cell_type = nn_type[, 1],
    top1_source_coarse = coarse_cite(nn_type[, 1]),
    top1_source_score = nn_score[, 1],
    top1_distance = d[, 1],
    no_celltype_used_in_pipeline = TRUE
  )
  cell_rows[[trait]] <- cell_dt

  edge_dt <- data.table(
    trait = trait,
    atac_barcode = rep(atac_bc, times = k),
    rank = rep(seq_len(k), each = length(atac_bc)),
    cite_barcode = as.vector(nn_bc),
    neighbor_weight = as.vector(w),
    neighbor_distance = as.vector(d),
    cite_cell_type = as.vector(nn_type),
    cite_coarse = coarse_cite(as.vector(nn_type)),
    cite_trait_score = as.vector(nn_score)
  )
  edge_dt[, atac_cell_type := as.character(atac_meta[atac_barcode, "cell_type"])]
  edge_dt[, `:=`(
    cite_is_trait_target = cite_cell_type %in% cite_target_types,
    cite_is_myeloid = cite_cell_type %in% CITE_MYELOID,
    is_atac_trait_target = atac_cell_type %in% atac_target_types,
    is_seed = seed_vec[match(atac_barcode, atac_bc)]
  )]
  edge_rows[[trait]] <- edge_dt

  source_s <- source_score_summary(cite_meta, trait)
  top_seed_sources <- edge_dt[is_seed == TRUE, .(
    total_weight = sum(neighbor_weight)
  ), by = cite_cell_type][order(-total_weight)]
  top_seed_sources_label <- paste(
    paste0(head(top_seed_sources$cite_cell_type, 10L), ":",
           round(head(top_seed_sources$total_weight, 10L), 1)),
    collapse = "; "
  )

  summary_rows[[trait]] <- cbind(
    source_s,
    data.frame(
      atac_n_cells = nrow(cell_dt),
      atac_n_target = sum(atac_target),
      atac_target_rate = mean(atac_target),
      atac_n_seed = sum(seed_vec),
      pseudo_z_target_auroc = compute_auc(pseudo_z, as.integer(atac_target)),
      source_target_weight_auroc = compute_auc(source_target_mass, as.integer(atac_target)),
      mean_source_target_weight_all = mean(source_target_mass),
      mean_source_target_weight_atac_target = mean(source_target_mass[atac_target]),
      mean_source_target_weight_atac_rest = mean(source_target_mass[!atac_target]),
      mean_source_target_weight_seed = mean(source_target_mass[seed_vec]),
      mean_source_myeloid_weight_all = mean(source_myeloid_mass),
      mean_source_myeloid_weight_seed = mean(source_myeloid_mass[seed_vec]),
      mean_source_score_seed = mean(source_score_weighted[seed_vec]),
      top_seed_source_celltypes_by_weight = top_seed_sources_label,
      stringsAsFactors = FALSE
    )
  )

  by_ct <- cell_dt[, .(
    n_cells = .N,
    n_seed = sum(is_seed),
    seed_rate = mean(is_seed),
    mean_pseudo_z = mean(pseudo_z_rna, na.rm = TRUE),
    mean_source_target_weight = mean(source_target_weight, na.rm = TRUE),
    mean_source_myeloid_weight = mean(source_myeloid_weight, na.rm = TRUE),
    common_predicted_source_type = top_table_label(predicted_source_cell_type, 5L),
    is_trait_target_celltype = unique(is_atac_trait_target)[1]
  ), by = .(trait, atac_cell_type)][order(trait, -seed_rate, -mean_pseudo_z)]
  celltype_rows[[trait]] <- by_ct

  plot_rows[[trait]] <- cell_dt[, .(
    trait,
    group = fifelse(is_seed, "seed",
                    fifelse(is_atac_trait_target, "ATAC_target", "ATAC_rest")),
    source_target_weight,
    source_myeloid_weight
  )]

  cat(sprintf(
    "  CITE score AUROC=%.4f; pseudo_z ATAC AUROC=%.4f; source-target-mass AUROC=%.4f\n",
    source_s$cite_score_target_auroc,
    compute_auc(pseudo_z, as.integer(atac_target)),
    compute_auc(source_target_mass, as.integer(atac_target))
  ))
  cat(sprintf(
    "  seed mean source target weight=%.3f; seed mean source myeloid weight=%.3f\n",
    mean(source_target_mass[seed_vec]), mean(source_myeloid_mass[seed_vec])
  ))
  cat("  top seed RNA sources by weight: ", top_seed_sources_label, "\n\n", sep = "")
}

summary_dt <- rbindlist(summary_rows, use.names = TRUE, fill = TRUE)
cell_dt <- rbindlist(cell_rows, use.names = TRUE, fill = TRUE)
celltype_dt <- rbindlist(celltype_rows, use.names = TRUE, fill = TRUE)
edge_dt <- rbindlist(edge_rows, use.names = TRUE, fill = TRUE)
plot_dt <- rbindlist(plot_rows, use.names = TRUE, fill = TRUE)

fwrite(summary_dt, file.path(OUT_DIR, "rna_source_diagnostics_summary.csv"))
fwrite(cell_dt, file.path(OUT_DIR, "rna_source_cell_level.csv.gz"))
fwrite(celltype_dt, file.path(OUT_DIR, "rna_source_by_atac_celltype.csv"))
fwrite(edge_dt, file.path(OUT_DIR, "rna_source_neighbor_edges.csv.gz"))

p1 <- ggplot(summary_dt, aes(x = trait, y = cite_score_target_auroc, fill = trait)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0.5, linetype = 2, color = "gray40") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "CITE source scDRS score quality", x = NULL, y = "AUROC") +
  theme_bw() +
  theme(legend.position = "none")
ggsave(file.path(FIG_DIR, "cite_source_score_auroc.png"), p1,
       width = 5.5, height = 4, dpi = 150)

p2 <- ggplot(summary_dt, aes(x = trait, y = mean_source_target_weight_seed, fill = trait)) +
  geom_col(width = 0.7) +
  geom_point(aes(y = mean_source_target_weight_all), color = "black", size = 2) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "RNA target-source weight in ATAC seeds",
       subtitle = "black point = average among all ATAC cells",
       x = NULL, y = "weighted fraction") +
  theme_bw() +
  theme(legend.position = "none")
ggsave(file.path(FIG_DIR, "seed_source_target_weight.png"), p2,
       width = 6, height = 4, dpi = 150)

p3 <- melt(
  plot_dt,
  id.vars = c("trait", "group"),
  measure.vars = c("source_target_weight", "source_myeloid_weight"),
  variable.name = "source_fraction",
  value.name = "weighted_fraction"
)
p3 <- ggplot(p3, aes(x = group, y = weighted_fraction, fill = group)) +
  geom_boxplot(outlier.size = 0.2) +
  facet_grid(source_fraction ~ trait) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "RNA source composition after trait-aware mapping",
       x = NULL, y = "weighted fraction") +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 25, hjust = 1))
ggsave(file.path(FIG_DIR, "source_weight_by_group.png"), p3,
       width = 9, height = 6, dpi = 150)

cat("\n=== RNA source diagnostics finished ===\n")
print(summary_dt)
cat("\nOutputs are in: ", OUT_DIR, "\n", sep = "")
