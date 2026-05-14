# module4 diagnostic: scDRS FDR-based seed candidates
# -----------------------------------------------------------------------------
# Purpose:
#   Replace the crude pseudo-zRNA top 1% seed diagnostic with a scDRS
#   significance-based diagnostic.
#
# Method boundary:
#   - Seed candidates are selected only from scDRS p-values / scores and mapping
#     weights.
#   - Cell type is used only for post hoc diagnostics.
#   - This script does not replace the main pipeline yet.
#
# Run:
#   setwd("D:/RStudio/work/module4_traitaware_pseudopair")
#   source("04_evaluation/diagnose_scdrs_fdr_seed_candidates.R")
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
ATAC_META_FILE <- file.path(ROOT, "04_evaluation/inputs/granja_cell_meta.rds")
CITE_SCORE_RDS <- file.path(ROOT, "01_data_prep_external_rna/output/cite_rna_scdrs_scores.rds")
OUT_DIR <- file.path(ROOT, "04_evaluation/module4_v2_diagnostics/fdr_seed_candidates")
FIG_DIR <- file.path(OUT_DIR, "figures")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

TRAITS <- c("MONO", "PLT", "LYM")
FDR_THRESHOLDS <- c(0.05, 0.10, 0.20)
TOP_N <- 350L

STRATEGIES <- data.table(
  strategy = c("weight_gt0", "weight_ge0p25", "weight_ge0p50",
               "weight_ge0p80", "top350_positive"),
  source_weight_cutoff = c(0, 0.25, 0.50, 0.80, NA_real_),
  top_n = c(NA_integer_, NA_integer_, NA_integer_, NA_integer_, TOP_N)
)

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

stop_missing <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0L) {
    stop("Missing required file(s):\n  ", paste(missing, collapse = "\n  "),
         call. = FALSE)
  }
}

compute_auc <- function(score, label) {
  keep <- is.finite(score) & !is.na(label)
  score <- score[keep]
  label <- label[keep]
  if (sum(label == 1L) < 10L || sum(label == 0L) < 10L) return(NA_real_)
  if (length(unique(score)) < 2L) return(NA_real_)
  as.numeric(pROC::auc(pROC::roc(
    label, score, quiet = TRUE, direction = "<", levels = c(0L, 1L)
  )))
}

safe_div <- function(num, den) {
  if (!is.finite(den) || den <= 0) return(NA_real_)
  num / den
}

top_table_label <- function(x, n = 10L) {
  tab <- sort(table(x), decreasing = TRUE)
  if (length(tab) == 0L) return("")
  paste(paste0(names(head(tab, n)), ":", as.integer(head(tab, n))), collapse = "; ")
}

top_weight_label <- function(dt, type_col = "cite_cell_type",
                             weight_col = "neighbor_weight", n = 10L) {
  if (nrow(dt) == 0L) return("")
  s <- dt[, .(total_weight = sum(get(weight_col), na.rm = TRUE)),
          by = type_col][order(-total_weight)]
  if (nrow(s) == 0L) return("")
  paste(paste0(head(s[[type_col]], n), ":", round(head(s$total_weight, n), 1)),
        collapse = "; ")
}

select_candidate <- function(sig_source_weight, strategy, cutoff, top_n) {
  x <- sig_source_weight
  out <- rep(FALSE, length(x))
  if (strategy == "top350_positive") {
    ok <- which(is.finite(x) & x > 0)
    if (length(ok) == 0L) return(out)
    ord <- ok[order(x[ok], decreasing = TRUE)]
    keep <- head(ord, min(top_n, length(ord)))
    out[keep] <- TRUE
    return(out)
  }
  if (identical(cutoff, 0) || isTRUE(cutoff == 0)) {
    return(is.finite(x) & x > 0)
  }
  is.finite(x) & x >= cutoff
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

cat(">>> Diagnose scDRS FDR-based seed candidates\n\n")

stop_missing(c(ATAC_META_FILE, CITE_SCORE_RDS))

atac_meta <- load_atac_meta(ATAC_META_FILE)
cite_meta <- readRDS(CITE_SCORE_RDS)
cite_meta <- as.data.frame(cite_meta)
if (!all(c("cell", "cell_type") %in% colnames(cite_meta))) {
  stop("CITE score RDS must contain cell and cell_type columns.", call. = FALSE)
}
rownames(cite_meta) <- cite_meta$cell

summary_rows <- list()
cell_rows <- list()
atac_by_celltype_rows <- list()
rna_by_celltype_rows <- list()
source_rows <- list()

for (trait in TRAITS) {
  cat("=== Trait ", trait, " ===\n", sep = "")

  score_col <- paste0(trait, "_norm_score")
  pval_col <- paste0(trait, "_pval")
  if (!all(c(score_col, pval_col) %in% colnames(cite_meta))) {
    stop("Missing CITE score columns for ", trait, ": ",
         paste(setdiff(c(score_col, pval_col), colnames(cite_meta)), collapse = ", "),
         call. = FALSE)
  }

  map_file <- file.path(MAP_DIR, paste0("granja_cite_mapping_", trait, ".rds"))
  stop_missing(map_file)
  mapping <- readRDS(map_file)

  atac_bc <- mapping$atac_barcodes
  nn_bc <- mapping$rna_neighbor_barcode
  w <- mapping$neighbor_weight
  k <- ncol(nn_bc)

  if (!all(atac_bc %in% rownames(atac_meta))) {
    stop("Some ATAC barcodes are missing from ATAC meta for ", trait, call. = FALSE)
  }
  if (!all(as.vector(nn_bc) %in% rownames(cite_meta))) {
    stop("Some mapped CITE neighbors are missing from CITE meta for ", trait, call. = FALSE)
  }

  score <- as.numeric(cite_meta[[score_col]])
  pval <- as.numeric(cite_meta[[pval_col]])
  qval <- rep(NA_real_, length(pval))
  ok_p <- is.finite(pval) & !is.na(pval)
  qval[ok_p] <- p.adjust(pval[ok_p], method = "BH")

  atac_cell_type <- as.character(atac_meta[atac_bc, "cell_type"])
  atac_target <- atac_cell_type %in% ATAC_TARGETS[[trait]]
  atac_target_rate <- mean(atac_target)
  cite_target <- cite_meta$cell_type %in% CITE_TARGETS[[trait]]
  cite_myeloid <- cite_meta$cell_type %in% CITE_MYELOID
  cite_target_rate <- mean(cite_target, na.rm = TRUE)

  nn_ix <- match(as.vector(nn_bc), cite_meta$cell)
  nn_type <- matrix(cite_meta$cell_type[nn_ix], nrow = nrow(nn_bc), ncol = k)
  nn_target <- matrix(cite_target[nn_ix], nrow = nrow(nn_bc), ncol = k)
  nn_myeloid <- matrix(cite_myeloid[nn_ix], nrow = nrow(nn_bc), ncol = k)
  source_target_weight <- rowSums(w * nn_target, na.rm = TRUE)
  source_myeloid_weight <- rowSums(w * nn_myeloid, na.rm = TRUE)

  edge_dt_base <- data.table(
    trait = trait,
    atac_barcode = rep(atac_bc, times = k),
    rank = rep(seq_len(k), each = length(atac_bc)),
    cite_barcode = as.vector(nn_bc),
    neighbor_weight = as.vector(w),
    cite_cell_type = as.vector(nn_type)
  )

  for (fdr_threshold in FDR_THRESHOLDS) {
    sig_rna <- is.finite(qval) & qval <= fdr_threshold & is.finite(score) & score > 0

    nn_sig <- matrix(sig_rna[nn_ix], nrow = nrow(nn_bc), ncol = k)
    sig_source_weight <- rowSums(w * nn_sig, na.rm = TRUE)

    sig_source_auc <- compute_auc(sig_source_weight, as.integer(atac_target))
    sig_rna_target_precision <- safe_div(sum(sig_rna & cite_target, na.rm = TRUE),
                                         sum(sig_rna, na.rm = TRUE))
    sig_rna_myeloid_rate <- safe_div(sum(sig_rna & cite_myeloid, na.rm = TRUE),
                                     sum(sig_rna, na.rm = TRUE))

    rna_dt <- data.table(
      trait = trait,
      fdr_threshold = fdr_threshold,
      cell = cite_meta$cell,
      cell_type = cite_meta$cell_type,
      norm_score = score,
      pval = pval,
      qval = qval,
      is_sig_rna = sig_rna,
      is_target_posthoc = cite_target,
      is_myeloid_posthoc = cite_myeloid
    )

    rna_by_ct <- rna_dt[, .(
      n_cells = .N,
      n_sig_rna = sum(is_sig_rna, na.rm = TRUE),
      sig_rna_rate = mean(is_sig_rna, na.rm = TRUE),
      mean_norm_score = mean(norm_score, na.rm = TRUE),
      min_qval = suppressWarnings(min(qval, na.rm = TRUE)),
      is_target_celltype_posthoc = unique(is_target_posthoc)[1],
      is_myeloid_celltype_posthoc = unique(is_myeloid_posthoc)[1]
    ), by = .(trait, fdr_threshold, cell_type)][order(trait, fdr_threshold, -sig_rna_rate)]
    rna_by_celltype_rows[[paste(trait, fdr_threshold, sep = "_")]] <- rna_by_ct

    cell_dt <- data.table(
      trait = trait,
      fdr_threshold = fdr_threshold,
      barcode = atac_bc,
      atac_cell_type_posthoc = atac_cell_type,
      is_atac_target_posthoc = atac_target,
      sig_source_weight = sig_source_weight,
      source_target_weight_posthoc = source_target_weight,
      source_myeloid_weight_posthoc = source_myeloid_weight,
      no_celltype_used_to_select = TRUE
    )

    for (si in seq_len(nrow(STRATEGIES))) {
      strategy <- STRATEGIES$strategy[si]
      cutoff <- STRATEGIES$source_weight_cutoff[si]
      top_n <- STRATEGIES$top_n[si]
      candidate <- select_candidate(sig_source_weight, strategy, cutoff, top_n)
      cell_dt[[paste0("candidate_", strategy)]] <- candidate

      n_candidate <- sum(candidate)
      n_candidate_target <- sum(candidate & atac_target)
      candidate_precision <- safe_div(n_candidate_target, n_candidate)
      candidate_recall <- safe_div(n_candidate_target, sum(atac_target))
      candidate_enrichment <- candidate_precision / max(atac_target_rate, 1e-12)

      edge_dt <- copy(edge_dt_base)
      edge_dt[, is_candidate_atac := candidate[match(atac_barcode, atac_bc)]]
      top_seed_sources <- top_weight_label(edge_dt[is_candidate_atac == TRUE])

      summary_rows[[paste(trait, fdr_threshold, strategy, sep = "_")]] <- data.frame(
        trait = trait,
        fdr_threshold = fdr_threshold,
        strategy = strategy,
        source_weight_cutoff = cutoff,
        top_n = top_n,
        n_cite_cells = nrow(cite_meta),
        cite_target_rate = cite_target_rate,
        n_sig_rna = sum(sig_rna, na.rm = TRUE),
        sig_rna_rate = mean(sig_rna, na.rm = TRUE),
        sig_rna_target_precision_posthoc = sig_rna_target_precision,
        sig_rna_target_enrichment_posthoc =
          sig_rna_target_precision / max(cite_target_rate, 1e-12),
        sig_rna_myeloid_rate_posthoc = sig_rna_myeloid_rate,
        sig_rna_top_celltypes_posthoc =
          top_table_label(cite_meta$cell_type[sig_rna], n = 10L),
        n_atac_cells = length(atac_bc),
        atac_target_rate = atac_target_rate,
        sig_source_weight_auc_posthoc = sig_source_auc,
        mean_sig_source_weight_all = mean(sig_source_weight, na.rm = TRUE),
        mean_sig_source_weight_atac_target = mean(sig_source_weight[atac_target], na.rm = TRUE),
        mean_sig_source_weight_atac_rest = mean(sig_source_weight[!atac_target], na.rm = TRUE),
        n_candidate_atac = n_candidate,
        candidate_rate = n_candidate / length(atac_bc),
        n_candidate_target_posthoc = n_candidate_target,
        candidate_precision_posthoc = candidate_precision,
        candidate_recall_posthoc = candidate_recall,
        candidate_enrichment_posthoc = candidate_enrichment,
        candidate_mean_source_target_weight_posthoc =
          mean(source_target_weight[candidate], na.rm = TRUE),
        candidate_mean_source_myeloid_weight_posthoc =
          mean(source_myeloid_weight[candidate], na.rm = TRUE),
        candidate_top_atac_celltypes_posthoc =
          top_table_label(atac_cell_type[candidate], n = 10L),
        candidate_top_rna_sources_by_weight_posthoc = top_seed_sources,
        no_celltype_used_to_select = TRUE,
        stringsAsFactors = FALSE
      )

      by_ct <- data.table(
        trait = trait,
        fdr_threshold = fdr_threshold,
        strategy = strategy,
        atac_cell_type = atac_cell_type,
        is_target_celltype_posthoc = atac_target,
        is_candidate = candidate,
        sig_source_weight = sig_source_weight,
        source_target_weight_posthoc = source_target_weight,
        source_myeloid_weight_posthoc = source_myeloid_weight
      )[, .(
        n_cells = .N,
        n_candidate = sum(is_candidate),
        candidate_rate = mean(is_candidate),
        mean_sig_source_weight = mean(sig_source_weight, na.rm = TRUE),
        mean_source_target_weight_posthoc = mean(source_target_weight_posthoc, na.rm = TRUE),
        mean_source_myeloid_weight_posthoc = mean(source_myeloid_weight_posthoc, na.rm = TRUE),
        is_target_celltype_posthoc = unique(is_target_celltype_posthoc)[1]
      ), by = .(trait, fdr_threshold, strategy, atac_cell_type)][
        order(trait, fdr_threshold, strategy, -candidate_rate, -n_candidate)
      ]
      atac_by_celltype_rows[[paste(trait, fdr_threshold, strategy, sep = "_")]] <- by_ct

      src <- edge_dt[is_candidate_atac == TRUE, .(
        total_weight = sum(neighbor_weight, na.rm = TRUE),
        n_edges = .N
      ), by = .(cite_cell_type)][order(-total_weight)]
      src[, `:=`(
        trait = trait,
        fdr_threshold = fdr_threshold,
        strategy = strategy
      )]
      setcolorder(src, c("trait", "fdr_threshold", "strategy",
                         "cite_cell_type", "total_weight", "n_edges"))
      source_rows[[paste(trait, fdr_threshold, strategy, sep = "_")]] <- src
    }

    cell_rows[[paste(trait, fdr_threshold, sep = "_")]] <- cell_dt
  }

  cat("  min qval: ", signif(min(qval, na.rm = TRUE), 4),
      "; q<=0.10 & score>0: ", sum(qval <= 0.10 & score > 0, na.rm = TRUE),
      "; source-weight AUROC at FDR 0.10 will be in summary\n", sep = "")
}

summary_dt <- rbindlist(summary_rows, use.names = TRUE, fill = TRUE)
cell_dt <- rbindlist(cell_rows, use.names = TRUE, fill = TRUE)
atac_by_ct_dt <- rbindlist(atac_by_celltype_rows, use.names = TRUE, fill = TRUE)
rna_by_ct_dt <- rbindlist(rna_by_celltype_rows, use.names = TRUE, fill = TRUE)
source_dt <- rbindlist(source_rows, use.names = TRUE, fill = TRUE)

fwrite(summary_dt, file.path(OUT_DIR, "fdr_seed_summary.csv"))
fwrite(cell_dt, file.path(OUT_DIR, "fdr_seed_cell_level.csv.gz"))
fwrite(atac_by_ct_dt, file.path(OUT_DIR, "fdr_seed_by_atac_celltype_posthoc.csv"))
fwrite(rna_by_ct_dt, file.path(OUT_DIR, "fdr_sig_rna_by_celltype_posthoc.csv"))
fwrite(source_dt, file.path(OUT_DIR, "fdr_seed_source_by_cite_celltype_posthoc.csv"))

plot_dt <- summary_dt[strategy %in% c("weight_gt0", "weight_ge0p50", "top350_positive")]

p1_dt <- plot_dt[is.finite(candidate_enrichment_posthoc)]
p1_line_dt <- p1_dt[, if (.N >= 2L) .SD, by = .(trait, strategy)]
p1 <- ggplot(p1_dt, aes(x = factor(fdr_threshold), y = candidate_enrichment_posthoc,
                          color = strategy, group = strategy)) +
  geom_hline(yintercept = 1, linetype = 2, color = "gray40") +
  geom_point(size = 2) +
  geom_line(data = p1_line_dt, linewidth = 0.7) +
  facet_wrap(~ trait, scales = "free_y") +
  labs(title = "Post hoc enrichment of FDR-based ATAC seed candidates",
       x = "RNA scDRS FDR threshold", y = "seed enrichment vs ATAC target rate",
       color = "strategy") +
  theme_bw()
ggsave(file.path(FIG_DIR, "fdr_seed_enrichment.png"), p1,
       width = 8, height = 4.5, dpi = 150)

p2_dt <- unique(summary_dt[, .(trait, fdr_threshold, sig_source_weight_auc_posthoc)])
p2_dt <- p2_dt[is.finite(sig_source_weight_auc_posthoc)]
p2_line_dt <- p2_dt[, if (.N >= 2L) .SD, by = trait]
p2 <- ggplot(p2_dt,
             aes(x = factor(fdr_threshold), y = sig_source_weight_auc_posthoc,
                 color = trait, group = trait)) +
  geom_hline(yintercept = 0.5, linetype = 2, color = "gray40") +
  geom_point(size = 2) +
  geom_line(data = p2_line_dt, linewidth = 0.7) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "ATAC AUROC of transferred significant-RNA source weight",
       x = "RNA scDRS FDR threshold", y = "AUROC", color = "trait") +
  theme_bw()
ggsave(file.path(FIG_DIR, "fdr_sig_source_weight_auroc.png"), p2,
       width = 6.5, height = 4, dpi = 150)

p3_dt <- plot_dt[is.finite(candidate_mean_source_myeloid_weight_posthoc)]
p3_line_dt <- p3_dt[, if (.N >= 2L) .SD, by = .(trait, strategy)]
p3 <- ggplot(p3_dt, aes(x = factor(fdr_threshold),
                          y = candidate_mean_source_myeloid_weight_posthoc,
                          color = strategy, group = strategy)) +
  geom_point(size = 2) +
  geom_line(data = p3_line_dt, linewidth = 0.7) +
  facet_wrap(~ trait, scales = "free_y") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Post hoc myeloid RNA-source weight in seed candidates",
       x = "RNA scDRS FDR threshold", y = "mean myeloid source weight",
       color = "strategy") +
  theme_bw()
ggsave(file.path(FIG_DIR, "fdr_seed_myeloid_source_weight.png"), p3,
       width = 8, height = 4.5, dpi = 150)

cat("\n=== FDR seed candidate diagnostics finished ===\n")
print(summary_dt[strategy == "top350_positive" & fdr_threshold == 0.10,
                 .(trait, fdr_threshold, strategy, n_sig_rna,
                   sig_rna_target_precision_posthoc,
                   sig_source_weight_auc_posthoc,
                   n_candidate_atac,
                   candidate_precision_posthoc,
                   candidate_enrichment_posthoc,
                   candidate_mean_source_myeloid_weight_posthoc)])
cat("\nOutputs are in: ", OUT_DIR, "\n", sep = "")
