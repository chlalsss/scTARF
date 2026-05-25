# Diagnose whether v3 top5-significant completion scores carry target signal.
# This does not build a graph or rerun RWR. It only evaluates cell-level
# diagnostics already written by 02_trait_mapping/07_build_top5_significant_completion.R.

ROOT <- "D:/RStudio/work/module4_traitaware_pseudopair"

if (!exists("V3_SCORE_TAGS")) {
  output_root <- file.path(ROOT, "02_trait_mapping/output")
  tags <- basename(list.dirs(output_root, recursive = FALSE, full.names = TRUE))
  V3_SCORE_TAGS <- sort(tags[grepl("^v3_.*top5.*significant", tags)])
}

if (!exists("V3_SCORE_TRAITS")) V3_SCORE_TRAITS <- c("MONO")
if (!exists("V3_SCORE_OUT_DIR")) {
  V3_SCORE_OUT_DIR <- file.path(ROOT, "04_evaluation/v3_top5_significant_diagnostics")
}

META_FILE <- file.path(ROOT, "04_evaluation/inputs/granja_cell_meta.rds")

TARGETS <- list(
  MONO = c("11_CD14.Mono.1", "12_CD14.Mono.2"),
  PLT  = c("02_Early.Eryth", "03_Late.Eryth", "04_Early.Baso"),
  LYM  = c("16_Pre.B", "17_B", "18_Plasma",
           "19_CD8.N", "23_CD8.EM", "24_CD8.CM",
           "20_CD4.N1", "21_CD4.N2", "22_CD4.M",
           "25_NK")
)

SCORE_COLUMNS <- c(
  "retained_gene_count",
  "retained_gene_fraction",
  "top5_real_expression_sum",
  "top5_significant_expression_sum",
  "retained_expression_fraction"
)

dir.create(V3_SCORE_OUT_DIR, recursive = TRUE, showWarnings = FALSE)

compute_auc_rank <- function(score, label) {
  keep <- is.finite(score) & !is.na(label)
  score <- score[keep]
  label <- as.integer(label[keep])
  n_pos <- sum(label == 1L)
  n_neg <- sum(label == 0L)
  if (n_pos < 10L || n_neg < 10L) return(NA_real_)
  if (!is.finite(stats::var(score)) || stats::var(score) < 1e-12) return(NA_real_)
  r <- rank(score, ties.method = "average")
  (sum(r[label == 1L]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

safe_div <- function(a, b) ifelse(b == 0, NA_real_, a / b)

top_metrics <- function(score, label, top_fraction) {
  keep <- is.finite(score) & !is.na(label)
  score <- score[keep]
  label <- as.integer(label[keep])
  n <- length(score)
  n_top <- max(1L, ceiling(n * top_fraction))
  idx <- order(score, decreasing = TRUE)[seq_len(n_top)]
  precision <- mean(label[idx] == 1L)
  target_rate <- mean(label == 1L)
  c(
    n_top = n_top,
    precision = precision,
    enrichment = precision / max(target_rate, 1e-12)
  )
}

meta <- readRDS(META_FILE)
if (is.null(rownames(meta)) && "barcode" %in% names(meta)) rownames(meta) <- meta$barcode
if (!"cell_type" %in% names(meta)) stop("meta missing cell_type")

rows <- list()

for (tag in V3_SCORE_TAGS) {
  for (trait_name in V3_SCORE_TRAITS) {
    diag_file <- file.path(
      ROOT,
      "02_trait_mapping/output",
      tag,
      paste0(tolower(trait_name), "_completion_diagnostics"),
      paste0("top5_significant_cell_diagnostics_", trait_name, ".csv.gz")
    )
    if (!file.exists(diag_file)) next

    dt <- read.csv(gzfile(diag_file), stringsAsFactors = FALSE, check.names = FALSE)
    common <- intersect(dt$barcode, rownames(meta))
    if (length(common) == 0L) next

    dt <- dt[match(common, dt$barcode), , drop = FALSE]
    cell_type <- meta[common, "cell_type"]
    label <- as.integer(cell_type %in% TARGETS[[trait_name]])
    target_rate <- mean(label == 1L)

    for (score_col in intersect(SCORE_COLUMNS, names(dt))) {
      score <- dt[[score_col]]
      top1 <- top_metrics(score, label, 0.01)
      top5 <- top_metrics(score, label, 0.05)
      top10 <- top_metrics(score, label, 0.10)

      rows[[length(rows) + 1L]] <- data.frame(
        output_tag = tag,
        trait = trait_name,
        score = score_col,
        n_cells = length(label),
        n_positive = sum(label == 1L),
        target_rate = target_rate,
        auroc = compute_auc_rank(score, label),
        mean_target = mean(score[label == 1L], na.rm = TRUE),
        mean_rest = mean(score[label == 0L], na.rm = TRUE),
        delta_target_minus_rest =
          mean(score[label == 1L], na.rm = TRUE) - mean(score[label == 0L], na.rm = TRUE),
        top1pct_n = unname(top1["n_top"]),
        top1pct_precision = unname(top1["precision"]),
        top1pct_enrichment = unname(top1["enrichment"]),
        top5pct_n = unname(top5["n_top"]),
        top5pct_precision = unname(top5["precision"]),
        top5pct_enrichment = unname(top5["enrichment"]),
        top10pct_n = unname(top10["n_top"]),
        top10pct_precision = unname(top10["precision"]),
        top10pct_enrichment = unname(top10["enrichment"]),
        source_file = diag_file,
        stringsAsFactors = FALSE
      )
    }
  }
}

score_summary <- if (length(rows) > 0L) do.call(rbind, rows) else data.frame()

out_csv <- file.path(V3_SCORE_OUT_DIR, "v3_completion_score_diagnostics.csv")
write.csv(score_summary, out_csv, row.names = FALSE)

cat("\n=== v3 completion score diagnostics ===\n")
cat("Tags scanned: ", paste(V3_SCORE_TAGS, collapse = ", "), "\n", sep = "")
cat("Output: ", out_csv, "\n", sep = "")

if (nrow(score_summary) > 0L) {
  print(score_summary[order(score_summary$output_tag, -score_summary$auroc),
                      c("output_tag", "trait", "score", "auroc",
                        "top5pct_precision", "top5pct_enrichment",
                        "mean_target", "mean_rest")],
        row.names = FALSE)
}
