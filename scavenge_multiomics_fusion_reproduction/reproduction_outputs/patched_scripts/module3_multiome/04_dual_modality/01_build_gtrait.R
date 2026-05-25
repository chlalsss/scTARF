# 04_dual_modality/01_build_gtrait.R
#
# §2.2 G_trait 构建. variant -> peak -> gene 完整链路, 每 trait 输出
# top-N (kneedle on bracket [100, 1000]) 个基因带 [0,1] 归一化权重 w_g.
# 完全自包含: 不复用 module2_gwas_gene_risk 旧产物 (旧值 d0=50kb 与本方案
# d0=5kb 不一致, 重写而非复用).
#
# Inputs:
#   paths$gwas_tsv                                   fine-mapped variants (variant/trait/PP)
#   <out>/03_scavenge_baseline/SE_atac_gcbiased.rds  peak GRanges (rowRanges)
#   <out>/01_data_qc/BMMC_rna_features.rds           gene_symbol whitelist
#   TxDb.Hsapiens.UCSC.hg38.knownGene + org.Hs.eg.db gene TSS (strand-aware)
#
# Outputs:
#   <out>/04_dual_modality/Gtrait_<trait>.rds        per-trait dt(gene_symbol, gene_score, w_g, rank)
#   <out>/04_dual_modality/Gtrait_summary.csv        n_genes / kneedle / top genes
#   <out>/04_dual_modality/figures/Gtrait_score_curves.png  normalized score curves
#
# Red line (§2.10): 本脚本输入只用 GWAS variants + peak GRanges + gene TSS +
# RNA gene 白名单, 不读 BMMC_cell_meta.rds 或任何 BioClassification 字段.
# 手动审计:  grep -nE "BMMC_cell_meta|cell_type|BioClassification" 01_build_gtrait.R
#           只能命中 header 的红线说明 + 此处审计提示, 不能命中执行行.

source("D:/RStudio/work/scavenge_multiomics_fusion_reproduction/reproduction_outputs/patched_scripts/module3_multiome/_bootstrap.R")
required_pkgs(c(
  "data.table", "GenomicRanges", "IRanges", "S4Vectors",
  "SummarizedExperiment", "GenomicFeatures", "AnnotationDbi",
  "GenomeInfoDb",
  "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db",
  "ggplot2"
))

paths   <- m3_paths()
out_dir <- file.path(paths$out_dir, "04_dual_modality")
fig_dir <- file.path(out_dir, "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# --- locked params --------------------------------------------------------
D0        <- 5000L            # peak->gene exp decay (SCAVENGE / ABC default)
CUTOFF    <- 100000L          # peak->gene distance hard cutoff
N_BRACKET <- c(100L, 1000L)   # kneedle clamp on top-N
HEAD_CAP  <- 5000L            # only consider top HEAD_CAP genes for kneedle
TRAITS    <- c("MONO", "PLT", "LYM")

# --- 1. GWAS variants -> GRanges -----------------------------------------

if (!file.exists(paths$gwas_tsv)) stop("Missing GWAS TSV: ", paths$gwas_tsv)
gwas <- data.table::fread(paths$gwas_tsv)
if ("var" %in% colnames(gwas) && !"variant" %in% colnames(gwas)) {
  data.table::setnames(gwas, "var", "variant")
}
gwas <- gwas[trait %in% TRAITS]
message("GWAS rows after MONO/PLT/LYM filter: ", nrow(gwas))

bed_dt <- parse_var_to_bed(gwas$variant)
gwas   <- gwas[bed_dt$ok]
bed_dt <- bed_dt[bed_dt$ok]
gwas[, `:=`(chr = bed_dt$chr, start = bed_dt$start, end = bed_dt$end)]

# parse_var_to_bed -> 0-based BED; GRanges is 1-based (start = pos).
variant_gr <- GenomicRanges::GRanges(
  seqnames = gwas$chr,
  ranges   = IRanges::IRanges(start = gwas$start + 1L, width = 1L),
  trait    = gwas$trait,
  PP       = as.numeric(gwas$PP)
)
message("variant_gr: ", length(variant_gr), " variants across ",
        length(unique(as.character(GenomicRanges::seqnames(variant_gr)))),
        " chroms")

# --- 2. peak GRanges + RNA whitelist -------------------------------------

SE <- readRDS(file.path(paths$out_dir, "03_scavenge_baseline",
                        "SE_atac_gcbiased.rds"))
peak_gr <- SummarizedExperiment::rowRanges(SE)
message("peak_gr: ", length(peak_gr), " peaks")

rna_feat <- readRDS(file.path(paths$out_dir, "01_data_qc",
                              "BMMC_rna_features.rds"))
rna_symbols <- unique(as.character(rna_feat$gene_symbol))
message("RNA gene whitelist: ", length(rna_symbols), " symbols")

# --- 3. variant -> peak risk ---------------------------------------------

ov <- GenomicRanges::findOverlaps(variant_gr, peak_gr)
hit_dt <- data.table::data.table(
  peak_idx = S4Vectors::subjectHits(ov),
  trait    = variant_gr$trait[S4Vectors::queryHits(ov)],
  PP       = variant_gr$PP[S4Vectors::queryHits(ov)]
)
peak_risk_dt <- hit_dt[, .(risk = sum(PP)), by = .(peak_idx, trait)]
peak_risk_dt <- peak_risk_dt[risk > 0]

for (tr in TRAITS) {
  rs <- peak_risk_dt[trait == tr, risk]
  if (length(rs) == 0) {
    message("  ", tr, ": 0 peaks with risk"); next
  }
  message(sprintf(
    "  %s: %d peaks with risk > 0 | risk q50=%.3f q90=%.3f max=%.3f",
    tr, length(rs), stats::median(rs),
    as.numeric(stats::quantile(rs, 0.9)), max(rs)
  ))
}

peak_idx_used <- sort(unique(peak_risk_dt$peak_idx))
peak_used_gr  <- peak_gr[peak_idx_used]
peak_center   <- as.integer(
  (GenomicRanges::start(peak_used_gr) + GenomicRanges::end(peak_used_gr)) %/% 2L
)
peak_center_gr <- GenomicRanges::GRanges(
  seqnames = GenomicRanges::seqnames(peak_used_gr),
  ranges   = IRanges::IRanges(peak_center, peak_center)
)
S4Vectors::mcols(peak_center_gr)$peak_idx <- peak_idx_used

# --- 4. gene TSS GRanges (whitelist-restricted) --------------------------

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
g <- GenomicFeatures::genes(txdb, single.strand.genes.only = TRUE)

sym <- AnnotationDbi::mapIds(
  org.Hs.eg.db::org.Hs.eg.db,
  keys = names(g), keytype = "ENTREZID", column = "SYMBOL",
  multiVals = "first"
)
keep <- !is.na(sym) & sym %in% rna_symbols
g    <- g[keep]
sym  <- sym[keep]

# strand-aware TSS
tss_pos <- ifelse(
  as.character(GenomicRanges::strand(g)) == "+",
  GenomicRanges::start(g), GenomicRanges::end(g)
)
tss_gr <- GenomicRanges::GRanges(
  seqnames    = GenomicRanges::seqnames(g),
  ranges      = IRanges::IRanges(tss_pos, tss_pos),
  strand      = GenomicRanges::strand(g),
  gene_symbol = unname(sym)
)
tss_gr <- GenomeInfoDb::keepStandardChromosomes(tss_gr,
                                                pruning.mode = "coarse")
# 1 symbol = 1 TSS; mapIds 'first' already picked one entrez per symbol,
# but defensive dedup in case of edge cases.
tss_gr <- tss_gr[!duplicated(S4Vectors::mcols(tss_gr)$gene_symbol)]
message("tss_gr after symbol-whitelist + dedup: ", length(tss_gr), " genes")

# --- 5. peak -> gene link (geometry only, trait-agnostic) ----------------

peak_window <- GenomicRanges::resize(
  peak_center_gr, width = 2L * CUTOFF + 1L, fix = "center"
)
ov_pg <- GenomicRanges::findOverlaps(peak_window, tss_gr)

link_dt <- data.table::data.table(
  peak_idx    = S4Vectors::mcols(peak_center_gr)$peak_idx[
                  S4Vectors::queryHits(ov_pg)],
  gene_symbol = S4Vectors::mcols(tss_gr)$gene_symbol[
                  S4Vectors::subjectHits(ov_pg)],
  d           = abs(
    GenomicRanges::start(peak_center_gr)[S4Vectors::queryHits(ov_pg)] -
    GenomicRanges::start(tss_gr)[S4Vectors::subjectHits(ov_pg)]
  )
)
link_dt <- link_dt[d <= CUTOFF]
link_dt[, weight := exp(-d / D0)]
message("peak-gene links within cutoff: ", nrow(link_dt))

# --- 6. gene_score per trait ---------------------------------------------

score_dt <- merge(link_dt, peak_risk_dt, by = "peak_idx",
                  allow.cartesian = TRUE)
gene_score_dt <- score_dt[, .(gene_score = sum(risk * weight)),
                          by = .(trait, gene_symbol)]
data.table::setorder(gene_score_dt, trait, -gene_score)

# --- 7. kneedle top-N + normalize ---------------------------------------

# Standard kneedle (Satopaa 2011) for monotone-decreasing L-curve mapped to
# the unit square [0,1]^2 with endpoints anchored at (0,1) and (1,0). The
# curve lies below the chord y = 1 - x; the elbow is the rank with maximum
# vertical drop from the chord, argmax_i (1 - x_i - y_i). y is normalized
# end-to-end so the curve actually touches both anchors, otherwise tail
# offset can pull which.max to the last index.
kneedle_pick <- function(scores, bracket = N_BRACKET, head_cap = HEAD_CAP) {
  s <- sort(scores, decreasing = TRUE)
  s <- s[seq_len(min(length(s), head_cap))]
  n <- length(s)
  if (n < bracket[1]) {
    return(list(idx_raw = n, idx_clamped = n, s = s))
  }
  x <- (seq_len(n) - 1) / (n - 1)
  s_min <- min(s); s_max <- max(s)
  y <- if (s_max > s_min) (s - s_min) / (s_max - s_min) else rep(1, n)
  knee <- which.max(1 - x - y)
  idx_clamped <- max(bracket[1], min(bracket[2], knee))
  list(idx_raw = knee, idx_clamped = idx_clamped, s = s)
}

summary_rows <- list()
curve_rows   <- list()

for (tr in TRAITS) {
  st <- gene_score_dt[trait == tr]
  if (nrow(st) == 0) {
    warning("No gene_score rows for trait: ", tr); next
  }
  k       <- kneedle_pick(st$gene_score)
  N_pick  <- k$idx_clamped
  pick    <- st[seq_len(N_pick)]
  pick[, w_g  := gene_score / max(gene_score)]
  pick[, rank := seq_len(.N)]

  out_file <- file.path(out_dir, paste0("Gtrait_", tr, ".rds"))
  saveRDS(pick, out_file)
  message(sprintf(
    "  %s: kneedle_raw=%d, picked=%d, score_max=%.4g, top1=%s",
    tr, k$idx_raw, N_pick, max(pick$gene_score), pick$gene_symbol[1]
  ))

  summary_rows[[tr]] <- data.table::data.table(
    trait              = tr,
    n_genes_picked     = N_pick,
    score_max          = max(pick$gene_score),
    kneedle_x_raw      = k$idx_raw,
    kneedle_x_clamped  = N_pick,
    top1_gene          = pick$gene_symbol[1],
    top5_genes         = paste(
      pick$gene_symbol[seq_len(min(5L, N_pick))], collapse = ","
    )
  )

  curve_n <- length(k$s)
  curve_rows[[tr]] <- data.table::data.table(
    trait      = tr,
    rank       = seq_len(curve_n),
    score      = k$s,
    score_norm = k$s / max(k$s),
    knee_at    = N_pick
  )
}

summary_dt <- data.table::rbindlist(summary_rows)
data.table::fwrite(summary_dt, file.path(out_dir, "Gtrait_summary.csv"))
print(summary_dt)

curve_dt <- data.table::rbindlist(curve_rows)
p <- ggplot2::ggplot(curve_dt,
                     ggplot2::aes(x = rank, y = score_norm)) +
  ggplot2::geom_line() +
  ggplot2::geom_vline(
    data = unique(curve_dt[, .(trait, knee_at)]),
    ggplot2::aes(xintercept = knee_at),
    color = "red", linetype = 2
  ) +
  ggplot2::facet_wrap(~ trait, scales = "free_x") +
  ggplot2::labs(
    title = "G_trait gene_score (normalized) with kneedle pick",
    x = "Rank", y = "score / max(score)"
  ) +
  ggplot2::theme_bw()
ggplot2::ggsave(
  file.path(fig_dir, "Gtrait_score_curves.png"),
  p, width = 9, height = 3, dpi = 150
)

message("\nG_trait built. Outputs in: ", out_dir)
