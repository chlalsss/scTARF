# 04b_seed_dual.R
# §2.D' 双端并集种子: ATAC z_score top 1% (来自 baseline) ∪ SCDRS norm_score top 1%.
#
# motivation: §2.D 当前 MultiOmicsFusion 种子只用 SCDRS norm_score, ATAC 端 GWAS 变异信号没进种子.
# 实测两端 top 1% Jaccard 重合度 < 1% (4-8 个交集细胞), 说明两端互补不冗余, 并集能利用两端信息.
# 交集太少 (4-8 个) RWR 无法跑, 必须并集.
#
# 输入:
#   - SCDRS norm_score: <SCORE_DIR>/<trait>.score.gz
#   - ATAC z_score:     <BASE_DIR>/SCAVENGE_baseline_<trait>.tsv.gz
#
# 输出:
#   - <OUT_DIR>/seed_dual_<trait>.rds (named numeric, 1370 cells/trait 量级, 1/0)
#   - <OUT_DIR>/seed_dual_summary.csv (n_atac, n_rna, n_inter, n_union)

suppressPackageStartupMessages({
  library(data.table)
})

OUT_DIR   <- "D:/SCAVENGE-reproducibility-main/BMMC/04_dual_modality"
BASE_DIR  <- "D:/SCAVENGE-reproducibility-main/BMMC/03_scavenge_baseline"
SCORE_DIR <- "D:/RStudio/work/module3_multiome/04_dual_modality/01_scdrs_python/output"

TRAITS <- c(
  MONO = "UKB_460K.blood_MONOCYTE_COUNT",
  PLT  = "UKB_460K.blood_PLATELET_COUNT",
  LYM  = "UKB_460K.blood_LYMPHOCYTE_COUNT"
)
SEED_TOP_PCT <- 0.01

# common cells = rna_z_block 的 cells (SCDRS-aligned 68303), 跟 04_seed_scdrs.R 一致
rna_z_ref <- readRDS(file.path(OUT_DIR, "rna_z_block_MultiOmicsFusion_MONO.rds"))
common <- rownames(rna_z_ref)
cat("common cells:", length(common), "\n\n")

seed_summary <- data.frame()

for (tname in names(TRAITS)) {
  trait <- TRAITS[[tname]]
  cat("=== ", tname, " ===\n", sep = "")

  # --- RNA 端: SCDRS norm_score top 1% ---
  dt <- fread(file.path(SCORE_DIR, paste0(trait, ".score.gz")))
  setnames(dt, 1, "barcode")
  setkey(dt, barcode)
  norm_score <- dt[common, norm_score]
  names(norm_score) <- common
  stopifnot(!any(is.na(norm_score)))

  thr_rna  <- as.numeric(quantile(norm_score, 1 - SEED_TOP_PCT))
  rna_seed <- common[norm_score >= thr_rna]
  cat("  SCDRS top 1% (RNA): ", length(rna_seed), " cells (thr=",
      sprintf("%.3f", thr_rna), ")\n", sep = "")

  # --- ATAC 端: baseline z_score top 1% ---
  bf <- fread(file.path(BASE_DIR, paste0("SCAVENGE_baseline_", tname, ".tsv.gz")))
  # baseline 的 cell 集合可能跟 SCDRS-aligned common 不一致, 取交集后再算阈值
  bf <- bf[cell %in% common]
  thr_atac <- as.numeric(quantile(bf$z_score, 1 - SEED_TOP_PCT))
  atac_seed <- bf$cell[bf$z_score >= thr_atac]
  cat("  ATAC top 1% (z_score): ", length(atac_seed),
      " cells (thr=", sprintf("%.3f", thr_atac), ")\n", sep = "")

  # --- 并集 ---
  inter <- intersect(rna_seed, atac_seed)
  uni   <- union(rna_seed, atac_seed)
  cat("  intersect (AND): ", length(inter), " cells\n", sep = "")
  cat("  union (OR):      ", length(uni), " cells\n", sep = "")

  # 写 seed 向量 (1/0)
  seed_vec <- numeric(length(common))
  names(seed_vec) <- common
  seed_vec[uni] <- 1
  saveRDS(seed_vec, file.path(OUT_DIR, paste0("seed_dual_", tname, ".rds")))
  cat("  saved seed_dual_", tname, ".rds (n_seed=", sum(seed_vec), ")\n", sep = "")

  seed_summary <- rbind(seed_summary, data.frame(
    trait      = tname,
    n_atac_top1pct = length(atac_seed),
    n_rna_top1pct  = length(rna_seed),
    n_intersect    = length(inter),
    n_union        = length(uni),
    jaccard        = round(length(inter) / length(uni), 4),
    thr_rna_norm_score = round(thr_rna, 3),
    thr_atac_z_score   = round(thr_atac, 3)
  ))
}

write.csv(seed_summary,
          file.path(OUT_DIR, "seed_dual_summary.csv"), row.names = FALSE)
cat("\n")
print(seed_summary)
cat("\n04b_seed_dual.R OK.\n")
