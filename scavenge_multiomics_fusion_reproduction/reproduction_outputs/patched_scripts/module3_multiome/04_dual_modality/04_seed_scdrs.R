# 04_seed_scdrs.R
# §2.D: SCDRS norm_score top 1% as seed cells (per trait, atac_w-INDEPENDENT).
# 直接 source(), 跑一次即可 (atac_w sweep 也用同一份 seed).
#
# 脚本功能：
# 读取每个 trait 的 scDRS 分数表，在共同细胞集合中按 norm_score 排序。
# 对每个 trait 取 top 1% 细胞，保存为 0/1 seed 向量。
# seed 向量的 names 是细胞 barcode，后续 RWR 会依赖这个顺序和 KNN 矩阵对齐。

suppressPackageStartupMessages({
  library(data.table)
})

OUT_DIR   <- "D:/RStudio/work/scavenge_multiomics_fusion_reproduction/reproduction_outputs/BMMC/04_dual_modality"
SCORE_DIR <- "D:/RStudio/work/scavenge_multiomics_fusion_reproduction/reproduction_outputs/patched_scripts/module3_multiome/04_dual_modality/01_scdrs_python/output"

# TRAITS 是本实验关注的三个血液性状。
# 左边 MONO/PLT/LYM 是代码内部简写，右边是 scDRS 输出文件使用的性状名。
TRAITS <- c(
  MONO = "UKB_460K.blood_MONOCYTE_COUNT",
  PLT  = "UKB_460K.blood_PLATELET_COUNT",
  LYM  = "UKB_460K.blood_LYMPHOCYTE_COUNT"
)

# 每个性状只取 scDRS norm_score 排名前 1% 的细胞作为 seed。
# 这样 seed 数量固定、解释简单，也避免把过多弱相关细胞放进传播起点。
SEED_TOP_PCT <- 0.01

# common cells = rna_z_block 的 cells (already aligned to SCDRS-filtered set, 68303)
# atac_z_block_MultiOmicsFusion.rds 是 align 前的 69249, 不能直接当 common
# 用 RNA block 的 rownames 定义共同细胞集合，是因为该对象已经和 scDRS 结果对齐。
# 后续 seed、KNN、RWR 都以这个共同细胞集合为准，避免不同矩阵细胞顺序不一致。
rna_z_ref <- readRDS(file.path(OUT_DIR, "rna_z_block_MultiOmicsFusion_MONO.rds"))
common <- rownames(rna_z_ref)
cat("common cells:", length(common), "\n\n")

seed_summary <- data.frame()
for (tname in names(TRAITS)) {
  trait <- TRAITS[[tname]]
  cat("=== ", tname, " ===\n", sep = "")

  # 读取当前 trait 的 scDRS 分数表。
  # 第一列重命名为 barcode，便于按细胞条形码和 common cells 对齐。
  dt <- fread(file.path(SCORE_DIR, paste0(trait, ".score.gz")))
  setnames(dt, 1, "barcode")
  setkey(dt, barcode)

  # 只取共同细胞的 norm_score，并把细胞名作为向量 names。
  # stopifnot 用来保证没有对齐失败产生的 NA。
  norm_score <- dt[common, norm_score]
  names(norm_score) <- common
  stopifnot(!any(is.na(norm_score)))

  # 计算 top 1% 阈值，并选择大于等于该阈值的细胞。
  # 因为分数可能有并列，实际 n_seed 可能和 1% 有极小差异。
  thr <- as.numeric(quantile(norm_score, 1 - SEED_TOP_PCT))
  seed_idx <- which(norm_score >= thr)
  cat("  norm_score top 1% threshold: ", round(thr, 3), "\n")
  cat("  n_seed: ", length(seed_idx), "\n")

  # 保存 0/1 seed 向量：1 表示 seed cell，0 表示非 seed cell。
  # RWR 脚本会读取这个向量，并把它归一化成随机游走的初始分布。
  seed_vec <- numeric(length(common))
  names(seed_vec) <- common
  seed_vec[seed_idx] <- 1
  saveRDS(seed_vec, file.path(OUT_DIR, paste0("seed_MultiOmicsFusion_", tname, ".rds")))
  cat("  saved seed_MultiOmicsFusion_", tname, ".rds\n", sep = "")

  # 记录每个 trait 的 seed 数量和阈值，供后续检查和复现使用。
  seed_summary <- rbind(seed_summary, data.frame(
    trait = tname, n_seed = length(seed_idx), thr_norm_score = round(thr, 3)
  ))
}

# 写出 seed 选择摘要，避免只保存 RDS 后无法快速查看阈值和数量。
write.csv(seed_summary,
          file.path(OUT_DIR, "seed_summary_MultiOmicsFusion.csv"), row.names = FALSE)
cat("\n")
print(seed_summary)
cat("\n04_seed_scdrs.R OK.\n")
