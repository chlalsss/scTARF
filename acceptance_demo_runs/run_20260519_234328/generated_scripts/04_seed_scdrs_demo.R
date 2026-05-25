# 04_seed_scdrs.R
# §2.D: SCDRS norm_score top 1% as seed cells (per trait, atac_w-INDEPENDENT).
# 直接 source(), 跑一次即可 (atac_w sweep 也用同一份 seed).

suppressPackageStartupMessages({
  library(data.table)
})

OUT_DIR   <- "D:/RStudio/work/acceptance_demo_runs/run_20260519_234328/outputs/algorithm_core_mono"
SCORE_DIR <- "D:/RStudio/work/module3_multiome/04_dual_modality/01_scdrs_python/output"

TRAITS <- c(
  MONO = "UKB_460K.blood_MONOCYTE_COUNT",
  PLT  = "UKB_460K.blood_PLATELET_COUNT",
  LYM  = "UKB_460K.blood_LYMPHOCYTE_COUNT"
)
SEED_TOP_PCT <- 0.01

# common cells = rna_z_block 的 cells (already aligned to SCDRS-filtered set, 68303)
# atac_z_block_v2.rds 是 align 前的 69249, 不能直接当 common
rna_z_ref <- readRDS(file.path(OUT_DIR, "rna_z_block_v2_MONO.rds"))
common <- rownames(rna_z_ref)
cat("common cells:", length(common), "\n\n")

seed_summary <- data.frame()
for (tname in names(TRAITS)) {
  trait <- TRAITS[[tname]]
  cat("=== ", tname, " ===\n", sep = "")
  dt <- fread(file.path(SCORE_DIR, paste0(trait, ".score.gz")))
  setnames(dt, 1, "barcode")
  setkey(dt, barcode)
  norm_score <- dt[common, norm_score]
  names(norm_score) <- common
  stopifnot(!any(is.na(norm_score)))

  thr <- as.numeric(quantile(norm_score, 1 - SEED_TOP_PCT))
  seed_idx <- which(norm_score >= thr)
  cat("  norm_score top 1% threshold: ", round(thr, 3), "\n")
  cat("  n_seed: ", length(seed_idx), "\n")

  seed_vec <- numeric(length(common))
  names(seed_vec) <- common
  seed_vec[seed_idx] <- 1
  saveRDS(seed_vec, file.path(OUT_DIR, paste0("seed_v2_", tname, ".rds")))
  cat("  saved seed_v2_", tname, ".rds\n", sep = "")

  seed_summary <- rbind(seed_summary, data.frame(
    trait = tname, n_seed = length(seed_idx), thr_norm_score = round(thr, 3)
  ))
}

write.csv(seed_summary,
          file.path(OUT_DIR, "seed_summary_v2.csv"), row.names = FALSE)
cat("\n")
print(seed_summary)
cat("\n04_seed_scdrs.R OK.\n")
