# module1_method_B_finalize_outputs.R
# 模块 1 · 方法 B：把 Python scSHEFT CSV 结果整理为统一接口契约
#
# 输入：
#   method_B_scSHEFT/shared_embedding.csv
#   method_B_scSHEFT/atac_topk_rna_neighbors.csv
#   method_B_scSHEFT/atac_label_transfer.csv
#
# 输出：
#   method_B_scSHEFT/shared_embedding.rds
#   method_B_scSHEFT/atac_topk_rna_neighbors.rds
#   method_B_scSHEFT/label_transfer_agreement.csv
#   method_B_scSHEFT/label_transfer_agreement_highconf.csv
#   method_B_scSHEFT/label_transfer_agreement_strict.csv

required_pkgs <- c("data.table", "jsonlite")
missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("缺少包: ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

out_dir <- "D:/SCAVENGE-reproducibility-main/blood_data/multiomics_module1_pairing/method_B_scSHEFT"

cat("========== 方法 B 输出整理 ==========\n")
cat("  目录:", out_dir, "\n")

label_path <- file.path(out_dir, "atac_label_transfer.csv")
embed_path <- file.path(out_dir, "shared_embedding.csv")
topk_path  <- file.path(out_dir, "atac_topk_rna_neighbors.csv")

stopifnot(file.exists(label_path), file.exists(embed_path), file.exists(topk_path))

cat("\n========== shared_embedding.rds ==========\n")
shared_embed <- fread(embed_path)
saveRDS(as.data.frame(shared_embed),
        file.path(out_dir, "shared_embedding.rds"),
        compress = FALSE)
cat("  已保存: shared_embedding.rds (", nrow(shared_embed), "x", ncol(shared_embed), ")\n", sep = "")

cat("\n========== atac_topk_rna_neighbors.rds ==========\n")
topk <- fread(topk_path)
atac_topk_rna <- split(
  topk[, .(rna_cell, distance)],
  topk$atac_cell
)
saveRDS(atac_topk_rna,
        file.path(out_dir, "atac_topk_rna_neighbors.rds"),
        compress = FALSE)
cat("  已保存: atac_topk_rna_neighbors.rds (", length(atac_topk_rna), " ATAC cells)\n", sep = "")

cat("\n========== agreement ==========\n")
label_df <- fread(label_path)
label_df_ok <- label_df[
  !is.na(original.label) & !is.na(predicted.id)
]

make_agreement <- function(dt, suffix = NULL) {
  if (nrow(dt) == 0) return(NA_real_)

  conf <- table(
    Original = dt$original.label,
    Transferred = dt$predicted.id
  )
  orig_types  <- rownames(conf)
  trans_types <- colnames(conf)
  shared_lab  <- intersect(orig_types, trans_types)
  agreement_per_type <- sapply(orig_types, function(ct) {
    if (ct %in% trans_types) conf[ct, ct] / sum(conf[ct, ]) else 0
  })
  overall <- sum(sapply(shared_lab, function(ct) conf[ct, ct])) / sum(conf)

  agreement_df <- data.frame(
    celltype = orig_types,
    n_cells = as.integer(rowSums(conf)),
    agreement = round(agreement_per_type * 100, 1),
    stringsAsFactors = FALSE
  )
  agreement_df <- agreement_df[order(-agreement_df$agreement), ]

  fname <- if (is.null(suffix)) {
    "label_transfer_agreement.csv"
  } else {
    paste0("label_transfer_agreement_", suffix, ".csv")
  }
  fwrite(agreement_df, file.path(out_dir, fname))
  cat("  已保存:", fname, "\n")
  overall
}

overall <- make_agreement(label_df_ok)
cat("  总体一致率:", round(overall * 100, 2), "%\n")

highconf_mask <- !is.na(label_df_ok$score) & label_df_ok$score >= 0.5
cat("  highconf score >= 0.5:", sum(highconf_mask), "/", nrow(label_df_ok), "\n")
hc_overall <- make_agreement(label_df_ok[highconf_mask], "highconf")
if (!is.na(hc_overall)) cat("  highconf 总体一致率:", round(hc_overall * 100, 2), "%\n")

strict_mask <- !is.na(label_df_ok$score) & label_df_ok$score >= 0.95
cat("  strict score >= 0.95:", sum(strict_mask), "/", nrow(label_df_ok), "\n")
strict_overall <- make_agreement(label_df_ok[strict_mask], "strict")
if (!is.na(strict_overall)) cat("  strict 总体一致率:", round(strict_overall * 100, 2), "%\n")

summary <- list(
  method = "scSHEFT official CSUBioGroup implementation",
  n_cells = nrow(label_df),
  overall_agreement = overall,
  highconf_n = sum(highconf_mask),
  highconf_agreement = hc_overall,
  strict_n = sum(strict_mask),
  strict_agreement = strict_overall,
  note = "anchor_weight=0 run; official full anchor disabled because it requires a dense RNA x ATAC distance matrix."
)
jsonlite::write_json(summary, file.path(out_dir, "finalize_summary.json"),
                     pretty = TRUE, auto_unbox = TRUE)

cat("\n======================================================\n")
cat("  方法 B 输出整理完成\n")
cat("======================================================\n")
cat("  shared_embedding.rds:", file.exists(file.path(out_dir, "shared_embedding.rds")), "\n")
cat("  atac_topk_rna_neighbors.rds:", file.exists(file.path(out_dir, "atac_topk_rna_neighbors.rds")), "\n")
cat("  overall agreement:", round(overall * 100, 2), "%\n")
cat("======================================================\n")
