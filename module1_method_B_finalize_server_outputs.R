# module1_method_B_finalize_server_outputs.R
# 模块 1 · 方法 B：整理服务器版全量官方 scSHEFT 结果
#
# 输入目录：
#   method_B_scSHEFT/server_full_official_anchor_1000steps_20260424/
#     - shared_embedding.csv
#     - atac_label_transfer.csv
#     - timing.json
#
# 输出目录（同输入目录）：
#   - atac_topk_rna_neighbors.csv
#   - atac_topk_rna_neighbors.rds
#   - shared_embedding.rds
#   - label_transfer_agreement.csv
#   - label_transfer_agreement_highconf.csv
#   - label_transfer_agreement_strict.csv
#   - finalize_summary.json

required_pkgs <- c("data.table", "RANN", "jsonlite")
missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("缺少包: ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(data.table)
  library(RANN)
  library(jsonlite)
})

base_dir <- "D:/SCAVENGE-reproducibility-main/blood_data/multiomics_module1_pairing/method_B_scSHEFT"
candidate_dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
server_dirs <- candidate_dirs[grepl("server_full_official", basename(candidate_dirs), fixed = TRUE)]

if (length(server_dirs) == 0) {
  stop("未找到服务器结果目录。请确认 method_B_scSHEFT 下存在 server_full_official* 目录。")
}
if (length(server_dirs) > 1) {
  stop("找到多个服务器结果目录，请保留一个或手动在脚本中指定：\n",
       paste(server_dirs, collapse = "\n"))
}

out_dir <- server_dirs[1]
label_path <- file.path(out_dir, "atac_label_transfer.csv")
embed_path <- file.path(out_dir, "shared_embedding.csv")
timing_path <- file.path(out_dir, "timing.json")

stopifnot(file.exists(label_path), file.exists(embed_path), file.exists(timing_path))

cat("========== 方法 B 服务器结果整理 ==========\n")
cat("  目录:", out_dir, "\n")

cat("\n========== 读取 shared embedding ==========\n")
shared_embed <- fread(embed_path)
if (!all(c("cell", "modality") %in% colnames(shared_embed))) {
  stop("shared_embedding.csv 缺少 cell/modality 列")
}

latent_cols <- setdiff(colnames(shared_embed), c("cell", "modality"))
if (length(latent_cols) == 0) stop("未找到 latent 列")

rna_dt  <- shared_embed[modality == "rna"]
atac_dt <- shared_embed[modality == "atac"]
if (nrow(rna_dt) == 0 || nrow(atac_dt) == 0) stop("RNA/ATAC 行为空")

rna_mat  <- as.matrix(rna_dt[, ..latent_cols])
atac_mat <- as.matrix(atac_dt[, ..latent_cols])
storage.mode(rna_mat)  <- "double"
storage.mode(atac_mat) <- "double"

cat("  RNA :", nrow(rna_mat),  "x", ncol(rna_mat),  "\n")
cat("  ATAC:", nrow(atac_mat), "x", ncol(atac_mat), "\n")

cat("\n========== 保存 shared_embedding.rds ==========\n")
saveRDS(as.data.frame(shared_embed),
        file.path(out_dir, "shared_embedding.rds"),
        compress = FALSE)
cat("  已保存: shared_embedding.rds\n")

cat("\n========== 计算 top-k RNA 邻居 ==========\n")
K_TOPK <- 30
nn <- RANN::nn2(data = rna_mat, query = atac_mat, k = K_TOPK)

topk_rows <- vector("list", nrow(atac_mat))
atac_topk_rna <- vector("list", nrow(atac_mat))
rna_cells <- rna_dt$cell
atac_cells <- atac_dt$cell

for (i in seq_len(nrow(atac_mat))) {
  df_i <- data.frame(
    rna_cell = rna_cells[nn$nn.idx[i, ]],
    distance = nn$nn.dists[i, ],
    stringsAsFactors = FALSE
  )
  atac_topk_rna[[i]] <- df_i
  topk_rows[[i]] <- data.frame(
    atac_cell = atac_cells[i],
    rank = seq_len(K_TOPK),
    rna_cell = df_i$rna_cell,
    distance = df_i$distance,
    stringsAsFactors = FALSE
  )
}
names(atac_topk_rna) <- atac_cells

topk_df <- rbindlist(topk_rows)
fwrite(topk_df, file.path(out_dir, "atac_topk_rna_neighbors.csv"))
saveRDS(atac_topk_rna,
        file.path(out_dir, "atac_topk_rna_neighbors.rds"),
        compress = FALSE)
cat("  已保存: atac_topk_rna_neighbors.csv\n")
cat("  已保存: atac_topk_rna_neighbors.rds\n")

cat("\n========== 生成 agreement 文件 ==========\n")
label_df <- fread(label_path)
label_df_ok <- label_df[!is.na(original.label) & !is.na(predicted.id)]

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

  out_df <- data.frame(
    celltype = orig_types,
    n_cells = as.integer(rowSums(conf)),
    agreement = round(agreement_per_type * 100, 1),
    stringsAsFactors = FALSE
  )
  out_df <- out_df[order(-out_df$agreement), ]

  fname <- if (is.null(suffix)) "label_transfer_agreement.csv"
  else paste0("label_transfer_agreement_", suffix, ".csv")
  fwrite(out_df, file.path(out_dir, fname))
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

timing <- jsonlite::read_json(timing_path, simplifyVector = TRUE)
summary <- list(
  source_timing = timing,
  overall_agreement = overall,
  highconf_n = sum(highconf_mask),
  highconf_agreement = hc_overall,
  strict_n = sum(strict_mask),
  strict_agreement = strict_overall,
  note = "Server full official anchor run; top-k neighbors computed post hoc from shared_embedding.csv"
)
jsonlite::write_json(summary, file.path(out_dir, "finalize_summary.json"),
                     pretty = TRUE, auto_unbox = TRUE)

cat("\n======================================================\n")
cat("  方法 B 服务器结果整理完成\n")
cat("======================================================\n")
cat("  overall agreement:", round(overall * 100, 2), "%\n")
cat("  输出目录:", out_dir, "\n")
cat("======================================================\n")
