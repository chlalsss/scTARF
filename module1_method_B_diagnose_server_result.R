# module1_method_B_diagnose_server_result.R
# 模块 1 · 方法 B：诊断服务器版全量官方 scSHEFT 结果
#
# 输出：
#   - score_summary.txt
#   - predicted_class_distribution.csv
#   - original_class_distribution.csv
#   - confusion_matrix_top.csv
#   - class_pair_mismatch_top.csv
#   - diagnosis_summary.json

required_pkgs <- c("data.table", "jsonlite")
missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("缺少包: ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

base_dir <- "D:/SCAVENGE-reproducibility-main/blood_data/multiomics_module1_pairing/method_B_scSHEFT"
candidate_dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
server_dirs <- candidate_dirs[grepl("server_full_official", basename(candidate_dirs), fixed = TRUE)]

if (length(server_dirs) == 0) {
  stop("未找到服务器结果目录。")
}
if (length(server_dirs) > 1) {
  stop("找到多个服务器结果目录，请保留一个或手动指定。")
}

out_dir <- server_dirs[1]
label_path <- file.path(out_dir, "atac_label_transfer.csv")
stopifnot(file.exists(label_path))

cat("========== 方法 B 服务器结果诊断 ==========\n")
cat("  目录:", out_dir, "\n")

dt <- fread(label_path)
dt_ok <- dt[!is.na(original.label) & !is.na(predicted.id)]
if (nrow(dt_ok) == 0) stop("无有效 label transfer 记录")

cat("\n========== score 分布 ==========\n")
score_summary <- summary(dt_ok$score)
score_quant <- quantile(dt_ok$score, probs = c(0, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1),
                        na.rm = TRUE)
score_lines <- c(
  "summary(score):",
  capture.output(print(score_summary)),
  "",
  "quantile(score):",
  capture.output(print(score_quant)),
  "",
  paste0("n(score>=0.5): ", sum(dt_ok$score >= 0.5, na.rm = TRUE)),
  paste0("n(score>=0.95): ", sum(dt_ok$score >= 0.95, na.rm = TRUE))
)
writeLines(score_lines, file.path(out_dir, "score_summary.txt"))
cat(paste(score_lines, collapse = "\n"), "\n")

cat("\n========== 类别分布 ==========\n")
pred_dist <- dt_ok[, .N, by = predicted.id][order(-N)]
orig_dist <- dt_ok[, .N, by = original.label][order(-N)]
setnames(pred_dist, c("predicted.id", "N"), c("class", "n_cells"))
setnames(orig_dist, c("original.label", "N"), c("class", "n_cells"))
pred_dist[, proportion := n_cells / sum(n_cells)]
orig_dist[, proportion := n_cells / sum(n_cells)]
fwrite(pred_dist, file.path(out_dir, "predicted_class_distribution.csv"))
fwrite(orig_dist, file.path(out_dir, "original_class_distribution.csv"))
cat("  已保存: predicted_class_distribution.csv\n")
cat("  已保存: original_class_distribution.csv\n")
cat("  预测前10类:\n")
print(pred_dist[1:min(10, .N)])

cat("\n========== 混淆与偏移 ==========\n")
conf_dt <- dt_ok[, .N, by = .(original.label, predicted.id)][order(-N)]
fwrite(conf_dt, file.path(out_dir, "confusion_matrix_top.csv"))
cat("  已保存: confusion_matrix_top.csv\n")

wrong_dt <- dt_ok[original.label != predicted.id, .N, by = .(original.label, predicted.id)][order(-N)]
fwrite(wrong_dt, file.path(out_dir, "class_pair_mismatch_top.csv"))
cat("  已保存: class_pair_mismatch_top.csv\n")
cat("  误配前15对:\n")
print(wrong_dt[1:min(15, .N)])

overall <- mean(dt_ok$original.label == dt_ok$predicted.id)

summary_json <- list(
  out_dir = out_dir,
  n_cells = nrow(dt_ok),
  overall_agreement = overall,
  score_summary = as.list(score_summary),
  score_quantile = as.list(score_quant),
  n_score_ge_0_5 = sum(dt_ok$score >= 0.5, na.rm = TRUE),
  n_score_ge_0_95 = sum(dt_ok$score >= 0.95, na.rm = TRUE),
  top_predicted_classes = pred_dist[1:min(10, .N)],
  top_original_classes = orig_dist[1:min(10, .N)],
  top_mismatches = wrong_dt[1:min(15, .N)]
)
jsonlite::write_json(summary_json, file.path(out_dir, "diagnosis_summary.json"),
                     pretty = TRUE, auto_unbox = TRUE)

cat("\n======================================================\n")
cat("  方法 B 服务器结果诊断完成\n")
cat("======================================================\n")
cat("  overall agreement:", round(overall * 100, 2), "%\n")
cat("  输出目录:", out_dir, "\n")
cat("======================================================\n")
