# module1_confidence_sweep.R
# 跨方法公平对比的小工具：对任意一份 atac_label_transfer.csv
#   - 扫描置信度阈值
#   - 画 "保留比例 vs 一致率" 曲线
#   - 在最严阈值下给出逐类型一致率
#
# 为什么需要：不同方法的 score 语义不同
#   - Seurat anchor: prediction.score.max (anchor softmax) → 天生两极
#   - LIGER: top-k 投票比例 → 分布偏高 (中位 0.83)
#   - bindSC: TBD（2026-04-23 从 GLUE 切换为 bindSC）
# 直接比较 "score ≥ 0.5" 的一致率是不公平的。用阈值曲线或 "保留相同比例细胞"
# 才是对齐口径。
#
# 运行方式：
#   source("D:/RStudio/work/module1_confidence_sweep.R")
#   sweep_confidence("method_C_liger")           # 只看 LIGER
#   sweep_confidence("method_A_seurat", "predicted.id", "score")  # 等 A 迁移过来后

library(data.table)
library(dplyr)
library(ggplot2)

data_dir <- "D:/SCAVENGE-reproducibility-main/blood_data"
base_dir <- file.path(data_dir, "multiomics_module1_pairing")

# ---------------------------------------------------------------
# 主函数：对任一 method 目录下的 atac_label_transfer.csv 做扫描
# ---------------------------------------------------------------
sweep_confidence <- function(method_subdir,
                              pred_col   = "predicted.id",
                              score_col  = "score",
                              orig_col   = "original.label",
                              out_prefix = NULL) {
  method_dir <- file.path(base_dir, method_subdir)
  csv_path   <- file.path(method_dir, "atac_label_transfer.csv")
  if (!file.exists(csv_path)) {
    stop("找不到 ", csv_path)
  }

  df <- fread(csv_path)
  need <- c(pred_col, score_col, orig_col)
  miss <- setdiff(need, names(df))
  if (length(miss) > 0) stop("csv 缺列: ", paste(miss, collapse = ", "))

  ok <- !is.na(df[[pred_col]]) & !is.na(df[[orig_col]]) & !is.na(df[[score_col]])
  df <- df[ok, ]
  n_total <- nrow(df)

  # 1) 阈值扫描
  thresholds <- seq(0, 1, by = 0.05)
  sweep_res <- lapply(thresholds, function(th) {
    keep <- df[[score_col]] >= th
    if (sum(keep) == 0) {
      return(data.frame(threshold = th, frac_kept = 0, agreement = NA_real_,
                         n_kept = 0))
    }
    agreement <- mean(df[[pred_col]][keep] == df[[orig_col]][keep])
    data.frame(threshold = th,
               frac_kept = sum(keep) / n_total,
               agreement = agreement,
               n_kept    = sum(keep))
  }) |> rbindlist()

  cat("\n========== 阈值扫描:", method_subdir, "==========\n")
  print(sweep_res, row.names = FALSE)

  # 2) 曲线图：阈值 x 一致率（同时显示保留比例）
  if (is.null(out_prefix)) out_prefix <- gsub("[^A-Za-z0-9]", "_", method_subdir)

  p_sweep <- ggplot(sweep_res, aes(x = frac_kept)) +
    geom_line(aes(y = agreement * 100), color = "#3C5488", size = 1) +
    geom_point(aes(y = agreement * 100), color = "#3C5488", size = 1.5) +
    scale_x_reverse(labels = scales::percent) +
    labs(title = paste0(method_subdir, ": Coverage–Agreement Trade-off"),
         subtitle = "越往左保留越少（更严），越往右保留越多",
         x = "保留细胞比例", y = "一致率 (%)") +
    theme_classic()

  fig_path <- file.path(method_dir, paste0("Fig_", out_prefix, "_confidence_sweep.png"))
  ggsave(fig_path, plot = p_sweep, width = 7, height = 5, dpi = 300)
  cat("已保存:", fig_path, "\n")

  # 3) 最严阈值（全票同类，score == 1.0）下的逐类型一致率
  strict_mask <- df[[score_col]] >= 1.0 - 1e-9
  cat("\n[最严阈值 score >= 1.0] 保留",
      sum(strict_mask), "/", n_total,
      sprintf("(%.1f%%)\n", 100 * sum(strict_mask) / n_total))

  if (sum(strict_mask) > 0) {
    df_s <- df[strict_mask, ]
    conf <- table(Original    = df_s[[orig_col]],
                  Transferred = df_s[[pred_col]])
    orig_types  <- rownames(conf)
    trans_types <- colnames(conf)
    agreement_per_type <- sapply(orig_types, function(ct) {
      if (ct %in% trans_types) conf[ct, ct] / sum(conf[ct, ]) else 0
    })
    strict_df <- data.frame(
      celltype        = orig_types,
      n_cells_strict  = as.integer(rowSums(conf)),
      agreement_strict = round(agreement_per_type * 100, 1)
    ) |> arrange(desc(agreement_strict))
    print(strict_df, row.names = FALSE)

    out_csv <- file.path(method_dir, "label_transfer_agreement_strict.csv")
    fwrite(strict_df, out_csv)
    cat("已保存:", out_csv, "\n")
  }

  invisible(sweep_res)
}

# ---------------------------------------------------------------
# 便捷入口：默认跑 LIGER
# ---------------------------------------------------------------
if (interactive() || !exists("SKIP_AUTO_RUN")) {
  sweep_confidence("method_C_liger")
}
