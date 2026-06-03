# Build a fully refreshed Chinese final-delivery package under final_outputs/最终输出2.
# This script scans D:/RStudio/work, collects process tables, supplements the
# formal AUROC evidence table, and redraws all report figures with Chinese text.

source("D:/RStudio/work/multiomics_fusion_final_delivery_plan_20260528/scripts/_final_delivery_common.R")
need_packages(c("data.table", "ggplot2"))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

OUT_ROOT <- file.path(FINAL_ROOT, "最终输出2")
out_root_norm <- normalizePath(OUT_ROOT, winslash = "/", mustWork = FALSE)
final_root_norm <- normalizePath(FINAL_ROOT, winslash = "/", mustWork = TRUE)
if (dir.exists(OUT_ROOT)) {
  if (!identical(dirname(out_root_norm), final_root_norm) || !identical(basename(out_root_norm), "最终输出2")) {
    stop("Refusing to clean unexpected output path: ", OUT_ROOT, call. = FALSE)
  }
  unlink(OUT_ROOT, recursive = TRUE, force = TRUE)
}
OUT_DIRS <- c(
  "00_全盘扫描",
  "01_主结果表",
  "02_主结果图",
  "03_消融与参数",
  "04_效率",
  "05_RNA参考依赖",
  "06_案例研究/MONO",
  "06_案例研究/MONO/图片",
  "07_过程数据",
  "07_过程数据/候选文件镜像",
  "99_来源记录"
)
invisible(lapply(file.path(OUT_ROOT, OUT_DIRS), dir_ensure))

write_csv2 <- function(x, path) {
  dir_ensure(dirname(path))
  data.table::fwrite(x, path, bom = TRUE)
  invisible(path)
}

write_text <- function(x, path) {
  dir_ensure(dirname(path))
  writeLines(enc2utf8(x), path, useBytes = TRUE)
  invisible(path)
}

source_exists <- function(path) {
  parts <- strsplit(as.character(path), ";", fixed = TRUE)
  vapply(parts, function(x) {
    x <- trimws(x)
    length(x) > 0L && all(file.exists(x))
  }, logical(1))
}

cn_trait <- function(x) {
  map <- c(
    "MONO" = "单核细胞",
    "PLT" = "血小板",
    "LYM" = "淋巴细胞",
    "LYM_wide" = "淋巴细胞宽定义",
    "EOS" = "嗜酸性粒细胞",
    "MCH" = "平均红细胞血红蛋白量",
    "RBC" = "红细胞",
    "RDW" = "红细胞分布宽度",
    "WBC" = "白细胞"
  )
  out <- unname(map[as.character(x)])
  out[is.na(out)] <- as.character(x)[is.na(out)]
  out
}

cn_pairing <- function(x) {
  map <- c(
    "paired" = "配对",
    "non_paired" = "非配对",
    "rna_reference_only" = "仅核糖核酸参考"
  )
  out <- unname(map[as.character(x)])
  out[is.na(out)] <- as.character(x)[is.na(out)]
  out
}

cn_dataset <- function(x) {
  map <- c(
    "BMMC paired multiome" = "BMMC配对多组学",
    "Granja non-paired ATAC + RNA reference" = "Granja非配对染色质+核糖核酸参考",
    "RNA reference scDRS" = "核糖核酸参考scDRS"
  )
  out <- unname(map[as.character(x)])
  out[is.na(out)] <- as.character(x)[is.na(out)]
  out
}

cn_version <- function(x) {
  map <- c(
    "Baseline" = "原始SCAVENGE",
    "scDRS-only" = "纯scDRS核糖核酸基线",
    "Paired-Fusion" = "配对核糖核酸-染色质融合",
    "v2.1" = "v2.1 基础伪配对融合",
    "v2.2" = "v2.2 性状感知FDR融合",
    "v3.1" = "v3.1 单参考细胞消融",
    "v3.2" = "v3.2 前5显著核糖核酸补全",
    "v3.3" = "v3.3 稀疏映射前5融合",
    "Final-best" = "最佳融合汇总"
  )
  out <- unname(map[as.character(x)])
  out[is.na(out)] <- as.character(x)[is.na(out)]
  out
}

cn_version_label <- function(x) cn_version(x)

cn_metric <- function(x) {
  map <- c(
    "Runtime (s)" = "运行时间（秒）",
    "Memory (GB)" = "峰值内存（GB）",
    "AUROC" = "曲线下面积",
    "ACC" = "前列命中率"
  )
  out <- unname(map[as.character(x)])
  out[is.na(out)] <- as.character(x)[is.na(out)]
  out
}

cn_method_label <- function(x) {
  y <- as.character(x)
  y <- gsub("Paired-Fusion: paired RNA-ATAC fusion", "配对核糖核酸-染色质融合", y, fixed = TRUE)
  y <- gsub("v3.3: sparse mapping + top5 significant fusion", "v3.3 稀疏映射前5融合", y, fixed = TRUE)
  y <- gsub("Non-paired v3.3", "非配对v3.3", y, fixed = TRUE)
  y <- gsub("Paired-Fusion", "配对融合", y, fixed = TRUE)
  y
}

plot_theme <- function(base_size = 12) {
  theme_bw(base_size = base_size, base_family = "sans") +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      legend.position = "bottom"
    )
}

save_plot <- function(plot, path, width, height, dpi = 300) {
  dir_ensure(dirname(path))
  if (tolower(tools::file_ext(path)) == "pdf" && isTRUE(capabilities("cairo"))) {
    ggsave(path, plot, width = width, height = height, dpi = dpi, device = grDevices::cairo_pdf)
  } else {
    ggsave(path, plot, width = width, height = height, dpi = dpi)
  }
  invisible(path)
}

safe_copy_name <- function(path, index) {
  name <- gsub("[:*?\"<>|/\\\\]", "_", basename(path))
  sprintf("%04d_%s", index, name)
}

cat(">>> 1. 全盘扫描候选过程文件\n")
all_files <- list.files(WORK_ROOT, recursive = TRUE, full.names = TRUE, all.files = FALSE, no.. = TRUE)
all_files <- all_files[!grepl("/final_outputs/最终输出2/", normalizePath(all_files, winslash = "/", mustWork = FALSE), fixed = TRUE)]
candidate_re <- "(auroc|auc|roc|scdrs|target_enrichment|result|summary|sweep|ablation|runtime|memory|efficiency|manifest)\\.(csv|tsv|txt|rds|xlsx)$"
candidates <- all_files[grepl(candidate_re, basename(all_files), ignore.case = TRUE)]

file_info <- file.info(candidates)
scan <- data.frame(
  源文件 = normalizePath(candidates, winslash = "/", mustWork = FALSE),
  文件名 = basename(candidates),
  扩展名 = tolower(tools::file_ext(candidates)),
  文件大小字节 = as.numeric(file_info$size),
  修改时间 = as.character(file_info$mtime),
  stringsAsFactors = FALSE
)
scan$类别 <- "其他候选"
scan$类别[grepl("scdrs|target_enrichment", scan$源文件, ignore.case = TRUE)] <- "scDRS或核糖核酸参考"
scan$类别[grepl("auroc|auc|roc", scan$源文件, ignore.case = TRUE)] <- "曲线下面积结果"
scan$类别[grepl("sweep|ablation", scan$源文件, ignore.case = TRUE)] <- "消融或参数扫描"
scan$类别[grepl("runtime|memory|efficiency|elapsed", scan$源文件, ignore.case = TRUE)] <- "效率或资源"
scan$类别[grepl("manifest", scan$源文件, ignore.case = TRUE)] <- "来源记录"

scan$行数 <- NA_integer_
scan$列数 <- NA_integer_
scan$列名 <- NA_character_
is_text_table <- scan$扩展名 %in% c("csv", "tsv", "txt")
for (i in which(is_text_table)) {
  path <- scan$源文件[i]
  sep <- if (tolower(tools::file_ext(path)) == "tsv") "\t" else "auto"
  x <- tryCatch(data.table::fread(path, nrows = 1000000, sep = sep, data.table = FALSE, showProgress = FALSE),
                error = function(e) NULL)
  if (!is.null(x)) {
    scan$行数[i] <- nrow(x)
    scan$列数[i] <- ncol(x)
    scan$列名[i] <- paste(names(x), collapse = ";")
  }
}
write_csv2(scan, file.path(OUT_ROOT, "00_全盘扫描", "全盘候选过程文件索引.csv"))

mirror_map <- scan[scan$扩展名 %in% c("csv", "tsv", "txt"), c("源文件", "文件名", "类别"), drop = FALSE]
mirror_map$镜像文件 <- NA_character_
for (i in seq_len(nrow(mirror_map))) {
  to <- file.path(OUT_ROOT, "07_过程数据", "候选文件镜像", safe_copy_name(mirror_map$源文件[i], i))
  dir_ensure(dirname(to))
  ok <- file.copy(mirror_map$源文件[i], to, overwrite = TRUE)
  if (isTRUE(ok)) mirror_map$镜像文件[i] <- normalizePath(to, winslash = "/", mustWork = FALSE)
}
write_csv2(mirror_map, file.path(OUT_ROOT, "07_过程数据", "候选文件镜像清单.csv"))

cat("候选文件数: ", nrow(scan), "；已镜像文本表: ", sum(!is.na(mirror_map$镜像文件)), "\n", sep = "")

cat(">>> 2. 抽取全盘曲线下面积候选证据\n")
extract_auc_rows <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (!ext %in% c("csv", "tsv", "txt")) return(NULL)
  sep <- if (ext == "tsv") "\t" else "auto"
  x <- tryCatch(data.table::fread(path, sep = sep, data.table = FALSE, showProgress = FALSE),
                error = function(e) NULL)
  if (is.null(x) || nrow(x) == 0L) return(NULL)
  value_cols <- intersect(
    c("auroc", "AUROC", "auc", "AUC", "target_auroc", "cite_score_target_auroc",
      "pseudo_z_target_auroc", "source_target_weight_auroc", "value"),
    names(x)
  )
  value_cols <- value_cols[grepl("auc|auroc|value", value_cols, ignore.case = TRUE)]
  if (length(value_cols) == 0L) return(NULL)
  trait_col <- intersect(c("trait", "target", "target_definition"), names(x))[1]
  method_col <- intersect(c("method", "method_version", "experiment_id", "method_label"), names(x))[1]
  feature_col <- intersect(c("feature_mode", "experiment_version", "table_version"), names(x))[1]
  atac_col <- intersect(c("atac_w", "atac_weight", "final_atac_w"), names(x))[1]
  npos_col <- intersect(c("n_positive", "n_positive_target", "atac_n_target"), names(x))[1]
  ntot_col <- intersect(c("n_total", "n_cells", "cite_n_cells", "atac_n_cells"), names(x))[1]

  rows <- lapply(value_cols, function(vc) {
    data.frame(
      源文件 = normalizePath(path, winslash = "/", mustWork = FALSE),
      指标列 = vc,
      性状 = if (!is.na(trait_col)) as.character(x[[trait_col]]) else NA_character_,
      方法 = if (!is.na(method_col)) as.character(x[[method_col]]) else NA_character_,
      特征或实验 = if (!is.na(feature_col)) as.character(x[[feature_col]]) else NA_character_,
      染色质权重 = if (!is.na(atac_col)) suppressWarnings(as.numeric(x[[atac_col]])) else NA_real_,
      数值 = suppressWarnings(as.numeric(x[[vc]])),
      阳性数 = if (!is.na(npos_col)) suppressWarnings(as.numeric(x[[npos_col]])) else NA_real_,
      总数 = if (!is.na(ntot_col)) suppressWarnings(as.numeric(x[[ntot_col]])) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  out <- bind_rows_fill(rows)
  out <- out[is.finite(out$数值), , drop = FALSE]
  if (nrow(out) == 0L) return(NULL)
  out
}

auc_evidence <- bind_rows_fill(lapply(scan$源文件, extract_auc_rows))
if (nrow(auc_evidence) > 0L) {
  write_csv2(auc_evidence, file.path(OUT_ROOT, "07_过程数据", "全盘曲线下面积候选证据.csv"))
}

cat("抽取到候选证据行: ", nrow(auc_evidence), "\n", sep = "")

cat(">>> 3. 补全正式主结果长表\n")
main_file <- file.path(FINAL_ROOT, "01_master_tables/main_result_comparison_table.csv")
if (!file.exists(main_file)) stop("缺少主结果长表: ", main_file, call. = FALSE)
main <- read_csv_if_exists(main_file)

extra_specs <- data.frame(
  source_file = c(
    file.path(MODULE4_ROOT, "04_evaluation/v3_dense_top1_missing_traits/neighborhood_fdr_ge0p80/v3_dense_top1_missing_traits_neighborhood_fdr_ge0p80_auroc.csv"),
    file.path(MODULE4_ROOT, "04_evaluation/v3_top5_significant/trait_screening_ge0p80/v3_top5_significant_trait_screening_ge0p80_auroc.csv"),
    file.path(MODULE4_ROOT, "04_evaluation/v3_top5_significant_p005/trait_screening_ge0p80/v3_top5_significant_p005_trait_screening_ge0p80_auroc.csv")
  ),
  dataset = "Granja non-paired ATAC + RNA reference",
  pairing_setting = "non_paired",
  method_version = c(
    "v3_dense_top1_missing_traits",
    "v3_top5_significant_trait_screening",
    "v3_top5_significant_p005_trait_screening"
  ),
  rna_reference = "CITE / Granja RNA reference",
  top_k = c(1, 5, 5),
  feature_mode = c(
    "top1_dense_missing_traits",
    "top5_significant_trait_screening",
    "top5_significant_p005_trait_screening"
  ),
  note = c(
    "全盘扫描补入：v3.1 top1 missing-traits 正式运行",
    "全盘扫描补入：v3.2 top5 significant trait-screening 正式运行",
    "全盘扫描补入：v3.2 top5 significant p<=0.05 trait-screening 正式运行"
  ),
  stringsAsFactors = FALSE
)

extra_rows <- list()
for (i in seq_len(nrow(extra_specs))) {
  spec <- extra_specs[i, ]
  x <- read_csv_if_exists(spec$source_file)
  if (is.null(x)) next
  extra_rows[[length(extra_rows) + 1L]] <- normalize_auroc_table(
    x = x,
    source_file = spec$source_file,
    dataset = spec$dataset,
    pairing_setting = spec$pairing_setting,
    method_version = spec$method_version,
    rna_reference = spec$rna_reference,
    top_k = spec$top_k,
    feature_mode = spec$feature_mode,
    note = spec$note
  )
}

main2 <- bind_rows_fill(list(main, bind_rows_fill(extra_rows)))
main2 <- main2[main2$metric == "AUROC" & is.finite(main2$value), , drop = FALSE]
main2$final_version <- classify_final_version(
  method_version = main2$method_version,
  method = main2$method,
  feature_mode = main2$feature_mode,
  method_group = main2$method_group,
  pairing_setting = main2$pairing_setting
)
main2$final_version_label <- final_version_label(main2$final_version)
main2 <- unique(main2)

best_method <- best_by_key(
  main2,
  keys = c("dataset", "pairing_setting", "trait", "method_group", "method_version", "feature_mode")
)
best_final <- best_by_key(
  main2,
  keys = c("dataset", "pairing_setting", "trait", "final_version")
)

main2_cn <- main2
name_map <- c(
  dataset = "数据集", pairing_setting = "数据类型", trait = "性状",
  target_definition = "目标定义", method_group = "方法组",
  method_version = "方法版本", method = "原始方法名",
  rna_reference = "核糖核酸参考", top_k = "参考细胞数",
  feature_mode = "特征模式", atac_w = "染色质权重",
  metric = "指标", value = "数值", ci_low = "置信区间下限",
  ci_high = "置信区间上限", n_positive = "阳性数",
  n_total = "总数", source_file = "来源文件",
  note = "备注", final_version = "最终版本",
  final_version_label = "最终版本标签"
)
mapped_names <- unname(name_map[names(main2_cn)])
mapped_names[is.na(mapped_names)] <- names(main2_cn)[is.na(mapped_names)]
names(main2_cn) <- mapped_names
write_csv2(main2_cn, file.path(OUT_ROOT, "01_主结果表", "正式主结果补全长表.csv"))

best_final_cn <- best_final
best_final_cn$中文数据集 <- cn_dataset(best_final_cn$dataset)
best_final_cn$中文数据类型 <- cn_pairing(best_final_cn$pairing_setting)
best_final_cn$中文性状 <- cn_trait(best_final_cn$trait)
best_final_cn$中文最终版本 <- cn_version(best_final_cn$final_version)
write_csv2(best_final_cn, file.path(OUT_ROOT, "01_主结果表", "正式主结果按最终版本最佳表.csv"))

write_csv2(best_method, file.path(OUT_ROOT, "01_主结果表", "正式主结果按方法版本最佳表.csv"))

make_wide <- function(x) {
  scdrs_lookup <- x[x$final_version == "scDRS-only", , drop = FALSE]
  scdrs_by_trait <- lapply(split(seq_len(nrow(scdrs_lookup)), scdrs_lookup$trait), function(ii) {
    best_by_key(scdrs_lookup[ii, , drop = FALSE], keys = c("trait"))
  })
  main_x <- x[x$pairing_setting != "rna_reference_only", , drop = FALSE]
  group_key <- paste(main_x$dataset, main_x$pairing_setting, main_x$trait, sep = "\r")
  rows <- lapply(split(seq_len(nrow(main_x)), group_key), function(ii) {
    y <- main_x[ii, , drop = FALSE]
    pick <- function(ver) {
      z <- y[y$final_version %in% ver, , drop = FALSE]
      if (nrow(z) == 0L) return(NULL)
      z[which.max(z$value), , drop = FALSE]
    }
    original <- pick("Baseline")
    scdrs <- pick("scDRS-only")
    if (is.null(scdrs) && identical(y$pairing_setting[1], "non_paired")) scdrs <- scdrs_by_trait[[y$trait[1]]]
    fusion <- pick(c("Paired-Fusion", "v2.1", "v2.2", "v3.1", "v3.2", "v3.3"))
    get_val <- function(z) if (is.null(z)) NA_real_ else z$value[1]
    get_src <- function(z) if (is.null(z)) NA_character_ else z$source_file[1]
    get_ver <- function(z) if (is.null(z)) NA_character_ else cn_version(z$final_version[1])
    data.frame(
      数据集 = cn_dataset(y$dataset[1]),
      数据类型 = cn_pairing(y$pairing_setting[1]),
      性状 = cn_trait(y$trait[1]),
      原始SCAVENGE曲线下面积 = get_val(original),
      纯scDRS曲线下面积 = get_val(scdrs),
      配对融合曲线下面积 = get_val(pick("Paired-Fusion")),
      v2.1基础伪配对曲线下面积 = get_val(pick("v2.1")),
      v2.2性状感知FDR曲线下面积 = get_val(pick("v2.2")),
      v3.1单参考细胞曲线下面积 = get_val(pick("v3.1")),
      v3.2前5显著补全曲线下面积 = get_val(pick("v3.2")),
      v3.3稀疏映射前5融合曲线下面积 = get_val(pick("v3.3")),
      最佳融合曲线下面积 = get_val(fusion),
      最佳融合版本 = get_ver(fusion),
      融合相对原始方法变化 = get_val(fusion) - get_val(original),
      纯scDRS来源 = get_src(scdrs),
      原始方法来源 = get_src(original),
      最佳融合来源 = get_src(fusion),
      stringsAsFactors = FALSE
    )
  })
  out <- bind_rows_fill(rows)
  out[order(out$数据集, out$性状), , drop = FALSE]
}

wide <- make_wide(main2)
write_csv2(wide, file.path(OUT_ROOT, "01_主结果表", "最终主结果宽表.csv"))

cat("正式主结果行数: ", nrow(main2), "；最终版本最佳行数: ", nrow(best_final), "\n", sep = "")

cat(">>> 4. 重画中文主结果图\n")
keep_versions <- c("Baseline", "scDRS-only", "Paired-Fusion", "v2.1", "v2.2", "v3.1", "v3.2", "v3.3")
plot_df <- best_final[best_final$final_version %in% keep_versions, , drop = FALSE]
plot_df$中文性状 <- factor(cn_trait(plot_df$trait), levels = cn_trait(sort(unique(plot_df$trait))))
plot_df$中文方法 <- factor(cn_version(plot_df$final_version), levels = cn_version(keep_versions))
plot_df$中文分面 <- paste(cn_pairing(plot_df$pairing_setting), cn_dataset(plot_df$dataset), sep = " / ")

p_main <- ggplot(plot_df, aes(x = 中文性状, y = value, fill = 中文方法)) +
  geom_col(position = position_dodge(width = 0.82), width = 0.75) +
  facet_wrap(~ 中文分面, scales = "free_x") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(x = NULL, y = "曲线下面积", fill = "方法", title = "主结果：不同方法在各性状上的表现") +
  plot_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p_main, file.path(OUT_ROOT, "02_主结果图", "主结果_各性状方法曲线下面积对比.png"), 12, 7)
save_plot(p_main, file.path(OUT_ROOT, "02_主结果图", "主结果_各性状方法曲线下面积对比.pdf"), 12, 7)

baseline <- best_final[best_final$final_version == "Baseline", c("dataset", "pairing_setting", "trait", "value")]
names(baseline)[names(baseline) == "value"] <- "baseline_auroc"
delta_df <- merge(plot_df[plot_df$method_group != "original_scavenge", ], baseline,
                  by = c("dataset", "pairing_setting", "trait"), all.x = TRUE)
delta_df$delta_vs_baseline <- delta_df$value - delta_df$baseline_auroc
delta_df <- delta_df[is.finite(delta_df$delta_vs_baseline), , drop = FALSE]

p_delta <- ggplot(delta_df, aes(x = 中文性状, y = delta_vs_baseline, fill = 中文方法)) +
  geom_hline(yintercept = 0, color = "grey35", linewidth = 0.4) +
  geom_col(position = position_dodge(width = 0.82), width = 0.75) +
  facet_wrap(~ 中文分面, scales = "free_x") +
  labs(x = NULL, y = "相对原始SCAVENGE的变化", fill = "方法", title = "主结果：相对原始方法的变化") +
  plot_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p_delta, file.path(OUT_ROOT, "02_主结果图", "主结果_相对原始方法变化.png"), 12, 7)
save_plot(p_delta, file.path(OUT_ROOT, "02_主结果图", "主结果_相对原始方法变化.pdf"), 12, 7)
write_csv2(delta_df, file.path(OUT_ROOT, "02_主结果图", "主结果_相对原始方法变化绘图数据.csv"))

heat0 <- plot_df
heat0$面板性状 <- paste(cn_pairing(heat0$pairing_setting), cn_trait(heat0$trait), sep = " / ")
panel_levels <- unique(heat0$面板性状[order(heat0$pairing_setting, heat0$trait)])
method_levels <- cn_version(keep_versions)
heat_grid <- expand.grid(中文方法 = method_levels, 面板性状 = panel_levels, stringsAsFactors = FALSE)
heat <- merge(
  heat_grid,
  heat0[c("中文方法", "面板性状", "value", "dataset", "pairing_setting", "trait", "source_file", "note", "final_version")],
  by = c("中文方法", "面板性状"),
  all.x = TRUE,
  sort = FALSE
)
heat$有结果 <- is.finite(heat$value)
heat$中文方法 <- factor(heat$中文方法, levels = method_levels)
heat$面板性状 <- factor(heat$面板性状, levels = rev(panel_levels))
write_csv2(heat, file.path(OUT_ROOT, "02_主结果图", "主结果_方法性状热图绘图数据.csv"))

p_heat <- ggplot(heat, aes(x = 中文方法, y = 面板性状, fill = value)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(有结果, sprintf("%.3f", value), "")), size = 3) +
  scale_fill_gradient(low = "#F7FBFF", high = "#08519C", limits = c(0, 1), na.value = "grey94", name = "曲线下面积") +
  labs(x = NULL, y = NULL, title = "主结果：方法与性状热图") +
  plot_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
save_plot(p_heat, file.path(OUT_ROOT, "02_主结果图", "主结果_方法性状热图.png"), 12, 8)
save_plot(p_heat, file.path(OUT_ROOT, "02_主结果图", "主结果_方法性状热图.pdf"), 12, 8)

cat(">>> 5. 重画中文消融与参数图\n")
abl <- main2
base <- main2[main2$final_version == "Baseline", c("dataset", "pairing_setting", "trait", "value")]
names(base)[names(base) == "value"] <- "baseline_auroc"
abl <- merge(abl, base, by = c("dataset", "pairing_setting", "trait"), all.x = TRUE)
abl$delta_vs_baseline <- abl$value - abl$baseline_auroc
abl$中文性状 <- cn_trait(abl$trait)
abl$中文版本 <- cn_version(abl$final_version)
abl$中文分面 <- paste(cn_pairing(abl$pairing_setting), cn_dataset(abl$dataset), sep = " / ")
best_abl <- best_by_key(abl, keys = c("dataset", "pairing_setting", "trait", "final_version"))
write_csv2(best_abl, file.path(OUT_ROOT, "03_消融与参数", "消融_按性状最佳表.csv"))

p_ab <- ggplot(best_abl, aes(x = 中文性状, y = value, fill = 中文版本)) +
  geom_col(position = position_dodge(width = 0.82), width = 0.75) +
  facet_wrap(~ 中文分面, scales = "free_x") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(x = NULL, y = "曲线下面积", fill = "版本", title = "消融分析：各版本在不同性状上的表现") +
  plot_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p_ab, file.path(OUT_ROOT, "03_消融与参数", "消融_各性状版本对比.png"), 12, 7)
save_plot(p_ab, file.path(OUT_ROOT, "03_消融与参数", "消融_各性状版本对比.pdf"), 12, 7)

weight <- abl[abl$method_group == "fusion" & is.finite(abl$atac_w), , drop = FALSE]
write_csv2(weight, file.path(OUT_ROOT, "03_消融与参数", "参数_染色质权重扫描表.csv"))
if (nrow(weight) > 0L) {
  weight$中文方法版本 <- cn_version(weight$final_version)
  p_weight <- ggplot(weight, aes(x = atac_w, y = value, color = 中文方法版本, group = interaction(中文方法版本, method_version))) +
    geom_point(size = 2) +
    geom_line(linewidth = 0.75) +
    facet_wrap(~ 中文性状 + 中文分面, scales = "free_y") +
    labs(x = "染色质权重", y = "曲线下面积", color = "方法版本", title = "参数分析：染色质与核糖核酸融合权重扫描") +
    plot_theme()
  save_plot(p_weight, file.path(OUT_ROOT, "03_消融与参数", "参数_融合权重扫描.png"), 12, 7)
  save_plot(p_weight, file.path(OUT_ROOT, "03_消融与参数", "参数_融合权重扫描.pdf"), 12, 7)
}

cat(">>> 6. 重画中文效率图\n")
eff_file <- file.path(FINAL_ROOT, "04_efficiency/efficiency_runtime_memory_table.csv")
if (!file.exists(eff_file)) eff_file <- file.path(FINAL_ROOT, "04_robustness_efficiency/efficiency/efficiency_runtime_memory_table.csv")
eff <- read_csv_if_exists(eff_file)
if (!is.null(eff) && nrow(eff) > 0L) {
  eff <- eff[eff$status == "OK", , drop = FALSE]
  eff$中文方法 <- ifelse(eff$pairing_setting == "paired", "配对融合", "非配对v3.3")
  eff$细胞规模 <- factor(eff$n_cells, levels = c(1000, 10000, 50000, 100000), labels = c("1千", "1万", "5万", "10万"))
  write_csv2(eff, file.path(OUT_ROOT, "04_效率", "效率_运行时间内存曲线下面积表.csv"))

  p_rt <- ggplot(eff, aes(x = n_cells, y = runtime_seconds, color = 中文方法)) +
    geom_line(linewidth = 0.8) + geom_point(size = 2.4) +
    scale_x_continuous(breaks = c(1000, 10000, 50000, 100000), labels = c("1千", "1万", "5万", "10万")) +
    labs(x = "细胞数", y = "运行时间（秒）", color = "方法", title = "效率：运行时间随数据规模变化") +
    plot_theme()
  save_plot(p_rt, file.path(OUT_ROOT, "04_效率", "效率_运行时间.png"), 8, 5)

  p_mem <- ggplot(eff, aes(x = n_cells, y = peak_memory_gb, color = 中文方法)) +
    geom_line(linewidth = 0.8) + geom_point(size = 2.4) +
    scale_x_continuous(breaks = c(1000, 10000, 50000, 100000), labels = c("1千", "1万", "5万", "10万")) +
    labs(x = "细胞数", y = "峰值内存（GB）", color = "方法", title = "效率：峰值内存随数据规模变化") +
    plot_theme()
  save_plot(p_mem, file.path(OUT_ROOT, "04_效率", "效率_峰值内存.png"), 8, 5)

  long_eff <- rbind(
    data.frame(中文方法 = eff$中文方法, 细胞规模 = eff$细胞规模, 指标 = "运行时间（秒）", 数值 = eff$runtime_seconds, 标签 = sprintf("%.1f", eff$runtime_seconds)),
    data.frame(中文方法 = eff$中文方法, 细胞规模 = eff$细胞规模, 指标 = "峰值内存（GB）", 数值 = eff$peak_memory_gb, 标签 = sprintf("%.2f", eff$peak_memory_gb)),
    data.frame(中文方法 = eff$中文方法, 细胞规模 = eff$细胞规模, 指标 = "曲线下面积", 数值 = eff$auroc, 标签 = sprintf("%.3f", eff$auroc)),
    data.frame(中文方法 = eff$中文方法, 细胞规模 = eff$细胞规模, 指标 = "前列命中率", 数值 = eff$acc_top_prevalence, 标签 = sprintf("%.3f", eff$acc_top_prevalence))
  )
  long_eff$指标 <- factor(long_eff$指标, levels = c("运行时间（秒）", "峰值内存（GB）", "曲线下面积", "前列命中率"))
  p_eff <- ggplot(long_eff, aes(x = 细胞规模, y = 数值, color = 中文方法, group = 中文方法)) +
    geom_line(linewidth = 0.8) + geom_point(size = 2.2) +
    geom_text(aes(label = 标签), vjust = -0.65, size = 3, show.legend = FALSE) +
    facet_wrap(~ 指标, scales = "free_y", ncol = 2) +
    labs(x = "细胞数", y = NULL, color = "方法", title = "效率：最终方法规模测试总览") +
    plot_theme()
  save_plot(p_eff, file.path(OUT_ROOT, "04_效率", "效率_总览.png"), 9, 6)
}

cat(">>> 7. 重画中文核糖核酸参考依赖图\n")
ref_file <- file.path(FINAL_ROOT, "05_rna_reference_dependency/rna_reference_dependency_table.csv")
ref <- read_csv_if_exists(ref_file)
if (!is.null(ref) && nrow(ref) > 0L) {
  ref$中文性状 <- cn_trait(ref$trait)
  ref$中文参考 <- ref$rna_reference
  ref$中文参考 <- gsub("Granja RNA", "Granja核糖核酸", ref$中文参考, fixed = TRUE)
  ref$中文参考 <- gsub("NeurIPS CITE RNA source diagnostics", "NeurIPS CITE核糖核酸来源诊断", ref$中文参考, fixed = TRUE)
  ref$中文参考 <- gsub("NeurIPS CITE RNA", "NeurIPS CITE核糖核酸", ref$中文参考, fixed = TRUE)
  ref$中文指标 <- ref$score_type
  ref$中文指标 <- gsub("target_auroc", "目标细胞曲线下面积", ref$中文指标, fixed = TRUE)
  ref$中文指标 <- gsub("cite_score_target_auroc", "CITE得分曲线下面积", ref$中文指标, fixed = TRUE)
  write_csv2(ref, file.path(OUT_ROOT, "05_RNA参考依赖", "核糖核酸参考依赖表.csv"))
  p_ref <- ggplot(ref, aes(x = 中文性状, y = auroc, fill = 中文参考)) +
    geom_col(position = position_dodge(width = 0.82), width = 0.75) +
    facet_wrap(~ 中文指标) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(x = NULL, y = "曲线下面积", fill = "核糖核酸参考", title = "核糖核酸参考依赖：不同参考来源的表现") +
    plot_theme() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot(p_ref, file.path(OUT_ROOT, "05_RNA参考依赖", "核糖核酸参考依赖对比.png"), 10, 6)
}

diag_file <- file.path(FINAL_ROOT, "05_rna_reference_dependency/rna_source_diagnostics_compact.csv")
diag <- read_csv_if_exists(diag_file)
if (!is.null(diag) && nrow(diag) > 0L) {
  diag_long <- data.frame()
  for (metric in intersect(c("cite_score_target_auroc", "pseudo_z_target_auroc", "source_target_weight_auroc"), names(diag))) {
    diag_long <- rbind(diag_long, data.frame(
      中文性状 = cn_trait(diag$trait),
      指标 = metric,
      数值 = suppressWarnings(as.numeric(diag[[metric]])),
      stringsAsFactors = FALSE
    ))
  }
  diag_long$指标 <- factor(
    diag_long$指标,
    levels = c("cite_score_target_auroc", "pseudo_z_target_auroc", "source_target_weight_auroc"),
    labels = c("CITE得分曲线下面积", "伪配对得分曲线下面积", "来源权重曲线下面积")
  )
  write_csv2(diag_long, file.path(OUT_ROOT, "05_RNA参考依赖", "核糖核酸来源诊断长表.csv"))
  p_diag <- ggplot(diag_long, aes(x = 中文性状, y = 数值, fill = 指标)) +
    geom_col(position = position_dodge(width = 0.82), width = 0.75) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(x = NULL, y = "曲线下面积", fill = "诊断指标", title = "核糖核酸来源诊断：目标定位能力") +
    plot_theme() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot(p_diag, file.path(OUT_ROOT, "05_RNA参考依赖", "核糖核酸来源诊断曲线下面积.png"), 9, 5)
}

cat(">>> 8. 重画中文案例研究图\n")
case_dir <- file.path(FINAL_ROOT, "06_case_study/MONO")
case <- read_csv_if_exists(file.path(case_dir, "case_cell_level_scores.csv"))
by_ct <- read_csv_if_exists(file.path(case_dir, "case_high_risk_celltype_table.csv"))
if (!is.null(case) && nrow(case) > 0L) {
  write_csv2(case, file.path(OUT_ROOT, "06_案例研究/MONO", "案例_细胞级风险分数表.csv"))
  case$融合前5 <- factor(ifelse(case$fusion_top, "融合前5%高风险", "其他细胞"), levels = c("其他细胞", "融合前5%高风险"))
  p_case_scatter <- ggplot(case, aes(x = baseline_trs, y = fusion_trs, color = 融合前5)) +
    geom_point(alpha = 0.35, size = 0.7) +
    scale_color_manual(values = c("其他细胞" = "grey70", "融合前5%高风险" = "#C81D25")) +
    labs(x = "原始SCAVENGE风险分数", y = "融合风险分数", color = "细胞分组",
         title = "案例研究：单核细胞原始方法与融合方法风险分数") +
    plot_theme()
  save_plot(p_case_scatter, file.path(OUT_ROOT, "06_案例研究/MONO/图片", "案例_原始与融合风险分数对比.png"), 6, 5)
}
if (!is.null(by_ct) && nrow(by_ct) > 0L) {
  write_csv2(by_ct, file.path(OUT_ROOT, "06_案例研究/MONO", "案例_高风险细胞类型富集表.csv"))
  top_ct <- head(by_ct[order(-by_ct$fusion_top_rate, -by_ct$fusion_top), ], 15)
  p_case_bar <- ggplot(top_ct, aes(x = reorder(cell_type, fusion_top_rate), y = fusion_top_rate)) +
    geom_col(fill = "#C81D25") +
    coord_flip() +
    labs(x = NULL, y = "融合风险分数前5%细胞占比", title = "案例研究：高风险细胞类型富集") +
    plot_theme()
  save_plot(p_case_bar, file.path(OUT_ROOT, "06_案例研究/MONO/图片", "案例_高风险细胞类型富集.png"), 8, 5)
}

meta_file <- file.path(MODULE4_ROOT, "04_evaluation/inputs/granja_cell_meta.rds")
if (!is.null(case) && file.exists(meta_file)) {
  meta <- readRDS(meta_file)
  if (inherits(meta, "data.table")) meta <- as.data.frame(meta)
  if (is.null(rownames(meta)) && "barcode" %in% names(meta)) rownames(meta) <- meta$barcode
  umap_cols <- intersect(c("UMAP_1", "UMAP_2", "umap_1", "umap_2", "UMAP1", "UMAP2"), names(meta))
  if (length(umap_cols) >= 2L) {
    case$umap_1 <- as.numeric(meta[case$barcode, umap_cols[1]])
    case$umap_2 <- as.numeric(meta[case$barcode, umap_cols[2]])
    p_umap <- ggplot(case, aes(x = umap_1, y = umap_2, color = fusion_trs)) +
      geom_point(size = 0.45, alpha = 0.8) +
      scale_color_viridis_c(option = "magma", name = "融合风险分数") +
      labs(x = "UMAP 1", y = "UMAP 2", title = "案例研究：单核细胞融合风险分数在UMAP上的分布") +
      theme_void(base_size = 12) +
      theme(legend.position = "right", plot.title = element_text(face = "bold"))
    save_plot(p_umap, file.path(OUT_ROOT, "06_案例研究/MONO/图片", "案例_融合风险分数UMAP分布.png"), 6, 5)
  }
}

cat(">>> 9. 来源记录与说明\n")
manifest <- unique(main2[c("source_file", "dataset", "pairing_setting", "trait", "method_group", "method_version", "method", "feature_mode", "final_version", "note")])
manifest$状态 <- ifelse(source_exists(manifest$source_file), "已找到", "缺失")
manifest$中文数据集 <- cn_dataset(manifest$dataset)
manifest$中文数据类型 <- cn_pairing(manifest$pairing_setting)
manifest$中文性状 <- cn_trait(manifest$trait)
manifest$中文最终版本 <- cn_version(manifest$final_version)
write_csv2(manifest, file.path(OUT_ROOT, "99_来源记录", "主结果来源记录.csv"))

summary_lines <- c(
  "# 最终输出2说明",
  "",
  paste0("- 生成时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("- 全盘扫描根目录：", WORK_ROOT),
  paste0("- 候选过程文件数：", nrow(scan)),
  paste0("- 已镜像文本过程表数：", sum(!is.na(mirror_map$镜像文件))),
  paste0("- 全盘曲线下面积候选证据行数：", nrow(auc_evidence)),
  paste0("- 补全后正式主结果长表行数：", nrow(main2)),
  paste0("- 按最终版本取最佳后的行数：", nrow(best_final)),
  "",
  "## 目录说明",
  "",
  "- `00_全盘扫描`：全 work 目录候选过程文件索引。",
  "- `01_主结果表`：补全后的正式主结果长表、最佳表和最终宽表。",
  "- `02_主结果图`：中文主结果图和对应绘图数据。",
  "- `03_消融与参数`：中文消融和权重扫描图。",
  "- `04_效率`：中文效率图和效率表。",
  "- `05_RNA参考依赖`：中文核糖核酸参考依赖与来源诊断图。",
  "- `06_案例研究`：中文案例研究表和图。",
  "- `07_过程数据`：镜像的文本过程表和全盘曲线下面积候选证据。",
  "- `99_来源记录`：主结果来源记录。"
)
write_text(summary_lines, file.path(OUT_ROOT, "最终输出2说明.md"))

cat("完成。最终输出2目录:\n  ", OUT_ROOT, "\n", sep = "")
