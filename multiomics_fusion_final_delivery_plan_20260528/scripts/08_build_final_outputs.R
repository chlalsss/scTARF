# Collect human-readable final files into final_outputs/最终输出.
# This script does not modify the original module3/module4 code.

source("D:/RStudio/work/multiomics_fusion_final_delivery_plan_20260528/scripts/_final_delivery_common.R")
ensure_final_dirs()

FINAL_READABLE_ROOT <- file.path(FINAL_ROOT, "最终输出")

final_dirs <- c(
  "00_方法流程",
  "01_主结果表",
  "02_主结果图",
  "03_消融与参数",
  "04_效率",
  "05_RNA参考依赖",
  "06_case_study",
  "07_论文汇报文字",
  "09_补充AUROC表",
  "99_来源记录"
)
invisible(lapply(file.path(FINAL_READABLE_ROOT, final_dirs), dir_ensure))

copy_if_exists <- function(from, to_dir, to_name = basename(from)) {
  if (!file.exists(from)) return(FALSE)
  dir_ensure(to_dir)
  file.copy(from, file.path(to_dir, to_name), overwrite = TRUE)
}

write_final_csv <- function(x, path) {
  dir_ensure(dirname(path))
  txt <- utils::capture.output(utils::write.csv(x, row.names = FALSE, na = ""))
  txt <- enc2utf8(txt)
  con <- tryCatch(file(path, open = "wb"), error = function(e) NULL)
  if (is.null(con)) {
    alt <- sub("\\.csv$", "_新生成_请关闭旧表后替换.csv", path)
    warning("Cannot write ", path, "; writing alternate file: ", alt)
    path <- alt
    con <- file(path, open = "wb")
  }
  on.exit(close(con), add = TRUE)
  writeBin(as.raw(c(0xEF, 0xBB, 0xBF)), con)
  writeBin(charToRaw(paste0(paste(txt, collapse = "\n"), "\n")), con)
  invisible(path)
}

copy_many <- function(paths, to_dir) {
  ok <- vapply(paths, copy_if_exists, logical(1), to_dir = to_dir)
  invisible(ok)
}

make_final_wide_table <- function() {
  main_file <- file.path(FINAL_ROOT, "01_master_tables/main_result_comparison_table.csv")
  if (!file.exists(main_file)) return(NULL)

  x <- read_csv_if_exists(main_file)
  x <- x[x$metric == "AUROC" & is.finite(x$value), , drop = FALSE]
  if (nrow(x) == 0L) return(NULL)
  if (!"final_version" %in% names(x)) {
    x$final_version <- classify_final_version(
      x$method_version, x$method, x$feature_mode, x$method_group, x$pairing_setting
    )
  }
  x$final_version_label <- final_version_label(x$final_version)

  pick_best <- function(y, idx) {
    z <- y[idx, , drop = FALSE]
    if (nrow(z) == 0L) return(NULL)
    z[which.max(z$value), , drop = FALSE]
  }

  scdrs_lookup <- x[x$final_version == "scDRS-only", , drop = FALSE]
  scdrs_by_trait <- lapply(split(seq_len(nrow(scdrs_lookup)), scdrs_lookup$trait), function(ii) {
    pick_best(scdrs_lookup, ii)
  })

  main_x <- x[x$pairing_setting != "rna_reference_only", , drop = FALSE]
  group_key <- paste(main_x$dataset, main_x$pairing_setting, main_x$trait, sep = "\r")
  rows <- lapply(split(seq_len(nrow(main_x)), group_key), function(ii) {
    y <- main_x[ii, , drop = FALSE]
    original <- pick_best(y, y$final_version == "Baseline")
    scdrs <- pick_best(y, y$final_version == "scDRS-only")
    if (is.null(scdrs) && identical(y$pairing_setting[1], "non_paired")) {
      scdrs <- scdrs_by_trait[[y$trait[1]]]
    }
    paired <- pick_best(y, y$final_version == "Paired-Fusion")
    v21 <- pick_best(y, y$final_version == "v2.1")
    v22 <- pick_best(y, y$final_version == "v2.2")
    v31 <- pick_best(y, y$final_version == "v3.1")
    v32 <- pick_best(y, y$final_version == "v3.2")
    v33 <- pick_best(y, y$final_version == "v3.3")
    fusion <- pick_best(y, y$final_version %in% c("Paired-Fusion", "v2.1", "v2.2", "v3.1", "v3.2", "v3.3"))

    get_val <- function(z) if (is.null(z)) NA_real_ else z$value[1]
    get_version <- function(z) if (is.null(z)) NA_character_ else z$final_version[1]
    get_method <- function(z) if (is.null(z)) NA_character_ else z$final_version_label[1]
    get_atac_w <- function(z) if (is.null(z)) NA_real_ else z$atac_w[1]
    get_src <- function(z) if (is.null(z)) NA_character_ else z$source_file[1]

    data.frame(
      数据集 = y$dataset[1],
      数据类型 = y$pairing_setting[1],
      性状 = y$trait[1],
      原始_SCAVENGE_AUROC = get_val(original),
      纯_scDRS_AUROC = get_val(scdrs),
      Paired_Fusion_AUROC = get_val(paired),
      `v2.1_基础伪配对融合_AUROC` = get_val(v21),
      `v2.2_trait_aware_FDR融合_AUROC` = get_val(v22),
      `v3.1_top1消融_AUROC` = get_val(v31),
      `v3.2_top5显著补全_AUROC` = get_val(v32),
      `v3.3_sparse_top5融合_AUROC` = get_val(v33),
      最佳融合_AUROC = get_val(fusion),
      最佳融合版本 = get_version(fusion),
      最佳融合方法 = get_method(fusion),
      最佳融合_ATAC权重 = get_atac_w(fusion),
      融合相对SCAVENGE变化 = get_val(fusion) - get_val(original),
      纯scDRS来源 = get_src(scdrs),
      原始SCAVENGE来源 = get_src(original),
      最佳融合来源 = get_src(fusion),
      stringsAsFactors = FALSE
    )
  })

  out <- bind_rows_fill(rows)
  out <- out[order(out$数据集, out$性状), , drop = FALSE]
  write_final_csv(out, file.path(FINAL_READABLE_ROOT, "01_主结果表", "最终主结果宽表.csv"))
  out
}

make_final_wide_table()

copy_many(c(
  file.path(FINAL_ROOT, "00_method_flow/method_protocol.md"),
  file.path(FINAL_ROOT, "00_method_flow/method_flowchart.mmd"),
  file.path(FINAL_ROOT, "00_method_flow/method_steps_table.csv"),
  file.path(FINAL_ROOT, "00_method_flow/method_version_definitions.csv"),
  file.path(FINAL_ROOT, "00_method_flow/method_version_definitions.md")
), file.path(FINAL_READABLE_ROOT, "00_方法流程"))

copy_many(c(
  file.path(FINAL_ROOT, "01_master_tables/main_result_comparison_table.csv"),
  file.path(FINAL_ROOT, "01_master_tables/main_result_best_by_method_version.csv"),
  file.path(FINAL_ROOT, "01_master_tables/main_result_best_by_final_version.csv"),
  file.path(FINAL_ROOT, "01_master_tables/method_version_definitions.csv")
), file.path(FINAL_READABLE_ROOT, "01_主结果表"))

copy_many(c(
  file.path(FINAL_ROOT, "02_main_figures/fig_main_auroc_by_trait.png"),
  file.path(FINAL_ROOT, "02_main_figures/fig_main_delta_vs_baseline.png"),
  file.path(FINAL_ROOT, "02_main_figures/fig_main_method_heatmap.png"),
  file.path(FINAL_ROOT, "02_main_figures/main_delta_vs_baseline_plot_data.csv"),
  file.path(FINAL_ROOT, "02_main_figures/main_method_heatmap_plot_data.csv")
), file.path(FINAL_READABLE_ROOT, "02_主结果图"))

copy_many(c(
  file.path(FINAL_ROOT, "09_missing_master_table/bmmc_paired_scdrs_only_auroc.csv"),
  file.path(FINAL_ROOT, "09_missing_master_table/nonpaired_trait_screening_scdrs_auroc.csv")
), file.path(FINAL_READABLE_ROOT, "09_补充AUROC表"))

copy_many(c(
  file.path(FINAL_ROOT, "03_ablation_parameter/ablation_result_table.csv"),
  file.path(FINAL_ROOT, "03_ablation_parameter/ablation_best_by_trait.csv"),
  file.path(FINAL_ROOT, "03_ablation_parameter/parameter_weight_sweep_table.csv"),
  file.path(FINAL_ROOT, "03_ablation_parameter/top5_significant_feature_filter_summary.csv"),
  file.path(FINAL_ROOT, "03_ablation_parameter/fig_ablation_by_trait.png"),
  file.path(FINAL_ROOT, "03_ablation_parameter/fig_weight_sweep.png"),
  file.path(FINAL_ROOT, "03_ablation_parameter/ablation_interpretation_template.md")
), file.path(FINAL_READABLE_ROOT, "03_消融与参数"))

copy_many(c(
  file.path(FINAL_ROOT, "04_robustness_efficiency/efficiency/efficiency_elapsed_step_times_from_demo_runs.csv"),
  file.path(FINAL_ROOT, "04_robustness_efficiency/efficiency/efficiency_total_runtime_from_demo_runs.csv"),
  file.path(FINAL_ROOT, "04_robustness_efficiency/efficiency/manual_efficiency_runs_template.csv"),
  file.path(FINAL_ROOT, "04_robustness_efficiency/efficiency/efficiency_runtime_memory_table.csv"),
  file.path(FINAL_ROOT, "04_robustness_efficiency/efficiency/fig_efficiency_step_time_demo.png"),
  file.path(FINAL_ROOT, "04_robustness_efficiency/efficiency/fig_efficiency_total_time_demo.png"),
  file.path(FINAL_ROOT, "04_robustness_efficiency/efficiency/fig_efficiency_runtime_scaling.png"),
  file.path(FINAL_ROOT, "04_robustness_efficiency/efficiency/fig_efficiency_memory_scaling.png"),
  file.path(FINAL_ROOT, "04_robustness_efficiency/efficiency/efficiency_notes.md")
), file.path(FINAL_READABLE_ROOT, "04_效率"))

copy_many(c(
  file.path(FINAL_ROOT, "05_rna_reference_dependency/rna_reference_dependency_table.csv"),
  file.path(FINAL_ROOT, "05_rna_reference_dependency/rna_source_diagnostics_compact.csv"),
  file.path(FINAL_ROOT, "05_rna_reference_dependency/fig_rna_reference_dependency.png"),
  file.path(FINAL_ROOT, "05_rna_reference_dependency/fig_rna_source_diagnostic_auroc.png"),
  file.path(FINAL_ROOT, "05_rna_reference_dependency/rna_reference_dependency_interpretation.md")
), file.path(FINAL_READABLE_ROOT, "05_RNA参考依赖"))

case_root <- file.path(FINAL_ROOT, "06_case_study")
if (dir.exists(case_root)) {
  case_files <- list.files(case_root, recursive = TRUE, full.names = TRUE)
  for (f in case_files[file.exists(case_files)]) {
    rel <- sub(paste0("^", gsub("\\\\", "/", normalizePath(case_root, winslash = "/", mustWork = TRUE)), "/?"),
               "", normalizePath(f, winslash = "/", mustWork = TRUE))
    to <- file.path(FINAL_READABLE_ROOT, "06_case_study", rel)
    dir_ensure(dirname(to))
    file.copy(f, to, overwrite = TRUE)
  }
}

copy_many(c(
  file.path(FINAL_ROOT, "07_manuscript_text/methods_final_workflow.md"),
  file.path(FINAL_ROOT, "07_manuscript_text/method_version_definitions.md"),
  file.path(FINAL_ROOT, "07_manuscript_text/results_main_comparison.md"),
  file.path(FINAL_ROOT, "07_manuscript_text/results_ablation.md"),
  file.path(FINAL_ROOT, "07_manuscript_text/results_efficiency.md"),
  file.path(FINAL_ROOT, "07_manuscript_text/results_rna_reference_dependency.md"),
  file.path(FINAL_ROOT, "07_manuscript_text/case_study.md")
), file.path(FINAL_READABLE_ROOT, "07_论文汇报文字"))

copy_many(c(
  file.path(FINAL_ROOT, "99_source_manifest/result_source_manifest.csv"),
  file.path(FINAL_ROOT, "99_source_manifest/source_manifest_master_results.csv")
), file.path(FINAL_READABLE_ROOT, "99_来源记录"))

readme <- c(
  "# 最终输出说明",
  "",
  "这个目录集中放可以给人看的最终整理结果；中间过程表仍保留在上一级各模块目录中。",
  "",
  "## 00_方法流程",
  "- `method_protocol.md`：最终方法文字说明。",
  "- `method_flowchart.mmd`：方法流程图源码，可以用 Mermaid 转成图片。",
  "- `method_steps_table.csv`：方法步骤表。",
  "",
  "## 01_主结果表",
  "- `最终主结果宽表.csv`：最推荐先看的主表。每个数据集/性状一行，列出原始 SCAVENGE、scDRS-only、Paired-Fusion、v2.1、v2.2、v3.1、v3.2、v3.3 和最佳融合结果。",
  "- `method_version_definitions.csv`：最终版本号定义表，说明 Baseline、scDRS-only、Paired-Fusion、v2.1、v2.2、v3.1、v3.2、v3.3 分别对应哪些原始结果。",
  "- `main_result_comparison_table.csv`：完整可追溯长表，保留所有方法版本、参数和来源。",
  "- `main_result_best_by_method_version.csv`：每个方法版本的最佳结果。",
  "- `main_result_best_by_final_version.csv`：按最终版本体系合并后的最佳结果，例如 v2.2 会合并 neighborhood-FDR 和 trait-screening 扩展运行。",
  "",
  "## 02_主结果图",
  "- `fig_main_auroc_by_trait.png`：不同方法的 AUROC 对比。",
  "- `fig_main_delta_vs_baseline.png`：融合方法相对原始 SCAVENGE 的变化。",
  "- `fig_main_method_heatmap.png`：方法 × 性状 AUROC 热图。",
  "",
  "## 03_消融与参数",
  "- 消融表、权重 sweep 表、消融图和解释模板。",
  "",
  "## 04_效率",
  "- demo 运行时间表和图。",
  "- `manual_efficiency_runs_template.csv`：正式 1K/10K/50K/100K 规模实验模板。",
  "",
  "## 05_RNA参考依赖",
  "- 不同 RNA reference 的 AUROC 和 source diagnostics 对比。",
  "",
  "## 06_case_study",
  "- 默认 MONO case 的细胞级分数表、高风险细胞类型表和图。",
  "",
  "## 07_论文汇报文字",
  "- methods/results/case 的文字模板。",
  "",
  "## 99_来源记录",
  "- 结果来源 manifest，用来追溯每个最终文件来自哪个原始结果。"
)
writeLines(readme, file.path(FINAL_READABLE_ROOT, "最终输出说明.md"), useBytes = TRUE)

cat("Final readable outputs collected under:\n  ", FINAL_READABLE_ROOT, "\n", sep = "")
