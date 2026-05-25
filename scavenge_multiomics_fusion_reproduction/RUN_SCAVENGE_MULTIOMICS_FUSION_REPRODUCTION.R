# SCAVENGE 多组学融合改进实验完整复现总入口
#
# 在 RStudio 中运行:
#   source("D:/RStudio/work/scavenge_multiomics_fusion_reproduction/RUN_SCAVENGE_MULTIOMICS_FUSION_REPRODUCTION.R")
#
# 输出位置:
#   D:/RStudio/work/scavenge_multiomics_fusion_reproduction/reproduction_outputs/
#
# 说明:
#   1. 本脚本会把本目录中的代码复制到 reproduction_outputs/patched_scripts/，
#      并把原始脚本里的硬编码输出路径改到 reproduction_outputs/ 下。
#   2. 大型原始数据不复制到本目录；需要提前按 DATASETS.md 准备。
#   3. scDRS 和 ATAC dropout sweep 耗时较长，完整运行可能需要数小时。

# =========================
# 0. 用户可修改配置
# =========================

# 大型外部数据根目录。只作为输入读取，不作为输出目录。
EXTERNAL_SCAVENGE_ROOT <- "D:/SCAVENGE-reproducibility-main"

# 健康造血 ATAC 输入。
BLOOD_ATAC_RDS <- file.path(
  EXTERNAL_SCAVENGE_ROOT,
  "blood_data",
  "scATAC-Healthy-Hematopoiesis-191120.rds"
)

# 已生成的单模态 SCAVENGE 结果目录。
# 如果没有该目录，脚本会跳过单模态图表重绘，只继续跑 BMMC 全流程。
BLOOD_SCAVENGE_OUT <- file.path(
  EXTERNAL_SCAVENGE_ROOT,
  "blood_data",
  "scavenge_out_small_all_traits"
)

# BMMC fine-mapped GWAS 输入。
GWAS_TSV <- file.path(
  EXTERNAL_SCAVENGE_ROOT,
  "blood_data",
  "CS.trans_all.PP01.strong_atac_overlap.tsv"
)

# 是否联网重新下载 scDRS gene set。
# 默认 FALSE：本提交目录已包含 bmmc_3traits.gs，避免验收机器无网络导致中断。
RUN_GENESET_DOWNLOAD <- FALSE

# 是否运行完整 ATAC dropout sweep。
# 完整配置为 8 个 p_drop x 3 seeds x 2 weights x 3 traits，耗时较长。
RUN_FULL_ROBUSTNESS_SWEEP <- TRUE

# =========================
# 1. 路径与工具函数
# =========================

get_this_file <- function() {
  frames <- sys.frames()
  for (i in rev(seq_along(frames))) {
    if (!is.null(frames[[i]]$ofile)) return(normalizePath(frames[[i]]$ofile, winslash = "/"))
  }
  normalizePath("RUN_SCAVENGE_MULTIOMICS_FUSION_REPRODUCTION.R", winslash = "/", mustWork = FALSE)
}

THIS_FILE <- get_this_file()
PROJECT_DIR <- dirname(THIS_FILE)
OUT_ROOT <- file.path(PROJECT_DIR, "reproduction_outputs")
BMMC_OUT <- file.path(OUT_ROOT, "BMMC")
BLOOD_OUT <- file.path(OUT_ROOT, "blood_atac_scavenge")
PATCH_ROOT <- file.path(OUT_ROOT, "patched_scripts")
PATCH_M3 <- file.path(PATCH_ROOT, "module3_multiome")
LOG_DIR <- file.path(OUT_ROOT, "logs")

dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)
dir.create(BMMC_OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(BMMC_OUT, "00_precheck"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(BMMC_OUT, "01_data_qc"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(BMMC_OUT, "03_scavenge_baseline"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(BMMC_OUT, "04_dual_modality"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(BMMC_OUT, "04_dual_modality", "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(BMMC_OUT, "05_robustness"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(BMMC_OUT, "05_robustness", "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(BLOOD_OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(BLOOD_OUT, "sparsity_output"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(BLOOD_OUT, "fig3_like"), recursive = TRUE, showWarnings = FALSE)
dir.create(PATCH_ROOT, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

r_path <- function(path) normalizePath(path, winslash = "/", mustWork = FALSE)

write_section <- function(...) {
  cat("\n", paste(rep("=", 78), collapse = ""), "\n", sep = "")
  cat(..., "\n", sep = "")
  cat(paste(rep("=", 78), collapse = ""), "\n", sep = "")
}

copy_dir <- function(from, to) {
  if (dir.exists(to)) unlink(to, recursive = TRUE, force = TRUE)
  dir.create(dirname(to), recursive = TRUE, showWarnings = FALSE)
  ok <- file.copy(from, dirname(to), recursive = TRUE, overwrite = TRUE)
  if (!ok) stop("Failed to copy directory: ", from, call. = FALSE)
}

patch_text_file <- function(path) {
  if (!grepl("\\.(R|r|py|md|txt)$", path)) return(invisible(FALSE))
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")

  replacements <- c(
    "D:/SCAVENGE-reproducibility-main/BMMC" = r_path(BMMC_OUT),
    "D:\\\\SCAVENGE-reproducibility-main\\\\BMMC" = r_path(BMMC_OUT),
    "D:/RStudio/work/module3_multiome" = r_path(PATCH_M3),
    "D:\\\\RStudio\\\\work\\\\module3_multiome" = r_path(PATCH_M3),
    "D:/SCAVENGE-reproducibility-main/blood_data/scavenge_out_small_all_traits" = r_path(file.path(BLOOD_OUT, "scavenge_out_small_all_traits")),
    "D:/SCAVENGE-reproducibility-main/sparsity_output" = r_path(file.path(BLOOD_OUT, "sparsity_output"))
  )
  for (old in names(replacements)) {
    lines <- gsub(old, replacements[[old]], lines, fixed = TRUE)
  }

  # 旧 module3 脚本常用 sys.frame(1)$ofile 推断 ../_bootstrap.R。
  # 被本总控脚本间接 source() 时，这个推断会指向总控脚本目录而不是子脚本目录。
  # 在运行副本中统一改成绝对 bootstrap 路径。
  lines <- gsub(
    'source(file.path(dirname(sys.frame(1)$ofile %||% "."), "..", "_bootstrap.R"))',
    paste0('source("', r_path(file.path(PATCH_M3, "_bootstrap.R")), '")'),
    lines,
    fixed = TRUE
  )

  # 多组学融合改进实验需要同时生成融合 w=0.7 和纯 ATAC 对照 w=1.0。
  lines <- gsub("ATAC_W_LIST <- c\\(1\\.0\\)", "ATAC_W_LIST <- c(0.7, 1.0)", lines)

  # 单模态绘图脚本从外部大型 SCAVENGE 结果目录读输入，
  # 但把新图片写到 reproduction_outputs 下。
  if (basename(path) == "plots.R") {
    lines <- gsub(
      '^out_dir <- .*$',
      paste0('out_dir <- "', r_path(BLOOD_SCAVENGE_OUT), '"'),
      lines
    )
    lines <- gsub(
      'fig_dir <- file.path\\(out_dir, "fig3_like"\\)',
      paste0('fig_dir <- "', r_path(file.path(BLOOD_OUT, "fig3_like")), '"'),
      lines
    )
  }

  if (basename(path) == "sparsity_scATAC.R") {
    lines <- gsub(
      "library\\(ggplot2\\)",
      paste0(
        "library(ggplot2)\n",
        "dir.create(\"", r_path(file.path(BLOOD_OUT, "sparsity_output")),
        "\", recursive = TRUE, showWarnings = FALSE)"
      ),
      lines
    )
  }

  # 对外提交不使用 MultiOmicsFusion 作为算法名称。运行副本里将历史内部文件名前缀
  # seed_MultiOmicsFusion / TRS_MultiOmicsFusion / eval_MultiOmicsFusion 等统一改为 MultiOmicsFusion。
  lines <- gsub("MultiOmicsFusion", "MultiOmicsFusion", lines, fixed = TRUE)
  lines <- gsub("MultiOmicsFusion", "MultiOmicsFusion", lines, fixed = TRUE)

  # 完整压力测试配置，和论文实验设置一致。
  if (basename(path) == "01_atac_sparsify.R" && RUN_FULL_ROBUSTNESS_SWEEP) {
    lines <- gsub(
      "P_DROP_LIST <- c\\([^\\n]+\\)",
      "P_DROP_LIST <- c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.8, 0.9)",
      lines
    )
    lines <- gsub(
      "SEED_LIST   <- c\\([^\\n]+\\)",
      "SEED_LIST   <- c(1, 2, 3)",
      lines
    )
  }

  writeLines(lines, path, useBytes = TRUE)
  invisible(TRUE)
}

run_r <- function(step_id, label, script) {
  write_section(step_id, " - ", label)
  script <- normalizePath(script, winslash = "/", mustWork = TRUE)
  log_file <- file.path(LOG_DIR, paste0(step_id, "_", tools::file_path_sans_ext(basename(script)), ".log"))
  con <- file(log_file, open = "wt", encoding = "UTF-8")
  sink(con, split = TRUE)
  sink(con, type = "message")
  start <- Sys.time()
  status <- "OK"
  err <- ""
  elapsed <- NA_real_
  res <- tryCatch(
    {
      step_env <- new.env(parent = globalenv())
      elapsed <- unname(system.time(source(script, encoding = "UTF-8", chdir = TRUE, local = step_env))[["elapsed"]])
      rm(step_env)
      invisible(gc())
      TRUE
    },
    error = function(e) {
      status <<- "ERROR"
      err <<- conditionMessage(e)
      FALSE
    }
  )
  end <- Sys.time()
  sink(type = "message")
  sink()
  close(con)

  rec <- data.frame(
    step_id = step_id,
    label = label,
    type = "R",
    script = script,
    status = status,
    elapsed_seconds = elapsed,
    start_time = format(start, "%Y-%m-%d %H:%M:%S"),
    end_time = format(end, "%Y-%m-%d %H:%M:%S"),
    log_file = log_file,
    error_message = err,
    stringsAsFactors = FALSE
  )
  if (!res) {
    print(rec)
    stop("Step failed: ", step_id, " - ", label, "\nSee log: ", log_file, call. = FALSE)
  }
  rec
}

run_skip_record <- function(step_id, label, reason) {
  write_section(step_id, " - ", label, " [skip]")
  data.frame(
    step_id = step_id,
    label = label,
    type = "R",
    script = "",
    status = "SKIPPED_EXISTING_OUTPUTS",
    elapsed_seconds = 0,
    start_time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    end_time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    log_file = "",
    error_message = reason,
    stringsAsFactors = FALSE
  )
}

all_files_exist <- function(paths) {
  all(file.exists(paths)) && all(file.info(paths)$size > 0)
}

run_py <- function(step_id, label, script) {
  write_section(step_id, " - ", label)
  script <- normalizePath(script, winslash = "/", mustWork = TRUE)
  py <- Sys.which("python")
  if (!nzchar(py)) stop("Cannot find python in PATH.", call. = FALSE)
  log_file <- file.path(LOG_DIR, paste0(step_id, "_", tools::file_path_sans_ext(basename(script)), ".log"))
  start <- Sys.time()
  t0 <- proc.time()[["elapsed"]]
  code <- system2(py, args = shQuote(script), stdout = log_file, stderr = log_file)
  elapsed <- proc.time()[["elapsed"]] - t0
  end <- Sys.time()
  rec <- data.frame(
    step_id = step_id,
    label = label,
    type = "Python",
    script = script,
    status = if (identical(code, 0L)) "OK" else "ERROR",
    elapsed_seconds = elapsed,
    start_time = format(start, "%Y-%m-%d %H:%M:%S"),
    end_time = format(end, "%Y-%m-%d %H:%M:%S"),
    log_file = log_file,
    error_message = if (identical(code, 0L)) "" else paste("exit code", code),
    stringsAsFactors = FALSE
  )
  if (!identical(code, 0L)) {
    print(rec)
    stop("Step failed: ", step_id, " - ", label, "\nSee log: ", log_file, call. = FALSE)
  }
  rec
}

require_file <- function(path, why) {
  if (!file.exists(path)) {
    stop("Missing required input: ", path, "\n", why, "\nSee DATASETS.md.", call. = FALSE)
  }
}

# =========================
# 2. 生成已改路径运行副本
# =========================

write_section("Preparing patched runnable copy")

if (dir.exists(PATCH_M3)) unlink(PATCH_M3, recursive = TRUE, force = TRUE)
dir.create(PATCH_M3, recursive = TRUE, showWarnings = FALSE)

# 组装 module3 运行副本。
file.copy(
  file.path(PROJECT_DIR, "02_bmmc_atac_only_baseline", c("_bootstrap.R", "_config.R", "install_deps.R", "README.md")),
  PATCH_M3,
  overwrite = TRUE
)
copy_dir(file.path(PROJECT_DIR, "02_bmmc_atac_only_baseline", "00_precheck"), file.path(PATCH_M3, "00_precheck"))
copy_dir(file.path(PROJECT_DIR, "02_bmmc_atac_only_baseline", "01_data_qc"), file.path(PATCH_M3, "01_data_qc"))
copy_dir(file.path(PROJECT_DIR, "02_bmmc_atac_only_baseline", "03_scavenge_baseline"), file.path(PATCH_M3, "03_scavenge_baseline"))
copy_dir(file.path(PROJECT_DIR, "03_multiomics_fusion", "04_dual_modality"), file.path(PATCH_M3, "04_dual_modality"))
copy_dir(file.path(PROJECT_DIR, "04_modality_degradation_robustness", "05_robustness"), file.path(PATCH_M3, "05_robustness"))
copy_dir(file.path(PROJECT_DIR, "05_paper_figures", "07_paper_figures"), file.path(PATCH_M3, "07_paper_figures"))

# 复制并修补ATAC 单模态复现实验辅助脚本。
PATCH_BLOOD <- file.path(PATCH_ROOT, "blood_atac_reproduction")
if (dir.exists(PATCH_BLOOD)) unlink(PATCH_BLOOD, recursive = TRUE, force = TRUE)
dir.create(PATCH_BLOOD, recursive = TRUE, showWarnings = FALSE)
file.copy(
  file.path(PROJECT_DIR, "01_scavenge_atac_reproduction", c("sparsity_scATAC.R", "plots.R", "plot_fig3f_heatmap.R")),
  PATCH_BLOOD,
  overwrite = TRUE
)

for (f in list.files(PATCH_ROOT, recursive = TRUE, full.names = TRUE)) patch_text_file(f)

Sys.setenv(
  M3_DATA_DIR = r_path(BMMC_OUT),
  M3_OUT_DIR = r_path(BMMC_OUT),
  M3_GWAS_TSV = r_path(GWAS_TSV),
  M3_AUTO_DOWNLOAD_LIMIT_MB = "5000"
)

cat("Project dir: ", PROJECT_DIR, "\n", sep = "")
cat("Output root: ", OUT_ROOT, "\n", sep = "")
cat("Patched scripts: ", PATCH_ROOT, "\n", sep = "")

require_file(GWAS_TSV, "BMMC precheck needs the fine-mapped GWAS variant table.")

# =========================
# 3. 健康造血 ATAC 单模态复现辅助输出
# =========================

records <- list()

write_section("ATAC single-modality reproduction")

if (file.exists(BLOOD_ATAC_RDS)) {
  # 生成稀疏度统计。该脚本只依赖原始 scATAC rds。
  records[[length(records) + 1]] <- run_r(
    "S1_1",
    "Healthy hematopoiesis scATAC sparsity statistics",
    file.path(PATCH_BLOOD, "sparsity_scATAC.R")
  )
} else {
  warning("Skip ATAC sparsity: missing BLOOD_ATAC_RDS: ", BLOOD_ATAC_RDS)
}

if (dir.exists(BLOOD_SCAVENGE_OUT)) {
  # 不复制大型结果；只把已有结果目录中的文件作为输入，重绘论文图。
  dir.create(file.path(BLOOD_OUT, "scavenge_out_small_all_traits"), recursive = TRUE, showWarnings = FALSE)
  # plots.R 默认画一个 trait；保留原脚本，可手动改 trait_to_plot 后重复 source。
  records[[length(records) + 1]] <- run_r(
    "S1_2",
    "Replot SCAVENGE cell-type/z-score/TRS figures from single-modality outputs",
    file.path(PATCH_BLOOD, "plots.R")
  )
} else {
  warning("Skip ATAC figures: missing BLOOD_SCAVENGE_OUT: ", BLOOD_SCAVENGE_OUT)
}

# =========================
# 4. BMMC ATAC-only baseline
# =========================

write_section("BMMC ATAC-only baseline")

records[[length(records) + 1]] <- run_r("S2_1", "Fetch/extract BMMC metadata", file.path(PATCH_M3, "00_precheck", "01_fetch_metadata.R"))
records[[length(records) + 1]] <- run_r("S2_2", "Extract BMMC peak coordinates", file.path(PATCH_M3, "00_precheck", "02b_extract_peaks_from_h5ad.R"))
records[[length(records) + 1]] <- run_r("S2_3", "Cell-type coverage precheck", file.path(PATCH_M3, "00_precheck", "02_celltype_coverage.R"))
records[[length(records) + 1]] <- run_r("S2_4", "Variant-peak overlap precheck", file.path(PATCH_M3, "00_precheck", "03_variant_peak_overlap.R"))
bmmc_matrix_outputs <- file.path(
  BMMC_OUT,
  "01_data_qc",
  c(
    "BMMC_atac_counts.rds",
    "BMMC_rna_counts.rds",
    "BMMC_cell_meta.rds",
    "BMMC_atac_features.rds",
    "BMMC_rna_features.rds"
  )
)
if (all_files_exist(bmmc_matrix_outputs)) {
  records[[length(records) + 1]] <- run_skip_record(
    "S2_5",
    "Extract BMMC matrices from h5ad",
    "Required BMMC matrix RDS files already exist under reproduction_outputs/BMMC/01_data_qc."
  )
} else {
  records[[length(records) + 1]] <- run_r("S2_5", "Extract BMMC matrices from h5ad", file.path(PATCH_M3, "01_data_qc", "01_extract_matrices_from_h5ad.R"))
}
records[[length(records) + 1]] <- run_r("S2_6", "Build ATAC SummarizedExperiment", file.path(PATCH_M3, "03_scavenge_baseline", "01_build_atac_se.R"))
records[[length(records) + 1]] <- run_r("S2_7", "Build background peaks and ATAC LSI", file.path(PATCH_M3, "03_scavenge_baseline", "02_bg_peaks_and_lsi.R"))
records[[length(records) + 1]] <- run_r("S2_8", "Run ATAC-only SCAVENGE baseline traits", file.path(PATCH_M3, "03_scavenge_baseline", "03_run_baseline_traits.R"))
records[[length(records) + 1]] <- run_r("S2_9", "Check ATAC-only baseline enrichment", file.path(PATCH_M3, "03_scavenge_baseline", "04_check_baseline.R"))

# =========================
# 5. Paired multi-omics fusion improvement
# =========================

write_section("Paired multi-omics fusion improvement")

if (RUN_GENESET_DOWNLOAD) {
  records[[length(records) + 1]] <- run_py("S3_1", "Download official scDRS gene set", file.path(PATCH_M3, "04_dual_modality", "01_scdrs_python", "download_gs.py"))
} else {
  cat("Skip gene-set download; using bundled bmmc_3traits.gs.\n")
}

records[[length(records) + 1]] <- run_r("S3_2", "Prepare BMMC RNA h5ad for scDRS", file.path(PATCH_M3, "04_dual_modality", "01_scdrs_python", "prepare_h5ad.R"))
records[[length(records) + 1]] <- run_py("S3_3", "Run official scDRS compute-score and downstream", file.path(PATCH_M3, "04_dual_modality", "01_scdrs_python", "run_scdrs.py"))
records[[length(records) + 1]] <- run_r("S3_4", "Build ATAC/RNA blocks and tune fusion inputs", file.path(PATCH_M3, "04_dual_modality", "02_fusion_autotune.R"))
records[[length(records) + 1]] <- run_r("S3_5", "Select scDRS seed cells", file.path(PATCH_M3, "04_dual_modality", "04_seed_scdrs.R"))
records[[length(records) + 1]] <- run_r("S3_6", "Build weighted KNN graphs for w=0.7 and w=1.0", file.path(PATCH_M3, "04_dual_modality", "03_build_weighted_knn.R"))
records[[length(records) + 1]] <- run_r("S3_7", "Run RWR/TRS for w=0.7 and w=1.0", file.path(PATCH_M3, "04_dual_modality", "05_rwr_trs.R"))
records[[length(records) + 1]] <- run_r("S3_8", "Evaluate fusion AUROC", file.path(PATCH_M3, "04_dual_modality", "06_evaluate_auroc.R"))
records[[length(records) + 1]] <- run_r("S3_9", "Plot TRS UMAP comparison figures", file.path(PATCH_M3, "07_paper_figures", "plot_trs_umap_comparison.R"))

# =========================
# 6. ATAC degradation robustness
# =========================

write_section("ATAC degradation robustness")

records[[length(records) + 1]] <- run_r("S4_1", "Run full ATAC dropout robustness sweep", file.path(PATCH_M3, "05_robustness", "01_atac_sparsify.R"))
records[[length(records) + 1]] <- run_r("S4_2", "Plot robustness curves and fusion advantage", file.path(PATCH_M3, "05_robustness", "04_plot_robustness.R"))

# =========================
# 7. 汇总输出索引
# =========================

write_section("Writing run summary")

summary_df <- do.call(rbind, records)
write.csv(summary_df, file.path(OUT_ROOT, "run_step_summary.csv"), row.names = FALSE)

sink(file.path(OUT_ROOT, "session_info.txt"))
print(sessionInfo())
sink()

all_outputs <- list.files(OUT_ROOT, recursive = TRUE, full.names = TRUE)
output_index <- data.frame(
  path = all_outputs,
  size_bytes = ifelse(file.exists(all_outputs), file.info(all_outputs)$size, NA_real_),
  stringsAsFactors = FALSE
)
write.csv(output_index, file.path(OUT_ROOT, "output_file_index.csv"), row.names = FALSE)

cat("\n完整复现流程结束。\n")
cat("输出目录: ", OUT_ROOT, "\n", sep = "")
cat("步骤日志: ", LOG_DIR, "\n", sep = "")
cat("步骤摘要: ", file.path(OUT_ROOT, "run_step_summary.csv"), "\n", sep = "")
cat("文件索引: ", file.path(OUT_ROOT, "output_file_index.csv"), "\n", sep = "")
