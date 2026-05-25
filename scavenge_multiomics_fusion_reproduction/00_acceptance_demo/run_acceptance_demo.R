# 毕业设计软件验收演示脚本
#
# 在 RStudio 中运行：
#   source("D:/RStudio/work/run_acceptance_demo.R")
#
# 设计目的：
#   1. 每次运行都创建一个带时间戳的独立结果文件夹；
#   2. 保留每一步控制台日志、运行耗时、R 环境信息；
#   3. 现场真正运行关键算法步骤，而不是只做结果可视化；
#   4. 为了控制验收演示时间，复用已经准备好的 ATAC/RNA 中间矩阵和 scDRS 结果，
#      但现场重新运行 MONO 的 seed 选择、带权互最近邻图、RWR/TRS 和 AUROC 评估。
#
# 每次运行会生成：
#   D:/RStudio/work/acceptance_demo_runs/run_YYYYMMDD_HHMMSS/

# =========================
# 1. 基础路径设置
# =========================

# 当前工作区根目录。
ROOT_DIR <- "D:/RStudio/work"

# 所有验收演示批次统一保存到该目录。
RUN_ROOT <- file.path(ROOT_DIR, "acceptance_demo_runs")

# BMMC 数据和正式实验结果所在目录。
# 本脚本只读取其中已有的输入和基准结果，不在正式目录里重跑核心算法。
BMMC_DIR <- "D:/SCAVENGE-reproducibility-main/BMMC"

# 当前批次时间戳，用于区分多次验收演示运行。
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

# 当前批次总目录。
batch_dir <- file.path(RUN_ROOT, paste0("run_", timestamp))

# 每一步控制台输出日志保存目录。
log_dir <- file.path(batch_dir, "logs")

# 每一步结果文件保存目录。
out_dir <- file.path(batch_dir, "outputs")

# 生成的临时演示脚本保存目录。
# 这些脚本由已有 module3 脚本自动改写路径得到，用于避免覆盖正式实验结果。
generated_script_dir <- file.path(batch_dir, "generated_scripts")

# 现场重跑核心算法的输出目录。
# 为了控制时间，只重跑 MONO + atac_w=0.7 这一条核心链路。
ALG_OUT_DIR <- file.path(out_dir, "algorithm_core_mono")

# Pressure-test plots are display materials for acceptance. Keep their generated
# files inside the current batch so the README paths always match the run.
ROBUST_OUT_DIR <- file.path(out_dir, "06_robustness_figures")

dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(generated_script_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ALG_OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(ROBUST_OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# 把 Windows 路径转成 R 脚本中更稳妥的正斜杠路径。
r_path <- function(path) {
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

# =========================
# 2. 工具函数
# =========================

# 复制单个文件到指定目录。
# 如果文件不存在，不中断流程，只给出 warning，便于验收时定位缺失文件。
copy_if_exists <- function(from, to_dir) {
  dir.create(to_dir, recursive = TRUE, showWarnings = FALSE)
  if (file.exists(from)) {
    file.copy(from, file.path(to_dir, basename(from)), overwrite = TRUE)
    return(TRUE)
  }
  warning("缺少结果文件：", from)
  FALSE
}

# 批量复制多个文件。
copy_many <- function(files, to_dir) {
  ok <- vapply(files, copy_if_exists, logical(1), to_dir = to_dir)
  invisible(ok)
}

# 判断脚本路径是否已经是绝对路径。
is_absolute_path <- function(path) {
  grepl("^[A-Za-z]:[/\\\\]", path) || startsWith(path, "/")
}

# 根据已有脚本生成一个“演示专用脚本”。
# 这样可以最大限度复用原有实验代码，同时把输出路径改到本次批次目录，
# 避免覆盖正式结果。
make_demo_script <- function(src_rel, new_name, patch_fun) {
  src <- file.path(ROOT_DIR, src_rel)
  lines <- readLines(src, encoding = "UTF-8", warn = FALSE)
  lines <- patch_fun(lines)
  out <- file.path(generated_script_dir, new_name)
  writeLines(lines, out, useBytes = TRUE)
  out
}

# 运行一个演示步骤：
#   1. 计时；
#   2. 保存控制台日志；
#   3. 复制关键结果；
#   4. 返回结构化运行记录，最后汇总成 elapsed_time_summary.csv。
run_step <- function(step_id, label, script, output_files = character()) {
  cat("\n==============================\n")
  cat(step_id, " - ", label, "\n", sep = "")
  cat("script: ", script, "\n", sep = "")
  cat("==============================\n")

  script_path <- if (is_absolute_path(script)) script else file.path(ROOT_DIR, script)
  log_file <- file.path(log_dir, paste0(step_id, "_console.log"))
  log_con <- file(log_file, open = "wt", encoding = "UTF-8")

  sink(log_con, split = TRUE)
  sink(log_con, type = "message")

  start_time <- Sys.time()
  elapsed <- NA_real_
  status <- "OK"
  error_message <- ""

# chdir = TRUE 表示 source() 某个脚本时，临时切换到该脚本所在目录。
# module3 的部分旧脚本会用 "../_bootstrap.R" 这类相对路径加载公共配置，
# 如果不切换工作目录，就会错误地去 D:/RStudio/work/../_bootstrap.R 查找。
  result_time <- tryCatch(
    {
      system.time(source(script_path, encoding = "UTF-8", chdir = TRUE))
    },
    error = function(e) {
      status <<- "ERROR"
      error_message <<- conditionMessage(e)
      NULL
    }
  )

  if (!is.null(result_time)) elapsed <- unname(result_time[["elapsed"]])
  end_time <- Sys.time()

  sink(type = "message")
  sink()
  close(log_con)

  step_out_dir <- file.path(out_dir, step_id)
  copy_many(output_files, step_out_dir)

  data.frame(
    step_id = step_id,
    label = label,
    script = script,
    status = status,
    elapsed_seconds = elapsed,
    start_time = format(start_time, "%Y-%m-%d %H:%M:%S"),
    end_time = format(end_time, "%Y-%m-%d %H:%M:%S"),
    error_message = error_message,
    stringsAsFactors = FALSE
  )
}

# =========================
# 3. 准备核心算法演示输入
# =========================

cat("本次验收演示批次目录：\n", batch_dir, "\n", sep = "")
cat("核心算法现场重跑输出目录：\n", ALG_OUT_DIR, "\n", sep = "")

# 复用已有 ATAC LSI z-score block 和 RNA PCA block。
# 这相当于验收现场跳过耗时较长的数据下载、scDRS 官方 CLI 和 RNA PCA 预处理，
# 但后面的 seed、KNN、RWR/TRS 会现场重新计算。
copy_if_exists(
  file.path(BMMC_DIR, "04_dual_modality", "atac_z_block_MultiOmicsFusion.rds"),
  ALG_OUT_DIR
)
copy_if_exists(
  file.path(BMMC_DIR, "04_dual_modality", "rna_z_block_MultiOmicsFusion_MONO.rds"),
  ALG_OUT_DIR
)

# 生成演示专用 baseline 检查脚本：
# 原始 04_check_baseline.R 使用调用栈推断 "../_bootstrap.R"，
# 从总入口脚本间接 source 时会误判路径。这里仅把 bootstrap 改成绝对路径，
# 其余 baseline 检查逻辑保持不变。
baseline_script <- make_demo_script(
  "module3_multiome/03_scavenge_baseline/04_check_baseline.R",
  "04_check_baseline_demo.R",
  function(lines) {
    lines <- gsub(
      'source(file.path(dirname(sys.frame(1)$ofile %||% "."), "..", "_bootstrap.R"))',
      'source("D:/RStudio/work/module3_multiome/_bootstrap.R")',
      lines,
      fixed = TRUE
    )
    lines
  }
)

# 生成演示专用 seed 脚本：
# 只改输出路径 OUT_DIR，让 seed 结果写入本次批次目录。
seed_script <- make_demo_script(
  "module3_multiome/04_dual_modality/04_seed_scdrs.R",
  "04_seed_scdrs_demo.R",
  function(lines) {
    lines <- gsub(
      'OUT_DIR   <- "D:/SCAVENGE-reproducibility-main/BMMC/04_dual_modality"',
      paste0('OUT_DIR   <- "', r_path(ALG_OUT_DIR), '"'),
      lines,
      fixed = TRUE
    )
    lines
  }
)

# 生成演示专用 KNN 脚本：
# 只跑 MONO 和 atac_w=0.7，现场重新构建带权 mutual KNN 图。
knn_script <- make_demo_script(
  "module3_multiome/04_dual_modality/03_build_weighted_knn.R",
  "03_build_weighted_knn_demo.R",
  function(lines) {
    lines <- gsub(
      'OUT_DIR <- "D:/SCAVENGE-reproducibility-main/BMMC/04_dual_modality"',
      paste0('OUT_DIR <- "', r_path(ALG_OUT_DIR), '"'),
      lines,
      fixed = TRUE
    )
    lines[grepl("^ATAC_W_LIST <-", lines)] <- "ATAC_W_LIST <- c(0.7)"
    lines[grepl("^TRAITS  <-", lines)] <- 'TRAITS  <- c("MONO")'
    lines
  }
)

# 生成演示专用 RWR/TRS 脚本：
# 只跑 MONO 和 atac_w=0.7，现场重新做带重启随机游走和 TRS 后处理。
rwr_script <- make_demo_script(
  "module3_multiome/04_dual_modality/05_rwr_trs.R",
  "05_rwr_trs_demo.R",
  function(lines) {
    lines <- gsub(
      'OUT_DIR <- "D:/SCAVENGE-reproducibility-main/BMMC/04_dual_modality"',
      paste0('OUT_DIR <- "', r_path(ALG_OUT_DIR), '"'),
      lines,
      fixed = TRUE
    )
    lines[grepl("^ATAC_W_LIST <-", lines)] <- "ATAC_W_LIST <- c(0.7)"
    lines[grepl("^TRAITS  <-", lines)] <- 'TRAITS  <- c("MONO")'
    lines
  }
)

# 生成演示专用 AUROC 脚本：
# 评估本次批次中刚刚生成的 TRS_MultiOmicsFusion_MONO_w0.7.tsv.gz。
auroc_script <- make_demo_script(
  "module3_multiome/04_dual_modality/06_evaluate_auroc.R",
  "06_evaluate_auroc_demo.R",
  function(lines) {
    lines <- gsub(
      'OUT_DIR      <- "D:/SCAVENGE-reproducibility-main/BMMC/04_dual_modality"',
      paste0('OUT_DIR      <- "', r_path(ALG_OUT_DIR), '"'),
      lines,
      fixed = TRUE
    )
    lines
  }
)

# Generate a batch-local robustness plotting script. The source script reads
# atac_sparsify_eval.csv from OUT_DIR and writes summary/figures under OUT_DIR,
# so copy the existing long-running robustness table into this acceptance batch.
copy_if_exists(
  file.path(BMMC_DIR, "05_robustness", "atac_sparsify_eval.csv"),
  ROBUST_OUT_DIR
)

robustness_script <- make_demo_script(
  "module3_multiome/05_robustness/04_plot_robustness.R",
  "07_plot_robustness_demo.R",
  function(lines) {
    lines <- gsub(
      'OUT_DIR  <- "D:/SCAVENGE-reproducibility-main/BMMC/05_robustness"',
      paste0('OUT_DIR  <- "', r_path(ROBUST_OUT_DIR), '"'),
      lines,
      fixed = TRUE
    )
    lines
  }
)

# =========================
# 4. 正式运行验收演示步骤
# =========================

steps <- list()

# 第 1 步：论文BMMC ATAC-only 基线实验，展示 ATAC-only SCAVENGE 基线。
# 这一步主要用于交代 baseline，不是核心改进算法。
steps[[length(steps) + 1]] <- run_step(
  step_id = "01_baseline_context",
  label = "ATAC-only SCAVENGE 基线检查",
  script = baseline_script,
  output_files = c(
    file.path(BMMC_DIR, "03_scavenge_baseline", "baseline_enrichment.csv"),
    file.path(BMMC_DIR, "03_scavenge_baseline", "figures", "baseline_top5_celltype.png")
  )
)

# 第 2 步：核心算法步骤 1。
# 从 scDRS norm_score 中选择 top 1% 作为 RWR 的种子细胞。
steps[[length(steps) + 1]] <- run_step(
  step_id = "02_seed_selection",
  label = "核心算法：scDRS 驱动的种子细胞选择",
  script = seed_script,
  output_files = c(
    file.path(ALG_OUT_DIR, "seed_summary_MultiOmicsFusion.csv")
  )
)

# 第 3 步：核心算法步骤 2。
# 在 ATAC/RNA 融合表示上现场重新构建带权 mutual KNN 图。
steps[[length(steps) + 1]] <- run_step(
  step_id = "03_weighted_knn",
  label = "核心算法：带权互最近邻图构建",
  script = knn_script,
  output_files = c(
    file.path(ALG_OUT_DIR, "knn_summary_MultiOmicsFusion.csv")
  )
)

# 第 4 步：核心算法步骤 3。
# 使用上一步构建的图和第 2 步的 seed，现场重新运行 RWR 并输出 TRS。
steps[[length(steps) + 1]] <- run_step(
  step_id = "04_rwr_trs",
  label = "核心算法：带重启随机游走与 TRS 计算",
  script = rwr_script,
  output_files = c(
    file.path(ALG_OUT_DIR, "trs_summary_MultiOmicsFusion.csv"),
    file.path(ALG_OUT_DIR, "TRS_MultiOmicsFusion_MONO_w0.7.tsv.gz")
  )
)

# 第 5 步：评估本次现场生成的 TRS。
# 这一步不是算法主体，但用于把第 4 步输出转成容易理解的 AUROC 数字。
steps[[length(steps) + 1]] <- run_step(
  step_id = "05_auroc_eval",
  label = "现场生成 TRS 的 AUROC 评估",
  script = auroc_script,
  output_files = c(
    file.path(ALG_OUT_DIR, "eval_MultiOmicsFusion_auroc.csv")
  )
)

# 第 6 步：论文ATAC 退化压力测试结果展示。
# 这一步只做已有压力测试结果的汇总和画图，放在核心算法之后作为展示材料。
steps[[length(steps) + 1]] <- run_step(
  step_id = "06_robustness_figures",
  label = "结果展示：ATAC 退化压力测试曲线",
  script = robustness_script
)

# =========================
# 5. 保存本次批次汇总
# =========================

summary_dt <- do.call(rbind, steps)

# 保存 R 运行环境信息，便于验收时说明复现实验环境。
writeLines(capture.output(sessionInfo()), file.path(batch_dir, "session_info.txt"))

# 给老师展示时，优先打开 README 中列出的文件。
readme_lines <- c(
  "毕业设计软件验收演示运行结果",
  paste0("本次批次目录：", batch_dir),
  paste0("运行时间戳：", timestamp),
  "",
  "本次现场真正重跑的核心算法链路：",
  "1. scDRS norm_score top 1% 种子细胞选择",
  "2. MONO + atac_w=0.7 的带权 mutual KNN 图构建",
  "3. 基于 seed 和 KNN 图的 RWR/TRS 计算",
  "4. 使用现场生成的 TRS 做 AUROC 评估",
  "",
  "建议验收时优先打开：",
  "1. logs/02_seed_selection_console.log",
  "2. logs/03_weighted_knn_console.log",
  "3. logs/04_rwr_trs_console.log",
  "4. outputs/algorithm_core_mono/seed_summary_MultiOmicsFusion.csv",
  "5. outputs/algorithm_core_mono/knn_summary_MultiOmicsFusion.csv",
  "6. outputs/algorithm_core_mono/trs_summary_MultiOmicsFusion.csv",
  "7. outputs/algorithm_core_mono/eval_MultiOmicsFusion_auroc.csv",
  "8. outputs/06_robustness_figures/*.png",
  "",
  "运行状态和耗时不再单独保存为 CSV，直接在 RStudio 控制台输出。"
)
writeLines(readme_lines, file.path(batch_dir, "README.txt"), useBytes = TRUE)

cat("\n\n====== 验收演示运行状态和耗时 ======\n")
for (i in seq_len(nrow(summary_dt))) {
  elapsed_txt <- if (is.na(summary_dt$elapsed_seconds[i])) {
    "NA"
  } else {
    sprintf("%.2f 秒", summary_dt$elapsed_seconds[i])
  }

  line <- sprintf(
    "%s | %s | %s",
    summary_dt$step_id[i],
    summary_dt$status[i],
    elapsed_txt
  )

  if (nzchar(summary_dt$error_message[i])) {
    line <- paste0(line, " | 错误：", summary_dt$error_message[i])
  }

  cat(line, "\n", sep = "")
}

core_steps <- c("02_seed_selection", "03_weighted_knn", "04_rwr_trs", "05_auroc_eval")
core_elapsed <- summary_dt$elapsed_seconds[summary_dt$step_id %in% core_steps]
core_elapsed <- core_elapsed[!is.na(core_elapsed)]
if (length(core_elapsed) > 0) {
  cat(sprintf("\n核心算法链路总耗时：%.2f 秒\n", sum(core_elapsed)))
}

cat("\n本次批次结果文件夹：\n", batch_dir, "\n", sep = "")
