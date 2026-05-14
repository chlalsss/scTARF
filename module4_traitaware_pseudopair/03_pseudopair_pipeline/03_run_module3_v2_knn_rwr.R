# module4 / Phase 1.5.3: run module3 v2 KNN + RWR/TRS on module4 inputs
# -----------------------------------------------------------------------------
# 目的:
#   最大限度复用 module3 的 v2 下游代码。
#
#   本脚本不会重新实现 KNN / RWR / TRS 算法。它会读取 module3 原始脚本:
#     - module3_multiome/04_dual_modality/03_build_weighted_knn.R
#     - module3_multiome/04_dual_modality/05_rwr_trs.R
#
#   然后只替换三类“运行配置”:
#     1. OUT_DIR       -> module4 的 module3_v2_compat 输出目录
#     2. ATAC_W_LIST   -> 默认 c(0.7, 1.0)
#     3. TRAITS        -> 默认 c("MONO", "PLT", "LYM")
#
# 输入目录:
#   03_pseudopair_pipeline/output/module3_v2_compat/
#     atac_z_block_v2.rds
#     rna_z_block_v2_<trait>.rds
#     seed_v2_<trait>.rds
#
# 输出同一目录:
#     knn_v2_<trait>_w0.7.rds
#     knn_v2_<trait>_w1.0.rds
#     TRS_v2_<trait>_w0.7.tsv.gz
#     TRS_v2_<trait>_w1.0.tsv.gz
#     knn_summary_v2.csv
#     trs_summary_v2.csv
#
# 运行:
#   setwd("D:/RStudio/work/module4_traitaware_pseudopair")
#   source("03_pseudopair_pipeline/03_run_module3_v2_knn_rwr.R")
#
# 可选:
#   只跑 MONO:
#     TRAITS_TO_RUN <- c("MONO")
#     source("03_pseudopair_pipeline/03_run_module3_v2_knn_rwr.R")
#
#   只跑主融合权重 0.7:
#     ATAC_W_LIST <- c(0.7)
#     source("03_pseudopair_pipeline/03_run_module3_v2_knn_rwr.R")
# -----------------------------------------------------------------------------

ROOT <- "D:/RStudio/work/module4_traitaware_pseudopair"
MODULE3_DIR <- "D:/RStudio/work/module3_multiome/04_dual_modality"
OUT_DIR <- file.path(ROOT, "03_pseudopair_pipeline/output/module3_v2_compat")
PATCHED_DIR <- file.path(OUT_DIR, "_patched_module3_scripts")
dir.create(PATCHED_DIR, recursive = TRUE, showWarnings = FALSE)

if (!exists("TRAITS_TO_RUN")) TRAITS_TO_RUN <- c("MONO", "PLT", "LYM")
if (!exists("ATAC_W_LIST")) ATAC_W_LIST <- c(0.7, 1.0)

MODULE3_KNN_SCRIPT <- file.path(MODULE3_DIR, "03_build_weighted_knn.R")
MODULE3_RWR_SCRIPT <- file.path(MODULE3_DIR, "05_rwr_trs.R")

stop_missing <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0L) {
    stop("Missing required file(s):\n  ", paste(missing, collapse = "\n  "),
         call. = FALSE)
  }
}

r_literal_chr <- function(x) {
  paste(sprintf('"%s"', x), collapse = ", ")
}

r_literal_num <- function(x) {
  paste(format(x, trim = TRUE, scientific = FALSE), collapse = ", ")
}

patch_module3_script <- function(script_path, patched_path) {
  txt <- readLines(script_path, warn = FALSE)

  txt <- sub(
    "^ATAC_W_LIST\\s*<-.*$",
    paste0("ATAC_W_LIST <- c(", r_literal_num(ATAC_W_LIST), ")"),
    txt
  )
  txt <- sub(
    "^OUT_DIR\\s*<-.*$",
    paste0('OUT_DIR <- "', normalizePath(OUT_DIR, winslash = "/", mustWork = TRUE), '"'),
    txt
  )
  txt <- sub(
    "^TRAITS\\s*<-.*$",
    paste0("TRAITS <- c(", r_literal_chr(TRAITS_TO_RUN), ")"),
    txt
  )

  writeLines(txt, patched_path, useBytes = TRUE)
  patched_path
}

cat(">>> Run module3 v2 KNN + RWR/TRS on module4-compatible inputs\n\n")
cat("OUT_DIR: ", OUT_DIR, "\n", sep = "")
cat("TRAITS: ", paste(TRAITS_TO_RUN, collapse = ", "), "\n", sep = "")
cat("ATAC_W_LIST: ", paste(ATAC_W_LIST, collapse = ", "), "\n\n", sep = "")

required_inputs <- c(
  file.path(OUT_DIR, "atac_z_block_v2.rds"),
  file.path(OUT_DIR, paste0("rna_z_block_v2_", TRAITS_TO_RUN, ".rds")),
  file.path(OUT_DIR, paste0("seed_v2_", TRAITS_TO_RUN, ".rds")),
  MODULE3_KNN_SCRIPT,
  MODULE3_RWR_SCRIPT
)
stop_missing(required_inputs)

patched_knn <- file.path(PATCHED_DIR, "module3_03_build_weighted_knn.module4_patched.R")
patched_rwr <- file.path(PATCHED_DIR, "module3_05_rwr_trs.module4_patched.R")

cat("[1/3] Creating patched module3 runner scripts ...\n")
patch_module3_script(MODULE3_KNN_SCRIPT, patched_knn)
patch_module3_script(MODULE3_RWR_SCRIPT, patched_rwr)
cat("  wrote: ", patched_knn, "\n", sep = "")
cat("  wrote: ", patched_rwr, "\n\n", sep = "")

cat("[2/3] Running module3 KNN script with module4 OUT_DIR ...\n")
source(patched_knn, local = new.env(parent = globalenv()))

cat("\n[3/3] Running module3 RWR/TRS script with module4 OUT_DIR ...\n")
source(patched_rwr, local = new.env(parent = globalenv()))

cat("\n=== module3 v2 KNN + RWR/TRS finished for module4 ===\n")
cat("Expected outputs in: ", OUT_DIR, "\n", sep = "")
cat("  knn_v2_<trait>_w<weight>.rds\n")
cat("  TRS_v2_<trait>_w<weight>.tsv.gz\n")
cat("  knn_summary_v2.csv\n")
cat("  trs_summary_v2.csv\n")
