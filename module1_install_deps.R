# module1_install_deps.R
# 模块 1 相关依赖一次性安装脚本。
#
# ⚠️ 必须在**全新的 R session**里运行（不要先 library Seurat/ggplot2 再跑本脚本），
#    否则会出现 "Package 'ggplot2' cannot be unloaded" 这类升级冲突。
#
# 推荐操作：
#   1. 关闭当前 RStudio 里所有已加载的包（或直接重启 R：Session → Restart R）
#   2. source("D:/RStudio/work/module1_install_deps.R")
#   3. 装完再运行 module1_method_C_liger.R

cat("========== 模块 1 依赖安装 ==========\n")

# 1) 先升级 ggplot2 等核心画图包到一致版本（避免 rliger 安装过程中的隐式升级）
core_pkgs <- c("ggplot2", "patchwork", "cowplot", "ggrepel", "Matrix")
for (p in core_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

# 2) 安装 rliger（LIGER 官方 R 包，CRAN 版本）
if (!requireNamespace("rliger", quietly = TRUE)) {
  cat("[安装] rliger (CRAN) ...\n")
  install.packages("rliger")
} else {
  cat("[已装] rliger", as.character(packageVersion("rliger")), "\n")
}

# 3) 安装 RcppPlanc（rliger 2.x iNMF 加速后端，来自 welch-lab 的 r-universe）
if (!requireNamespace("RcppPlanc", quietly = TRUE)) {
  cat("[安装] RcppPlanc (welch-lab r-universe) ...\n")
  tryCatch(
    install.packages(
      "RcppPlanc",
      repos = c("https://welch-lab.r-universe.dev",
                "https://cloud.r-project.org")
    ),
    error = function(e) {
      warning("RcppPlanc 安装失败，rliger 将回退到慢速 CPU 版本。错误: ",
              conditionMessage(e))
    }
  )
} else {
  cat("[已装] RcppPlanc", as.character(packageVersion("RcppPlanc")), "\n")
}

# 4) 其他工具包
utils_pkgs <- c("RANN", "data.table", "dplyr", "jsonlite")
for (p in utils_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

# 5) 汇报版本
cat("\n========== 依赖版本汇报 ==========\n")
for (p in c("rliger", "RcppPlanc", "Seurat", "Signac", "harmony",
            "RANN", "ggplot2", "Matrix", "data.table", "jsonlite")) {
  v <- tryCatch(as.character(packageVersion(p)), error = function(e) "NOT INSTALLED")
  cat(sprintf("  %-12s : %s\n", p, v))
}

cat("\n依赖安装完成。现在建议：\n")
cat("  1) Session → Restart R（清空 namespace）\n")
cat("  2) source(\"D:/RStudio/work/module1_method_C_liger.R\")\n")
