# install_deps.R
# 仅装包, 不做业务. 在 RStudio 里**单独一个 session**跑, 跑完关掉这个 session 再开新 session 跑 prepare_h5ad.R.
#
# 装 zellkonverter (R <-> AnnData 桥) + Matrix (rds 里的 sparse 矩阵会用到).
# zellkonverter 第一次 library() 时会用 basilisk 拉一个内嵌 conda env, 那一步在 prepare_h5ad.R 里发生.

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("zellkonverter", update = FALSE, ask = FALSE)

# 检查
cat("zellkonverter version: ",
    as.character(packageVersion("zellkonverter")), "\n")
cat("Matrix version:        ",
    as.character(packageVersion("Matrix")), "\n")
cat("install_deps.R OK. 关掉本 session, 新开 session 跑 prepare_h5ad.R.\n")
