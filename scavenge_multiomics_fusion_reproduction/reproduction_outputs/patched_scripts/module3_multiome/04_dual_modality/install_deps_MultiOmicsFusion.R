# install_deps_MultiOmicsFusion.R
# 仅装包. 在 RStudio 单独 session 跑, 跑完关掉, 重开 session 跑 02_fusion_autotune.R.
# §2.B 需要 bluster (Bioconductor) 算 70k cells 的 approximate silhouette.

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("bluster", update = FALSE, ask = FALSE)

cat("bluster: ", as.character(packageVersion("bluster")), "\n")
cat("install_deps_MultiOmicsFusion.R OK. 关 session, 新开 session 跑 02_fusion_autotune.R.\n")
