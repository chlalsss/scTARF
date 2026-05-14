source("D:/RStudio/work/module2_wnn_lite_gated/_bootstrap.R")

gated_version <- "celltype_gated"
gated_label <- "celltype_gated"
k <- env_integer("WNN_LITE_KNN_K", env_integer("MODULE2_KNN_K", 30))
gated_knn_file <- file.path(
  gated_paths(gated_version)$out_dir,
  sprintf("wnn_lite_celltype_gated_mutualknn_k%d.rds", k)
)

source(file.path(wnn_gated_code_dir, "_diagnostics_gated_common.R"))
