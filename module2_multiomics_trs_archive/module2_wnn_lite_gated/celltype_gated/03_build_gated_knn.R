source("D:/RStudio/work/module2_wnn_lite_gated/_bootstrap.R")

gated_version <- "celltype_gated"
gated_label <- "celltype_gated"
gated_weights_file <- file.path(
  gated_paths(gated_version)$out_dir,
  "wnn_lite_celltype_gated_weights.csv"
)

source(file.path(wnn_gated_code_dir, "_build_gated_knn_common.R"))
