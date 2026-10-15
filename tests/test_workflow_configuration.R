# Check configuration before package loading or dataset I/O.
read_config <- function(path, overrides = list()) {
  e <- list2env(overrides, parent = baseenv())
  for (expr in parse(path, encoding = "UTF-8")) {
    if (is.call(expr) && identical(expr[[1]], as.name("suppressPackageStartupMessages"))) break
    eval(expr, e)
  }
  e
}
drop_path <- "scripts/02_sparsification/01_paired_atac_sparsify.R"
default <- read_config(drop_path)
stopifnot(identical(default$P_DROP_LIST, seq(0, 0.9, by = 0.1)))
custom <- read_config(drop_path, list(P_DROP_LIST = 0, SEED_LIST = 7L, ATAC_W_LIST = 1))
stopifnot(identical(custom$P_DROP_LIST, 0), identical(custom$SEED_LIST, 7L),
          identical(custom$ATAC_W_LIST, 1))
for (p in c("15_build_weighted_knn.R", "17_rwr_trs.R")) {
  e <- read_config(file.path("scripts/01_paired_analysis", p), list(ATAC_W_LIST = c(0.5, 0.6, 0.7, 0.8, 0.9, 1)))
  stopifnot(identical(e$ATAC_W_LIST, c(0.5, 0.6, 0.7, 0.8, 0.9, 1)))
}
cat("PASS: complete dropout defaults and caller-supplied grids preserved\n")
