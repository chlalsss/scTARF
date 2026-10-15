# Exercise the real Python LSI adapter through the R launcher in a spaced path.
source("scripts/_paths.R")
scratch <- tempfile("scTARF path with spaces ")
dir.create(scratch)
python <- Sys.which("python")
stopifnot(nzchar(python))
npy <- file.path(scratch, "input matrix.npy")
cells <- file.path(scratch, "cell names.csv")
write.csv(data.frame(cell = c("cell_A", "cell_B", "cell_C")), cells, row.names = FALSE)
code <- sprintf("import numpy as np; np.save(%s, np.arange(6).reshape(3,2))", encodeString(npy, quote = '"'))
stopifnot(system2(python, c("-c", shQuote(code))) == 0L)
e <- new.env(parent = globalenv())
for (expr in parse("scripts/03_nonpaired_analysis/44_prepare_module3_v2_inputs.R")) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      is.call(expr[[3]]) && identical(expr[[3]][[1]], as.name("function"))) eval(expr, e)
}
e$PYTHON_EXE <- sctarf_python()
e$PY_HELPER <- normalizePath("scripts/03_nonpaired_analysis/42_export_lsi_npy.py", winslash = "/")
e$ATAC_LSI_NPY <- npy
e$ATAC_CELLS_CSV <- cells
e$ATAC_LSI_CSV <- file.path(scratch, "output matrix.csv")
e$FORCE_EXPORT_LSI <- TRUE
e$export_lsi_if_needed()
result <- read.csv(e$ATAC_LSI_CSV)
stopifnot(identical(result$barcode, c("cell_A", "cell_B", "cell_C")),
          identical(result$LSI_1, c(0L, 2L, 4L)), identical(result$LSI_2, c(1L, 3L, 5L)))
cat("PASS: PATH executable and spaced paths preserve LSI values/barcodes\n")
