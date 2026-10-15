# Synthetic input -> existing weighted graph -> existing RWR -> existing TRS.
# This does not calculate GWAS/scDRS scores or reproduce a manuscript AUROC.

load_sctarf_core <- function(repo_root) {
  files <- file.path(repo_root, "module3_multiome", "04_dual_modality",
                     c("03_build_weighted_knn.R", "05_rwr_trs.R"))
  if (!all(file.exists(files))) stop("Missing original paired core source scripts")
  for (pkg in c("Matrix", "BiocNeighbors")) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop("Install required package: ", pkg)
  }
  suppressPackageStartupMessages({library(Matrix); library(BiocNeighbors)})
  extract <- function(path, keep) {
    env <- new.env(parent = environment(load_sctarf_core))
    found <- character()
    for (expression in parse(path, encoding = "UTF-8")) {
      if (is.call(expression) && as.character(expression[[1]])[1] %in% c("<-", "=")) {
        name <- if (is.symbol(expression[[2]])) as.character(expression[[2]]) else ""
        if (name %in% keep) {
          eval(expression, envir = env)
          found <- c(found, name)
        }
      }
    }
    if (!all(keep %in% found)) stop("Original script no longer has the expected core definitions: ", path)
    env
  }
  # Separate environments preserve the graph's MAX_ITER=1 and RWR's MAX_ITER=200.
  graph <- extract(files[1], c("K_INIT", "K_MAX", "ISO_TOL", "MAX_ITER", "build_knn"))
  propagation <- extract(files[2], c("GAMMA", "TOL", "MAX_ITER", "CAP_Q", "SCALE_FAC", "rwr", "post_process"))
  list(build_knn = graph$build_knn, rwr = propagation$rwr,
       post_process = propagation$post_process, source_files = files)
}

run_core_demo <- function(repo_root, output_dir) {
  core <- load_sctarf_core(repo_root)
  set.seed(20261015)
  n <- 600L
  cells <- sprintf("synthetic_%04d", seq_len(n))
  # These continuous synthetic signals have no biological or disease meaning.
  shared <- rnorm(n)
  atac <- scale(cbind(shared + rnorm(n), matrix(rnorm(n * 29L), n)))
  rna <- scale(cbind(shared + rnorm(n, sd = 0.4), matrix(rnorm(n * 9L), n)))
  rownames(atac) <- rownames(rna) <- cells
  rna_score <- as.numeric(rna[, 1])
  selected <- order(rna_score, decreasing = TRUE)[seq_len(ceiling(0.01 * n))]
  seeds <- setNames(as.numeric(seq_len(n) %in% selected), cells)
  scores <- data.frame(cell = cells, synthetic_rna_score = rna_score, seed = seeds,
                       row.names = NULL)
  summaries <- list()
  started <- proc.time()[["elapsed"]]
  for (w in c(1.0, 0.7)) {
    set.seed(20261015)
    embedding <- cbind(sqrt(w) * atac, sqrt(1 - w) * rna)
    rownames(embedding) <- cells
    result <- core$build_knn(embedding)
    graph <- result$adj
    stopifnot(identical(rownames(graph), cells), identical(colnames(graph), cells))
    rowsums <- as.numeric(Matrix::rowSums(graph))
    stopifnot(all(abs(rowsums[rowsums > 0] - 1) < 1e-10))
    raw <- core$rwr(graph, seeds)
    trs <- core$post_process(raw)
    stopifnot(length(trs) == n, all(is.finite(trs)), min(trs) >= 0, max(trs) <= 1)
    column <- if (w == 1) "trs_w1" else "trs_w0.7"
    scores[[column]] <- trs
    summaries[[column]] <- data.frame(
      input = "synthetic standardized representations and synthetic RNA scores",
      n_cells = n, n_seeds = sum(seeds), atac_weight = w, k_final = result$k_final,
      isolated_fraction = result$iso_frac, trs_min = min(trs), trs_max = max(trs),
      graph_object_bytes = as.numeric(object.size(graph)),
      graph_source_md5 = unname(tools::md5sum(core$source_files[1])),
      propagation_source_md5 = unname(tools::md5sum(core$source_files[2])))
  }
  summary <- do.call(rbind, summaries)
  summary$elapsed_seconds_total <- proc.time()[["elapsed"]] - started
  rownames(summary) <- NULL
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(scores, file.path(output_dir, "scores.csv"), row.names = FALSE)
  write.csv(summary, file.path(output_dir, "run_summary.csv"), row.names = FALSE)
  writeLines(c("Synthetic core-computation example; not a manuscript reproduction.",
               capture.output(sessionInfo())), file.path(output_dir, "session_info.txt"))
  cat("Wrote synthetic core example to ", normalizePath(output_dir), "\n", sep = "")
  print(summary[c("n_cells", "n_seeds", "atac_weight", "k_final", "isolated_fraction", "trs_min", "trs_max")], row.names = FALSE)
  invisible(list(scores = scores, summary = summary))
}

if (sys.nframe() == 0L) {
  script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
  repo <- normalizePath(file.path(dirname(script), ".."), mustWork = TRUE)
  args <- commandArgs(trailingOnly = TRUE)
  out <- if (length(args)) args[1] else file.path(repo, "examples", "output", "core_demo")
  run_core_demo(repo, out)
}
