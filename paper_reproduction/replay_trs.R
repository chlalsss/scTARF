# Replay the manuscript's saved graphs/seeds using its original RWR/TRS code.
# Does not rebuild graphs, change seeds, or modify any historical output.
replay_trs <- function(repo_root, graph_file, seed_file, expected_file, tolerance = 1e-12) {
  source(file.path(repo_root, "examples", "run_core_demo.R"), local = TRUE)
  core <- load_sctarf_core(repo_root)
  graph <- readRDS(graph_file)
  seed <- readRDS(seed_file)
  cells <- rownames(graph)
  stopifnot(nrow(graph) == ncol(graph), length(seed) == nrow(graph),
            !is.null(cells), !anyDuplicated(cells), identical(colnames(graph), cells),
            all(is.finite(seed)), sum(seed) > 0)
  # Preserve original positional seed use; named inputs must confirm that order.
  if (!is.null(names(seed))) stopifnot(identical(names(seed), cells))
  trs <- core$post_process(core$rwr(graph, seed))
  expected <- utils::read.delim(expected_file, check.names = FALSE, stringsAsFactors = FALSE)
  stopifnot(all(c("barcode", "TRS") %in% names(expected)), !anyDuplicated(expected$barcode),
            setequal(expected$barcode, cells), all(is.finite(trs)))
  difference <- max(abs(trs - expected$TRS[match(cells, expected$barcode)]))
  list(scores = data.frame(barcode = cells, TRS = trs), n_cells = length(cells),
       max_abs_difference = difference, matched = difference <= tolerance,
       tolerance = tolerance)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 5L) stop("Usage: replay_trs.R <repo> <graph.rds> <seed.rds> <reference.tsv.gz> <new output.tsv>")
  if (file.exists(args[5])) stop("Choose a new output file; refusing to overwrite: ", args[5])
  result <- replay_trs(args[1], args[2], args[3], args[4])
  dir.create(dirname(args[5]), recursive = TRUE, showWarnings = FALSE)
  write.table(result$scores, args[5], sep = "\t", row.names = FALSE, quote = FALSE)
  print(result[c("n_cells", "max_abs_difference", "matched", "tolerance")])
  if (!result$matched) stop("Replayed TRS did not match the historical cell scores")
}
