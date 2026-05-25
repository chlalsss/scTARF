# Random walk with restart (RWR) supporting per-cell heterogeneous restart
# probability gamma_i.
#
# Iteration:
#   v_{s+1}[i] = (1 - gamma_i) * sum_j(T[i,j] * v_s[j]) + gamma_i * v_0[i]
# where T = column-normalized adjacency.
#
# When gamma is a scalar, this reduces exactly to SCAVENGE's randomWalk_sparse
# behavior (verified by C0/C1 sanity in 02_run_scavenge.R).
#
# Convergence is guaranteed for any gamma_i in (0, 1) by personalized
# PageRank theory (Avrachenkov, Litvak, Nemirovsky, Osipova 2007).

# Build column-stochastic transition matrix from a binary mutual-KNN graph.
# Drops dangling nodes (col sum 0); caller should ensure no isolates upstream.
column_normalize_adj <- function(intM) {
  intM <- methods::as(intM, "dgCMatrix")
  col_sums <- Matrix::colSums(intM)
  if (any(col_sums == 0)) {
    stop("intM has dangling nodes (col_sum == 0); remove isolated cells first.",
         call. = FALSE)
  }
  D_inv <- Matrix::Diagonal(n = ncol(intM), x = 1 / col_sums)
  intM %*% D_inv
}

# Main RWR. Returns a named numeric vector (stationary distribution).
random_walk_hetgamma <- function(intM,
                                 queryCells,
                                 gamma,
                                 eps = 1e-5,
                                 max_iter = 200,
                                 verbose = FALSE) {
  cells <- rownames(intM)
  if (is.null(cells)) stop("intM must have rownames", call. = FALSE)
  if (!identical(cells, colnames(intM))) {
    stop("intM rownames and colnames must match", call. = FALSE)
  }
  n <- length(cells)

  if (length(gamma) == 1) {
    gamma_vec <- rep(as.numeric(gamma), n)
  } else {
    if (is.null(names(gamma))) {
      stop("gamma must be a named vector when length > 1", call. = FALSE)
    }
    missing_cells <- setdiff(cells, names(gamma))
    if (length(missing_cells) > 0) {
      stop("gamma is missing values for ", length(missing_cells), " cells; e.g. ",
           paste(utils::head(missing_cells, 3), collapse = ", "), call. = FALSE)
    }
    gamma_vec <- as.numeric(gamma[cells])
  }
  if (any(!is.finite(gamma_vec)) || any(gamma_vec <= 0) || any(gamma_vec >= 1)) {
    stop("gamma_i must lie strictly in (0, 1)", call. = FALSE)
  }

  # Transition matrix
  Tmat <- column_normalize_adj(intM)

  # Initial / restart distribution v0
  seed_in_graph <- intersect(queryCells, cells)
  if (length(seed_in_graph) == 0) {
    stop("No queryCells present in intM rownames", call. = FALSE)
  }
  v0 <- numeric(n)
  names(v0) <- cells
  v0[seed_in_graph] <- 1 / length(seed_in_graph)

  v <- v0
  delta <- Inf
  iter_used <- 0L
  for (iter in seq_len(max_iter)) {
    Tv <- as.numeric(Tmat %*% v)
    v_new <- (1 - gamma_vec) * Tv + gamma_vec * v0
    delta <- max(abs(v_new - v))
    v <- v_new
    iter_used <- iter
    if (verbose && (iter <= 5 || iter %% 10 == 0)) {
      cat(sprintf("    iter %3d  delta = %.2e\n", iter, delta))
    }
    if (delta < eps) break
  }
  if (delta >= eps) {
    warning("RWR did not converge in ", max_iter,
            " iterations; final delta = ", signif(delta, 3))
  } else if (verbose) {
    cat(sprintf("    converged at iter %d, delta = %.2e\n", iter_used, delta))
  }

  names(v) <- cells
  attr(v, "n_iter") <- iter_used
  attr(v, "final_delta") <- delta
  attr(v, "n_seed_in_graph") <- length(seed_in_graph)
  v
}

# Helper: cap the top tail at a quantile, then min-max scale to [0, 1].
# Mirrors SCAVENGE's capOutlierQuantile() %>% max_min_scale().
cap_and_scale <- function(x, q_ceiling = 0.95) {
  cap <- as.numeric(stats::quantile(x, probs = q_ceiling, na.rm = TRUE))
  y <- pmin(x, cap)
  rng <- range(y, na.rm = TRUE)
  if (diff(rng) <= 0) return(rep(0, length(y)))
  (y - rng[1]) / (rng[2] - rng[1])
}

# Helper: scale factor = mean of top scale_pct of z_score (matches SCAVENGE).
compute_scale_factor <- function(z_score, scale_pct = 0.01) {
  z <- z_score[is.finite(z_score)]
  cap_n <- max(1L, ceiling(length(z) * scale_pct))
  z_top <- sort(z, decreasing = TRUE)[seq_len(cap_n)]
  mean(z_top, na.rm = TRUE)
}
