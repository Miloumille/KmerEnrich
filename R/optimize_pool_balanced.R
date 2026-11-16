#' Balanced maximum-coverage primer pool of exactly K primers
#'
#' Chooses K primers to maximise the coverage of the **least-covered** target
#' (max-min), breaking ties by total coverage. Unlike [optimize_pool()] there is no
#' hard floor, so a pool always exists for every K — coverage simply drops for small
#' K. This keeps every genome in focus at once instead of favouring the easiest one.
#'
#' @param cov A [coverage_matrix()] result.
#' @param k Pool size.
#' @return List with `patterns`, `k`, `coverage`, `per_target`, `min_cov`,
#'   `regularity`, `status`.
#' @export
optimize_pool_balanced <- function(cov, k) {
  m <- cov$matrix
  nb <- nrow(m); np <- ncol(m); bt <- cov$bin_target
  vs <- sort(unique(bt))
  if (np > 250L)
    warning(sprintf(paste("optimize_pool_balanced(): %d primers; the ILP grows very slow",
                          "beyond ~250. Keep the top candidates by enrichment first."), np),
            call. = FALSE)

  nvar <- np + nb + 1L                        # x (np primers), y (nb bins), z (min fraction)
  zi <- nvar
  obj <- c(rep(0, np), rep(1, nb), nb + 1)    # max (nb+1)*z + sum(y): z dominant, total tie-break

  rows <- vector("list", nb + 1L + length(vs) + 1L)
  dir <- character(length(rows)); rhs <- numeric(length(rows))
  for (b in seq_len(nb)) {                    # y_b <= sum of chosen primers covering b
    r <- numeric(nvar); r[np + b] <- 1; r[which(m[b, ] == 1L)] <- -1
    rows[[b]] <- r; dir[b] <- "<="; rhs[b] <- 0
  }
  r <- numeric(nvar); r[seq_len(np)] <- 1
  rows[[nb + 1L]] <- r; dir[nb + 1L] <- "=="; rhs[nb + 1L] <- k
  for (j in seq_along(vs)) {                  # sum_{b in v} y_b >= z * nbins_v
    idx <- which(bt == vs[j])
    r <- numeric(nvar); r[np + idx] <- 1; r[zi] <- -length(idx)
    rows[[nb + 1L + j]] <- r; dir[nb + 1L + j] <- ">="; rhs[nb + 1L + j] <- 0
  }
  r <- numeric(nvar); r[zi] <- 1
  rows[[length(rows)]] <- r; dir[length(rows)] <- "<="; rhs[length(rows)] <- 1

  sol <- lpSolve::lp("max", obj, do.call(rbind, rows), dir, rhs, binary.vec = seq_len(np + nb))
  chosen <- cov$patterns[sol$solution[seq_len(np)] > 0.5]
  covered <- sol$solution[np + seq_len(nb)] > 0.5
  depth <- rowSums(m[, chosen, drop = FALSE]); cd <- depth[depth > 0]
  per <- stats::setNames(vapply(vs, function(v) mean(covered[bt == v]), numeric(1)),
                         cov$targets[vs])
  list(patterns = chosen, k = k, coverage = mean(covered), per_target = per,
       min_cov = if (length(per)) min(per) else NA_real_,
       regularity = if (length(cd) <= 1) 1 else max(0, 1 - stats::sd(cd) / mean(cd)),
       status = sol$status)
}

#' Balanced saturation sweep: optimal pool for every size N
#'
#' Solves [optimize_pool_balanced()] for K = 1..k_max, returning the coverage curve,
#' the plateau N* and every pool so a report can switch between them.
#'
#' @param cov A [coverage_matrix()] result.
#' @param k_max Largest pool size to evaluate.
#' @param plateau Relative-gain threshold marking the plateau.
#' @return List with `curve` (k, coverage, min_cov), `n_star` and `pools`.
#' @export
saturation_balanced <- function(cov, k_max, plateau = 0.05) {
  k_max <- min(k_max, length(cov$patterns))
  pools <- lapply(seq_len(k_max), function(k) suppressWarnings(optimize_pool_balanced(cov, k)))
  cvg <- vapply(pools, function(p) p$coverage, numeric(1))
  curve <- data.frame(k = seq_len(k_max), coverage = cvg,
                      min_cov = vapply(pools, function(p) p$min_cov, numeric(1)))
  n_star <- k_max
  for (i in seq_len(k_max)) {
    nxt <- if (i < k_max) cvg[i + 1] else NA
    if (is.na(nxt) || cvg[i] <= 0 || (nxt - cvg[i]) / cvg[i] < plateau) { n_star <- i; break }
  }
  list(curve = curve, n_star = n_star, pools = pools)
}
