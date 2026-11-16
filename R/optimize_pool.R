#' Maximum-coverage primer pool of exactly K primers
#'
#' Exact ILP (re-solved per K): choose K primers maximising covered bins, subject
#' to a per-target coverage floor.
#'
#' @param cov A [coverage_matrix()] result.
#' @param k Pool size.
#' @param floor Minimum fraction of bins each target must cover.
#' @return List with `patterns`, `k`, `coverage`, `per_target`, `regularity`, `status`.
#' @export
optimize_pool <- function(cov, k, floor = 0.5) {
  m <- cov$matrix
  nb <- nrow(m); np <- ncol(m); bt <- cov$bin_target
  if (np > 250L)
    warning(sprintf(paste("optimize_pool(): %d primers to optimise; the max-coverage ILP",
                          "grows very slow beyond ~250 (minutes to hours). Keep only the top",
                          "candidates by enrichment before optimising."), np), call. = FALSE)
  obj <- c(rep(0, np), rep(1, nb))

  rows <- vector("list", nb + 1L + length(unique(bt)))
  dir <- character(length(rows)); rhs <- numeric(length(rows))
  for (i in seq_len(nb)) {
    r <- numeric(np + nb); r[np + i] <- 1; r[which(m[i, ] == 1L)] <- -1
    rows[[i]] <- r; dir[i] <- "<="; rhs[i] <- 0
  }
  rows[[nb + 1L]] <- c(rep(1, np), rep(0, nb)); dir[nb + 1L] <- "=="; rhs[nb + 1L] <- k
  for (j in seq_along(unique(bt))) {
    t <- unique(bt)[j]; idx <- which(bt == t)
    r <- numeric(np + nb); r[np + idx] <- 1
    rows[[nb + 1L + j]] <- r; dir[nb + 1L + j] <- ">="; rhs[nb + 1L + j] <- ceiling(floor * length(idx))
  }

  sol <- lpSolve::lp("max", obj, do.call(rbind, rows), dir, rhs, all.bin = TRUE)
  if (sol$status != 0) {
    return(list(patterns = NA_character_, k = k, coverage = NA_real_,
                per_target = NA_real_, regularity = NA_real_, status = sol$status))
  }
  chosen <- cov$patterns[sol$solution[seq_len(np)] > 0.5]
  covered <- sol$solution[np + seq_len(nb)] > 0.5
  depth <- rowSums(m[, chosen, drop = FALSE])
  cd <- depth[depth > 0]
  list(
    patterns = chosen, k = k,
    coverage = mean(covered),
    per_target = stats::setNames(
      vapply(sort(unique(bt)), function(t) mean(covered[bt == t]), numeric(1)),
      cov$targets[sort(unique(bt))]
    ),
    regularity = if (length(cd) <= 1) 1 else max(0, 1 - stats::sd(cd) / mean(cd)),
    status = sol$status
  )
}
