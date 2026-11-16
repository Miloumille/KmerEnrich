#' Coverage saturation curve and plateau size N*
#'
#' Sweeps K = 1..k_max, re-solving the pool each time, and picks N* where the
#' relative coverage gain from adding one more primer drops below `plateau`.
#'
#' @param cov A [coverage_matrix()] result.
#' @param k_max Largest pool size to evaluate (capped at the primer count).
#' @param floor Per-target coverage floor.
#' @param plateau Relative-gain threshold marking the plateau.
#' @return List with `curve` (data.frame) and `n_star`.
#' @export
saturation_curve <- function(cov, k_max, floor = 0.5, plateau = 0.05) {
  k_max <- min(k_max, length(cov$patterns))
  if (length(cov$patterns) > 250L)
    warning(sprintf(paste("saturation_curve(): %d primers x %d ILP solves; this can be",
                          "extremely slow (the ILP is intractable beyond ~250 primers).",
                          "Pre-filter to the top candidates by enrichment."),
                    length(cov$patterns), k_max), call. = FALSE)
  curve <- do.call(rbind, lapply(seq_len(k_max), function(k) {
    r <- suppressWarnings(optimize_pool(cov, k, floor))
    data.frame(k = k, coverage = r$coverage, status = r$status)
  }))
  curve$marginal_gain <- c(NA, diff(curve$coverage))

  n_star <- NA_integer_
  for (i in seq_len(nrow(curve))) {
    cur <- curve$coverage[i]
    if (is.na(cur)) next
    nxt <- if (i < nrow(curve)) curve$coverage[i + 1] else NA
    if (is.na(nxt) || cur <= 0 || (nxt - cur) / cur < plateau) {
      n_star <- curve$k[i]; break
    }
  }
  list(curve = curve, n_star = n_star)
}
