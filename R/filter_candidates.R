#' Admit specific, thermo-valid primers into the optimiser
#'
#' Keeps scored primers whose enrichment clears `min_enrichment` and whose
#' sequence passes the thermodynamic gate (spec: specificity is an entry
#' threshold, not part of the coverage objective).
#'
#' @param scored data.frame with `pattern` and `enrichment` columns.
#' @param min_enrichment Minimum enrichment fold to admit.
#' @param tm_min,tm_max,gc_min,gc_max,gc_clamp_len Passed to [passes_thermo()].
#' @return The admitted subset of `scored`.
#' @export
filter_candidates <- function(scored, min_enrichment = 2,
                              tm_min = 15, tm_max = 45,
                              gc_min = 0.2, gc_max = 0.8, gc_clamp_len = 2) {
  ok <- vapply(scored$pattern, passes_thermo, logical(1),
               tm_min = tm_min, tm_max = tm_max,
               gc_min = gc_min, gc_max = gc_max, gc_clamp_len = gc_clamp_len)
  keep <- !is.na(scored$enrichment) & scored$enrichment >= min_enrichment & ok
  out <- scored[keep, , drop = FALSE]
  rownames(out) <- NULL
  out
}
