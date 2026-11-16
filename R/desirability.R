#' Single desirability score (0-100) for a primer pool
#'
#' Weighted geometric mean of the three visible metrics; enrichment is first
#' saturated to 0-1. Geometric so a weak metric sinks the score.
#'
#' @param enrichment Enrichment fold (>= 0).
#' @param coverage Genome coverage fraction (0-1).
#' @param regularity Coverage evenness (0-1).
#' @param enrich_ref Enrichment giving a 0.5 desirability.
#' @param w_enrich Weight on enrichment; coverage and regularity split the rest.
#' @return Numeric score in 0-100.
#' @export
desirability <- function(enrichment, coverage, regularity, enrich_ref = 10, w_enrich = 0.4) {
  d <- c(enrichment / (enrichment + enrich_ref), coverage, regularity)
  w <- c(w_enrich, (1 - w_enrich) / 2, (1 - w_enrich) / 2)
  prod(d ^ w) * 100
}
