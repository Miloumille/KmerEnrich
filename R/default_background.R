#' Parameter-free background model from base frequencies
#'
#' A null model for when no background reads are available: each base occurs at
#' the given frequency (uniform 25\% by default), independently of its neighbours
#' (order-0). Drop-in replacement for a fitted [markov_background()] model.
#'
#' @param freqs Named A/C/G/T base frequencies (default: uniform).
#' @return A model list with `pi` and `trans`, usable by [score_enrichment()].
#' @export
default_background <- function(freqs = c(A = 0.25, C = 0.25, G = 0.25, T = 0.25)) {
  bases <- c("A", "C", "G", "T")
  p <- freqs[bases] / sum(freqs[bases])
  trans <- matrix(p, 4L, 4L, byrow = TRUE, dimnames = list(bases, bases))
  list(pi = p, trans = trans)
}
