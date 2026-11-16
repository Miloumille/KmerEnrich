#' Concrete bases behind an IUPAC letter
#'
#' @param letter Single IUPAC letter.
#' @return Character vector of the bases it stands for.
#' @export
iupac_bases <- function(letter) {
  c(A = "A", C = "C", G = "G", T = "T",
    R = "AG", Y = "CT", S = "CG", W = "AT", K = "GT", M = "AC",
    B = "CGT", D = "AGT", H = "ACT", V = "ACG", N = "ACGT")[[letter]] |>
    strsplit("", fixed = TRUE) |>
    unlist()
}

#' Expected frequency of a primer under an order-1 Markov model
#'
#' Forward recursion over the pattern's per-position IUPAC base sets, summing the
#' probability of every concrete instantiation.
#'
#' @param pattern Primer pattern (IUPAC allowed).
#' @param model Model from [markov_background()].
#' @return Numeric expected frequency.
#' @export
expected_freq <- function(pattern, model) {
  sets <- lapply(strsplit(pattern, "", fixed = TRUE)[[1]], iupac_bases)
  f <- stats::setNames(model$pi[sets[[1]]], sets[[1]])
  if (length(sets) == 1) return(sum(f))
  for (i in 2:length(sets)) {
    nxt <- stats::setNames(numeric(length(sets[[i]])), sets[[i]])
    for (b in sets[[i]]) nxt[b] <- sum(f * model$trans[names(f), b])
    f <- nxt
  }
  sum(f)
}
