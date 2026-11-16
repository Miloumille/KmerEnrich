#' Observed per-position frequency of a primer in a sequence set
#'
#' @param pattern Primer pattern (IUPAC allowed).
#' @param seqs A [Biostrings::DNAStringSet].
#' @return Numeric frequency (matches / start positions), 0 if no positions.
#' @export
observed_freq <- function(pattern, seqs) {
  k <- nchar(pattern)
  hits <- sum(Biostrings::vcountPattern(pattern, seqs, fixed = c(pattern = FALSE, subject = TRUE)))
  positions <- sum(pmax(0L, Biostrings::width(seqs) - k + 1L))
  if (positions == 0) 0 else hits / positions
}

#' Score viral enrichment of primers as a fold value
#'
#' Observed viral frequency divided by the frequency the pre-fitted order-1 Markov
#' background model predicts. Higher means more viral-specific.
#'
#' @param patterns Character vector of primer patterns.
#' @param target_seqs A [Biostrings::DNAStringSet] of one target's strains.
#' @param model A fitted model from [markov_background()] or [markov_background_files()].
#' @return data.frame with `pattern`, `freq_virus`, `exp_background`, `enrichment`.
#' @export
score_enrichment <- function(patterns, target_seqs, model) {
  obs <- vapply(patterns, observed_freq, numeric(1), seqs = target_seqs)
  exp <- vapply(patterns, expected_freq, numeric(1), model = model)
  data.frame(pattern = patterns, freq_virus = obs, exp_background = exp,
             enrichment = obs / exp, row.names = NULL, stringsAsFactors = FALSE)
}
