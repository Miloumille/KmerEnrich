#' Every concrete instantiation of an IUPAC pattern
#'
#' @param pattern Primer pattern (IUPAC allowed).
#' @return Character vector of concrete `A/C/G/T` primers.
#' @export
expand_iupac <- function(pattern) {
  sets <- lapply(strsplit(pattern, "", fixed = TRUE)[[1]], iupac_bases)
  do.call(paste0, expand.grid(sets, stringsAsFactors = FALSE))
}

#' Thermodynamic entry gate for a primer
#'
#' All IUPAC instantiations must sit within the Tm and GC windows, and the 3'
#' anchor (first `gc_clamp_len` bases) must carry a G or C.
#'
#' @param primer Primer pattern (IUPAC allowed).
#' @param tm_min,tm_max Melting-temperature window (degrees Celsius).
#' @param gc_min,gc_max GC-content window (fraction).
#' @param gc_clamp_len Length of the 3' anchor checked for a G/C.
#' @return Single logical.
#' @export
passes_thermo <- function(primer, tm_min = 15, tm_max = 45,
                          gc_min = 0.2, gc_max = 0.8, gc_clamp_len = 2) {
  variants <- expand_iupac(primer)
  tm <- vapply(variants, primer_tm, numeric(1))
  gc <- vapply(variants, gc_content, numeric(1))
  clamp <- grepl("[GC]", substr(primer, 1, gc_clamp_len))
  clamp && all(tm >= tm_min & tm <= tm_max) && all(gc >= gc_min & gc <= gc_max)
}
