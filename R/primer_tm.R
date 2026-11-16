#' G/C fraction of a primer
#'
#' @param primer DNA string.
#' @return Fraction of bases that are G or C.
#' @export
gc_content <- function(primer) {
  b <- strsplit(primer, "", fixed = TRUE)[[1]]
  mean(b %in% c("G", "C"))
}

#' Nearest-neighbour melting temperature (SantaLucia 1998)
#'
#' @param primer Concrete DNA oligo (`A/C/G/T`).
#' @param na Sodium concentration (M) for the salt correction.
#' @param conc Primer strand concentration (M).
#' @return Melting temperature in degrees Celsius.
#' @export
primer_tm <- function(primer, na = 0.05, conc = 5e-7) {
  h <- c(AA = -7.9, AT = -7.2, TA = -7.2, CA = -8.5, GT = -8.4,
         CT = -7.8, GA = -8.2, CG = -10.6, GC = -9.8, GG = -8.0)
  s <- c(AA = -22.2, AT = -20.4, TA = -21.3, CA = -22.7, GT = -22.4,
         CT = -21.0, GA = -22.2, CG = -27.2, GC = -24.4, GG = -19.9)
  comp <- c(A = "T", T = "A", G = "C", C = "G")
  key <- function(d) if (!is.na(h[d])) d else paste0(comp[[substr(d, 2, 2)]], comp[[substr(d, 1, 1)]])
  init <- function(base) if (base %in% c("G", "C")) c(0.1, -2.8) else c(2.3, 4.1)

  b <- strsplit(primer, "", fixed = TRUE)[[1]]
  n <- length(b)
  dh <- 0; ds <- 0
  for (i in seq_len(n - 1)) {
    k <- key(paste0(b[i], b[i + 1]))
    dh <- dh + h[[k]]; ds <- ds + s[[k]]
  }
  e5 <- init(b[1]); e3 <- init(b[n])
  dh <- dh + e5[1] + e3[1]
  ds <- ds + e5[2] + e3[2] + 0.368 * (n - 1) * log(na)
  1000 * dh / (ds + 1.987 * log(conc / 4)) - 273.15
}
