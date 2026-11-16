#' IUPAC ambiguity letter for a set of bases
#'
#' @param bases Character vector of `A`/`C`/`G`/`T`.
#' @return Single IUPAC letter covering the distinct bases in `bases`.
#' @export
iupac_code <- function(bases) {
  key <- paste(sort(unique(bases)), collapse = "")
  codes <- c(A = "A", C = "C", G = "G", T = "T",
             AG = "R", CT = "Y", CG = "S", AT = "W", GT = "K", AC = "M",
             CGT = "B", AGT = "D", ACT = "H", ACG = "V", ACGT = "N")
  codes[[key]]
}
