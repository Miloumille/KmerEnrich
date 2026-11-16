#' Merge Hamming-neighbour k-mers into degenerate IUPAC primers
#'
#' K-mers sharing the fixed 3' anchor are greedily grouped when their differing
#' positions stay within `max_degenerate`; each group becomes one IUPAC consensus.
#'
#' @param kmers Character vector of equal-length k-mers.
#' @param max_degenerate Maximum number of wobble positions per primer.
#' @param anchor Number of leading (3'-anchor) positions kept fixed.
#' @return Character vector of primers (IUPAC letters where merged).
#' @export
collapse_to_iupac <- function(kmers, max_degenerate = 3L, anchor = 2L) {
  if (length(kmers) == 0) return(character(0))
  mat <- do.call(rbind, strsplit(kmers, "", fixed = TRUE))
  free <- seq.int(anchor + 1L, ncol(mat))
  prefix <- apply(mat[, seq_len(anchor), drop = FALSE], 1, paste0, collapse = "")

  merge_bucket <- function(rows) {
    m <- mat[rows, , drop = FALSE]
    taken <- logical(nrow(m))
    primers <- character(0)
    for (i in which(!taken)) {
      if (taken[i]) next
      members <- i
      wobble <- integer(0)
      for (j in setdiff(which(!taken), i)) {
        d <- free[m[i, free] != m[j, free]]
        if (length(union(wobble, d)) <= max_degenerate) {
          members <- c(members, j)
          wobble <- union(wobble, d)
        }
      }
      taken[members] <- TRUE
      cons <- vapply(seq_len(ncol(m)), function(p) iupac_code(m[members, p]), character(1))
      primers <- c(primers, paste0(cons, collapse = ""))
    }
    primers
  }

  unlist(lapply(split(seq_along(kmers), prefix), merge_bucket), use.names = FALSE)
}
