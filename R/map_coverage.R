#' Map primers to the genome intervals they cover
#'
#' Finds IUPAC-aware binding sites of each primer in every strain and returns the
#' interval each read covers under the SMART geometry (extension toward the 5' end
#' by default).
#'
#' @param patterns Character vector of primer patterns (IUPAC allowed).
#' @param target_seqs A [Biostrings::DNAStringSet] of one target's strains.
#' @param read_len Effective read length extending from the binding site.
#' @param direction `"5prime"` (default) or `"3prime"`.
#' @return data.frame with `pattern`, `strain`, `start`, `cov_start`, `cov_end`.
#' @export
map_coverage <- function(patterns, target_seqs, read_len, direction = c("5prime", "3prime")) {
  direction <- match.arg(direction)
  empty <- data.frame(pattern = character(0), strain = character(0),
                      start = integer(0), cov_start = integer(0), cov_end = integer(0),
                      stringsAsFactors = FALSE)

  hits_for <- function(pat, seq, strain) {
    k <- nchar(pat)
    w <- nchar(as.character(seq))
    if (w < k) return(NULL)
    m <- Biostrings::matchPattern(pat, seq, fixed = c(pattern = FALSE, subject = TRUE))
    if (length(m) == 0) return(NULL)
    starts <- Biostrings::start(m)
    if (direction == "5prime") {
      cs <- pmax(1L, starts - as.integer(read_len)); ce <- starts + k - 1L
    } else {
      cs <- starts; ce <- pmin(w, starts + k - 1L + as.integer(read_len))
    }
    data.frame(pattern = pat, strain = strain, start = starts,
               cov_start = cs, cov_end = ce, stringsAsFactors = FALSE)
  }

  rows <- list()
  for (pat in patterns) {
    for (i in seq_along(target_seqs)) {
      rows[[length(rows) + 1L]] <- hits_for(pat, target_seqs[[i]], names(target_seqs)[i])
    }
  }
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0) return(empty)
  do.call(rbind, rows)
}
