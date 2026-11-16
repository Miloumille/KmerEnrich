#' Enumerate k-mers conserved across all strains of a target
#'
#' For each strain, extract the set of distinct k-mers (optionally only within
#' the first `window` bases from the 5' end), then intersect across strains so
#' only k-mers present in every strain remain. Whole-genome by default, matching
#' the SMART coverage geometry where priming sites are needed genome-wide.
#'
#' @param target_seqs A [Biostrings::DNAStringSet] of one target's strains.
#' @param k k-mer length (default 9, the SMART-9N anchor length).
#' @param window Optional integer; restrict to the first `window` bases of each
#'   strain. `NULL` (default) uses the whole genome.
#' @return A character vector of conserved k-mers.
#' @export
enumerate_conserved_kmers <- function(target_seqs, k = 9, window = NULL) {
  kmer_set <- function(seq) {
    s <- as.character(seq)
    if (!is.null(window)) s <- substr(s, 1, window)
    n <- nchar(s)
    if (n < k) return(character(0))
    unique(substring(s, 1:(n - k + 1), k:n))
  }

  sets <- lapply(seq_along(target_seqs), function(i) kmer_set(target_seqs[[i]]))
  Reduce(intersect, sets)
}
