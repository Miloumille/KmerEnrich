#' Binary bins x primers coverage matrix across targets
#'
#' @param patterns Admitted primer patterns.
#' @param targets Named list of [Biostrings::DNAStringSet] (from [load_targets()]).
#' @param read_len Effective read length (passed to [map_coverage()]).
#' @param bin_size Genome discretisation step.
#' @return List with `matrix` (bins x primers), `bin_target`, `patterns`, `targets`.
#' @export
coverage_matrix <- function(patterns, targets, read_len, bin_size = 100L) {
  bin_of <- function(pos) (pos - 1L) %/% bin_size + 1L
  blocks <- vector("list", length(targets))
  bin_target <- integer(0)

  for (t in seq_along(targets)) {
    seqs <- targets[[t]]
    n_bins <- as.integer(ceiling(max(nchar(as.character(seqs))) / bin_size))
    m <- matrix(0L, n_bins, length(patterns), dimnames = list(NULL, patterns))
    cov <- map_coverage(patterns, seqs, read_len)
    for (i in seq_len(nrow(cov))) {
      m[bin_of(cov$cov_start[i]):bin_of(cov$cov_end[i]), cov$pattern[i]] <- 1L
    }
    blocks[[t]] <- m
    bin_target <- c(bin_target, rep(t, n_bins))
  }

  list(matrix = do.call(rbind, blocks), bin_target = bin_target,
       patterns = patterns, targets = names(targets))
}
