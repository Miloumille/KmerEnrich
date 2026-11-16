#' Build an order-1 Markov model from raw mono/di counts
#'
#' @param mono Named length-4 vector of base counts.
#' @param di Named length-16 vector of dinucleotide counts.
#' @return List with `pi` and `trans` (Laplace-smoothed).
#' @export
markov_from_counts <- function(mono, di) {
  bases <- c("A", "C", "G", "T")
  mono <- mono + 1
  di <- di + 1
  pi <- mono[bases] / sum(mono[bases])
  trans <- matrix(0, 4L, 4L, dimnames = list(bases, bases))
  for (a in bases) for (b in bases) trans[a, b] <- di[[paste0(a, b)]]
  trans <- trans / rowSums(trans)
  list(pi = pi, trans = trans)
}

#' Fit an order-1 Markov background model from in-memory reads
#'
#' @param reads Character vector of background read sequences.
#' @return List with `pi` and `trans`.
#' @export
markov_background <- function(reads) {
  ss <- Biostrings::DNAStringSet(reads)
  markov_from_counts(colSums(Biostrings::oligonucleotideFrequency(ss, 1L)),
                     colSums(Biostrings::oligonucleotideFrequency(ss, 2L)))
}

#' Fit an order-1 Markov background by streaming FASTQ files
#'
#' Streams each FASTQ file in chunks, accumulating only mono/di counts, so memory
#' stays constant regardless of background size. FASTA files are read whole.
#'
#' @param paths Character vector of file paths (`.fastq/.fq[.gz]` or FASTA).
#' @param chunk_size Reads per FASTQ chunk.
#' @return List with `pi` and `trans`.
#' @export
markov_background_files <- function(paths, chunk_size = 5e5) {
  mono <- NULL; di <- NULL
  accumulate <- function(ss) {
    m <- colSums(Biostrings::oligonucleotideFrequency(ss, 1L))
    d <- colSums(Biostrings::oligonucleotideFrequency(ss, 2L))
    mono <<- if (is.null(mono)) m else mono + m
    di <<- if (is.null(di)) d else di + d
  }
  for (p in paths) {
    if (grepl("\\.(fastq|fq)(\\.gz)?$", p, ignore.case = TRUE)) {
      strm <- ShortRead::FastqStreamer(p, n = chunk_size)
      tryCatch(
        repeat {
          chunk <- ShortRead::yield(strm)
          if (length(chunk) == 0) break
          accumulate(ShortRead::sread(chunk))
        },
        finally = close(strm)
      )
    } else {
      accumulate(Biostrings::readDNAStringSet(p))
    }
  }
  markov_from_counts(mono, di)
}
