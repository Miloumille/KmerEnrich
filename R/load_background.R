#' Load background reads from a `background/` folder
#'
#' Parses the input structure `background/<source>/<reads>.fastq[.gz]`. Each
#' subfolder is one background source (vector, host, ...). FASTQ reads and any
#' tolerated FASTA sequences inside are pooled into one character vector per
#' source. Paired-end `_1`/`_2` files are both read and pooled (no pairing logic
#' is needed for k-mer counting).
#'
#' @param path Path to the `background/` folder.
#' @return A named list of character vectors of read sequences, one per source.
#' @importFrom ShortRead readFastq sread
#' @importFrom Biostrings readDNAStringSet
#' @export
load_background <- function(path) {
  if (!dir.exists(path)) stop("Background path does not exist: ", path)
  fastq_ext <- "\\.(fastq|fq)(\\.gz)?$"
  fasta_ext <- "\\.(fasta|fa|fna)$"

  subdirs <- list.dirs(path, recursive = FALSE, full.names = TRUE)
  dirs <- if (length(subdirs) > 0) subdirs else path

  background <- list()
  for (d in dirs) {
    reads <- character(0)
    for (f in list.files(d, pattern = fastq_ext, full.names = TRUE, ignore.case = TRUE)) {
      reads <- c(reads, as.character(ShortRead::sread(ShortRead::readFastq(f))))
    }
    for (f in list.files(d, pattern = fasta_ext, full.names = TRUE, ignore.case = TRUE)) {
      reads <- c(reads, unname(as.character(Biostrings::readDNAStringSet(f))))
    }
    if (length(reads) > 0) {
      name <- if (length(subdirs) > 0) basename(d) else "background"
      background[[name]] <- reads
    }
  }

  if (length(background) == 0) stop("No FASTQ/FASTA files found under: ", path)
  background
}
