#' Load viral target genomes from a `virus/` folder
#'
#' Parses the input structure `virus/<target>/<strains>.fasta`. Each subfolder is
#' one viral target; the FASTA files inside are its strains. FASTA files placed
#' directly in `path` (no subfolder) are treated as a single target named
#' `"target"`.
#'
#' @param path Path to the `virus/` folder.
#' @return A named list of [Biostrings::DNAStringSet] objects, one per target.
#' @importFrom Biostrings readDNAStringSet
#' @export
load_targets <- function(path) {
  if (!dir.exists(path)) stop("Target path does not exist: ", path)
  fasta_ext <- "\\.(fasta|fa|fna)$"

  read_all <- function(files) {
    do.call(c, lapply(files, Biostrings::readDNAStringSet))
  }

  targets <- list()
  subdirs <- list.dirs(path, recursive = FALSE, full.names = TRUE)
  for (d in subdirs) {
    files <- list.files(d, pattern = fasta_ext, full.names = TRUE, ignore.case = TRUE)
    if (length(files) > 0) targets[[basename(d)]] <- read_all(files)
  }

  bare <- list.files(path, pattern = fasta_ext, full.names = TRUE, ignore.case = TRUE)
  if (length(bare) > 0) targets[["target"]] <- read_all(bare)

  if (length(targets) == 0) stop("No FASTA files found under: ", path)
  targets
}
