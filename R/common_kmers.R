#' Generate common observed Kmers
#'
#' @description The common_kmers function generates a list of all kmers that have been found in common in all sequences of all fasta files between position 0 and x.
#'
#' @param virus_genome_path path to the folder containing the fasta files or to the fasta file. Can be a folder or a fasta file. (files must be in .fasta format).
#' @param k value indicating the length of the kmers that need to be searched for.
#' @param x value indicating the position until which the kmers need to be searched for. If x is larger than the length of the sequence, the function will take the whole sequence.
#'
#' @return The function will return a list of all kmers that have been found in common in all sequences of all fasta files between position 0 and x.
#'
#' @examples
#' # Example of how to use the function
#' k_mers <- gen_obs_kmers_one_fasta("data_processed/virus/reference/Mayaro_virus.fasta", 4, 20)
#' print(k_mers)
#' k_mers <- gen_obs_kmers_one_fasta("../data_processed/virus/all_dengue/", 4, 100)
#' print(k_mers)
#'
#' @export
#'
common_kmers <- function(virus_genome_path, k, x) {
  find_kmers <- function(sequence, k) {
    kmers <- sapply(1:(nchar(sequence) - k + 1), function(i) substr(sequence, i, i + k - 1))
    return(unique(kmers))
  }
  if (substr(virus_genome_path, nchar(virus_genome_path), nchar(virus_genome_path)) == "/") {
    fasta_files <- list.files(virus_genome_path, pattern = "\\.fasta$", recursive = TRUE, full.names = TRUE)
    kmer_lists <- list()
    for (file in fasta_files) {
      sequences <- readLines(file)
      concatenated_sequence <- paste(sequences[!grepl("^>", sequences)], collapse = "")
      truncated_sequence <- substr(concatenated_sequence, 1, x)
      kmer_lists[[file]] <- find_kmers(truncated_sequence, k)
    }
    all_kmers <- Reduce(intersect, kmer_lists)
  }
  else{
    sequences <- readLines(virus_genome_path)
    concatenated_sequence <- paste(sequences[!grepl("^>", sequences)], collapse = "")
    truncated_sequence <- substr(concatenated_sequence, 1, x)
    all_kmers <- find_kmers(truncated_sequence, k)
  }
  return(all_kmers)
}

