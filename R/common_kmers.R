#' Generate common observed Kmers
#'
#' @description The common_kmers function generates a list of all kmers that have been found in common in all sequences of all fasta files between position 0 and x.
#' Documentation has been partially AI generated
#'
#' @param virus_genome_path Path to the folder containing FASTA files or to a single FASTA file.
#'   If a folder, the function finds kmers common to all sequences across all FASTA files within that folder.
#'   If a single file, it finds unique kmers within that file's sequences.
#'   (Files must be in .fasta format).
#' @param k Value indicating the length of the kmers that need to be searched for.
#' @param x Value indicating the position until which the kmers need to be searched for.
#'   If `x` is larger than the length of the sequence, the function will take the whole sequence.
#'
#' @return A character vector of all kmers that have been found in common
#'   across all sequences (if `virus_genome_path` is a folder) or unique within
#'   the single sequence (if `virus_genome_path` is a file), extracted between position 0 and `x`.
#'
#' @details
#' This function concatenates all sequences within each FASTA file (ignoring header lines)
#' before extracting kmers. Kmers are extracted from position 0 up to (and including) position `x`.
#'
#' @examples
#' # Create dummy fasta files for demonstration
#' # (These files will be created and removed by the example itself)
#' # Use a temporary directory for examples to ensure clean execution
#' temp_example_dir <- file.path(tempdir(), "common_kmers_example")
#' dir.create(temp_example_dir, showWarnings = FALSE)
#'
#' # Create a subdirectory for the fasta files within the temp_example_dir
#' temp_fasta_files_dir <- file.path(temp_example_dir, "temp_fasta_dir")
#' dir.create(temp_fasta_files_dir, showWarnings = FALSE)
#'
#' # Write dummy fasta files to the temporary directory
#' writeLines(c(">seq1", "ATGCATGC", ">seq2", "GCATGCAT"),
#'            file.path(temp_fasta_files_dir, "genome1.fasta"))
#' writeLines(c(">seqA", "ATGCAATGCA", ">seqB", "TGCATGCA"),
#'            file.path(temp_fasta_files_dir, "genome2.fasta"))
#' writeLines(c(">seqX", "ATGCATGCGTAG", ">seqY", "TGCATGCGTAC"),
#'            file.path(temp_fasta_files_dir, "genome3.fasta"))
#'
#' # Example 1: Find common kmers in a folder
#' common_kmers_folder <- common_kmers(temp_fasta_files_dir, k = 3, x = 10)
#' print(common_kmers_folder)
#'
#' # Example 2: Find unique kmers in a single fasta file
#' common_kmers_single_file <- common_kmers(file.path(temp_fasta_files_dir, "genome1.fasta"), k = 3, x = 10)
#' print(common_kmers_single_file)
#'
#' # Clean up dummy files and directory (important for clean examples)
#' unlink(temp_example_dir, recursive = TRUE)
#' @export
common_kmers <- function(virus_genome_path, k, x) {
  find_kmers <- function(sequence, k) {
    if (nchar(sequence) < k) { # Handle case where sequence is shorter than k
      return(character(0))
    }
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
