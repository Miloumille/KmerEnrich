#' Generate statistics about the kmers included in the RNAseq files
#'
#' @description The get_vector_stats_kmers function generates a list of all kmers that have been found in common in all sequences of all fasta files between position 0 and x.
#' @param v_kmers_stats Dataframe containing the kmers that need to be searched for. The dataframe must contain a column "Kmer" with the kmers that need to be searched for.
#' @param vector_genome_path path to the folder containing the fastq files or to the fastq file for the mosquito RNA sequences. Can be a folder or a fastq file. (files must be in .fastq format). The fastq files must end with _1.fastq and _2.fastq in order to be accessed properly. The input path on the contrary must be the part just before the _1.fastq or _2.fastq.
#'
#' @return The function will complete the "v_kmers_stats" Dataframe by adding for each kmer the frequence at which it appears in the fastq files as well as the amount of time they appeared.
#'
#' @examples
#' # Example of how to use the function
#' k_mers <- gen_obs_kmers_one_fasta("data_processed/virus/reference/Mayaro_virus.fasta", 4, 20)
#' print(k_mers)
#'
#' @export
get_vector_stats_kmers <- function(v_kmers_stats, vector_genome_path) {
  count_kmer_in_fastq <- function(fq, specific_kmer) {
    kmer_count <- sum(str_count(fq, specific_kmer))
    total_size <- sum(nchar(fq))
    freq_m <- if (total_size > 0) kmer_count / total_size else 0
    return(list(kmer_count = kmer_count, freq_m = freq_m))
  }
  all_results_df <- list()
  for (path in vector_genome_path) {
    fq1 <- readFastq(paste0(path, "_1.fastq"))
    fq2 <- readFastq(paste0(path, "_2.fastq"))
    fq <- c(as.character(sread(fq1)), as.character(sread(fq2)))
    temp_results <- v_kmers_stats
    temp_results$kmer_count <- integer(length(temp_results$Kmer))
    temp_results$freq_m <- numeric(length(temp_results$Kmer))
    for (i in seq_along(temp_results$Kmer)) {
      kmer <- temp_results$Kmer[i]
      fastq_count <- count_kmer_in_fastq(fq, kmer)
      temp_results$kmer_count[i] <- fastq_count$kmer_count
      temp_results$freq_m[i] <- fastq_count$freq_m
    }
    temp_results$accession_id <- basename(path)
    all_results_df[[length(all_results_df) + 1]] <- temp_results
  }
  all_results_df <- rbindlist(all_results_df)
  return(all_results_df)
}
