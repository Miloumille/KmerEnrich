#' Generate statistics about the kmers included in the RNAseq files
#'
#' @description Calculates the count and frequency of specified k-mers within
#'   RNA-seq data from paired-end FASTQ files, enriching a provided k-mer
#'   statistics dataframe.
#'
#'   Documentation has been partially AI generated
#'
#' @param v_kmers_stats A `data.frame` (or `tibble`) containing k-mer statistics,
#'   which will be enriched by this function. It MUST include a column named
#'   `Kmer` (type `character`) with the k-mers to be analyzed. Additionally,
#'   it MUST contain `mean_freq_virus`, `mean_count`, and
#'   `max_max_log10_segm_size` columns (type `numeric`), which are used
#'   to calculate the `Enrichment_score` and `Distribution_score`.
#' @param vector_genome_path A `character` vector where each element is the
#'   BASE PATH (filename without `_1.fastq` or `_2.fastq` suffix) for
#'   a pair of paired-end FASTQ files. For example, if your files are
#'   `sampleA_1.fastq` and `sampleA_2.fastq`, the corresponding element in
#'   this vector should be `'path/to/sampleA'`. The function expects both
#'   `_1.fastq` and `_2.fastq` files to exist for each base path. These files
#'   contain the RNA sequences (e.g., mosquito RNA sequences).
#'
#' @return A `data.frame` (or `tibble`) that builds upon the input `v_kmers_stats`.
#'   It includes the following new columns for each k-mer and each input FASTQ file pair:
#'   \itemize{
#'     \item `kmer_count`: The total number of times the k-mer appeared across
#'       both `_1.fastq` and `_2.fastq` reads for a given sample.
#'     \item `freq_vector`: The frequency of the k-mer (k-mer count divided
#'       by total sequence length) in the FASTQ files for a given sample.
#'     \item `accession_id`: The base name of the FASTQ file pair (e.g., 'sampleA'
#'       from 'path/to/sampleA').
#'     \item `Enrichment_score`: Calculated as `mean_freq_virus / freq_vector`.
#'     \item `Distribution_score`: Calculated as `mean_count / max_max_log10_segm_size`.
#'   }
#'
#' @details
#' This function processes paired-end FASTQ files, expecting them to be named
#' with `_1.fastq` and `_2.fastq` suffixes (e.g., `sample_A_1.fastq`,
#' `sample_A_2.fastq`). The `vector_genome_path` parameter should provide
#' the base path *without* these suffixes.
#'
#' For each k-mer in `v_kmers_stats$Kmer` and each pair of FASTQ files,
#' the function concatenates all reads from both files and calculates the
#' k-mer's total count and its frequency (count / total sequence length).
#'
#' The `Enrichment_score` and `Distribution_score` are calculated using
#' pre-existing columns from the input `v_kmers_stats` dataframe:
#' `mean_freq_virus`, `mean_count`, and `max_max_log10_segm_size`.
#' Ensure these columns are present and contain relevant data from
#' your viral genome analysis.
#'
#' This function relies on `ShortRead::readFastq` and `stringr::str_count`.
#'
#' @export
get_vector_stats_kmers <- function(v_kmers_stats, vector_genome_path) {
  count_kmer_in_fastq <- function(fq, specific_kmer) {
    kmer_count <- sum(str_count(fq, specific_kmer))
    total_size <- sum(nchar(fq))
    freq_vector <- if (total_size > 0) kmer_count / total_size else 0
    return(list(kmer_count = kmer_count, freq_vector = freq_vector))
  }
  all_results_df <- list()
  for (path in vector_genome_path) {
    fq1 <- readFastq(paste0(path, "_1.fastq"))
    fq2 <- readFastq(paste0(path, "_2.fastq"))
    fq <- c(as.character(sread(fq1)), as.character(sread(fq2)))
    temp_results <- v_kmers_stats
    temp_results$kmer_count <- integer(length(temp_results$Kmer))
    temp_results$freq_vector <- numeric(length(temp_results$Kmer))
    for (i in seq_along(temp_results$Kmer)) {
      kmer <- temp_results$Kmer[i]
      fastq_count <- count_kmer_in_fastq(fq, kmer)
      temp_results$kmer_count[i] <- fastq_count$kmer_count
      temp_results$freq_vector[i] <- fastq_count$freq_vector
    }
    temp_results$accession_id <- basename(path)
    all_results_df[[length(all_results_df) + 1]] <- temp_results
  }
  all_results_df <- rbindlist(all_results_df)

  all_results_df$Enrichment_score <- all_results_df$mean_freq_virus / all_results_df$freq_vector
  all_results_df$Distribution_score <- all_results_df$mean_count / all_results_df$max_max_log10_segm_size
  return(all_results_df)
}
