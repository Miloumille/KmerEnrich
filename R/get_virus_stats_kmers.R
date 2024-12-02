#' Generate the statistics of each kmer
#'
#' @description The get_stats_kmers function will take the output of the kmers_pos_df and gather statistical information about the frequency of appearence of each kmer, the mean of the segments separating kmers, the mean count of each kmers in their virus, the mean standart deviation of the segments seperating kmers and the stadart deviation of the mean count in order to check how regularly this kmer appears in the virus genome.
#'
#' @param df_virus The input is a Dataframe that is typically the output of the "kmers_pos_df" function. A Dataframe containing a columns for each Kmer, a column containing the SequenceName, a column containing the Position of the Kmer, a column containing the Type of the virus (subtype), a column containing the segment size (distance between previous kmer and actual kmer) and a last column containing the log10 value of the segment size.
#'
#' @return This function will return a dataframe containing the statistical information about the frequency of appearence of each kmer, the mean of the segments separating kmers, the mean count of each kmers in their virus, the mean standart deviation of the segments seperating kmers and the stadart deviation of the mean count in order to check how regularly this kmer appears in the virus genome.
#'
#' @export
get_virus_stats_kmers <- function(df_virus) {
  # Drop lines that correspond to the distance between first nucleotide and first kmer
  df_virus <- df_virus %>%
    filter(!is.na(segm_size))

  virus_values <- df_virus %>%
    group_by(SequenceName, Kmer) %>%
    summarize(mean_segm_size = mean(log10(segm_size), na.rm = TRUE),
              sd_segm_size = ifelse(n() > 1, sd(log10(segm_size), na.rm = TRUE), 0),
              count = n(),
              freq = n()/max(Position),
              .groups = 'drop')

  results <- virus_values %>%
    group_by(Kmer) %>%
    summarize(
      mean_mean_log10_segments = mean(mean_segm_size, na.rm = TRUE),
      mean_sd_log10_segments = mean(sd_segm_size, na.rm = TRUE),
      mean_count = mean(count, na.rm = TRUE),
      sd_count = ifelse(n() > 1, sd(count, na.rm = TRUE), NA),
      mean_freq_virus = mean(freq, na.rm = TRUE),
      .groups = 'drop'
    )
  return(results)
}


