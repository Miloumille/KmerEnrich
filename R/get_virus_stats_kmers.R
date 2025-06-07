#' Generate statistics of each k-mer
#'
#' @description Calculates various statistical metrics for each k-mer based on
#'   their positions within virus genomes. It aggregates information such as
#'   appearance frequency, segment sizes between k-mers, and k-mer counts
#'   across different sequences.
#'
#'   Documentation has been partially AI generated
#'
#' @param df_virus A `data.frame` (or `tibble`), typically the output of the
#'   `kmers_pos_df` function. It MUST contain the following columns:
#'   \itemize{
#'     \item `Kmer` (`character`): The k-mer sequence.
#'     \item `SequenceName` (`character`): The name of the sequence where the k-mer was found.
#'     \item `Position` (`numeric`): The starting position of the k-mer in the sequence.
#'     \item `Type` (`character`): The subtype of the virus (if applicable).
#'     \item `segm_size` (`numeric`): The size of the segment between the current k-mer
#'       and the preceding k-mer. This column will have `NA` for the first k-mer
#'       of each sequence, which the function handles.
#'     \item `log10_segm_size` (`numeric`): The base-10 logarithm of (`segm_size` + 1).
#'   }
#'
#' @return A `data.frame` (or `tibble`) summarized by `Kmer`, containing the
#'   following statistical metrics:
#'   \itemize{
#'     \item `mean_mean_log10_segm_size`: The mean of the `log10_segm_size` values
#'       (segment sizes between k-mers) across all occurrences of a given k-mer,
#'       averaged per sequence and then across sequences.
#'     \item `max_max_log10_segm_size`: The maximum `log10_segm_size` observed
#'       for a k-mer across all sequences.
#'     \item `mean_sd_log10_segm_size`: The mean of the standard deviations of
#'       `log10_segm_size` for each k-mer within each sequence, averaged across
#'       all sequences.
#'     \item `mean_count`: The mean number of times a k-mer appeared in individual
#'       virus sequences.
#'     \item `sd_count`: The standard deviation of the k-mer counts across
#'       different virus sequences. Is `NA` if there's only one observation.
#'     \item `mean_freq_virus`: The mean frequency of the k-mer across different
#'       virus sequences (calculated as k-mer count / total sequence length in
#'       each sequence).
#'   }
#'
#' @details
#' This function performs a two-step aggregation:
#' \enumerate{
#'   \item First Level Aggregation (by `SequenceName`, `Kmer`):
#'     It calculates `mean_log10_segm_size`, `max_log10_segm_size`,
#'     `sd_log10_segm_size`, `count` (k-mer count within that sequence),
#'     and `freq` (k-mer count / sequence length) for each unique k-mer
#'     within each individual virus sequence. Rows where `segm_size` is `NA`
#'     (corresponding to the implied segment before the first k-mer, or end of sequence)
#'     are filtered out before this step, as `segm_size` is not defined for them.
#'   \item Second Level Aggregation (by `Kmer`):
#'     These per-sequence statistics are then further summarized across
#'     all sequences to yield the final, overall statistics for each k-mer.
#' }
#' The `sd_log10_segm_size` is set to `0` if `n() <= 1` to prevent `sd()` errors
#' when there's insufficient data for a standard deviation. Similarly, `sd_count`
#' is set to `NA` under the same condition.
#'
#' @export
get_virus_stats_kmers <- function(df_virus) {
  # Drop lines that correspond to the distance between first nucleotide and first kmer
  df_virus <- df_virus %>%
    filter(!is.na(segm_size))

  virus_values <- df_virus %>%
    group_by(SequenceName, Kmer) %>%
    summarize(
      mean_log10_segm_size = mean(log10(segm_size), na.rm = TRUE),
      max_log10_segm_size = max(log10(segm_size), na.rm = TRUE),
      sd_log10_segm_size = ifelse(n() > 1, sd(log10(segm_size), na.rm = TRUE), 0),
      count = n(),
      freq = n() / max(Position),
      .groups = 'drop'
    )

  results <- virus_values %>%
    group_by(Kmer) %>%
    summarize(
      mean_mean_log10_segm_size = mean(mean_log10_segm_size, na.rm = TRUE),
      max_max_log10_segm_size = max(max_log10_segm_size, na.rm = TRUE),
      mean_sd_log10_segm_size = mean(sd_log10_segm_size, na.rm = TRUE),
      mean_count = mean(count, na.rm = TRUE),
      sd_count = ifelse(n() > 1, sd(count, na.rm = TRUE), NA),
      mean_freq_virus = mean(freq, na.rm = TRUE),
      .groups = 'drop'
    )
  return(results)
}


