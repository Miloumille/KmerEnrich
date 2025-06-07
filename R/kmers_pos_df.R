#' Generate a dataframe of k-mer positions and segment sizes
#'#'
#' @description Generates a data frame detailing the absolute positions of
#'   specified k-mers within each sequence of input FASTA file(s). It also
#'   calculates the size of segments between consecutive k-mers, as well as the
#'   initial segment from the sequence start to the first k-mer and the final
#'   segment from the last k-mer to the sequence end.
#'
#'   Documentation has been partially AI generated
#'
#' @param virus_genome_path A `character` string specifying the path to either:
#'   \itemize{
#'     \item A single `.fasta` file. In this case, the `Type` column in the
#'       output will be extracted from the base name of the file
#'       (e.g., `'virusA'` from `'path/to/virusA.fasta'`).
#'     \item A folder containing one or more `.fasta` files. If a folder is
#'       provided, the function will process all `.fasta` files within it
#'       (including subdirectories due to `recursive = TRUE`). The `Type` column
#'       for each sequence will be derived from the name of the immediate parent
#'       directory of the FASTA file.
#'   }
#' @param kmers_list A `character` vector containing the k-mer sequences to
#'   search for within the genomes.
#'
#' @return A `data.frame` (or `tibble`) with one row per k-mer occurrence and
#'   per sequence boundary marker. It includes the following columns:
#'   \itemize{
#'     \item `SequenceName` (`character`): The name of the sequence from the FASTA file.
#'     \item `Kmer` (`character`): The k-mer sequence found.
#'     \item `Position` (`numeric`): The 1-based starting position of the k-mer in the
#'       sequence. For the 'segment to sequence end' marker, this is the
#'       `sequence_length`.
#'     \item `Type` (`character`): The virus subtype, extracted as described
#'       in `virus_genome_path`.
#'     \item `segm_size` (`numeric`): The size of the segment. This is calculated
#'       as `Position - dplyr::lag(Position)`. For the very first k-mer of a
#'       sequence (or the initial segment from position 0), this value will be `NA`.
#'     \item `log10_segm_size` (`numeric`): The base-10 logarithm of (`segm_size` + 1).
#'       This value will be `NA` where `segm_size` is `NA`.
#'   }
#'
#' @details
#' For each sequence and each k-mer in `kmers_list`, the function identifies
#' all occurrences.
#'
#' To facilitate segment size calculation, the function adds an artificial
#' 'final position entry' for each k-mer equal to the `sequence_length`.
#' This allows calculation of the segment from the last observed k-mer to
#' the end of the sequence.
#'
#' The `segm_size` column is calculated by subtracting the position of the
#' previous k-mer (or sequence start) from the current k-mer's position.
#' The first `segm_size` for each k-mer within each sequence will be `NA`
#' because there's no preceding k-mer to calculate the distance from.
#' These `NA` values typically represent the segment from the start of the
#' sequence to the first observed k-mer.
#'
#' The `log10_segm_size` is calculated as `log10(segm_size + 1)`. The `+ 1`
#' is added to handle potential `segm_size` values of 0 (though unlikely
#' for k-mer distances) and to ensure valid logarithmic calculations.
#'
#' If `virus_genome_path` points to a folder, the function recursively searches
#' for all `.fasta` files within that folder.
#'
#' @export
kmers_pos_df <- function(virus_genome_path, kmers_list) {
  # Small function to search kmers
  find_kmer_positions <- function(sequence, kmer) {
    positions <- gregexpr(kmer, sequence)[[1]]
    positions <- positions[positions > 0]

    return(positions)
  }

  # In case the virus_genome_path if a folder that contains fasta files
  if (substr(virus_genome_path, nchar(virus_genome_path), nchar(virus_genome_path)) == "/") {
    fasta_files <- list.files(virus_genome_path, pattern = "\\.fasta$", recursive = TRUE, full.names = TRUE)
    all_kmers_df <- list()
    for (fasta_file in fasta_files) {
      # Extract the virus type from the file name
      virus_type <- basename(dirname(fasta_file))

      genomes <- readDNAStringSet(fasta_file)
      kmers_df <- data.frame(SequenceName = character(), Kmer = character(), Position = integer(), Type = character())

      for (i in seq_along(genomes)) {
        sequence_name <- names(genomes)[i]
        sequence <- as.character(genomes[i])
        sequence_length <- nchar(sequence)

        for (kmer in kmers_list) {
          positions <- find_kmer_positions(sequence, kmer)
          if (length(positions) > 0) {
            temp_df <- data.frame(SequenceName = sequence_name, Kmer = kmer, Position = positions, Type = virus_type)
            kmers_df <- rbind(kmers_df, temp_df)
          }

          #Adding the final position entry (end of sequence)
          final_position_df <- data.frame(SequenceName = sequence_name, Kmer = kmer, Position = sequence_length, Type = virus_type)
          kmers_df <- rbind(kmers_df, final_position_df)
        }
      }
      # structuring dataframe to contain all information
      kmers_df <- kmers_df %>%
        filter(Kmer %in% kmers_list) %>%
        arrange(SequenceName, Kmer, Position) %>%
        group_by(SequenceName, Kmer) %>%
        mutate(
          segm_size = Position - dplyr::lag(Position),
          log10_segm_size = ifelse(is.na(segm_size), NA, log10(segm_size + 1))
        ) %>%
        ungroup()

      # Combine data of all fasta files
      all_kmers_df[[length(all_kmers_df) + 1]] <- kmers_df
    }

    # Transform list into dataframe containing all information
    final_kmers_df <- dplyr::bind_rows(all_kmers_df)
  }


  # In case the virus_genome_path is directly a fasta file
  else{
    virus_type<- sub(".*/(.*)\\.fasta$", "\\1", virus_genome_path)
    genomes <- readDNAStringSet(virus_genome_path)
    kmers_df <- data.frame(SequenceName = character(), Kmer = character(), Position = integer(), Type = character())

    for (i in seq_along(genomes)) {
      sequence_name <- names(genomes)[i]
      sequence <- as.character(genomes[i])
      sequence_length <- nchar(sequence)

      for (kmer in kmers_list) {
        positions <- find_kmer_positions(sequence, kmer)
        if (length(positions) > 0) {
          temp_df <- data.frame(SequenceName = sequence_name, Kmer = kmer, Position = positions, Type = virus_type)
          kmers_df <- rbind(kmers_df, temp_df)
        }
        final_position_df <- data.frame(SequenceName = sequence_name, Kmer = kmer, Position = sequence_length, Type = virus_type)
        kmers_df <- rbind(kmers_df, final_position_df)
      }
    }
    # structuring dataframe to contain all information
    final_kmers_df <- kmers_df %>%
      filter(Kmer %in% kmers_list) %>%
      arrange(SequenceName, Kmer, Position) %>%
      group_by(SequenceName, Kmer) %>%
      mutate(
        segm_size = Position - dplyr::lag(Position),
        log10_segm_size = ifelse(is.na(segm_size), NA, log10(segm_size + 1))
      ) %>%
      ungroup()
  }
  return(final_kmers_df)
}
