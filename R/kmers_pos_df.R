#' Generate a dataframe of position of all kmers, distance between kmers for the fasta file(s)
#'
#' @description The kmers_pos_df function is used to generate a dataframe containing the position of every kmer on every sequence of every fasta file. It also contains the size of the segments betweem 2 kmers or between the first nucleotide and the first kmer or the last kmer and the last nucleotide.
#'
#' @param virus_genome_path The input is a string that can be either a path to a fasta file or a path to a folder containing fasta files.
#' @param kmers_list The input is a list of strings containing the kmers to search for.
#'
#' @return This function will return a dataframe containing the position of every kmer on every sequence of every fasta file with specific characteristics as the sequence name, the position of the kmer, the type of the virus (subtype), the segment size (distance between previous kmer and actual kmer) and the log10 value of the segment size.
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
          segm_size = Position - lag(Position),
          log10_segm_size = ifelse(is.na(segm_size), NA, log10(segm_size + 1))
        ) %>%
        ungroup()

      # Combine data of all fasta files
      all_kmers_df[[length(all_kmers_df) + 1]] <- kmers_df
    }

    # Transform list into dataframe containing all information
    final_kmers_df <- bind_rows(all_kmers_df)
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
        segm_size = Position - lag(Position),
        log10_segm_size = ifelse(is.na(segm_size), NA, log10(segm_size + 1))
      ) %>%
      ungroup()
  }
  return(final_kmers_df)
}
