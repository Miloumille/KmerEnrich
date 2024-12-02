#' Generate an HTML report
#'
#' @description The KmerEnrichFullReport function generates an HTML report from a template included in the package. The report includes a summary of the results of the KmerEnrich pipeline. This summary will give the main characteristics of the kmers found in the virus genome and the vector genome. The report will also include a plot of the enchrichment fo the Kmer compared in the virus versus the mosquito. The x value is the ratio of the mean count of the Kmer in the virus genome to the mean of the standard deviation of the segments separating the Kmers. The y value is the log10 of the mean frequency of the Kmer in the virus genome minus the log10 of the frequency of the Kmer in the vector genome. The report will also include a plot of the position of the Kmers in the virus genome(s).
#'
#' @param virus_genome_path The input is a string that can be either a path to a fasta file or a path to a folder containing fasta files.
#' @param vector_genome_path The input is a string that can be either a path to a fasta file or a path to a folder containing fasta files.
#' @param k_list A list of integers indicating the Kmers to search for in the virus genome. The report will include the Kmers that appear in all the virus genomes.
#' @param x Value of x indicating the x first bases in the Virus gneome in which the Kmers will be searched. The report will only keep the Kmers appearing in all the virus genomes.
#' @param report_name The name of the report file to generate. Default is "report.html".
#'
#' @return This function will generate an HTML report in the working directory.
#' @export
KmerEnrichFullReport <- function(virus_genome_path, vector_genome_path, k_list, x, report_name = "report.html") {
  template_path <- system.file("rmd", "KmerEnrichFullReport.Rmd", package = "KmerEnrich")
  if (template_path == "") {
    stop("Template Rmd file not found in the package.")
  }

  ############# Pipeline ###############
  results = list()
  for (k in c(k_list)){
    best_kmers_list <- common_kmers(virus_genome_path,k,x)
    v_kmer_pos_df <- kmers_pos_df(virus_genome_path,best_kmers_list)
    v_kmers_stats <- get_virus_stats_kmers(v_kmer_pos_df)
    complete_kmer_stats <- get_vector_stats_kmers(v_kmers_stats, vector_genome_path)
    results[[paste0("k_", k)]] <- complete_kmer_stats
  }

  results <- lapply(results, function(table) {
    if (is.data.frame(table)) {
      required_columns <- c("mean_count", "mean_sd_log10_segments", "mean_freq_virus", "freq_m", "accession_id", "Kmer")
      if (all(required_columns %in% names(table))) {
        # Calculate x_value and y_value
        table$x_value <- table$mean_count / mean(c(mean(table$mean_sd_log10_segments, na.rm = TRUE),
                                                   table$mean_sd_log10_segments), na.rm = TRUE)
        table$y_value <- log10(table$mean_freq_virus) - log10(table$freq_m)
      }
    }
    return(table)
  })

  ########### G position ###############

  fasta_files <- list.files(virus_genome_path, pattern = "\\.fasta$", recursive = TRUE, full.names = TRUE)

  find_kmer_positions <- function(virus_genome_path, kmer) {
    genome <- readDNAStringSet(virus_genome_path)[[1]]
    kmer_pos <- matchPattern(kmer, genome)
    starts <- Biostrings::start(kmer_pos)
    return(data.frame(Position = starts, Genome = basename(virus_genome_path),kmer = kmer))
  }

  plot_data <- list()
  for (result in results){
    sorted_table <- result[order(-result$x_value), ]
    top_5_kmers <- head(sorted_table$Kmer, 5)
    for (kmer in top_5_kmers){
      kmer_data_list <- list()
      for (i in 1:length(fasta_files)) {
        kmer_data_list[[i]] <- find_kmer_positions(fasta_files[i], kmer)
      }
      G_positions <- do.call(rbind, kmer_data_list)
      plot_data[[paste0(kmer)]] <- list(
        result = result,
        kmer = kmer,
        kmer_positions = G_positions,
        Kmer_size = nchar(kmer)
      )
    }
  }



  ############# Render #################

  # Render the Rmd file
  rmarkdown::render(
    input = template_path,
    output_file = report_name,
    output_dir = "reports",
    params = list(virus_genome_path=virus_genome_path,
                  vector_genome_path=vector_genome_path,
                  k_list=k_list,
                  x=x,
                  results=results,
                  plot_data=plot_data),
    envir = new.env()
  )

  message("Report generated at: reports/", report_name)
}
