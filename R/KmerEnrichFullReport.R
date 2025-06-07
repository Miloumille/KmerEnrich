#' Generate an HTML report
#'
#' @description The `KmerEnrichFullReport` function executes the KmerEnrich pipeline and generates a comprehensive HTML report from an R Markdown template included within the package. This report summarizes the characteristics of identified k-mers within the virus and vector genomes, visualizes k-mer enrichment, and plots k-mer positions across virus genomes.
#'
#'   Documentation has been partially AI generated
#'
#' @param virus_genome_path A string specifying the path to a FASTA file or a folder containing multiple FASTA files representing the virus genome(s). K-mers will be searched within these sequences.
#' @param vector_genome_path A string specifying the path to a FASTA file or a folder containing multiple FASTA files representing the vector genome(s). K-mer frequencies in the vector genome are used for enrichment calculations.
#' @param k_list A numeric vector or list of integers, where each integer represents a k-mer length. The report will process and display results for k-mers of these specified lengths.
#' @param x An integer value indicating the number of first bases (prefix length) in each virus genome within which k-mers will be initially searched for the `common_kmers` step. Only k-mers found in this prefix across *all* virus genomes (if multiple are provided) are considered.
#' @param report_name A string specifying the base name for the generated HTML report file and the report's containing directory. The default is "report".
#'
#' @details This function orchestrates several steps of the KmerEnrich pipeline:
#' \enumerate{
#'   \item \strong{K-mer Identification and Statistics:} For each `k` in `k_list`:
#'     \itemize{
#'       \item `common_kmers`: Identifies k-mers of length `k` that are present in the first `x` bases of all virus genomes specified by `virus_genome_path`.
#'       \item `kmers_pos_df`: Determines the positions of these common k-mers within all virus genomes.
#'       \item `get_virus_stats_kmers`: Calculates statistical properties for these k-mers based on their distribution and counts within the virus genomes.
#'       \item `get_vector_stats_kmers`: Incorporates k-mer frequency data from the `vector_genome_path` to compute enrichment scores. The results for each `k` are stored in the `results` list, accessible as `results[[paste0("k_", k)]]`. Each element of this list is a data frame containing columns like `Kmer`, `mean_count`, `sd_count`, `Distribution_score`, and `Enrichment_score`.
#'     }
#'   \item \strong{K-mer Position Data for Plotting:} For each `k` in `k_list`, the top 5 k-mers (based on `Distribution_score`) from the `results` data frame for that `k` are selected. For each of these top k-mers, their exact genomic positions across *all* virus genomes are identified. This data is structured and stored in the `plot_data` list. Each entry in `plot_data` (keyed by the k-mer sequence, e.g., "CATGG") contains:
#'     \itemize{
#'       \item `result`: The full statistical table (`complete_kmer_stats`) for the k-mer's length.
#'       \item `kmer`: The k-mer sequence itself.
#'       \item `kmer_positions`: A data frame with `Position` and `Genome` columns indicating where the k-mer was found.
#'       \item `Kmer_size`: The length of the k-mer.
#'     }
#'   \item \strong{Report Generation and Output:}
#'     \itemize{
#'       \item A dated subdirectory is created within a "reports" folder (e.g., `reports/report_20250607_183000`) to store all output files.
#'       \item For each `k` in `k_list`, the complete k-mer results (`results[[paste0("k_", k)]]`) are saved as a CSV file (e.g., `kmer_results_k5.csv`) within this directory.
#'       \item An HTML report (`report_name.html`) is rendered into the same dated directory using the internal R Markdown template. This report includes:
#'         \itemize{
#'           \item A table of the top 5 k-mers by `Enrichment_score` for each `k`.
#'           \item A scatter plot visualizing the `Distribution_score` (x-axis: "count per virus / sd segm virus") against the `Enrichment_score` (y-axis: "log10(frequency virus) - log10(frequency mos)") for all k-mers for each `k`. K-mers are connected by dashed grey lines, and points are colored by `accession_id`.
#'           \item Individual plots showing the genomic positions of the top 5 k-mers (from `plot_data`) across different virus genomes. The subtitle for these plots attempts to include an "Experiment" name, which would ideally be a descriptive string related to the analysis run.
#'         }
#'     }
#' }
#'
#' @return This function primarily generates an HTML report and associated CSV files in a newly created dated directory within the working directory. It does not return an R object directly.
#' @export
KmerEnrichFullReport <- function(virus_genome_path, vector_genome_path, k_list, x, report_name = "report") {
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
    sorted_table <- result[order(-result$Distribution_score), ]
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

  ############# Create dated report directory ###############
  # Get current date and time for unique folder name
  current_datetime <- format(Sys.time(), "%Y%m%d_%H%M%S")
  report_dir <- file.path("reports", paste0(report_name," ",current_datetime))

  # Create the directory if it doesn't exist
  if (!dir.exists(report_dir)) {
    dir.create(report_dir, recursive = TRUE)
  }

  ############# Save results to CSV files ###############
  for (k_size in k_list) {
    kmer_data <- results[[paste0("k_", k_size)]]
    if (!is.null(kmer_data) && is.data.frame(kmer_data)) {
      csv_file_name <- paste0("kmer_results_k", k_size, ".csv")
      write.csv(kmer_data, file = file.path(report_dir, csv_file_name), row.names = FALSE)
      message("Saved Kmer results for k=", k_size, " to: ", file.path(report_dir, csv_file_name))
    } else {
      warning("No data found for k-mer size: ", k_size, ". Skipping CSV generation.")
    }
  }

  ############# Render HTML report #################
  # Render the Rmd file into the newly created directory
  rmarkdown::render(
    input = template_path,
    output_file = paste0(report_name,".html"),
    output_dir = report_dir,
    params = list(
      virus_genome_path = virus_genome_path,
      vector_genome_path = vector_genome_path,
      k_list = k_list,
      x = x,
      results = results,
      plot_data = plot_data
    ),
    envir = new.env()
  )

  message("Report and CSV files generated in: ", report_dir)
}


