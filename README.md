---

# `KmerEnrich` R Package

---

## Overview

The `KmerEnrich` R package provides a robust pipeline for analyzing k-mer enrichment in virus genomes compared to a vector (e.g., mosquito) genome. It helps identify significant k-mers, visualize their distribution within viral sequences, and compare their frequencies between viral and vector hosts. The package culminates in a comprehensive HTML report summarizing all findings.

Parts of this documentation has been AI generated.

## Features

* **Common K-mer Identification**: Find k-mers shared across multiple virus sequences or unique to a single one.
* **K-mer Position Mapping**: Map the exact genomic locations of k-mers and calculate segment sizes between them.
* **Virus K-mer Statistics**: Compute statistical metrics for k-mers within virus genomes, including counts, frequencies, and segment distribution.
* **Vector K-mer Statistics**: Calculate k-mer counts and frequencies in vector RNA-seq data, and derive enrichment scores.
* **Automated HTML Reporting**: Generate a detailed, dynamic HTML report with tables and plots for easy interpretation of results.

---

## Installation

You can install the `KmerEnrich` package directly from GitHub using the `devtools` package. If you don't have `devtools` installed, you'll need to install it first.

```R
# Install devtools if you haven't already
install.packages("devtools")

# Install Biostrings (a dependency)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Biostrings")
BiocManager::install("ShortRead") # Also a dependency

# Install stringr and data.table (other dependencies)
install.packages("stringr")
install.packages("data.table") # For rbindlist
install.packages("dplyr")     # For data manipulation
install.packages("ggplot2")    # For plotting
install.packages("rmarkdown")  # For report generation
install.packages("knitr")      # For tables in report

# Install KmerEnrich from GitHub
devtools::install_github("your-github-username/KmerEnrich") # Replace "your-github-username" with your actual GitHub username or organization
```

---

## Usage

The primary function to use this package is `KmerEnrichFullReport`, which orchestrates the entire analysis pipeline and generates the final HTML report.

### Core Function: `KmerEnrichFullReport`

This function runs the full pipeline, from k-mer identification to generating a detailed HTML report and CSV files of the results.

```R
KmerEnrichFullReport(
  virus_genome_path,
  vector_genome_path,
  k_list,
  x,
  report_name = "report"
)
```

#### Arguments

* **`virus_genome_path`**: `character`. Path to a FASTA file or a folder containing multiple FASTA files representing the virus genome(s). K-mers will be searched within these sequences.
* **`vector_genome_path`**: `character`. Path to a FASTA file or a folder containing multiple FASTA files representing the vector genome(s) (e.g., mosquito RNA-seq data). K-mer frequencies in the vector genome are used for enrichment calculations. For RNA-seq data, this should be the base path for paired-end FASTQ files (e.g., `'path/to/sampleA'` for `sampleA_1.fastq` and `sampleA_2.fastq`).
* **`k_list`**: `numeric vector` or `list`. Integers representing the k-mer lengths to be analyzed (e.g., `c(3, 4, 5)`).
* **`x`**: `integer`. The number of first bases (prefix length) in each virus genome within which k-mers will be initially searched for the `common_kmers` step. Only k-mers found in this prefix across *all* virus genomes (if multiple are provided) are considered.
* **`report_name`**: `character`. The base name for the generated HTML report file and its containing dated directory (e.g., setting to `"MyAnalysis"` will create `reports/MyAnalysis_YYYYMMDD_HHMMSS/MyAnalysis.html`). Default is `"report"`.

#### Output

The function generates:
* A dated subdirectory within a `reports/` folder (e.g., `reports/report_20250607_183000`).
* CSV files (e.g., `kmer_results_k5.csv`) containing detailed k-mer statistics for each `k` in `k_list`, saved within the dated directory.
* An HTML report (`report_name.html`) within the dated directory, providing:
    * Tables of top k-mers by enrichment score.
    * Scatter plots showing k-mer distribution vs. enrichment.
    * Plots illustrating the genomic positions of top k-mers within virus genomes.

---

### Structured image of the Pipeline

<img src="inst/img/Pipeline_image.png" alt="Logo" width="800" />

---

### Underlying Functions (for advanced users)

While `KmerEnrichFullReport` streamlines the process, you can also use the individual functions that make up the pipeline for more granular control or custom analyses.

#### `common_kmers`

Generates a list of k-mers commonly observed across all sequences in specified FASTA files (or unique k-mers within a single file's sequences) up to a given position `x`.

```R
common_kmers(virus_genome_path, k, x)
```

* **`virus_genome_path`**: Path to FASTA file(s).
* **`k`**: Length of k-mers to search for.
* **`x`**: Position until which k-mers should be searched (prefix length).

**Returns**: A character vector of common k-mers.

#### `kmers_pos_df`

Generates a data frame detailing the absolute positions of specified k-mers within each sequence of input FASTA file(s), including segment sizes between k-mers.

```R
kmers_pos_df(virus_genome_path, kmers_list)
```

* **`virus_genome_path`**: Path to FASTA file(s).
* **`kmers_list`**: A character vector of k-mer sequences to search for.

**Returns**: A data frame with `SequenceName`, `Kmer`, `Position`, `Type`, `segm_size`, and `log10_segm_size` columns.

#### `get_virus_stats_kmers`

Calculates various statistical metrics for each k-mer based on their positions and counts within virus genomes.

```R
get_virus_stats_kmers(df_virus)
```

* **`df_virus`**: A data frame (typically from `kmers_pos_df`) containing `Kmer`, `SequenceName`, `Position`, `Type`, `segm_size`, and `log10_segm_size` columns.

**Returns**: A data frame summarized by `Kmer`, including `mean_mean_log10_segm_size`, `max_max_log10_segm_size`, `mean_sd_log10_segm_size`, `mean_count`, `sd_count`, and `mean_freq_virus`.

#### `get_vector_stats_kmers`

Calculates k-mer counts and frequencies within RNA-seq data from paired-end FASTQ files, enriching a provided k-mer statistics data frame with enrichment and distribution scores.

```R
get_vector_stats_kmers(v_kmers_stats, vector_genome_path)
```

* **`v_kmers_stats`**: A data frame containing k-mer statistics (must have `Kmer`, `mean_freq_virus`, `mean_count`, and `max_max_log10_segm_size` columns).
* **`vector_genome_path`**: A character vector where each element is the base path for a pair of paired-end FASTQ files (e.g., `'path/to/sampleA'`).

**Returns**: An enriched data frame including `kmer_count`, `freq_vector`, `accession_id`, `Enrichment_score`, and `Distribution_score`.

---

## Examples

Here's how you might set up and run a basic analysis using `KmerEnrichFullReport`.

```R
# 1. Create dummy virus FASTA files for demonstration
#    (These will be created in a temporary directory and removed afterwards)
temp_virus_dir <- file.path(tempdir(), "virus_genomes")
dir.create(temp_virus_dir, recursive = TRUE, showWarnings = FALSE)

writeLines(c(">VirusA_Seq1", "ATGCATGCATGCATGC",
             ">VirusA_Seq2", "GCATGCATGCATGCAT"),
           file.path(temp_virus_dir, "virusA.fasta"))

writeLines(c(">VirusB_Seq1", "CCGATGCATGCATGCC",
             ">VirusB_Seq2", "TGCATGCATGCATGCA"),
           file.path(temp_virus_dir, "virusB.fasta"))

# 2. Create dummy vector (RNA-seq) FASTQ files
temp_vector_dir <- file.path(tempdir(), "vector_rnas")
dir.create(temp_vector_dir, recursive = TRUE, showWarnings = FALSE)

# For ShortRead::readFastq to work with empty files, ensure content
writeLines(c("@read1", "ATGCATGCATGC", "+", "###########",
             "@read2", "GCATGCATGCAT", "+", "###########"),
           file.path(temp_vector_dir, "mosquito_sample_1.fastq"))
writeLines(c("@read1", "ATGCATGCATGC", "+", "###########",
             "@read2", "TGCATGCATGCAT", "+", "###########"),
           file.path(temp_vector_dir, "mosquito_sample_2.fastq"))

# Set the base path for vector RNA-seq files
vector_base_path <- file.path(temp_vector_dir, "mosquito_sample")

# 3. Define analysis parameters
k_lengths <- c(3, 4) # Analyze 3-mers and 4-mers
prefix_length <- 15 # Search for common kmers in the first 15 bases

# 4. Run the full KmerEnrich pipeline
#    This will create a dated folder inside a 'reports' directory
#    e.g., 'reports/my_enrichment_report_20250607_183000/'
KmerEnrichFullReport(
  virus_genome_path = temp_virus_dir,
  vector_genome_path = vector_base_path,
  k_list = k_lengths,
  x = prefix_length,
  report_name = "my_enrichment_report"
)

# 5. Clean up dummy files and directories (important for reproducibility)
unlink(temp_virus_dir, recursive = TRUE)
unlink(temp_vector_dir, recursive = TRUE)
# Note: The 'reports' directory will remain with the generated report and CSVs.
```

---

## Contributing

We welcome contributions! If you have suggestions for improvements, bug fixes, or new features, please open an issue or submit a pull request on the GitHub repository.

---

## License

This project is licensed under the MIT License - see the `LICENSE` file for details.

---

## Contact

For questions or support, please open an issue on the GitHub repository.

---
