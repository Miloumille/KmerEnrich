kmer_palette <- function(n) {
  pal <- c("#2a78d6", "#1baf7a", "#eda100", "#008300",
           "#4a3aa7", "#e34948", "#e87ba4", "#eb6834")
  if (n <= length(pal)) pal[seq_len(n)] else grDevices::colorRampPalette(pal)(n)
}

#' Plot where primers bind along a viral genome
#'
#' Maps the primers onto the target's first strain (reference genome), drawn as a
#' single horizontal bar: each primer is a coloured covered-region segment plus a
#' binding-site marker, labelled with its sequence.
#'
#' @param patterns Primer patterns (IUPAC allowed).
#' @param targets Named list of [Biostrings::DNAStringSet].
#' @param read_len Effective read length (passed to [map_coverage()]).
#' @param target Name of the target to plot (default: first).
#' @return A ggplot object.
#' @export
plot_primer_map <- function(patterns, targets, read_len, target = names(targets)[1]) {
  ref <- targets[[target]][1]
  glen <- as.integer(Biostrings::width(ref)[1])
  cov <- map_coverage(patterns, ref, read_len)
  cov$pattern <- factor(cov$pattern, levels = patterns)
  cov$mid <- (cov$cov_start + cov$cov_end) / 2

  ggplot2::ggplot(cov) +
    ggplot2::annotate("rect", xmin = 1, xmax = glen, ymin = 0.9, ymax = 1.1,
                      fill = "grey90", colour = "grey75", linewidth = 0.3) +
    ggplot2::geom_segment(
      ggplot2::aes(x = .data$cov_start, xend = .data$cov_end, y = 1, yend = 1,
                   colour = .data$pattern),
      linewidth = 5, alpha = 0.9, lineend = "butt") +
    ggplot2::geom_point(ggplot2::aes(x = .data$start, y = 1), colour = "grey15", size = 1.7) +
    ggplot2::geom_text(ggplot2::aes(x = .data$mid, label = .data$pattern),
                       y = 1.22, size = 2.6, colour = "grey25") +
    ggplot2::scale_colour_manual(values = kmer_palette(length(patterns)), drop = FALSE) +
    ggplot2::scale_y_continuous(limits = c(0.7, 1.45), breaks = NULL) +
    ggplot2::labs(x = paste0("Position sur ", target, " (bp)"), y = NULL, colour = "Primer") +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), legend.position = "right")
}
