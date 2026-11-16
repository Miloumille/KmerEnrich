#' Plot the coverage saturation curve
#'
#' Coverage vs pool size K, with the plateau N* marked. Single series, no legend.
#'
#' @param sat A [saturation_curve()] result.
#' @return A ggplot object.
#' @export
plot_saturation <- function(sat) {
  ggplot2::ggplot(sat$curve, ggplot2::aes(.data$k, .data$coverage)) +
    ggplot2::geom_vline(xintercept = sat$n_star, linetype = "dashed", colour = "#52514e") +
    ggplot2::geom_line(linewidth = 1, colour = "#2a78d6") +
    ggplot2::geom_point(size = 2.5, colour = "#2a78d6") +
    ggplot2::annotate("text", x = sat$n_star, y = min(sat$curve$coverage, na.rm = TRUE),
                      label = paste0("N* = ", sat$n_star), hjust = -0.15, colour = "#52514e") +
    ggplot2::scale_y_continuous(labels = function(x) paste0(round(x * 100), "%")) +
    ggplot2::labs(x = "Nombre de primers (K)", y = "Couverture du génome") +
    ggplot2::theme_minimal()
}
