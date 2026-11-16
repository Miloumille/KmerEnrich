test_that("plot_saturation builds a ggplot marking N*", {
  sat <- list(
    curve = data.frame(k = 1:3, coverage = c(0.5, 0.8, 0.82),
                       status = 0L, marginal_gain = c(NA, 0.3, 0.02)),
    n_star = 2L
  )
  p <- plot_saturation(sat)

  expect_s3_class(p, "ggplot")
  expect_silent(ggplot2::ggplot_build(p))
  expect_true(any(vapply(p$layers, function(l) inherits(l$geom, "GeomVline"), logical(1))))
})
