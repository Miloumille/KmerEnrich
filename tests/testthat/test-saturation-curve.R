make_cov <- function() list(
  matrix = matrix(c(1, 1, 0,  0, 0, 1), 3, 2, dimnames = list(NULL, c("A", "B"))),
  bin_target = c(1L, 1L, 1L), patterns = c("A", "B"), targets = "t1"
)

test_that("saturation_curve reports coverage for each K", {
  sc <- saturation_curve(make_cov(), k_max = 2, floor = 0, plateau = 0.05)
  expect_equal(nrow(sc$curve), 2L)
  expect_equal(sc$curve$coverage, c(2 / 3, 1))
})

test_that("saturation_curve stops at the plateau", {
  # gain from K=1 (2/3) to K=2 (1) is 50%
  expect_equal(saturation_curve(make_cov(), 2, floor = 0, plateau = 0.60)$n_star, 1L)  # 50% < 60% -> plateau at 1
  expect_equal(saturation_curve(make_cov(), 2, floor = 0, plateau = 0.05)$n_star, 2L)  # 50% > 5%  -> go to max
})
