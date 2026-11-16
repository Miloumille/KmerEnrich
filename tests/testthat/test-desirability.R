test_that("desirability is a weighted geometric mean scaled to 100", {
  # enrichment 10 -> d_e = 10/20 = 0.5; coverage 1; regularity 1
  # score = (0.5^0.4 * 1^0.3 * 1^0.3) * 100
  expect_equal(round(desirability(10, 1, 1), 2), 75.79)
})

test_that("desirability collapses to zero if any metric is zero", {
  expect_equal(desirability(10, 0, 1), 0)   # coverage 0 sinks the score
})

test_that("desirability rises with coverage", {
  expect_gt(desirability(10, 0.8, 1), desirability(10, 0.4, 1))
})
