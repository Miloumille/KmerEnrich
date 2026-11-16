test_that("default_background is a uniform independence model", {
  m <- default_background()
  expect_equal(unname(m$pi), rep(0.25, 4))
  expect_equal(unname(rowSums(m$trans)), rep(1, 4))
  expect_equal(expected_freq("AC", m), 0.25^2)  # independent bases
})

test_that("default_background honours custom base frequencies", {
  m <- default_background(c(A = 0.1, C = 0.4, G = 0.4, T = 0.1))
  expect_equal(unname(m$pi[c("C", "T")]), c(0.4, 0.1))
  expect_equal(expected_freq("CC", m), 0.4 * 0.4)
})
