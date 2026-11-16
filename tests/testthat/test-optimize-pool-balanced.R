test_that("optimize_pool_balanced balances coverage across targets", {
  m <- matrix(c(1, 1, 0, 0,   0, 0, 1, 1,   1, 0, 0, 0), nrow = 4)
  colnames(m) <- c("p1", "p2", "p3")
  cov <- list(matrix = m, bin_target = c(1, 1, 2, 2),
              patterns = c("p1", "p2", "p3"), targets = c("vA", "vB"))

  r <- optimize_pool_balanced(cov, 2)

  expect_setequal(r$patterns, c("p1", "p2"))   # one primer per virus, not two on one
  expect_equal(unname(r$per_target), c(1, 1))
  expect_equal(r$min_cov, 1)
})

test_that("optimize_pool_balanced always returns a pool (no infeasibility)", {
  m <- matrix(c(1, 1, 0, 0,   0, 0, 1, 1), nrow = 4)
  colnames(m) <- c("p1", "p2")
  cov <- list(matrix = m, bin_target = c(1, 1, 2, 2),
              patterns = c("p1", "p2"), targets = c("vA", "vB"))

  r <- optimize_pool_balanced(cov, 1)          # 1 primer can't cover both -> min_cov 0, still a pool
  expect_length(r$patterns, 1L)
  expect_equal(r$status, 0)
  expect_equal(r$min_cov, 0)
})
