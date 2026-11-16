test_that("optimize_pool maximises covered bins for a given K", {
  cov <- list(
    matrix = matrix(c(1, 1, 0,  0, 0, 1), 3, 2, dimnames = list(NULL, c("A", "B"))),
    bin_target = c(1L, 1L, 1L), patterns = c("A", "B"), targets = "t1"
  )
  expect_equal(optimize_pool(cov, 1, floor = 0)$patterns, "A")   # A covers 2 of 3 bins
  expect_equal(optimize_pool(cov, 1, floor = 0)$coverage, 2 / 3)
  expect_equal(optimize_pool(cov, 2, floor = 0)$coverage, 1)
})

test_that("optimize_pool enforces the per-target floor", {
  cov <- list(
    matrix = matrix(c(1, 1, 0, 0,  0, 0, 1, 1), 4, 2, dimnames = list(NULL, c("A", "B"))),
    bin_target = c(1L, 1L, 2L, 2L), patterns = c("A", "B"), targets = c("t1", "t2")
  )
  expect_true(optimize_pool(cov, 1, floor = 0.5)$status != 0)     # one primer can't cover both targets
  r2 <- optimize_pool(cov, 2, floor = 0.5)
  expect_equal(r2$coverage, 1)
  expect_setequal(r2$patterns, c("A", "B"))
})

test_that("optimize_pool returns per_target named and ordered by target index", {
  cov <- list(
    matrix = matrix(c(0, 0, 1, 1,  1, 1, 0, 0), 4, 2, dimnames = list(NULL, c("A", "B"))),
    bin_target = c(2L, 2L, 1L, 1L), patterns = c("A", "B"), targets = c("t1", "t2")
  )
  r <- optimize_pool(cov, 2, floor = 0)
  expect_equal(names(r$per_target), c("t1", "t2"))
  expect_equal(unname(r$per_target["t1"]), 1)   # B covers t1's bins (rows 3-4)
  expect_equal(unname(r$per_target["t2"]), 1)   # A covers t2's bins (rows 1-2)
})

test_that("optimize_pool and saturation_curve warn when given too many primers", {
  np <- 260L; nb <- 3L
  set.seed(1)
  m <- matrix(0L, nb, np); m[cbind(sample(nb, np, TRUE), seq_len(np))] <- 1L
  colnames(m) <- paste0("p", seq_len(np))
  cov <- list(matrix = m, bin_target = rep(1L, nb),
              patterns = paste0("p", seq_len(np)), targets = c("v"))
  expect_warning(optimize_pool(cov, 2L, floor = 0), "beyond ~250")
  expect_warning(saturation_curve(cov, k_max = 2L, floor = 0), "intractable beyond ~250")
})
