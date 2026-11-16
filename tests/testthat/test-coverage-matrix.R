test_that("coverage_matrix marks the bins a primer covers", {
  targets <- list(t1 = Biostrings::DNAStringSet(c(s1 = "AAACGTTTTT")))
  cm <- coverage_matrix("CGT", targets, read_len = 2, bin_size = 3L)

  # match at start 4 -> covered interval [2, 6] -> bins 1 and 2 (of ceil(10/3)=4)
  expect_equal(dim(cm$matrix), c(4L, 1L))
  expect_equal(as.integer(cm$matrix[, 1]), c(1L, 1L, 0L, 0L))
  expect_equal(cm$bin_target, c(1L, 1L, 1L, 1L))
})

test_that("coverage_matrix stacks targets and tags each bin", {
  targets <- list(
    t1 = Biostrings::DNAStringSet(c(s1 = "CGTAAA")),
    t2 = Biostrings::DNAStringSet(c(s1 = "AAACGT"))
  )
  cm <- coverage_matrix("CGT", targets, read_len = 1, bin_size = 100L)

  expect_equal(nrow(cm$matrix), 2L)              # one bin per target
  expect_equal(cm$bin_target, c(1L, 2L))
  expect_equal(as.integer(cm$matrix[, 1]), c(1L, 1L))  # CGT present in both
})
