test_that("collapse_to_iupac merges a single mid-position variant", {
  expect_equal(collapse_to_iupac(c("ACGAT", "ACGGT"), max_degenerate = 1L),
               "ACGRT")
})

test_that("collapse_to_iupac never wobbles the 3' anchor", {
  out <- collapse_to_iupac(c("ACGTT", "GCGTT"), anchor = 2L, max_degenerate = 1L)
  expect_setequal(out, c("ACGTT", "GCGTT"))
})

test_that("collapse_to_iupac respects the degenerate budget", {
  out <- collapse_to_iupac(c("ACGAT", "ACGGT", "ACTAT"), max_degenerate = 1L)
  expect_setequal(out, c("ACGRT", "ACTAT"))
})

test_that("collapse_to_iupac returns empty on empty input", {
  expect_length(collapse_to_iupac(character(0)), 0)
})
