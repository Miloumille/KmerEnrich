test_that("gc_content returns the G/C fraction", {
  expect_equal(gc_content("GGCC"), 1)
  expect_equal(gc_content("ATAT"), 0)
  expect_equal(gc_content("ACGT"), 0.5)
})

test_that("primer_tm ranks GC-rich above AT-rich of equal length", {
  expect_gt(primer_tm("GCGCGCGCG"), primer_tm("ATATATATA"))
})

test_that("primer_tm increases with length", {
  expect_gt(primer_tm("GCGCGCGCGCGC"), primer_tm("GCGCGCGCG"))
})

test_that("primer_tm returns a finite number", {
  expect_true(is.finite(primer_tm("ACGTACGTA")))
})
