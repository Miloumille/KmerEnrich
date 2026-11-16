test_that("iupac_code maps base sets to ambiguity letters", {
  expect_equal(iupac_code(c("A", "A")), "A")
  expect_equal(iupac_code(c("A", "G")), "R")
  expect_equal(iupac_code(c("G", "C")), "S")
  expect_equal(iupac_code(c("A", "C", "G", "T")), "N")
  expect_equal(iupac_code(c("C", "G", "T")), "B")
})
