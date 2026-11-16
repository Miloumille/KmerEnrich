test_that("expand_iupac enumerates instantiations", {
  expect_setequal(expand_iupac("ACR"), c("ACA", "ACG"))
  expect_equal(expand_iupac("ACG"), "ACG")
})

test_that("passes_thermo accepts a valid primer in a wide window", {
  expect_true(passes_thermo("GCGCGCGCG", tm_min = -100, tm_max = 100, gc_min = 0, gc_max = 1))
})

test_that("passes_thermo rejects Tm outside the window", {
  expect_false(passes_thermo("GCGCGCGCG", tm_min = 90, tm_max = 100, gc_min = 0, gc_max = 1))
})

test_that("passes_thermo rejects GC content outside the window", {
  expect_false(passes_thermo("GCGCGCGCG", tm_min = -100, tm_max = 100, gc_min = 0, gc_max = 0.5))
})

test_that("passes_thermo requires a GC clamp at the 3' anchor", {
  expect_false(passes_thermo("ATATATATA", tm_min = -100, tm_max = 100, gc_min = 0, gc_max = 1))
})

test_that("passes_thermo evaluates every IUPAC instantiation", {
  expect_true(passes_thermo("GCGCGCGCS", tm_min = -100, tm_max = 100, gc_min = 0, gc_max = 1))
})
