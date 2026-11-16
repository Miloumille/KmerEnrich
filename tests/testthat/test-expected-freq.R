uniform <- list(
  pi = c(A = .25, C = .25, G = .25, T = .25),
  trans = matrix(.25, 4, 4, dimnames = list(c("A", "C", "G", "T"), c("A", "C", "G", "T")))
)

test_that("iupac_bases expands ambiguity letters", {
  expect_equal(iupac_bases("A"), "A")
  expect_equal(iupac_bases("R"), c("A", "G"))
  expect_equal(iupac_bases("N"), c("A", "C", "G", "T"))
})

test_that("expected_freq matches hand-computed values under a uniform model", {
  expect_equal(expected_freq("A", uniform), 0.25)
  expect_equal(expected_freq("AC", uniform), 0.0625)   # pi[A] * T[A,C]
  expect_equal(expected_freq("N", uniform), 1)         # any single base
  expect_equal(expected_freq("R", uniform), 0.5)       # A or G
})

test_that("expected_freq honours a directional model", {
  m <- list(
    pi = c(A = 1, C = 0, G = 0, T = 0),
    trans = matrix(c(0, 1, 0, 0,  .25, .25, .25, .25,  .25, .25, .25, .25,  .25, .25, .25, .25),
                   4, 4, byrow = TRUE,
                   dimnames = list(c("A", "C", "G", "T"), c("A", "C", "G", "T")))
  )
  expect_equal(expected_freq("AC", m), 1)   # pi[A]=1, T[A,C]=1
  expect_equal(expected_freq("AG", m), 0)   # T[A,G]=0
})
