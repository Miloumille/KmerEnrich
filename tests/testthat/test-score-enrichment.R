test_that("observed_freq counts IUPAC matches per start position", {
  seqs <- Biostrings::DNAStringSet(c(s1 = "ACGTACGT"))
  expect_equal(observed_freq("ACG", seqs), 2 / 6)   # 2 hits, 8-3+1=6 positions
})

test_that("score_enrichment ranks a viral-abundant motif above an absent one", {
  target <- Biostrings::DNAStringSet(c(s1 = "ACGACGACGACG"))
  model <- markov_background(c("TTTTTTTT", "TTTTTTTT"))

  res <- score_enrichment(c("ACG", "TTT"), target, model)

  acg <- res$enrichment[res$pattern == "ACG"]
  ttt <- res$enrichment[res$pattern == "TTT"]
  expect_equal(ttt, 0)
  expect_gt(acg, ttt)
  expect_true(is.finite(acg))
})

test_that("observed_freq counts degenerate IUPAC patterns", {
  seqs <- Biostrings::DNAStringSet(c(s1 = "ACAACG"))
  expect_equal(observed_freq("ACR", seqs), 2 / 4)   # ACA (pos 1) + ACG (pos 4), R = A/G
})

test_that("score_enrichment scores a degenerate primer with a finite fold", {
  target <- Biostrings::DNAStringSet(c(s1 = "ACGACAACGACA"))
  model <- markov_background(c("TTTTTTTT", "TTTTTTTT"))

  res <- score_enrichment("ACR", target, model)

  expect_equal(nrow(res), 1L)
  expect_gt(res$enrichment, 0)
  expect_true(is.finite(res$enrichment))
})
