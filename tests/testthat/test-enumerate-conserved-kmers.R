test_that("enumerate_conserved_kmers keeps only k-mers present in all strains", {
  seqs <- Biostrings::DNAStringSet(c(s1 = "ACGTACGT", s2 = "ACGTACGA"))

  kmers <- enumerate_conserved_kmers(seqs, k = 4)

  expect_true("ACGT" %in% kmers)   # present in both strains
  expect_false("ACGA" %in% kmers)  # only in s2
})

test_that("enumerate_conserved_kmers respects the optional 5' window", {
  seqs <- Biostrings::DNAStringSet(c(s1 = "ACGTACGT", s2 = "ACGTTTTT"))

  kmers <- enumerate_conserved_kmers(seqs, k = 4, window = 4)

  expect_equal(kmers, "ACGT")
})

test_that("enumerate_conserved_kmers returns empty when nothing is shared", {
  seqs <- Biostrings::DNAStringSet(c(s1 = "AAAAAAAA", s2 = "CCCCCCCC"))

  expect_length(enumerate_conserved_kmers(seqs, k = 4), 0)
})
