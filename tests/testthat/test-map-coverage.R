test_that("map_coverage returns 5' covered intervals for a match", {
  seqs <- Biostrings::DNAStringSet(c(s1 = "AAACGTTTTT"))
  cov <- map_coverage("CGT", seqs, read_len = 2, direction = "5prime")

  expect_equal(cov$start, 4L)
  expect_equal(cov$cov_start, 2L)   # start - read_len
  expect_equal(cov$cov_end, 6L)     # start + k - 1
})

test_that("map_coverage matches IUPAC patterns against the genome", {
  seqs <- Biostrings::DNAStringSet(c(s1 = "AAACGTTTTT"))
  cov <- map_coverage("CST", seqs, read_len = 2)  # S = C or G -> CGT matches

  expect_equal(nrow(cov), 1L)
  expect_equal(cov$start, 4L)
})

test_that("map_coverage clamps coverage to the genome bounds", {
  seqs <- Biostrings::DNAStringSet(c(s1 = "CGTAAAAAAA"))
  cov <- map_coverage("CGT", seqs, read_len = 5, direction = "5prime")

  expect_equal(cov$cov_start, 1L)   # clamped, not -1
})

test_that("map_coverage returns no rows when nothing matches", {
  seqs <- Biostrings::DNAStringSet(c(s1 = "AAAAAAAAAA"))
  expect_equal(nrow(map_coverage("CGT", seqs, read_len = 2)), 0L)
})
