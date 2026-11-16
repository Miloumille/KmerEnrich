test_that("markov_background returns normalised pi and transition rows", {
  m <- markov_background(c("ACGT", "ACGT", "ACGT"))

  expect_equal(sort(names(m$pi)), c("A", "C", "G", "T"))
  expect_equal(sum(m$pi), 1)
  expect_equal(unname(rowSums(m$trans)), rep(1, 4))
  expect_equal(dim(m$trans), c(4L, 4L))
})

test_that("markov_background reflects observed transitions", {
  m <- markov_background(c("ACACACAC", "ACACACAC"))

  # 'C' almost always follows 'A' here, so P(C|A) should dominate P(A|A)
  expect_gt(m$trans["A", "C"], m$trans["A", "A"])
})

test_that("markov_background_files matches the in-memory fit across chunks", {
  root <- withr::local_tempdir()
  reads <- c("ACGTACGT", "TTTTGGGG", "CACACACA", "GATTACAG")
  writeLines(
    unlist(lapply(seq_along(reads),
                  function(i) c(paste0("@r", i), reads[i], "+", strrep("I", nchar(reads[i]))))),
    file.path(root, "bg.fastq")
  )

  m_stream <- markov_background_files(file.path(root, "bg.fastq"), chunk_size = 2)  # forces 2 chunks
  m_mem <- markov_background(reads)

  expect_equal(m_stream$pi, m_mem$pi)
  expect_equal(m_stream$trans, m_mem$trans)
})

test_that("markov_background_files tolerates a FASTA background", {
  root <- withr::local_tempdir()
  writeLines(c(">h1", "ACGTACGT"), file.path(root, "ref.fasta"))

  m_file <- markov_background_files(file.path(root, "ref.fasta"))
  m_mem <- markov_background("ACGTACGT")

  expect_equal(m_file$pi, m_mem$pi)
})
