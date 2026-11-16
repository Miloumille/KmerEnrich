test_that("kmer_report writes a self-contained HTML file with the data injected", {
  seq <- Biostrings::DNAStringSet(c(
    s1 = paste0(strrep("A", 5), "GCAAGTTCA", strrep("A", 10), "CATGCACTG", strrep("A", 5))
  ))
  targets <- list(dengue = seq)
  admit <- data.frame(pattern = c("GCAAGTTCA", "CATGCACTG"), enrichment = c(5, 3))
  cov <- coverage_matrix(admit$pattern, targets, read_len = 8, bin_size = 5L)

  f <- tempfile(fileext = ".html")
  out <- kmer_report(targets, admit, cov = cov, read_len = 8, file = f, min_sites = 1L, k_max = 2L)

  expect_true(file.exists(out))
  html <- paste(readLines(out), collapse = "\n")
  expect_match(html, "KMerEnrich")
  expect_match(html, "GCAAGTTCA")                              # data injected
  expect_false(grepl("__DATA_JSON__", html, fixed = TRUE))     # placeholder replaced
})

test_that("kmer_report drops single-site primers by default", {
  seq <- Biostrings::DNAStringSet(c(
    s1 = paste0(strrep("A", 5), "GCAAGTTCA", strrep("A", 10), "GCAAGTTCA", strrep("A", 10), "CATGCACTG", strrep("A", 5))
  ))
  targets <- list(dengue = seq)                                # GCAAGTTCA binds twice, CATGCACTG once
  admit <- data.frame(pattern = c("GCAAGTTCA", "CATGCACTG"), enrichment = c(5, 3))
  cov <- coverage_matrix(admit$pattern, targets, read_len = 8, bin_size = 5L)

  html <- paste(readLines(kmer_report(targets, admit, cov = cov, read_len = 8, k_max = 2L,
                                      file = tempfile(fileext = ".html"))), collapse = "\n")
  expect_match(html, "GCAAGTTCA")        # two sites -> kept
  expect_no_match(html, "CATGCACTG")     # single site -> dropped by default
})
