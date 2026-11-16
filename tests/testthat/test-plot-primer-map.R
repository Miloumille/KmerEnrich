test_that("kmer_palette returns fixed-order distinct colours", {
  expect_equal(kmer_palette(1), "#2a78d6")
  expect_equal(kmer_palette(3), c("#2a78d6", "#1baf7a", "#eda100"))
  expect_equal(length(kmer_palette(8)), 8L)
})

test_that("plot_primer_map maps primers onto one reference-genome lane", {
  seq <- Biostrings::DNAStringSet(c(
    s1 = paste0(strrep("A", 5), "GCAAGTTCA", strrep("A", 10), "CATGCACTG", strrep("A", 5)),
    s2 = strrep("A", 40)  # second strain, ignored: only the first is the reference
  ))
  p <- plot_primer_map(c("GCAAGTTCA", "CATGCACTG"), list(dengue = seq), read_len = 5)

  expect_s3_class(p, "ggplot")
  built <- ggplot2::ggplot_build(p)
  is_seg <- vapply(p$layers, function(l) inherits(l$geom, "GeomSegment"), logical(1))
  seg <- built$data[[which(is_seg)]]

  expect_equal(nrow(seg), 2L)                  # one segment per primer on the reference
  expect_equal(length(unique(seg$colour)), 2L) # coloured by primer
  expect_equal(length(unique(seg$y)), 1L)      # single genome lane
})
