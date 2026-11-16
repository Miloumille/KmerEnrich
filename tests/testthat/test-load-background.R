test_that("load_background collects reads per source and pools pairs", {
  root <- withr::local_tempdir()
  dir.create(file.path(root, "aedes"))
  writeLines(c("@r1", "ACGTACGT", "+", "IIIIIIII",
               "@r2", "TTTTGGGG", "+", "IIIIIIII"),
             file.path(root, "aedes", "sample_1.fastq"))
  writeLines(c("@r3", "CCCCAAAA", "+", "IIIIIIII"),
             file.path(root, "aedes", "sample_2.fastq"))

  bg <- load_background(root)

  expect_named(bg, "aedes")
  expect_length(bg$aedes, 3)
  expect_true(all(nchar(bg$aedes) == 8))
})

test_that("load_background tolerates FASTA background", {
  root <- withr::local_tempdir()
  dir.create(file.path(root, "host"))
  writeLines(c(">h1", "GGGGCCCC"), file.path(root, "host", "ref.fasta"))

  bg <- load_background(root)

  expect_named(bg, "host")
  expect_equal(bg$host, "GGGGCCCC")
})

test_that("load_background errors on empty input", {
  root <- withr::local_tempdir()
  expect_error(load_background(root), "No FASTQ/FASTA")
})
