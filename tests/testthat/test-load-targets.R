test_that("load_targets groups strains by subfolder", {
  root <- withr::local_tempdir()
  dir.create(file.path(root, "dengue2"))
  dir.create(file.path(root, "west_nile"))
  writeLines(c(">s1", "ACGTACGTACGT"), file.path(root, "dengue2", "a.fasta"))
  writeLines(c(">s2", "ACGTACGTAAAA"), file.path(root, "dengue2", "b.fasta"))
  writeLines(c(">s3", "TTTTACGTACGT"), file.path(root, "west_nile", "c.fasta"))

  targets <- load_targets(root)

  expect_setequal(names(targets), c("dengue2", "west_nile"))
  expect_length(targets$dengue2, 2)
  expect_length(targets$west_nile, 1)
})

test_that("load_targets treats bare FASTA as a single target", {
  root <- withr::local_tempdir()
  writeLines(c(">s1", "ACGTACGT"), file.path(root, "only.fasta"))

  targets <- load_targets(root)

  expect_named(targets, "target")
  expect_length(targets$target, 1)
})

test_that("load_targets errors on empty input", {
  root <- withr::local_tempdir()
  expect_error(load_targets(root), "No FASTA")
})
