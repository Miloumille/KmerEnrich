test_that("package loads and exposes its namespace", {
  expect_true("KMerEnrich" %in% loadedNamespaces())
})
