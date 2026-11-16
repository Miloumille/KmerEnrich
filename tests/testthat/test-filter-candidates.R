test_that("filter_candidates keeps enriched, thermo-valid primers only", {
  scored <- data.frame(
    pattern = c("GCGCGCGCG", "CGCGCGCGC", "ATATATATA"),
    enrichment = c(10, 1, 10),
    stringsAsFactors = FALSE
  )

  out <- filter_candidates(scored, min_enrichment = 2,
                           tm_min = -100, tm_max = 100, gc_min = 0, gc_max = 1)

  # GCGCGCGCG: enriched + GC clamp -> kept
  # CGCGCGCGC: enrichment 1 < 2 -> dropped
  # ATATATATA: no GC clamp -> dropped
  expect_equal(out$pattern, "GCGCGCGCG")
})

test_that("filter_candidates returns zero rows when nothing qualifies", {
  scored <- data.frame(pattern = "ATATATATA", enrichment = 10, stringsAsFactors = FALSE)
  out <- filter_candidates(scored, min_enrichment = 2, gc_min = 0, gc_max = 1)
  expect_equal(nrow(out), 0L)
})

test_that("filter_candidates drops rows with NA/NaN enrichment", {
  scored <- data.frame(
    pattern = c("GCGCGCGCG", "GCGCGCGCG"),
    enrichment = c(10, NaN),
    stringsAsFactors = FALSE
  )
  out <- filter_candidates(scored, min_enrichment = 2,
                           tm_min = -100, tm_max = 100, gc_min = 0, gc_max = 1)
  expect_equal(nrow(out), 1L)
  expect_false(anyNA(out$pattern))
})
