data("mouseData")


test_that("compute distance matrix", {
  aggdat <- aggFeatures(mouseData, level = "order")
  distmat <- computeDistMat(aggdat, dist_method = "bray")
  expect_equal(round(sum(distmat)), 1740)
  distmat <- computeDistMat(aggdat, dist_method = "euclidean")
  expect_equal(round(sum(distmat)), 89262)
})