data("mouseData")

test_that("normalize data methods", {
  proportionMR <- normalizeData(mouseData, norm_method = "Proportion")
  expect_equal(sum(normFactors(proportionMR)),315593)
  cssMR <- normalizeData(mouseData, norm_method = "CSS")
  expect_equal(sum(normFactors(cssMR)),33914)
})

test_that("filter and normalize", {
  initialNormFactors <- normFactors(mouseData)
  filteredMR <- filterMEData(mouseData, minpresence = 2, minfeats = 10, minreads = 1500)
  filteredNormFactors <- normFactors(normalizeData(filteredMR, norm_method = "CSS"))
  expect_failure(expect_equal(filteredNormFactors, normFactors(filteredMR)))
  expect_equal(sum(filteredNormFactors),26488)
  filteredNormFactors <- normFactors(normalizeData(filteredMR, norm_method = "Proportion"))
  expect_failure(expect_equal(sum(filteredNormFactors),33914))
})

test_that("aggregate features", {
  genusMD <- aggFeatures(mouseData, level = "genus")
  expect_equal(ncol(genusMD), ncol(mouseData))
  expect_failure(expect_equal(nrow(genusMD),nrow(mouseData)))
  expect_equal(as.numeric(nrow(genusMD)), 61)
})