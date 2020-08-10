data("mouseData")
aggdat <- aggFeatures(mouseData, level = "genus")
aggmat <- MRcounts(aggdat, norm = TRUE)
aggmat <- log2(aggmat + 1)
feat1 <- "Bacteroides"
feat2 <- "Prevotella"
x <- aggmat[which(rownames(aggmat) == feat1), ]
y <- aggmat[which(rownames(aggmat) == feat2), ]
df <- data.frame(colnames(aggmat), x, y)
colnames(df) <- c("samname", "feat1", "feat2")

test_that("compute confidence interval - spearman", {
  mS <- suppressWarnings(stats::cor.test(df$feat2,df$feat1, method = "spearman"))
  ci_interval <- computeCI_Interval(num = nrow(df), mS = mS, method = "spearman")
  expect_equal(as.character(ci_interval),c("0.171106874573396", "0.46912194979915"))
})

test_that("compute confidence interval - pearson", {
  mS <- suppressWarnings(stats::cor.test(df$feat2,df$feat1, method = "pearson"))
  ci_interval <- computeCI_Interval(num = nrow(df), mS = mS, method = "pearson")
  expect_equal(as.character(ci_interval),c("0.265416496698058", "0.542758284278827"))
})

test_that("compute confidence interval - kendall", {
  mS <- suppressWarnings(stats::cor.test(df$feat2,df$feat1, method = "kendall"))
  ci_interval <- computeCI_Interval(num = nrow(df), mS = mS, method = "kendall")
  expect_equal(as.character(ci_interval),c("0.0930203803066753", "0.316044546486128"))
})

test_that("check error for small numbers", {
  mS <- suppressWarnings(stats::cor.test(df$feat2,df$feat1, method = "spearman"))
  ci_interval <- computeCI_Interval(num = 4, mS = mS, method = "spearman")
  expect_true(is.na(ci_interval[["lower"]]))
  expect_true(is.na(ci_interval[["upper"]]))
})
