data("mouseData")
pData(mouseData) <- data.frame(SAMPLE_ID = rownames(pData(mouseData)),
                           pData(mouseData))
md_counts <- system.file("extdata", "md_counts.tsv", package = "microbiomeExplorer")
md_feats <- system.file("extdata", "md_feats.tsv", package = "microbiomeExplorer")
md_pheno <- system.file("extdata", "md_pheno.tsv", package = "microbiomeExplorer")
pheno_add <- system.file("extdata", "add_pheno.tsv", package = "microbiomeExplorer")
mr_base <- readData(md_counts, type = "TAB")
df <- reshape2::melt(MRcounts(mouseData), varnames = c("OTU", "samname"), 
                     value.name = "Reads")

test_that("adding pheno data", {
  mr_pheno <- addPhenoData(mouseData,pheno_add)
  expect_equal(ncol(pData(mr_pheno)),2)
  expect_equal(nrow(pData(mr_pheno)),139)
  expect_equal(names(pData(mr_pheno)),c("SAMPLE_ID","rand_num"))
})

test_that("extending pheno data", {
  mr_pheno <- extendPhenoData(mouseData,pheno_add)
  expect_equal(ncol(pData(mr_pheno)),7)
  expect_equal(nrow(pData(mr_pheno)),139)
  expect_equal(names(pData(mr_pheno)),c("SAMPLE_ID","mouseID", "date", "diet", "relativeTime", "status", "rand_num"))
})

test_that("reading data", {
  expect_equal(ncol(pData(mr_base)),1)
  expect_equal(nrow(pData(mr_base)),139)
  expect_equal(names(pData(mr_base)),"SAMPLE_ID")
  expect_equal(ncol(fData(mr_base)),0)
  expect_equal(nrow(fData(mr_base)),10172)
  expect_equal(sum(MRcounts(mr_base)),315593)
  expect_equal(row.names(MRcounts(mr_base))[1:5],c("Prevotellaceae:1", "Lachnospiraceae:1", "Unclassified-Screened:1",
                                                   "Clostridiales:1", "Clostridiales:2"))
})

test_that("adding feat data", {
  mr_feats <- addFeatData(mr_base, featdata = md_feats)
  expect_equal(ncol(pData(mr_feats)),1)
  expect_equal(nrow(pData(mr_feats)),139)
  expect_equal(names(pData(mr_feats)),"SAMPLE_ID")
  expect_equal(ncol(fData(mr_feats)),7)
  expect_equal(nrow(fData(mr_feats)),10172)
  expect_equal(sum(MRcounts(mr_feats)),315593)
  expect_equal(names(fData(mr_feats)), c("superkingdom", "phylum", "class", "order", "family", "genus", "OTU" ))
  expect_equal(row.names(fData(mr_feats))[1:5],c("Prevotellaceae:1", "Lachnospiraceae:1", "Unclassified-Screened:1",
                                                   "Clostridiales:1", "Clostridiales:2"))
})

test_that("filtering min presence", {
  expect_equal(dim(filterMEData(mouseData, minpresence = 4))[[1]], 2248)
  expect_equal(dim(filterMEData(mouseData, minpresence = 10))[[1]], 1063)
  expect_equal(dim(filterMEData(mouseData, minpresence = 10))[[2]], 139)
})

test_that("filtering min feats", {
  expect_equal(dim(filterMEData(mouseData, minfeats = 500))[[1]],3309)
  expect_equal(dim(filterMEData(mouseData, minfeats = 500))[[2]],13)
})

test_that("filtering min reads", {
  expect_equal(dim(filterMEData(mouseData, minreads = 1500))[[1]],9628)
  expect_equal(dim(filterMEData(mouseData, minreads = 1500))[[2]],129)
})

test_that("filtering combined", {
  test <- filterMEData(mouseData, minpresence = 2, minfeats = 200, minreads = 1500)
  expect_false(any(rowSums(MRcounts(test) > 0) < 2))
  expect_false(any(colSums(MRcounts(test) > 0) < 200))
  expect_false(any(colSums(MRcounts(test)) < 1500))
  expect_equal(dim(test)[[1]], 4031)
  expect_equal(dim(test)[[2]], 127)
})

test_that("filter by pheno", {
  expect_equal(dim(filterByPheno(mouseData, rm_phenovalues = list("diet" = "Western")))[[1]], 10172)
  expect_equal(dim(filterByPheno(mouseData, rm_phenovalues = list("diet" = "Western")))[[2]], 85)
})





