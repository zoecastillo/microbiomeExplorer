data("mouseData")
aggdat <- aggFeatures(mouseData, level = "genus")

test_that("run differential test - limma", {
  limMD <- runDiffTest(aggdat = aggdat,
              level = "genus",
              phenotype = "diet",
              method = "limma")
  expect_equal(names(limMD),c("genus", "Western-BK", "AveExpr", "t", "P.Value", "adj.P.Val", "B"))
  limMD <- limMD[order(limMD$`Western-BK`),]
  expect_equal(as.character(limMD[1,"genus"]),"Prevotella")
  expect_equal(as.character(limMD[nrow(limMD),"genus"]),"Enterococcus")
})

test_that("run differential test - Kruskal-Wallis", {
  kwMD <- runDiffTest(aggdat = aggdat,
                       level = "genus",
                       phenotype = "diet",
                       method = "Kruskal-Wallis")
  expect_equal(names(kwMD),c("genus", "Western-BK", "AveExpr", "KW-Statistic", "P.Value", "adj.P.Val", "B"))
  kwMD <- kwMD[order(kwMD$`Western-BK`),]
  expect_equal(as.character(kwMD[1,"genus"]),"Prevotella")
  expect_equal(as.character(kwMD[nrow(kwMD),"genus"]),"Enterococcus")
  expect_equal(as.character(kwMD[1,"KW-Statistic"]),"92.951")
})

#ZILN IS NOT AVAILABLE ATM DUE TO ISSUE WITH LIMMA UPDATE IN METAGENOMESEQ
# test_that("run differential test - ZILN", {
#   expect_warning(runDiffTest(aggdat = aggdat,
#                              level = "genus",
#                              phenotype = "diet",
#                              method = "ZILN"))
#   zilnMD <- suppressWarnings(runDiffTest(aggdat = aggdat,
#                       level = "genus",
#                       phenotype = "diet",
#                       method = "ZILN"))
#   expect_equal(names(zilnMD),c("genus", "+samples in group 0", "+samples in group 1", "counts in group 0", "counts in group 1",
#                                "logFC", "se", "pvalues", "adjPvalues"))
#   zilnMD <- zilnMD[order(zilnMD$logFC),]
#   expect_equal(as.character(zilnMD[1,"genus"]),"Prevotella")
#   expect_equal(as.character(zilnMD[2,"genus"]),"Ruminococcus")
#   expect_equal(as.character(zilnMD[nrow(zilnMD),"genus"]),"Enterococcus")
#   expect_equal(as.character(zilnMD[1,"logFC"]),"-4.916")
# })

test_that("run differential test - DeSeq2", {
  skip_on_cran()
  deseq <- suppressWarnings(runDiffTest(aggdat = aggdat,
                      level = "genus",
                      phenotype = "diet",
                      method = "DESeq2"))
  expect_equal(names(deseq),c("genus", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"))
  deseq <- deseq[order(deseq$log2FoldChange),]
  expect_equal(as.character(deseq[1,"genus"]),"Prevotella")
  expect_equal(as.character(deseq[2,"genus"]),"Dorea")
  expect_equal(as.character(deseq[nrow(deseq),"genus"]),"Enterococcus")
  expect_equal(as.character(deseq[3,"log2FoldChange"]),"-4.664")
})

test_that("run differential test - phenolevels", {
  leveledDiff <- suppressWarnings(runDiffTest(aggdat = aggdat,
                                        level = "genus",
                                        phenotype = "mouseID",
                                        phenolevels = c("PM2","PM5"),
                                        method = "limma"))
  leveledDiff <- leveledDiff[order(leveledDiff$`PM5-PM2`),]
  expect_equal(as.character(leveledDiff[1,"genus"]),"Prevotella")
  expect_equal(as.character(leveledDiff[2,"genus"]),"Faecalibacterium")
  expect_equal(as.character(leveledDiff[nrow(leveledDiff),"genus"]),"Enterococcus")
  expect_equal(as.character(leveledDiff[5,"PM5-PM2"]),"-2.467")
})
