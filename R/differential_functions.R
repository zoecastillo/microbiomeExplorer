## DIFFERENTIAL ANALYSIS


#' Produce design matrix of pairwise comparisons
#'
#' This function takes in the levels of a factor phenotype and outputs a
#' design matrix of all pairwise comparisons.
#'
#' @param levels Character vector of the levels of a factor phenotype
#'
#' @return A model matrix
designPairs <- function(levels) {
  n <- length(levels)
  design <- matrix(0, n, choose(n, 2))
  rownames(design) <- levels
  colnames(design) <- seq_len(choose(n, 2))
  k <- 0
  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      k <- k + 1
      design[j, k] <- 1
      design[i, k] <- -1
      colnames(design)[k] <- paste(levels[j], "-", levels[i], sep = "")
    }
  }
  design
}


#' Performs differential abundance testing
#'
#' This function performs differential abundance testing between groups of a
#' specified phenotype. Four methods are available: limma, Kruskal-Wallis, 
#' ZILN and DESeq2 (see details).
#'
#' limma is a differential expression tool for microarray data using
#' linear models. It can also be applied to microbiome data.
#'
#' The Kruskal-Wallis test is a non-parametric rank test which examines if 
#' groups come from the same distribution. A significant result indicates at 
#' least one group is distributionally different than another group.
#'
#' ZILN is a zero-inflated log-normal model implemented in
#' \code{\link[metagenomeSeq]{fitFeatureModel}()} of the \code{metagenomeSeq}
#' package.
#' 
#' DeSeq2 performs differential gene expression analysis based on the negative 
#' binomial distribution
#'
#' @param aggdat aggregated MRExperiment
#' @param level Feature level.
#' @param phenotype Phenotype to test.
#' @param phenolevels levels of the phenotype to restrict the comparison to
#' @param log Log2 transform data. Default is TRUE.
#' @param method Differential testing method. One of "limma" (default),
#' "Kruskal-Wallis", "ZILN", or "DESeq2".
#' @param coef Numeric which indicates which pairwise comparison to analyze
#' when there are more than two groups. Corresponds to the column number of the
#' model matrix produced by \code{\link{designPairs}()}. If NULL, a test of any
#' difference between all groups is performed.
#'
#' @importFrom metagenomeSeq fitFeatureModel MRtable MRcounts
#' @importFrom Biobase pData
#' 
#' @return data.frame holding results of the differential analysis
#' 
#' @examples 
#' data("mouseData", package = "metagenomeSeq")
#' aggdat <- aggFeatures(mouseData, level = "genus")
#' runDiffTest(aggdat = aggdat,level = "genus", 
#'             phenotype = "diet", method = "Kruskal-Wallis")
#'
#' @export
runDiffTest <- function(
  aggdat, level, phenotype, phenolevels = NULL,
                        log = TRUE, coef = NULL,
                        method = c("limma", 
                                   "Kruskal-Wallis", 
                                   "ZILN", 
                                   "DESeq2")) {
  phenoTable <- pData(aggdat)
  method <- match.arg(method)
  
  phenoTypeDF <- phenoTable %>% dplyr::select(phenotype)
  pd <- forcats::fct_explicit_na(
    factor(phenoTypeDF[,phenotype]),
    na_level = "NA")
  chosenlevels <- pd
  na <- integer(0)
  if (!is.null(phenolevels)) {
    if (length(phenolevels) < 2) {
      stop("Only one phenolevel specified.")
    }
    chosenlevels <- pd %in% phenolevels
    na <- which(is.na(pd) | !chosenlevels)
    pd <- factor(pd[pd %in% phenolevels])
  } else if (length(levels(pd)) < 2) {
    stop("Phenotype only contains one group.")
  }
  
  if(length(levels(pd)) > 2)
    stop("Please restrict to two phenolevels for differential testing")
  
  norm <- (method != "DESeq2")
  aggmat <- MRcounts(aggdat, norm = norm)
  
  if (length(na) > 0) {
    aggmat <- aggmat[, -na]
    phenoTypeDF <- phenoTypeDF %>%
      dplyr::filter(!dplyr::row_number() %in% na)
  }
  
  if (log && !(method %in% "DESeq2")) {
    aggmat <- log2(aggmat + 1)
  }
  
  design <- stats::model.matrix(~ 0 + pd)
  colnames(design) <- levels(pd)
  
  if(method == "DESeq2"){
    mod <- stats::as.formula(paste("~",phenotype,sep = ""))
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = aggmat, 
                                          colData = phenoTypeDF, 
                                          design = mod)
    dds <- DESeq2::estimateSizeFactors(object = dds, type = "poscounts")
    dds <- DESeq2::DESeq(
      dds, 
      fitType = "parametric", 
      test = "Wald", 
      quiet = TRUE)
    out <- DESeq2::results(
      dds, 
      independentFiltering = FALSE, 
      cooksCutoff = FALSE)
  } else if (method == "ZILN") {
    good.ind <- which(rowSums(design[, c(1,2)]) == 1)
    ## app requires selection of specific phenolevels instead
    if (length(levels(pd)) > 2) {
      warning("Only using first two levels of phenotype")
    }
    design2 <- stats::model.matrix(~pd)
    fit <- fitFeatureModel(aggdat[, good.ind],
                           mod = design2[good.ind, c(1,2)],
                           spos = FALSE
    )
    out <- MRtable(fit, group = 3, number = Inf)
  } else if (method == "limma") {
    fit <- limma::lmFit(aggmat, design)
    contrast.matrix <- designPairs(colnames(design))
    fit2 <- limma::contrasts.fit(fit, contrast.matrix)
    fit2 <- limma::eBayes(fit2)
    out <- limma::topTable(fit2, coef = coef, adjust = "BH", number = Inf)
    
    if (is.null(coef)) {
      colnames(out)[seq_len(ncol(contrast.matrix))] <- 
        colnames(contrast.matrix)
    } else {
      colnames(out)[1] <- colnames(contrast.matrix)[coef]
    }
    
  } else if (method == "Kruskal-Wallis"){ 
    fit <- limma::lmFit(aggmat, design)
    contrast.matrix <- designPairs(colnames(design))
    fit2 <- limma::contrasts.fit(fit, contrast.matrix)
    fit2 <- limma::eBayes(fit2)
    out <- limma::topTable(fit2, coef = coef, adjust = "BH", number = Inf)
    
    colnames(out)[seq_len(ncol(contrast.matrix))] <- colnames(contrast.matrix)
    kw <- base::apply(
      aggmat, 1, function(x) broom::tidy(stats::kruskal.test(x, pd)))
    ## bind together and add row names
    kw.df <- dplyr::bind_rows(kw)
    kw.df$names <- names(kw)
    kw.df <- kw.df %>% tibble::column_to_rownames("names")
    matched <- match(rownames(out), rownames(kw.df))
    if ("F" %in% colnames(out)) {
      out$F <- kw.df$statistic[matched]
      colnames(out)[colnames(out) == "F"] <- "KW-Statistic"
    } else {
      out$t <- kw.df$statistic[matched]
      colnames(out)[colnames(out) == "t"] <- "KW-Statistic"
    }
    out$P.Value <- kw.df$p.value[matched]
    out$adj.P.Val <- stats::p.adjust(out$P.Value, method = "BH")
    out <- out[order(out$adj.P.Val), ]
  }
  out2 <- data.frame(rownames(out), out, check.names = FALSE) %>%
    dplyr::mutate_if(is.numeric, round, digits = getOption("me.round_digits"))
  colnames(out2)[1] <- level
  return(out2)
}
