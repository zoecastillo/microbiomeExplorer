


#' Calls appropriate normalization functions depending on input parameter
#' The two available methods included in the package are based on either calculating proportions 
#' or by using cumulative sum scaling (CSS), Paulson, et al. Nat Meth 2013.
#'
#' @param MRobj the MRexperiment
#' @param norm_method method to use for normalization; CSS or Proportional
#'
#' @importFrom metagenomeSeq cumNorm cumNormStatFast normFactors
#'
#' @return the normalized MRobj
#' @export
normalizeData <- function(MRobj, norm_method) {
  if (norm_method == "CSS") {
    MRobj <- cumNorm(MRobj, suppressMessages(cumNormStatFast(MRobj)))
  } else if (norm_method == "Proportion") {
    normFactors(MRobj) <- colSums(MRcounts(MRobj))
  }
  MRobj
}


## AGGREGATION 


#' Aggregates counts by level
#'
#' This function aggregates counts by a level specified in the featureData slot
#' of the MRexperiment object.
#'
#' @param MRobj An MRexperiment object.
#' @param level Level to aggregate over. If NULL, no aggregation occurs.
#' @param sort boolean determining if resulting aggregated MRexperiment should 
#' be sorted based on rowSums; default is TRUE
#'
#' @return Aggregated MRexperiment object or matrix depending on \code{out}.
#'
#' @importFrom metagenomeSeq aggTax MRcounts normFactors normFactors<-
#' @importFrom Biobase pData<- fData<-
#'
#' @export
aggFeatures <- function(MRobj, level = NULL, sort = TRUE) {
  if (!is.null(level)) {
    ## account for absent feature table
    if (level == "unavailable") {
      aggobj <- MRobj
    } else {
      aggobj <- aggTax(MRobj, lvl = level)
      if(class(normFactors(MRobj)) == "numeric")
        normFactors(aggobj) <- normFactors(MRobj)
    }
  } else {
    aggobj <- MRobj
  }
  aggcts <- MRcounts(aggobj)
  ## account for absent feature table
  if(ncol(fData(aggobj)) > 0){
    ## replace empty string with unnamed
    rownames(aggcts)[rownames(aggcts) %in% ""] <- "unnamed"
    rownames(fData(aggobj))[rownames(fData(aggobj)) %in% ""] <- "unnamed"
    levels(fData(aggobj)[, level]) <- gsub("^$", "unnamed", 
                                           levels(fData(aggobj)[, level]))
  }
  aggobj@assayData$counts <- aggcts
  if (sort) {
    aggobj <- aggobj[order(rowSums(MRcounts(aggobj)), decreasing = TRUE), ]
  }
  return(aggobj)
}
