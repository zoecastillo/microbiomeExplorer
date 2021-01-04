## AVERAGE FEATURE ABUNDANCE ANALYSIS 

#' Plot average relative abundance
#'
#' This function plots the average relative abundance of the top abundant features.
#'
#' @param aggdat aggregated MRExperiment object
#' @param level Feature level.
#' @param ind Indices of top abundant features to plot. Rest of features are
#' aggregated and displayed as "other".
#' @param plotTitle Plot title. Default shows no title.
#' @param ylab Y-axis label. Default is "Reads"
#' @param facet1 Phenotype for facet 1.
#' @param facet2 Phenotype for facet 2.
#' @param source name of the plot (needed for event handling); default is "A"
#' @param pwidth overall plot width; default is 500
#' @param pheight overall plot height; default is 150
#'
#' @author Janina Reeder
#'
#' @importFrom metagenomeSeq MRcounts
#' @importFrom Biobase pData
#'
#' @return plotly plot
#' 
#' @examples 
#' data("mouseData", package = "metagenomeSeq")
#' aggdat <- aggFeatures(mouseData, level = "genus")
#' plotAvgAbundance(aggdat, level = "genus")
#'
#' @export
plotAvgAbundance <- function(aggdat, level, ind = seq_len(10), plotTitle = "", 
                          ylab = "Reads", facet1 = NULL, facet2 = NULL, 
                          source = "A", pwidth = 500, pheight = 150) {
  facets <- NULL
  facet2s <- NULL 
  yval <- NULL
  Reads <- NULL
  norm <- (ylab == "Percentage")
  aggmat <- MRcounts(aggdat, norm = norm)
  phenoTable <- pData(aggdat)
  ordmat <- aggmat
  
  ## combine all other features as "other"
  if (nrow(aggmat) > (max(ind) + 1)) {
    ordmat <- rbind(aggmat[ind, ], 
                    other = colSums(aggmat[seq((max(ind) + 1),nrow(aggmat))
                                           , ]))
  }
  
  
  df <- as.data.frame(t(ordmat))
  
  df$samname <- rownames(df)
  df2 <- buildPlottingDF(df, phenoTable, facet1 = facet1, facet2 = facet2)
  
  ## find facet levels or add "nofacets"
  if (!is.null(facet1)) {
    facetvals <- levels(df2[, facet1])
  } else {
    facet1 <- "nofacet"
    df2$nofacet <- "nofacets"
    facetvals <- "nofacets"
  }

  ## find facet levels or add "nofacets:
  if (!is.null(facet2)) {
    facetvals2 <- levels(df2[, facet2])
  } else {
    facet2 <- "nofacet2"
    df2$nofacet2 <- "nofacets"
    facetvals2 <- "nofacets"
  }
  
  ## define color palette
  pal <- grDevices::colorRampPalette(
    c(RColorBrewer::brewer.pal(min(length(ind), 12), "Paired")))
  colvalues <- c(pal(max(ind)), "gray")
  
  # percentage 
  if (ylab == "Percentage") {
    rs <- rowSums(df2[,seq(1,ncol(df2)-3)])
    ## we should not reach this. Remove samples that are all 0 to avoid NaN.
    df2 <- df2[rs != 0,]
    df2[,seq(1,ncol(df2)-3)] <- df2[,seq(1,ncol(df2)-3)]/rs * 100
  }

  
  meanvals <- stats::aggregate(
    df2[,seq(1,ncol(df2)-3)],by = list(df2[[facet1]], df2[[facet2]]), mean)
  sdvals <- stats::aggregate(
    df2[,seq(1,ncol(df2)-3)],by = list(df2[[facet1]], df2[[facet2]]), 
    function(x) stats::sd(x)/sqrt(length(x)))
  meanvals <- reshape2::melt(
    meanvals, id.vars = c("Group.1","Group.2"), value.name = "mean")
  sdvals <- reshape2::melt(
    sdvals, id.vars = c("Group.1","Group.2"), value.name = "sd")
  plotData <- dplyr::full_join(meanvals, sdvals)
  names(plotData) <- c("facets", "facet2s", "feature", "mean", "sd")

  maxj <- length(facetvals2)
  maxi <- length(facetvals)
  totalwidth <- pwidth + 200 * length(facetvals)
  totalheight <- pheight + 300 * length(facetvals2)
  xaxis_text <- ""
  yaxis_text <- ""

  ## iterate over facet2 around facet1
  ## for each facetvalue in facet2 we need a row
  ## for each facetvalue in facet1 we need a column
  plotlist <- lapply(facetvals2, function(fv2) {
    # group and filter for each facetvalue in facet2
    j <- which(facetvals2 %in% fv2)
    sdf <- plotData %>%
      dplyr::group_by(facet2s) %>%
      dplyr::filter(facet2s %in% fv2) 

    if (facet2 != "nofacet2") {
      yaxis_text <- 'if'(is.na(fv2),"",
                   paste0(facet2, " ", fv2)
      )
    }
    
    ## for each facetvalue in facet1, we need a column
    p <- lapply(facetvals, function(f) {
      ## group and filter for each facetvalue in facet1 given a specific 
      ## facet2
      i <- which(facetvals %in% f)
      sp <- sdf %>%
        dplyr::group_by(facets) %>%
        dplyr::filter(facets %in% f) 
      
      showL <- FALSE
      ## We need to add xaxis labels in the last row only
      if (j == maxj) {
        xaxis_text <- f
        if(i == maxi){
          showL <- TRUE
        }
      }
  
      ## no values available for this combination of facet2 and facet1:
      ## draw empty plot
      if (nrow(sp) == 0) {
        return(buildEmptyPlotlyPlot(xaxis_text, ylab))
      }
      
      sp <- sp %>% plotly::plot_ly(
                           x = ~mean,
                           y = ~feature,
                           color = ~feature,
                           colors = colvalues,
                           width = totalwidth,
                           height = totalheight,
                           legendgroup = ~feature) %>%
        plotly::add_bars(error_x = ~list(array = sd,
                                         color = '#000000'),
                         showlegend = showL) %>%
        add_plotly_layout(plotTitle = plotTitle, 
                          xaxis_text = xaxis_text,
                          ylab = yaxis_text)
      sp
    })
    ## combine all columns in one row
    plotly::subplot(p, nrows = 1, shareY = TRUE, 
                    shareX = TRUE, titleX = TRUE)
  })
  
  p <- plotly::subplot(plotlist, nrows = length(plotlist), 
                       shareY = FALSE, titleX = TRUE, titleY = TRUE,
                       margin = c(0.02,0.02,0.05,0.05)) %>%
      add_plotly_config()
  
  return(p)
}



