## INTRA SAMPLE ANALYSIS 

#' Plot relative abundance
#'
#' This function plots the relative abundance of the top abundant features.
#'
#' @param aggdat aggregated MRExperiment object
#' @param level Feature level.
#' @param x_var Phenotype to aggregate over on X-axis. Default by "SAMPLE_ID".
#' @param ind Indices of top abundant features to plot. Rest of features are
#' aggregated and displayed as "other".
#' @param plotTitle Plot title. Default shows no title.
#' @param ylab Y-axis label. Default is "Reads"
#' @param facet1 Phenotype for facet 1.
#' @param facet2 Phenotype for facet 2.
#' @param source name of the plot (needed for event handling); default is "A"
#' @param pwidth overall plot width; default is 650 
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
#' plotAbundance(aggdat, level = "genus", x_var = "diet")
#'
#' @export
plotAbundance <- function(aggdat, level, x_var = "SAMPLE_ID",
                          ind = seq_len(10), plotTitle = "", ylab = "Reads", 
                          facet1 = NULL, facet2 = NULL, source = "A",
                          pwidth = 650, pheight = 150) {
  facets <- NULL
  facet2s <- NULL 
  yval <- NULL
  Reads <- NULL
  norm <- (ylab == "Percentage")
  aggmat <- MRcounts(aggdat, norm = norm)
  phenoTable <- pData(aggdat)
  ordmat <- aggmat
  xlab <- as.character(x_var)
  
  ## combine all other features as "other"
  if (nrow(aggmat) > max(ind)) {
    ordmat <- rbind(aggmat[ind, ], 
                    other = colSums(aggmat[seq((max(ind) + 1),nrow(aggmat))
                                           , ]))
  }
  
  ## prepare datastructure for plotting: join with relevant phenodata
  df <- reshape2::melt(ordmat, varnames = c(level, "samname"), 
                       value.name = "Reads")
  df[, level] <- as.factor(df[, level])
  df[, "samname"] <- as.factor(df[, "samname"])
  
  df2 <- buildPlottingDF(df = df, 
                         phenoTable = phenoTable,
                         x_var = x_var, 
                         facet1 = facet1,
                         facet2 = facet2)
  
  ## find facet levels or add "nofacets"
  if (!is.null(facet1)) {
    facetvals <- levels(df2[, facet1])
    xlab <- paste0("", "<br>", xlab, "<br>", 
                   paste(facet1, facetvals, sep = " "))
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
  
  ## get yval for every group to be shown (by feature, x_var and both 
  ## facets as needed)
  ## mean (of x_var group) used for numbers
  df2 <- df2 %>%
    dplyr::group_by(Family = get(level), x_var = get(x_var), 
                    facets = get(facet1), facet2s = get(facet2)) %>%
    dplyr::summarise(yval = mean(Reads, na.rm = TRUE))
  
  collevels <- levels(df2[,"Family"])
  
  ## percentage (of x_var group) for percentage
  if (ylab == "Percentage") {
    dftotal <- df2 %>%
      dplyr::group_by(x_var, facets, facet2s) %>%
      dplyr::summarise(total = sum(yval))
    df2 <- suppressMessages(dplyr::left_join(df2, dftotal))
    df2$yval <- 'if'(df2$total == 0, 0, (df2$yval / df2$total) * 100)
  }
  
  
  ## get number of values for each facet element
  ## this determines how much of the total width each facet gets
  ## determine number of elements in each facets to set widths
  widths <- getWidths(df2, facets = "facets", x_var = "x_var")
  
  df2$text <- paste0(
    paste0("<b>", df2$Family, "</b>"),
    paste0("<br>", x_var, ": ", df2$x_var),
    paste0(
      "<br>", ylab, ": ",
      round(df2$yval, digits = getOption("me.round_digits")),
      'if'(ylab == "Percentage", "%", "")
    )
  )
  
  xaxis_text <- ""
  maxj <- length(facetvals2)
  totalheight <- pheight + 300 * length(facetvals2)
  
  legendlevel <- getLegendLevel(df2, facets = "facets", facet2s = "facet2s")
  
  
  ## iterate over facet2 around facet1
  ## for each facetvalue in facet2 we need a row
  ## for each facetvalue in facet1 we need a column
  plotlist <- lapply(facetvals2, function(fv2) {
    # group and filter for each facetvalue in facet2
    j <- which(facetvals2 %in% fv2)
    sdf <- df2 %>%
      dplyr::group_by(facet2s) %>%
      dplyr::filter(facet2s %in% fv2) %>%
      droplevels()
    
    if (facet2 != "nofacet2") {
      ylab <- 'if'(is.na(fv2),
                   paste0("NA<br>", ylab),
                   paste0(facet2, " ", fv2, "<br>", ylab)
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
        xaxis_text <- xlab[i]
        if (f == legendlevel) {
          ##need to find first non-empty facet
          showL <- TRUE
        }
      }
      
      ## no values available for this combination of facet2 and facet1:
      ## draw empty plot
      if (nrow(sp) == 0) {
        return(buildEmptyPlotlyPlot(xaxis_text, ylab))
      }
      
      ## something to plot is available; drop all other levels
      sp <- sp %>% droplevels()
      ## make stacked barplot for this specific combination of facet2 
      ## and facet1
      sp <- sp %>%
        plotly::plot_ly(
          x = ~x_var, y = ~yval, hoverinfo = "text", 
          height = totalheight, width = pwidth,
          color = ~Family, showlegend = showL, 
          legendgroup = ~Family
        ) %>%
        plotly::add_bars(
          colors = colvalues,
          text = ~text,
          hoverlabel = list(
            font = list(
              color = "black"
            ),
            bordercolor = "black",
            borderwidth = 2
          )
        ) %>%
        plotly::layout(barmode = "stack") %>%
        add_plotly_layout(plotTitle = plotTitle, 
                          xaxis_text = xaxis_text, 
                          ylab = ylab)
      sp
    })
    ## combine all columns in one row
    plotly::subplot(p, nrows = 1, shareY = TRUE, 
                    titleX = TRUE, widths = widths)
  })
  
  p <- plotly::subplot(plotlist, nrows = length(plotlist), 
                       shareY = FALSE, titleX = TRUE, titleY = TRUE) %>%
    add_plotly_config()
  
  return(p)
}




#' Plot features
#'
#' This function plots the reads of a particular feature or set of features.
#'
#' @param aggdat aggregated MRExperiment
#' @param feature Feature to plot.
#' @param x_var Phenotype to aggregate over on X-axis. Default by "SAMPLE_ID".
#' @param ind Indices of top abundant features to plot. Needed to determine 
#' appropriate color
#' @param plotTitle Plot title. Default shows no title.
#' @param ylab Y-axis label. Default is "Reads"
#' @param xlab X-axis label. If NULL, x_var will be used as label.
#' @param facet1 Phenotype for facet 1.
#' @param facet2 Phenotype for facet 2.
#' @param log Log2 transform data. Default is FALSE.
#' @param showPoints add points for each sample on plot
#' @param fixedHeight sets a specific plot height (differential analysis)
#' @param x_levels restrict to specific levels of x_var (differential analysis)
#' @param pwidth overall plot width; default is 650 
#' 
#' @return plotly plot object
#'
#' @author Janina Reeder
#' 
#' @importFrom metagenomeSeq MRcounts
#' @importFrom Biobase pData
#' 
#' @examples 
#' data("mouseData", package = "metagenomeSeq")
#' aggdat <- aggFeatures(mouseData, level = "genus")
#' plotSingleFeature(aggdat, feature = "Prevotella", x_var = "diet")
#'
#' @export
plotSingleFeature <- function(aggdat, feature = "other", x_var = "SAMPLE_ID",
                              ind = seq_len(10), plotTitle = NULL, 
                              ylab = "Reads",
                              xlab = NULL,
                              facet1 = NULL, facet2 = NULL,
                              log = FALSE, showPoints = FALSE,
                              fixedHeight = NULL, x_levels = NULL,
                              pwidth = 500) {
  facets <- NULL
  facet2s <- NULL
  norm <- (ylab == "Percentage")
  aggmat <- MRcounts(aggdat, norm = norm)
  phenoTable <- pData(aggdat)
  if (log == TRUE) {
    aggmat <- log2(aggmat + 1)
    ylab <- paste0("Log ", ylab)
  }
  
  if (feature == "other") {
    aggmat <- rbind(aggmat[ind, ], 
                    other = colSums(aggmat[seq((max(ind) + 1),
                                               nrow(aggmat))
                                           , ]))
  }
  feat_pos <- match(feature, rownames(aggmat))
  ## feature not found; this happens when data is reaggregated in App
  if (is.na(feat_pos)) {
    return(NULL)
  }
  
  df <- data.frame(samname = colnames(aggmat), 
                   yval = aggmat[feat_pos, ],
                   total = colSums(aggmat))
  
  df2 <- buildPlottingDF(df = df, 
                         phenoTable = phenoTable,
                         x_var = x_var, 
                         facet1 = facet1,
                         facet2 = facet2)
  
  if(!is.null(x_levels)){
    df2 <- df2 %>%
      dplyr::filter(get(x_var) %in% x_levels)
  }
  
  if(is.null(xlab))
    xlab <- as.character(x_var)
  
  ## get facetvalues or set as "nofacets"
  if (!is.null(facet1)) {
    facetvals <- levels(df2[, facet1])
    xlab <- paste0("", "<br>", xlab, "<br>", 
                   paste(facet1, facetvals, sep = " "))
  } else {
    facet1 <- "nofacet"
    df2$nofacet <- "nofacets"
    facetvals <- "nofacets"
  }
  
  ## get facetvalues or set as "nofacets"
  if (!is.null(facet2)) {
    facetvals2 <- levels(df2[, facet2])
  } else {
    facet2 <- "nofacet2"
    df2$nofacet2 <- "nofacets"
    facetvals2 <- "nofacets"
  }
  
  ## set up color palette (same as relative abundance for color consistency)
  pal <- grDevices::colorRampPalette(
    c(RColorBrewer::brewer.pal(min(length(ind), 12), "Paired")))
  colvalues <- c(pal(max(ind)), "gray")
  
  if (feat_pos <= max(ind)) {
    colvalues <- colvalues[feat_pos]
  } else {
    colvalues <- "gray"
  }
  
  ## group data by given settings: x_var, facet1 and facet2
  df2 <- df2 %>%
    dplyr::group_by(x_var = get(x_var), 
                    facets = get(facet1), facet2s = get(facet2))
  
  ## percentage (of x_var group) for percentage
  if (ylab == "Percentage") {
    df2$yval <- 'if'(df2$total == 0, 0, (df2$yval / df2$total) * 100)
  }
  
  ## add name of feature dataframe columns
  df2$Feature <- feature
  df2$text <- paste0("<b>",df2$samname,"</b><br />",
                     df2$Feature,": ",
                     round(df2$yval, digits = getOption("me.round_digits")))
  
  ## determine number of elements in each facets to set widths
  drop <- TRUE
  if(!is.null(x_levels))
    drop <- FALSE
  widths <- getWidths(df2, facets = "facets", x_var = "x_var", drop = drop)
  
  if (showPoints) {
    showPoints <- "all"
  }
  
  xaxis_text <- ""
  maxj <- length(facetvals2)
  totalheight <- 120 + 300 * length(facetvals2)
  if(!is.null(fixedHeight)){
    totalheight <- fixedHeight
  }
  
  ## iterate over facet2 around facet1
  ## for each facetvalue in facet2 we need a row
  plotlist <- lapply(facetvals2, function(fv2) {
    j <- which(facetvals2 %in% fv2)
    sdf <- df2 %>%
      dplyr::group_by(facet2s) %>%
      dplyr::filter(facet2s %in% fv2) %>%
      droplevels()
    
    if (facet2 != "nofacet2") {
      ylab <- 'if'(is.na(fv2),
                   paste0("NA<br>", ylab),
                   paste0(facet2, " ", fv2, "<br>", ylab)
      )
    }
    ## iterating over facet1: a column for each
    p <- lapply(facetvals, function(f) {
      i <- which(facetvals %in% f)
      sp <- sdf %>%
        dplyr::group_by(facets) %>%
        dplyr::filter(facets %in% f)
      ## add xaxis labels in last row
      if (j == maxj) {
        xaxis_text <- xlab[i]
      }
      ## no entries: return empty plot
      if (nrow(sp) == 0) {
        return(buildEmptyPlotlyPlot(xaxis_text, ylab))
      }
      sp <- sp %>% droplevels()
      sp <- sp %>%
        plotly::plot_ly(
          x = ~x_var, y = ~yval, color = ~Feature, 
          hoverinfo = "text", height = totalheight, width = pwidth,
          showlegend = FALSE, legendgroup = ~Feature
        ) %>%
        plotly::add_boxplot(
          colors = colvalues,
          boxpoints = showPoints,
          text = ~text,
          pointpos = 0
        ) %>%
        add_plotly_layout(plotTitle = plotTitle, 
                          xaxis_text = xaxis_text, 
                          ylab = ylab)
      sp
    })
    plotly::subplot(p, nrows = 1, shareY = TRUE, 
                    titleX = TRUE, widths = widths)
  })
  
  ## nrows will be determined by levels of facet2
  p <- plotly::subplot(plotlist, nrows = length(plotlist), 
                       shareY = FALSE, titleX = TRUE, titleY = TRUE) %>%
    plotly::layout(showlegend = TRUE) %>%
    add_plotly_config()
  
  return(p)
}


#' Plot alpha diversity
#'
#' This function plots the alpha diversity. See ?vegan::diversity for details
#' on the available index
#'
#' @param aggdat aggregated MRExperiment
#' @param level Feature level
#' @param index Diversity index, one of "shannon", "simpson", "invsimpson" or 
#' "richness" (=number of features). Default is "shannon".
#' @param x_var Phenotype to aggregate over on X-axis. Default by "SAMPLE_ID".
#' @param ylab Y-axis label. Default is "Reads".
#' @param col_by Phenotype for coloring.
#' @param facet1 Phenotype for facet 1.
#' @param facet2 Phenotype for facet 2. 
#' @param plotTitle Plot title. By default, no title is used.
#' @param pwidth overall plot width; default is 650 
#' @param pheight overall plot height; default is 150
#'
#' @return plotly plot object
#' 
#' @importFrom metagenomeSeq MRcounts
#' @importFrom Biobase pData
#' 
#' @examples 
#' data("mouseData", package = "metagenomeSeq")
#' aggdat <- aggFeatures(mouseData, level = "genus")
#' plotAlpha(aggdat, level = "genus", index = "shannon", x_var = "diet")
#'
#' @export
plotAlpha <- function(aggdat, level,
                      index = c("shannon", "simpson", 
                                "invsimpson", "richness"),
                      x_var = "SAMPLE_ID", ylab = index,
                      col_by = NULL, facet1 = NULL, facet2 = NULL,
                      plotTitle = "", pwidth = 500, pheight = 150) {
  facets <- NULL
  facet2s <- NULL
  aggmat <- MRcounts(aggdat)
  phenoTable <- pData(aggdat)
  
  ## DETERMINE ALPHA DIVERSITY 
  if (index == "richness") {
    mat <- colSums(aggmat > 0)
    Index <- "Num of Features Observed"
  } else {
    mat <- vegan::diversity(aggmat, index = index, MARGIN = 2)
    Index <- paste0(stringr::str_to_title(index), "_Diversity")
  }
  
  df <- data.frame(names(mat), mat)
  
  colnames(df) <- c("samname", "Diversity")
  legendtitle <- ""
  
  if (!is.null(col_by)) {
    legendtitle <- list(yref = "paper",xref = "paper",
                        y=1.06,x=1.2, 
                        text = col_by,showarrow = FALSE)
    col_name <- paste0(col_by, "_color")
    pwidth <- pwidth + 150
  } else {
    col_name <- NULL
  }
  xlab <- as.character(x_var)
  
  df2 <- buildPlottingDF(df = df, 
                         phenoTable = phenoTable,
                         x_var = x_var, 
                         facet1 = facet1,
                         facet2 = facet2,
                         col_by = col_by,
                         col_name = col_name)
  col_by <- col_name
  
  ## set up for facets and color
  if (!is.null(facet1)) {
    facetvals <- levels(df2[, facet1])
    xlab <- paste0("", "<br>", xlab, "<br>", 
                   paste(facet1, facetvals, sep = " "))
  } else {
    facet1 <- "nofacet"
    df2$nofacet <- "nofacets"
    facetvals <- "nofacets"
  }
  
  if (!is.null(facet2)) {
    facetvals2 <- levels(df2[, facet2])
  } else {
    facet2 <- "nofacet2"
    df2$nofacet2 <- "nofacets"
    facetvals2 <- "nofacets"
  }
  
  collevels <- levels(df2[,col_by])
  ## will be zero if col_by is NULL
  numofcols <- length(collevels)
  colvalues <- getOption("me.colorvalues")[seq(1, numofcols)]
  if (is.null(col_by)) {
    col_by <- "nocolor"
    collevels <- ""
    df2$nocolor <- index
    colvalues <- "#a5a39f"
  } 
  
  df2$text <- paste0(
    "<b>",df2$samname,"</b><br />",
    round(df2$Diversity, digits = getOption("me.round_digits")))
  df2 <- df2 %>%
    dplyr::group_by(x_var = get(x_var), 
                    facets = get(facet1), facet2s = get(facet2))
  
  ## determine number of elements in each facets to set widths
  widths <- getWidths(df2, facets = "facets", x_var = "x_var")
  
  xaxis_text <- ""
  yaxis_text <- Index
  maxj <- length(facetvals2)
  totalheight <- pheight + 300 * length(facetvals2)
  
  ## iterate over facet2 around facet1
  plotlist <- lapply(facetvals2, function(fv2) {
    j <- which(facetvals2 %in% fv2)
    sdf <- df2 %>%
      dplyr::group_by(facet2s) %>%
      dplyr::filter(facet2s %in% fv2)
    
    
    if (facet2 != "nofacet2") {
      yaxis_text <- 'if'(is.na(fv2),
                         paste0("NA<br>", Index),
                         paste0(facet2, " ", fv2, "<br>", Index)
      )
    }
    ## iterating over facet1
    p <- lapply(facetvals, function(f) {
      i <- which(facetvals %in% f)
      sp <- sdf %>%
        dplyr::group_by(facets) %>%
        dplyr::filter(facets %in% f)
      if (j == maxj) {
        xaxis_text <- xlab[i]
      }
      if (nrow(sp) == 0) {
        return(buildEmptyPlotlyPlot(xaxis_text, ylab))
      }
      sp <- sp %>%
        dplyr::mutate(colorfact = as.factor(get(col_by)))
      present_levels <- unique(sp$colorfact)
      showL <- FALSE
      if(any(present_levels %in% collevels))
        showL <- TRUE
      sp <- sp %>%
        plotly::plot_ly(
          x = ~x_var, y = ~Diversity, color = ~colorfact, 
          hoverinfo = "text", height = totalheight, width = pwidth,
          showlegend = showL, legendgroup = ~colorfact
        ) %>%
        plotly::add_boxplot(colors = colvalues, 
                            boxpoints = "all", 
                            text = ~text,
                            pointpos = 0) %>%
        add_plotly_layout(plotTitle = plotTitle, 
                          xaxis_text = xaxis_text, 
                          ylab = yaxis_text)
      ## using global environment assignment as we are substituting in place 
      collevels <<- collevels[!collevels %in% present_levels]
      sp
    })
    plotly::subplot(p, nrows = 1, shareY = TRUE, 
                    titleX = TRUE, widths = widths)
  })
  
  p <- plotly::subplot(plotlist, nrows = length(plotlist), 
                       shareY = FALSE, titleX = TRUE, titleY = TRUE) %>%
    plotly::layout(showlegend = TRUE,
                   colorway = getOption("me.colorvalues")) %>%
    add_plotly_config()
  
  return(p)
}
