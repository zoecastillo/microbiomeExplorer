## INTER SAMPLE ANALYSIS 


#' Function to compute the distance matrix using vegdist from the vegan package
#'
#' @param aggdat aggregated MRExperiment
#' @param dist_method distance method from vegan package (See ?vegan::vegdist for details)
#' @param log transform count matrix to log2; default is TRUE
#' @param nfeatures number of features to use; default is all
#'
#' @return distance as dist
#' @export
computeDistMat <- function(aggdat, dist_method, 
                           log = TRUE, nfeatures = nrow(aggmat)){
  aggmat <- MRcounts(aggdat, norm = TRUE)
  mat <- aggmat[order(matrixStats::rowSds(aggmat), decreasing = TRUE), ]
  mat <- mat[seq(1,nfeatures), ]
  if (log == TRUE) {
    mat <- log2(mat + 1)
  }
  
  mat <- t(mat)
  vegan::vegdist(mat, method = tolower(dist_method))
}


#' Function to compute the PCAs for a given distance matrix
#'
#' @param distmat the distance matrix
#' @param pcas 2-element vector of PCAs to include in results
#'
#' @return the x slot limited to pcas after calling stats::prcomp on distmat
#' @export
calculatePCAs <- function(distmat, pcas){
  pcaRes <- stats::prcomp(distmat)
  pcaRes$x[,pcas]
}


#' Plot beta diversity
#'
#' This functions plots the beta diversity as a PCoA plot.
#'
#' @param aggdat aggregated MRExperiment
#' @param dim Vector of length 2 specifying which dimensions to plot.
#' @param log Log2 transform data. Default is TRUE.
#' @param dist_method Which distance method to use. See ?vegan::vegdist for details
#' \code{\link[vegan]{vegdist}()} for options. Default is "bray".
#' @param pcas precalculated pcas to avoid recalculation via CalcPCs
#' @param nfeatures Number of top features in terms of standard deviation.
#' Default is all.
#' @param col_by Phenotype for coloring.
#' @param shape_by Phenotype for shape. 
#' @param plotTitle Plot title. By default, becomes PCoA (code{dist.method}).
#' @param xlab X-axis label. By default, shows dimension and percent variance
#' explained.
#' @param ylab Y-axis label. By default, shows dimension and percent variance
#' explained.
#' @param pt_size the size of the markers
#' @param plotText adonis text to be added to plot
#' @param confInterval numeric value indicating confidence level for ellipses
#' @param allowWebGL boolean indicating if WebGL should be used for large dataadded
#' @param pwidth overall plot width; default is 550 (125 are added for legend)
#' @param pheight overall plot height; default is 550
#' 
#' @importFrom metagenomeSeq MRcounts
#' @importFrom Biobase pData
#'
#' @return plotly plot object
#'
#' @export
plotBeta <- function(aggdat, dim = 1:2,
                     log = TRUE, dist_method = "bray", pcas = NULL,
                     nfeatures = nrow(aggdat), col_by = NULL, shape_by = NULL,
                     plotTitle = "", xlab = NULL, ylab = NULL, pt_size = 8,
                     plotText = NULL, confInterval = NULL, allowWebGL = TRUE,
                     pwidth = 550, pheight = 550) {
  color <- NULL
  aggmat <- MRcounts(aggdat, norm = TRUE)
  phenoTable <- pData(aggdat)
  
  legendtitle <- ""
  
  if(is.null(pcas)){
    distmat <- computeDistMat(aggdat,
                              dist_method = dist_method,
                              log = log,
                              nfeatures = nfeatures)
    pcaRes <- stats::prcomp(distmat)
    pcas <- pcaRes$x[, dim]
  }
  
  legend <- FALSE
  testellipses <- NULL
  if (!is.null(col_by)) {
    colgroup <- forcats::fct_explicit_na(
      factor(phenoTable[, which(colnames(phenoTable) == col_by)]),
      na_level = "NA")
    if(!is.null(confInterval)){
      phenoTable[, which(colnames(phenoTable) == col_by)] <- colgroup
      ## need to through out single levels
      dup <- unique(colgroup[duplicated(colgroup)])
      elligroup <- factor(colgroup[colgroup %in% dup])
      ellipseData <- as.data.frame(pcas)
      ellipseData$SAMPLE_ID <- rownames(ellipseData)
      ellipseData <- dplyr::left_join(ellipseData,
                                      phenoTable[,c("SAMPLE_ID",col_by)])
      ellipseData <- ellipseData[ellipseData[[col_by]] %in% dup,]
      ## dataEllipse requires a plot even if draw = FALSE
      graphics::plot.new()
      testellipses <- car::dataEllipse(x = ellipseData[,dim[1]],
                                       y = ellipseData[,dim[2]],
                                       groups = factor(ellipseData[[col_by]]),
                                       draw = FALSE, 
                                       levels = confInterval,
                                       col = getOption("me.colorvalues")[seq(1, length(levels(elligroup)))])
      grDevices::dev.off()
      keepindices <- sapply(levels(elligroup), function(c){
        subEllipse <- ellipseData[ellipseData[[col_by]] == c,]
        if(sum(abs(stats::cor(subEllipse[,1:2]))) == 4){
          return(FALSE)
        }
        return(TRUE)
      })
      testellipses <- testellipses[keepindices]
    }
    
    pwidth <- pwidth + 125
    legendtitle <- list(yref = "paper",xref = "paper",
                        y=1.06,x=1.15, 
                        text = col_by, showarrow = FALSE)
    legend <- TRUE
  } else {
    colgroup <- I("black")
  }
  
  if (!is.null(shape_by)) {
    shapegroup <- forcats::fct_explicit_na(
      factor(phenoTable[, which(colnames(phenoTable) == shape_by)]),
      na_level = "NA")
    if(legend){
      legendtitle <- list(yref = "paper",xref = "paper",
                          y=1.06, x=1.15, 
                          text = paste0(col_by," - ",shape_by), 
                          showarrow = FALSE)
    } else {
      legendtitle <- list(yref = "paper",xref = "paper",
                          y=1.06,x=1.15, 
                          text = shape_by, showarrow = FALSE)
    }
    legend <- TRUE
    if(is.null(col_by))
      pwidth <- pwidth + 125
  } else {
    shapegroup <- I("circle")
  }
  
  
  plotdata <- tibble::tibble(
    text = paste0(
      paste0("<b>", colnames(aggmat), "</b>"),
      paste0("<br>", dim[1], ": ", 
             round(pcas[, 1], digits = getOption("me.round_digits"))),
      paste0("<br>", dim[2], ": ", 
             round(pcas[, 2], digits = getOption("me.round_digits")))
    ),
    SAMPLE_ID = colnames(aggmat), 
    x = pcas[, 1], 
    y = pcas[, 2], 
    color = colgroup,
    shape = shapegroup,
    adonis = "TEST"
  ) %>% dplyr::arrange(color)
  
  useWebGL <- FALSE
  if(nrow(plotdata) > getOption("me.minwebgl") && allowWebGL)
    useWebGL <- TRUE
  
  max_y <- max(plotdata$y)
  max_y <- max_y + 0.05 * max_y
  min_x <- min(plotdata$x) - 0.05 * min(plotdata$x)
  max_x <- max(plotdata$x)
  max_x <- max_x + 0.05 * max_x
  center_x <- (max_x + min_x)/2
  if(is.null(plotText))
    plotText <- ""
  
  plot_xlab <- 'if'(is.null(xlab), colnames(pcas)[1], xlab)
  plot_ylab <- 'if'(is.null(ylab), colnames(pcas)[2], ylab)
  
  p <- plotly::plot_ly(plotdata, hoverinfo = "text", 
                       height = pheight, width = pwidth,
                       colors = 
                         getOption("me.colorvalues")[seq(1, length(levels(colgroup)))]) 
  if(!is.null(testellipses)){
    res <- lapply(seq(1,length(testellipses)), function(t){
      p <<- p %>% plotly::add_polygons(x = testellipses[[t]][,1], 
                                       y = testellipses[[t]][,2],
                                       fillcolor = 'transparent',
                                       name = names(testellipses)[[t]],
                                       line = list(color = getOption("me.colorvalues")[t]))
    })
  }
  p <- p %>%
    plotly::add_text(
      x = center_x,
      y = max_y,
      text = plotText,
      showlegend = FALSE
    ) %>%
    plotly::add_markers(
      x = ~x,
      y = ~y,
      color = ~color,
      symbol = ~shape,
      text = ~text,
      marker = list(size = pt_size),
      customdata = ~SAMPLE_ID,
      showlegend = legend
    ) 
  if(useWebGL)
    p <- p %>% plotly::toWebGL()
  p <- p %>% 
    add_plotly_config() %>%
    plotly::layout(
      title = list(
        text = plotTitle,
        font = list(
          size = 14
        )),
      xaxis = list(title = paste0("&nbsp;\n", plot_xlab),
                   tickangle = 45),
      yaxis = list(title = paste0(plot_ylab, "\n&nbsp;")),
      margin = list(l = 75, b = 75, t = 50, r = 50),
      annotations = legendtitle
    )
  
  return(p)
}



#' Plot heatmap
#'
#' This function plots a heatmap of feature abundance.
#'
#' @param aggdat aggregated MRExperiment
#' @param features Vector of features to plot. If NULL, the top `nfeat`
#' features in terms of `sort_by` will be plotted.
#' @param log Log2 transform data. Default is TRUE.
#' @param sort_by Dispersion measure to sort features, one of "Fano", "MAD",
#' and "Variance"
#' @param nfeat Number of features to display. Default is 50.
#' @param col_by Vector of phenotypes for coloring.
#' @param row_by Name of feature level for coloring.
#' @param plotTitle Plot title. By default, no title.
#'
#' @importFrom metagenomeSeq MRcounts
#' @importFrom Biobase pData
#'
#' @return plotly heatmap
#' @export
plotHeatmap <- function(aggdat, features = NULL,
                        log = TRUE, sort_by = c("Fano", "MAD", "Variance"),
                        nfeat = 50, col_by = NULL, 
                        row_by = NULL, plotTitle = "") {
  aggmat <- MRcounts(aggdat, norm = TRUE)
  phenoTable <- pData(aggdat)
  
  if (is.null(nfeat)) {
    nfeat <- 50
  }
  
  ## show only specified features or return NULL if not found
  if (!is.null(features)) {
    aggmat <- aggmat[rownames(aggmat) %in% features, ]
    if (nrow(aggmat) == 0) {
      return(NULL)
    }
    aggmat <- aggmat[seq(1,min(length(features), nfeat)), ]
  } else {
    ind <- switch(sort_by,
                  "Fano" = matrixStats::rowVars(sqrt(aggmat)) / 
                    matrixStats::rowMeans2(sqrt(aggmat)), 
                  "MAD" = matrixStats::rowMads(log2(aggmat + 1)), 
                  "Variance" = matrixStats::rowVars(log2(aggmat + 1))
    ) 
    ind <- order(ind, decreasing = TRUE)[seq(1,min(length(ind), nfeat))]
    aggmat <- aggmat[ind, ]
  }
  
  if (log) {
    aggmat <- log2(aggmat + 1)
  }
  
  col_group <- NULL
  if (!is.null(col_by)) {
    col_group <- dplyr::select(phenoTable, col_by)
  }
  
  row_group <- NULL
  if (!is.null(row_by)) { 
    if (!row_by == "") {
      ## assign each feature to its group: rownames of aggregated fData 
      ## always correspond to aggregation level
      featData <- dplyr::select(fData(aggdat), row_by)
      rownames(featData)[rownames(featData) == ""] <- "unnamed"
      ## select the correct higher feature level based on rowname order
      row_group <- data.frame(featData[rownames(aggmat), ])
      names(row_group) <- row_by
    }
  }
  
  aggmat <- round(aggmat, digits = getOption("me.round_digits"))
  
  pal <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdYlBu")))
  colvalues <- pal(100)
  minval <- min(aggmat)
  maxval <- max(aggmat)
  
  if (is.null(col_group)) {
    if (is.null(row_group)) {
      hm <- heatmaply::heatmaply(aggmat,
                                 cexRow = 0.9,
                                 dendrogram = "column",
                                 col = colvalues,
                                 limits = c(minval, maxval),
                                 Rowv = FALSE,
                                 show_legend = FALSE,
                                 branches_lwd = 0.3,
                                 plot_method = "plotly"
      ) %>%
        plotly::colorbar(
          tickfont = list(size = 10), titlefont = list(size = 10),
          lenmode = "fraction", len = 0.2, which = 1
        )
    } else {
      hm <- heatmaply::heatmaply(aggmat,
                                 cexRow = 0.9,
                                 dendrogram = "column",
                                 col = colvalues,
                                 limits = c(minval, maxval),
                                 Rowv = FALSE,
                                 row_side_colors = row_group,
                                 show_legend = FALSE,
                                 branches_lwd = 0.3,
                                 plot_method = "plotly"
      ) %>%
        plotly::colorbar(
          tickfont = list(size = 10), titlefont = list(size = 10),
          lenmode = "fraction", len = 0.2, which = 1
        ) %>%
        plotly::colorbar(
          tickfont = list(size = 10), titlefont = list(size = 10),
          lenmode = "fraction", len = 0.3, which = 2
        )
    }
  } else {
    if (is.null(row_group)) {
      hm <- heatmaply::heatmaply(aggmat,
                                 cexRow = 0.9,
                                 dendrogram = "column",
                                 col = colvalues,
                                 limits = c(minval, maxval),
                                 Rowv = FALSE,
                                 col_side_colors = col_group,
                                 show_legend = FALSE,
                                 branches_lwd = 0.3,
                                 plot_method = "plotly"
      ) %>%
        plotly::colorbar(
          tickfont = list(size = 10), titlefont = list(size = 10),
          lenmode = "fraction", len = 0.3, which = 1
        ) %>%
        plotly::colorbar(
          tickfont = list(size = 10), titlefont = list(size = 10),
          lenmode = "fraction", len = 0.2, which = 2
        )
    } else {
      hm <- heatmaply::heatmaply(aggmat,
                                 cexRow = 0.9,
                                 dendrogram = "column",
                                 col = colvalues,
                                 limits = c(minval, maxval),
                                 Rowv = FALSE,
                                 col_side_colors = col_group,
                                 row_side_colors = row_group,
                                 show_legend = FALSE,
                                 branches_lwd = 0.3,
                                 plot_method = "plotly"
      ) %>%
        plotly::colorbar(
          tickfont = list(size = 10), titlefont = list(size = 10),
          lenmode = "fraction", len = 0.3, which = 1
        ) %>%
        plotly::colorbar(
          tickfont = list(size = 10), titlefont = list(size = 10),
          lenmode = "fraction", len = 0.2, which = 2
        ) %>%
        plotly::colorbar(
          tickfont = list(size = 10), titlefont = list(size = 10),
          lenmode = "fraction", len = 0.3, y = 0.5, which = 3
        )
    }
  }
  
  hm <- hm %>%
    plotly::layout(
      title = list(
        text = plotTitle,
        font = list(
          size = 14
        )),
      margin = list(l = 75, b = 75, t = 75, r = 50),
      showlegend = FALSE) %>%
  add_plotly_config()
  
  return(hm)
}
