


## QC PLOTS ###################################################################


#' Plots sequencing statistics scatterplot
#'
#' This function makes a scatterplot of read and feature counts
#' for each sample. It was adjusted based on original work by Mo Huang
#'
#' @param MRobj metagenomeSeq object to be plotted
#' @param col_by factor by which to color the points
#' @param log character indicating which (if any) axes should be shown as log
#' @param filter_feat Numeric Y-coordinate to draw horizontal dashed line to
#' indicate feature filtering. If 0 (default), no line is drawn.
#' @param filter_read Numeric X-coordinate to draw vertical dashed line to
#' indicate read count filtering. If 0 (default), no line is drawn.
#' @param allowWebGL boolean indicating if webGL should be added
#' @param pwidth overall plot width; default is 550 (125 are added for legend)
#' @param pheight overall plot height; default is 550
#'
#' @author Janina Reeder
#'
#' @importFrom metagenomeSeq MRcounts
#'
#' @return the plotly QC plot
#' 
#' @examples
#' data("mouseData", package = "metagenomeSeq")
#' makeQCPlot(mouseData)
#'
#' @export
makeQCPlot <- function(MRobj, col_by = NULL, log = "none",
                       filter_feat = 0, 
                       filter_read = 0,
                       allowWebGL = TRUE, 
                       pwidth = 550, pheight = 550) {
  cts <- MRcounts(MRobj)
  ## how many features are present
  nfeatures <- colSums(cts > 0)
  ## how many reads are present 
  nreads <- colSums(cts)
  legendtitle <- ""
  
  if (!is.null(col_by) && col_by != "No color") {
    group <- forcats::fct_explicit_na(
      factor(pData(MRobj)[, which(colnames(pData(MRobj)) == col_by)]),
      na_level = "NA")
    ## adjust overall plotted area if legend is shown
    pwidth <- pwidth + 125
    legendtitle <- list(yref = "paper",xref = "paper",
                        y=1.06,x=1.2, 
                        text = col_by,showarrow = FALSE)
  } else {
    group <- I("black")
    legend <- FALSE
  }
  
  if(length(log) > 0){
    if(log == "x axis"){
      nreads <- log(nreads + 1)
    } else if(log == "y axis"){
      nfeatures <- log(nfeatures + 1)
    } else if(log == "both"){
      nreads <- log(nreads + 1)
      nfeatures <- log(nfeatures + 1)
    }
  }
  
  
  plotdata <- tibble::tibble(
    text = paste0(paste0("<b>", colnames(cts), "</b>"), 
                  paste0("<br>Reads: ", nreads), 
                  paste0("<br>Features: ", nfeatures)),
    x = nreads, y = nfeatures, color = group
  )
  
  useWebGL <- FALSE
  if(nrow(plotdata) > getOption("me.minwebgl") && allowWebGL)
    useWebGL <- TRUE
  
  ## ensure viewport shows everything
  max_x <- max(plotdata$x)
  max_x <- max_x + 0.05 * max_x
  min_x <- 0 - 0.05 * max_x
  max_y <- max(plotdata$y)
  max_y <- max_y + 0.05 * max_y
  min_y <- 0 - 0.05 * max_y
  
  ## make lines invisible if no filters are applied
  featalpha <- 0
  featalphafill <- 0
  if (filter_feat != 0) {
    featalpha <- 1
    featalphafill <- 0.2
  }
  readalpha <- 0
  readalphafill <- 0
  if (filter_read != 0) {
    readalpha <- 1
    readalphafill <- 0.2
  }
  
  
  ## generate actual plot
  p <- plotly::plot_ly(
    plotdata, hoverinfo = "text", 
    height = pheight, width = pwidth,
    colors = getOption("me.colorvalues")[seq(1, length(levels(group)))]) %>%
    plotly::add_segments(
      x = 0, y = filter_feat, xend = max_x, yend = filter_feat,
      line = list(color = paste0("rgba(119, 0, 0,", featalpha, ")"), 
                  width = 1, dash = "dash"), showlegend = FALSE,
      fill = "tozeroy", 
      fillcolor = paste0("rgba(219,219,219,", featalphafill, ")")
    ) %>%
    plotly::add_segments(
      x = filter_read, y = 0, xend = filter_read, yend = max_y,
      line = list(color = paste0("rgba(119, 0, 0,", readalpha, ")"), 
                  width = 1, dash = "dash"), showlegend = FALSE,
      fill = "tozerox", 
      fillcolor = paste0("rgba(219,219,219,", readalphafill, ")")
    ) %>%
    plotly::add_markers(
      x = ~x,
      y = ~y,
      color = ~color,
      text = ~text
    ) 
  if(useWebGL)
    p <- p %>% plotly::toWebGL() 
  
  p <- p %>% plotly::config(
    toImageButtonOptions = list(
      format = "svg"
    ),
    displaylogo = FALSE,
    modeBarButtons = getOption("me.modebar")
  ) %>%
    plotly::layout(
      xaxis = list(title = "&nbsp;\nNumber of reads", 
                   tickangle = 45,
                   range = c(min_x, max_x)),
      yaxis = list(title = paste0("Number of features\n&nbsp;"), 
                   range = c(min_y, max_y)),
      margin = list(l = 75, b = 75, t = 50, r = 50),
      annotations = legendtitle
    ) 
  return(p)
}


#' Function plotting a plotly histogram on the given histvalue
#'
#' @param histvalue the value to plot as a histogram
#' @param plotTitle title of the plot
#' @param xaxisTitle name of xaxis; default is ""
#' @param yaxisTitle name of yaxis; default is ""
#' @param pwidth overall plot width; default is 200 
#' @param pheight overall plot height; default is 200
#'
#' @return plotly plot object
#' 
#' @examples 
#' data("mouseData", package = "metagenomeSeq")
#' plotlyHistogram(histvalue = colSums(MRcounts(mouseData) > 0),
#'   plotTitle = "Feature distribution",
#'   xaxisTitle = "features", yaxisTitle = "frequency")
#' 
#' @export
plotlyHistogram <- function(histvalue, plotTitle, 
                            xaxisTitle = "", yaxisTitle = "",
                            pwidth = 200, pheight = 200) {
  ## ensure viewport includes 0
  xminus <- 0 - 0.01 * max(histvalue)
  
  p <- plotly::plot_ly(x = ~histvalue,
                       width = pwidth,
                       height = pheight) %>%
    plotly::add_histogram(
      marker = list(
        color = "rgba(219, 219, 219,0.8)",
        line = list(
          color = "black",
          width = 1
        )
      ),
      hovertemplate = "Frequency of <br> (%{x}) <br><b> %{y}</b>"
    ) %>%
    plotly::add_segments(
      x = xminus, y = 0, xend = xminus, yend = 0,
      line = list(color = "#444", width = 1), 
      showlegend = FALSE, hoverinfo = "none"
    ) %>%
    plotly::layout(
      title = list(
        text = plotTitle,
        font = list(
          size = 14
        )
      ),
      xaxis = list(
        title = list(
          text = xaxisTitle,
          font = list(
            size = 12
          )
        ),
        showgrid = FALSE,
        tickfont = list(
          size = 9
        )
      ),
      yaxis = list(
        title = list(
          text = yaxisTitle,
          font = list(
            size = 12
          )
        ),
        visible = TRUE,
        showgrid = FALSE,
        tickfont = list(
          size = 9
        )
      ),
      margin = list(
        pad = 1, t = 20,b = 20,l = 20,r = 20
      )
    ) %>%
    plotly::config(
      toImageButtonOptions = list(
        format = "svg"
      ),
      displaylogo = FALSE,
      modeBarButtons = list(list("toImage"))
    )
  p
}

#' Function plotting a barplot showing number of OTUs per samples
#'
#' @param MRobj containing data to plot
#' @param col_by phenotype to color bars by; default is NULL
#' @param xaxisTitle name of xaxis; default is ""
#' @param yaxisTitle name of yaxis; default is ""
#' @param pwidth overall plot width; default is 600 
#' @param pheight overall plot height; default is 450
#' @param sortbyfreq boolean determining if bars should be sorted by frequency; 
#' default is FALSE
#' @param pheno_sort order of pheno levels to sort by; 
#' ignored if sortbyfreq is TRUE
#' @param x_levels character vector holding x values in order to be shown
#' 
#' @importFrom metagenomeSeq MRcounts
#'
#' @return plotly plot object
#' 
#' @examples 
#' data("mouseData", package = "metagenomeSeq")
#' plotlySampleBarplot(mouseData)
#' 
#' @export
plotlySampleBarplot <- function(MRobj, col_by = NULL,
                                xaxisTitle = "", yaxisTitle = "",
                                pwidth = 600, pheight = 450,
                                sortbyfreq = FALSE,
                                pheno_sort = NULL,
                                x_levels = NULL) {
  
  cts <- MRcounts(MRobj)
  numofotus <- colSums(cts > 0)
  showlegend <- FALSE
  
  if (!is.null(col_by) && col_by != "No color") {
    group <- forcats::fct_explicit_na(
      factor(pData(MRobj)[, which(colnames(pData(MRobj)) == col_by)]),
      na_level = "NA")
    showlegend <- TRUE
    pwidth <- pwidth + 100
  } else {
    group <- I("rgba(219, 219, 219,0.8)")
  }
  
  plotdata <- tibble::tibble(
    text = paste0(paste0("<b>", colnames(cts), "</b>"), 
                  paste0("<br>Number of unique features: ", numofotus)),
    x = colnames(cts), y = numofotus, color = group
  )
  if(isTRUE(sortbyfreq)){
    plotdata[["x"]] <- factor(
      plotdata[["x"]], 
      levels = unique(plotdata[["x"]])[order(plotdata[["y"]], 
                                             decreasing = TRUE)])
  } else if(!is.null(pheno_sort) && !is.null(x_levels)){
    x_factor <- forcats::fct_explicit_na(
      factor(pData(MRobj)[, which(colnames(pData(MRobj)) == pheno_sort)]),
      na_level = "NA")
    plotdata$x_factor <- forcats::fct_relevel(x_factor,as.character(x_levels))
    plotdata[["x"]] <- factor(
      plotdata[["x"]], 
      levels = plotdata[order(plotdata$x_factor), ][["x"]])
  }
  
  p <- plotly::plot_ly(
    plotdata, hoverinfo = "text",
    height = pheight, width = pwidth,
    colors = getOption("me.colorvalues")[seq(1, 
                                             length(levels(group)))]) %>%
    plotly::add_bars(
      x = ~x,
      y = ~y,
      color = ~color,
      text = ~text,
      showlegend = showlegend
    )  %>%
    plotly::layout(
      xaxis = list(title = "&nbsp;\nSamples",
                   tickfont = list(
                     size = 9
                   )),
      yaxis = list(title = "Number of features\n&nbsp;"),
      margin = list(l = 75, b = 120, t = 50, r = 50)
    ) %>%
    plotly::config(
      toImageButtonOptions = list(
        format = "svg"
      ),
      displaylogo = FALSE,
      modeBarButtons = getOption("me.modebar")
    )
  
  p
}

