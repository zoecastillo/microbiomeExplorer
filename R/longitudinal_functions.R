
#' Plot longitudinal features
#'
#' This function plots the reads of a particular feature over different time points.
#'
#' @param aggdat aggregated MRExperiment
#' @param feature Feature to plot.
#' @param x_var Phenotype to show along on X-axis.
#' @param id_var phenotype used to connect data points. Default is "SAMPLE_ID"
#' @param plotTitle Plot title. Default shows no title.
#' @param ylab Y-axis label. Default is "Reads"
#' @param log Log2 transform data. Default is FALSE.
#' @param showLines add lines between the points
#' @param fixedHeight sets a specific plot height (differential analysis)
#' @param x_levels restrict to specific levels of x_var (differential analysis)
#' @param pwidth overall plot width; default is 650 
#' 
#' @return plotly object holding long feature plot
#'
#' @author Janina Reeder, Mo Huang
#'
#' @importFrom metagenomeSeq MRcounts
#' @importFrom Biobase pData
#' 
#' @examples 
#' data("mouseData", package = "metagenomeSeq")
#' aggdat <- aggFeatures(mouseData, level = "genus")
#' plotLongFeature(aggdat, feature = "Prevotella", x_var = "diet", 
#'                 id_var = "mouseID")
#'
#' @export
plotLongFeature <- function(aggdat, feature, 
                            x_var, id_var = "SAMPLE_ID",
                            plotTitle = NULL, ylab = "Reads",
                            log = FALSE, showLines = TRUE,
                            fixedHeight = NULL, x_levels = NULL,
                            pwidth = 650) {
  yval <- NULL
  norm <- (ylab == "Percentage")
  aggmat <- MRcounts(aggdat, norm = norm)
  phenoTable <- pData(aggdat)
  if (log == TRUE) {
    aggmat <- log2(aggmat + 1)
    ylab <- paste0("Log ", ylab)
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
                         id_var = id_var)
  
  if(!is.null(x_levels)){
    df2 <- df2 %>%
      dplyr::filter(get(x_var) %in% x_levels)
  }
  
  ## set up color palette (same as relative abundance for color consistency)
  pal <- grDevices::colorRampPalette(c(RColorBrewer::brewer.pal(12, "Paired")))
  colvalues <- c(pal(12), "gray")
  
  if (feat_pos <= 12) {
    colvalues <- colvalues[feat_pos]
  } else {
    colvalues <- "gray"
  }
  
  ## group and summarize data by given settings: x_var, id_var  
  ## percentage (of x_var group) for percentage
  if (ylab == "Percentage") {
    df2$yval <- ifelse(df2$total == 0, 0, (df2$yval / df2$total) * 100)
  }
  
  df2 <- df2 %>%
    dplyr::group_by(x_var = get(x_var),
                    id_var = get(id_var)) %>%
    dplyr::mutate(id_mean = mean(yval, na.rm = TRUE))
  
  
  if(!is.null(x_levels)){
    df2$x_var <- forcats::fct_relevel(df2$x_var,as.character(x_levels))
  }
  
  ## add name of feature dataframe columns
  df2$Feature <- feature
  df2$text <- paste0("<b>",df2$samname,"</b><br />",
                     df2$Feature,": ",
                     round(df2$yval, digits = getOption("me.round_digits")))
  
  xaxis_text <- ""
  totalheight <- 420
  if(!is.null(fixedHeight)){
    totalheight <- fixedHeight
  }
  
  df2 <- df2 %>% 
    droplevels() %>%
    dplyr::group_by(id_var)
  
  p <- df2 %>%
    plotly::highlight_key(~id_var, group = id_var) %>%
    plotly::plot_ly(
      x = ~x_var, y = ~id_mean, color = ~Feature, 
      hoverinfo = "text", height = totalheight, width = pwidth
    ) %>%
    plotly::add_boxplot(
      colors = colvalues,
      pointpos = 0,
      showlegend = FALSE
    ) 
  if(showLines){
    p <- p %>% plotly::add_lines(
      color = I("black"),
      opacity = 0.2,
      name = id_var,
      text = ~id_var,
      showlegend = TRUE
    )
  }
  p <- p %>% plotly::add_markers(
    text = ~text,
    showlegend = FALSE
  ) %>%
    add_plotly_layout(plotTitle = plotTitle, 
                      xaxis_text = xaxis_text, 
                      ylab = ylab) %>%
    plotly::highlight(on = "plotly_click", 
                      off = "plotly_doubleclick",
                      selectize = TRUE,
                      dynamic = TRUE,
                      color = getOption("me.colorvalues"),
                      opacityDim = 0.4,
                      selected = plotly::attrs_selected(name = p$key))
  
  p <- p %>%
    plotly::layout(showlegend = TRUE)  %>%
    add_plotly_config()
  
  return(p)
}

