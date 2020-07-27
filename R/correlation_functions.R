## CORRELATION ANALYSIS 


#' Helper function to calculate the confidence interval for a cor.test
#'
#' @param num number of samples
#' @param mS results of cor.test
#' @param method statistical method used for cor.test
#'
#' @return named vector holding lower and upper thresholds
computeCI_Interval <- function(num, mS, method){
  if(num <= 4)
    return(c("lower" = NA, "upper" = NA))
  if(method == "spearman"){
    rho <- mS$estimate
    stderr <- 1.0 / sqrt(num - 3)
    delta <- 1.96 * stderr
    lower <- tanh(atanh(rho) - delta)
    upper <- tanh(atanh(rho) + delta)
  } else if(method == "pearson"){
    lower <- mS$conf.int[1]
    upper <- mS$conf.int[2]
  } else if(method == "kendall"){
    tau <- mS$estimate
    z_tau <- 0.5 * log((1+tau)/(1-tau))
    lower <- z_tau - 1.959964 * sqrt(0.437/(num - 4))
    upper <- z_tau + 1.959964 * sqrt(0.437/(num - 4))
  }
  ci_interval <- c(lower, upper)
  names(ci_interval) <- c("lower","upper")
  return(ci_interval)
}


#' Scatterplot of two features
#'
#' This function plots a scatterplot of two features along with sample
#' correlation statistics.
#'
#' @param aggdat aggregated MRExperiment
#' @param feat1 Feature 1.
#' @param feat2 Feature 2.
#' @param log Log2 transform data. Default is TRUE.
#' @param method Correlation coefficient. One of "spearman" (default),
#' "pearson", or "kendall".
#' @param addRegression boolean parameter indicating whether linear regression
#' line should be drawn; default: TRUE
#' @param col_by Phenotype for coloring.
#' @param facet1 Phenotype for facet 1.
#' @param facet2 Phenotype for facet 2. 
#' @param plotTitle Plot title. Default is no title.
#' @param xlab X-axis label. Default is \code{feat1}.
#' @param ylab Y-axis label. Default is \code{feat2}.
#' @param allowWebGL boolean indicating if WebGL should be used for large data
#' @param pwidth overall plot width; default is 550 
#' @param pheight overall plot height; default is 200
#'
#' @return list holding plotly plot and lm fit
#'
#' @importFrom metagenomeSeq MRcounts
#' @importFrom Biobase pData
#' @importFrom magrittr %>%
#'
#' @export
corrFeature <- function(aggdat, feat1 = NULL, feat2 = NULL,
                        log = TRUE, method = c("spearman", "pearson", "kendall"),
                        addRegression = TRUE,
                        col_by = NULL, facet1 = NULL, facet2 = NULL,
                        plotTitle = "", xlab = NULL, ylab = NULL, 
                        allowWebGL = TRUE, pwidth = 550, pheight = 200) {
  facet2s <- NULL
  facets <- NULL
  aggmat <- MRcounts(aggdat, norm = TRUE)
  phenoTable <- pData(aggdat)
  
  if (is.null(feat1) || is.null(feat2)) {
    return(NULL)
  }
  
  method <- match.arg(method)
  if (log) {
    aggmat <- log2(aggmat + 1)
  }
  x <- aggmat[which(rownames(aggmat) == feat1), ]
  y <- aggmat[which(rownames(aggmat) == feat2), ]
  
  ## feature not found: return immediately
  if (length(x) == 0 || length(y) == 0) {
    return(list(plot = NULL, stats = NULL))
  }
  
  df <- data.frame(colnames(aggmat), x, y)
  colnames(df) <- c("samname", "feat1", "feat2")
  showlegend <- TRUE
  
  legendtitle <- ""
  addwidth <- 0
  
  legendtitle <- ""
  
  if (!is.null(col_by)) {
    group <- forcats::fct_explicit_na(
      factor(phenoTable[, which(colnames(phenoTable) == col_by)]),
      na_level = "NA")
    addwidth <- 200
    legendtitle <- list(yref = "paper",xref = "paper",
                        y=1.06,x=1.1,
                        text = col_by,showarrow = FALSE)
    col_name <- paste0(col_by, "_color")
    
    pwidth <- pwidth + 150
  } else {
    col_name <- NULL
    group <- I("black")
    showlegend <- FALSE
  }
  
  df2 <- buildPlottingDF(df = df, 
                         phenoTable = phenoTable,
                         facet1 = facet1,
                         facet2 = facet2,
                         col_by = col_by,
                         col_name = col_name)
  col_by <- col_name
  
  logtext <- ""
  if (log == TRUE)
    logtext <- "Log "
  
  if(is.null(xlab))
    xlab <- feat1
  
  if (!is.null(facet1)) {
    facetvals <- levels(df2[, facet1])
    xlab <- paste0(xlab, "<br>", paste(facet1, facetvals, sep = " "))
    pwidth <- length(facetvals) * 300
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
  
  df2$group <- group
  
  df2 <- df2 %>%
    dplyr::group_by(facets = get(facet1), facet2s = get(facet2))
  
  useWebGL <- FALSE
  if(nrow(df2) > getOption("me.minwebgl") && allowWebGL)
    useWebGL <- TRUE
  
  xaxis_text <- ""
  
  if(is.null(ylab))
    ylab <- feat2
  yaxis_text <- ylab
  
  maxj <- length(facetvals2)
  totalheight <- pheight + 400 * length(facetvals2)
  pwidth <- pwidth + addwidth
  
  fitlist <- list()
  
  ## iterate over facet2 around facet1
  plotlist <- lapply(facetvals2, function(fv2) {
    j <- which(facetvals2 %in% fv2)
    sdf <- df2 %>%
      dplyr::group_by(facet2s) %>%
      dplyr::filter(facet2s %in% fv2)
    
    if (facet2 != "nofacet2") {
      yaxis_text <- 'if'(is.na(fv2),
                         paste0(ylab, "<br>NA<br> "),
                         paste0(ylab, "<br>", facet2, " ", fv2, "<br> ")
      )
    }
    
    ## iterating over facet1
    plist <- lapply(facetvals, function(f) {
      i <- which(facetvals %in% f)
      sp <- sdf %>%
        dplyr::group_by(facets) %>%
        dplyr::filter(facets %in% f)
      
      if (j == maxj) {
        xaxis_text <- paste0(logtext, xlab[i])
      }
      if (nrow(sp) == 0) {
        return(buildEmptyPlotlyPlot(xaxis_text, paste0(logtext, yaxis_text)))
      }
      
      facetgroup <- sp$group
      if (is.null(col_by)) {
        facetgroup <- I("black")
      }
      
      ## PERFORM CORRELATION 
      fit <- stats::lm(feat2 ~ feat1, data = sp, na.action = stats::na.exclude)
      mS_label <- tryCatch({
        mS <- stats::cor.test(sp$feat2,sp$feat1, method = method)
        ci_interval <- computeCI_Interval(num = nrow(sp), mS = mS, method = method)
        mS$lower <- round(ci_interval[["lower"]], getOption("me.round_digits"))
        mS$upper <- round(ci_interval[["upper"]], getOption("me.round_digits"))

        ## store fit in fitlist
        fitlist[[gsub("<br>", "_", 
                      paste0(yaxis_text, " <br/> ", xlab[i]))]] <<- mS
        paste0("\n", names(mS$estimate)," = ",
               round(mS$estimate,getOption("me.round_digits")),
               "; p = ",round(mS$p.value,getOption("me.round_digits")),
               "\nlower CI = ", mS$lower,
               "\nupper CI = ", mS$upper)
      }, error = function(e){
        ""
      })
      
      plotdata <- tibble::tibble(
        samname = sp$samname,
        x = sp$feat1,
        y = sp$feat2,
        color = facetgroup,
        text = paste0(paste0("<b>", sp$samname, "</b>"), 
                      paste0("<br>", feat1, ": ",
                             round(sp$feat1, digits = getOption("me.round_digits"))), 
                      paste0("<br>", feat2, ": ", 
                             round(sp$feat2, digits = getOption("me.round_digits"))))
      ) 
      
      fitlabel <- paste(
        "Adj R2: ", signif(summary(fit)$adj.r.squared, 5),
        "<br>Intercept: ", signif(fit$coef[[1]], 5),
        "<br>Slope: ", signif(fit$coef[[2]], 5),
        "<br>P: ", 'if'(nrow(summary(fit)$coef) < 2 || 
                          ncol(summary(fit)$coef < 4),
                        NA,
                        signif(summary(fit)$coef[2, 4], 5)
        )
      )
      
      p <- plotdata %>%
        plotly::plot_ly(x = ~x, hoverinfo = "text", 
                        width = pwidth, height = totalheight) %>%
        plotly::add_markers(
          y = ~y,
          color = ~color,
          colors = getOption("me.colorvalues")[seq(1, length(levels(group)))],
          text = ~text,
          showlegend = showlegend
        ) %>%
        plotly::add_lines(
          y = ~x,
          line = list(
            dash = "dash",
            color = "rgba(220,220,220,1)"
          ),
          name = "onetoone",
          text = "",
          showlegend = FALSE
        )
      if(addRegression){
        p <- p %>%
          plotly::add_lines(
            y = stats::predict(fit),
            line = list(
              color = "rgba(119,0,0, 1)",
              shape = "spline"
            ),
            name = "Linear Regression",
            text = fitlabel,
            showlegend = FALSE
          ) %>%
          plotly::add_ribbons(
            data = broom::augment(fit, se_fit = TRUE),
            x = ~feat1,
            ymin = ~ .fitted - 1.96 * .se.fit,
            ymax = ~ .fitted + 1.96 * .se.fit,
            line = list(color = "rgba(219,219,219, 0.35)"),
            fillcolor = "rgba(219,219,219, 0.3)",
            name = "Standard Error",
            showlegend = FALSE
          ) 
      }
      if(useWebGL)
        p <- p %>% plotly::toWebGL()
      p <- p %>% 
        add_plotly_layout(plotTitle = plotTitle, 
                          xaxis_text = paste0(xaxis_text, mS_label), 
                          ylab = paste0(logtext, yaxis_text)) %>%
        plotly::layout(margin = list(b = 150))
      p
    })
    plotly::subplot(plist, nrows = 1, shareY = TRUE, titleX = TRUE)
  })
  
  
  p <- plotly::subplot(plotlist, nrows = length(plotlist),   
               shareY = FALSE, titleX = TRUE, titleY = TRUE) %>%

  add_plotly_config()
  
  return(list(plot = p, stats = fitlist))
}

#' Scatterplot of a feature and a phenotype
#'
#' This function plots a scatterplot of a feature and a phenotype along with sample
#' correlation statistics.
#'
#' @param aggdat aggregated MRExperiment
#' @param feature Feature input.
#' @param phenotype Phenotype input (must be numeric)
#' @param log Log2 transform data. Default is TRUE.
#' @param method Correlation coefficient. One of "spearman" (default),
#' "pearson", or "kendall".
#' @param addRegression boolean parameter indicating whether linear regression
#' line should be drawn; default: TRUE
#' @param col_by Phenotype for coloring.
#' @param facet1 Phenotype for facet 1.
#' @param facet2 Phenotype for facet 2. (WIP/TODO)
#' @param plotTitle Plot title. Default is no title.
#' @param xlab X-axis label. Default is \code{feat1}.
#' @param ylab Y-axis label. Default is \code{feat2}.
#' @param allowWebGL boolean indicating if WebGL should be used for large data
#' @param pwidth overall plot width; default is 550 
#' @param pheight overall plot height; default is 200
#'
#' @return list holding plotly plot and lm fit
#'
#' @importFrom metagenomeSeq MRcounts
#' @importFrom Biobase pData
#'
#' @export
corrPhenotype <- function(aggdat, feature = NULL, phenotype = NULL, log = TRUE, 
                          method = c("spearman", "pearson", "kendall"),
                          addRegression = TRUE,
                          col_by = NULL, facet1 = NULL, facet2 = NULL,
                          plotTitle = "", xlab = NULL, ylab = NULL,
                          allowWebGL = TRUE, pwidth = 550, pheight = 200) {
  facet2s <- NULL
  facets <- NULL
  aggmat <- MRcounts(aggdat, norm = TRUE)
  phenoTable <- pData(aggdat)
  
  if (is.null(feature) | is.null(phenotype)) {
    return(NULL)
  }
  
  method <- match.arg(method)
  if (log) {
    aggmat <- log2(aggmat + 1)
  }
  
  phenoval <- phenoTable[, phenotype]
  
  featureval <- aggmat[which(rownames(aggmat) == feature), ]
  if (length(featureval) == 0) {
    return(list(plot = NULL, stats = NULL))
  }
  
  df <- data.frame(colnames(aggmat), featureval, phenoval)
  colnames(df) <- c("samname", "feature", "phenotype")
  showlegend <- TRUE
  
  legendtitle <- ""
  addwidth <- 0
  
  if (!is.null(col_by)) {
    group <- forcats::fct_explicit_na(
      factor(phenoTable[, which(colnames(phenoTable) == col_by)]),
      na_level = "NA")
    addwidth <- 200
    legendtitle <- list(yref = "paper",xref = "paper",
                        y=1.06,x=1.1,
                        text = col_by,showarrow = FALSE)
    col_name <- paste0(col_by, "_color")
    
    pwidth <- pwidth + 150
  } else {
    col_name <- NULL
    group <- I("black")
    showlegend <- FALSE
  }
  
  df2 <- buildPlottingDF(df = df, 
                         phenoTable = phenoTable,
                         facet1 = facet1,
                         facet2 = facet2,
                         col_by = col_by,
                         col_name = col_name)
  col_by <- col_name
  
  logtext <- ""
  if(log == TRUE)
    logtext <- "Log "
  
  if(is.null(xlab))
    xlab <- feature
  
  if (!is.null(facet1)) {
    facetvals <- levels(df2[, facet1])
    xlab <- paste0(xlab, "<br>", paste(facet1, facetvals, sep = " "))
    pwidth <- length(facetvals) * 300
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
  
  df2$group <- group
  df2 <- df2 %>%
    dplyr::group_by(facets = get(facet1), facet2s = get(facet2))
  
  useWebGL <- FALSE
  if(nrow(df2) > getOption("me.minwebgl") && allowWebGL)
    useWebGL <- TRUE
  
  xaxis_text <- ""
  if(is.null(ylab))
    ylab <- phenotype
  yaxis_text <- ylab
  
  maxj <- length(facetvals2)
  totalheight <- pheight + 400 * length(facetvals2)
  pwidth <- pwidth + addwidth
  
  fitlist <- c()
  ## iterate over facet2 around facet1
  plotlist <- lapply(facetvals2, function(fv2) {
    j <- which(facetvals2 %in% fv2)
    sdf <- df2 %>%
      dplyr::group_by(facet2s) %>%
      dplyr::filter(facet2s %in% fv2)
    
    if (facet2 != "nofacet2") {
      yaxis_text <- 'if'(is.na(fv2),
                         paste0(phenotype, "<br>NA<br> "),
                         paste0(phenotype, "<br>", facet2, " ", fv2, "<br> ")
      )
    }
    # iterating over facet1
    plist <- lapply(facetvals, function(f) {
      i <- which(facetvals %in% f)
      sp <- sdf %>%
        dplyr::group_by(facets) %>%
        dplyr::filter(facets %in% f)
      
      if (j == maxj) {
        xaxis_text <- paste0(logtext, xlab[i])
      }
      
      emptyplot <- FALSE
      if (nrow(sp) == 0) {
        emptyplot <- TRUE
      } else if (is.na(sp$phenotype)) {
        emptyplot <- TRUE
      }
      if (emptyplot) {
        return(buildEmptyPlotlyPlot(xaxis_text, paste0(logtext, yaxis_text)))
      }
      
      facetgroup <- sp$group
      if (is.null(col_by)) {
        facetgroup <- I("black")
      }
      
      ## PERFORM CORRELATION 
      fit <- stats::lm(phenotype ~ feature, data = sp, 
                       na.action = stats::na.exclude)       
      mS_label <- tryCatch({
        mS <- stats::cor.test(sp$phenotype, sp$feature, method = method)
        ci_interval <- computeCI_Interval(num = nrow(sp), mS = mS, method = method)
        mS$lower <- round(ci_interval[["lower"]], getOption("me.round_digits"))
        mS$upper <- round(ci_interval[["upper"]], getOption("me.round_digits"))
        ## store fit in fitlist
        fitlist[[gsub("<br>", "_", 
                      paste0(yaxis_text, " <br/> ", xlab[i]))]] <<- mS
        paste0("\n", names(mS$estimate)," = ",
               round(mS$estimate,getOption("me.round_digits")),
               "; p = ",round(mS$p.value,getOption("me.round_digits")),
               "\nlower CI = ", mS$lower,
               "\nupper CI = ", mS$upper)
      }, error = function(e){
        ""
      })
      
      pt <- f
      if (is.na(pt)) {
        pt <- "NA"
      } else if (pt == "nofacets") {
        pt <- ""
      }
      
      plotdata <- tibble::tibble(
        samname = sp$samname,
        x = sp$feature,
        y = sp$phenotype,
        color = facetgroup,
        text = paste0(paste0("<b>", sp$samname, "</b>"), 
                      paste0("<br>", feature, ": ", 
                             round(sp$feature, 
                                   digits = getOption("me.round_digits"))), 
                      paste0("<br>", phenotype, ": ", 
                             round(sp$phenotype, 
                                   digits = getOption("me.round_digits"))))
      )
      
      fitlabel <- paste(
        "Adj R2: ", signif(summary(fit)$adj.r.squared, 5),
        "<br>Intercept: ", signif(fit$coef[[1]], 5),
        "<br>Slope: ", signif(fit$coef[[2]], 5),
        "<br>P: ", 'if'(nrow(summary(fit)$coef) < 2 || 
                          ncol(summary(fit)$coef < 4),
                        NA,
                        signif(summary(fit)$coef[2, 4], 5)
        )
      )
      
      p <- plotdata %>%
        plotly::plot_ly(x = ~x, hoverinfo = "text", 
                        width = pwidth, height = totalheight) %>%
        plotly::add_markers(
          y = ~y,
          color = ~color,
          colors = getOption("me.colorvalues")[seq(1, length(levels(group)))],
          text = ~text,
          showlegend = showlegend
        ) %>%
        plotly::add_lines(
          y = ~x,
          line = list(
            dash = "dash",
            color = "rgba(220,220,220,1)"
          ),
          name = "onetoone",
          text = "",
          showlegend = FALSE
        ) 
      if(addRegression){
        p <- p %>%
          plotly::add_lines(
            y = stats::predict(fit),
            line = list(
              color = "rgba(119,0,0, 1)",
              shape = "spline"
            ),
            name = "Linear Regression",
            text = fitlabel,
            showlegend = FALSE
          ) %>%
          plotly::add_ribbons(
            data = broom::augment(fit, se_fit = TRUE),
            x = ~feature,
            ymin = ~ .fitted - 1.96 * .se.fit,
            ymax = ~ .fitted + 1.96 * .se.fit,
            line = list(color = "rgba(219,219,219, 0.35)"),
            fillcolor = "rgba(219,219,219, 0.3)",
            name = "Standard Error",
            showlegend = FALSE
          ) 
      }
      if(useWebGL)
        p <- p %>% plotly::toWebGL() 
      p <- p %>%
        add_plotly_layout(plotTitle = plotTitle, 
                          xaxis_text = paste0(xaxis_text, mS_label),  
                          ylab = paste0(logtext, yaxis_text)) %>%
        plotly::layout(margin = list(b = 150))
      p
    })
    plotly::subplot(plist, nrows = 1, shareY = TRUE, titleX = TRUE)
  })
  
  p <- plotly::subplot(plotlist, nrows = length(plotlist), 
               shareY = FALSE, titleX = TRUE, titleY = TRUE) %>%
    add_plotly_config()
  
  return(list(plot = p, stats = fitlist))
}


