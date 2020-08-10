

## Modules for inter sample analysis tab

#' Abundance Heatmap module - UI
#'
#' @param id namespace identifier
#'
#' @author Janina Reeder
#'
#' @return box holding the UI code
abundanceHeatmapUI <- function(id) {
  ns <- NS(id)
  
  box(
    title = "ABUNDANCE HEATMAP", solidHeader = TRUE, 
    collapsible = TRUE, width = 12,
    fluidRow(
      column(
        width = 11,
        shinycssloaders::withSpinner(
          plotly::plotlyOutput(ns("abhmplot"), 
                               width = "auto", 
                               height = "800px"),
          type = 3, color = "#424242", 
          color.background = "#fdfdfc"),
        br(),
        shinyjs::disabled(
          shinyWidgets::dropdownButton(
            tags$h3("Plot Options"),
            numericInput(
              inputId = ns("hmfeats"),
              label = "Number of features",
              value = 50, min = 1
            ),
            selectInput(
              inputId = ns("hmcol"),
              label = "Phenotype annotations",
              choices = "",
              multiple = TRUE
            ),
            selectInput(
              inputId = ns("rowcol"),
              label = "Feature annotation",
              choices = "",
              multiple = FALSE
            ),
            shinyWidgets::switchInput(
              inputId = ns("hmlog"),
              label = "Log Scale",
              size = "mini",
              value = TRUE,
              labelWidth = "80px"
            ),
            actionButton(
              ns("changeHeatmapSettings"), 
              label = "GO", width = "50px"),
            circle = FALSE, status = "danger", 
            icon = icon("gear"), width = "300px",
            label = "Plot Options",
            inputId = ns("optionbutton")
          ))
      )
    )
  )
}

#' Abundance Heatmap module - server
#'
#' @param input shiny input
#' @param output shiny output
#' @param session shiny session
#' @param aggDat aggregated MRExperiment
#' @param featLevel chosen feature level (aggregation level)
#' @param colorOptions reactive storing filters selected via data input
#' @param levelOpts all available level choices for this dataset
#' @param hmSort reactive storing sorting method for heatmap
#' @param hmFeatList reactive storing list of features to include in heatmap
#' @param reset boolean reactive which resets the module if TRUE
#'
#' @author Janina Reeder
#' 
#' @return R code needed to generate the heatmap
abundanceHeatmap <- function(input, output, session, 
                             aggDat, 
                             featLevel, 
                             colorOptions, 
                             levelOpts, 
                             hmSort, 
                             hmFeatList,
                             reset) {
  ns <- session$ns
  
  ## update with available phenotype annotation options
  observe({
    if(is.null(colorOptions())){
      updateSelectInput(session, "hmcol", choices = "")
    } else {
      updateSelectInput(session, "hmcol", choices = c("",colorOptions()))
    }
  })
  
  ## update with available feature level annotation options
  observe({
    req(levelOpts(), featLevel())
    uptoind <- which(levelOpts() == featLevel())
    updateSelectInput(session, "rowcol", 
                      choices = c("", levelOpts()[seq_len(uptoind)]))
  })
  
  ## reatives storing plot options and whether plot update is required
  hmColors <- reactiveVal(NULL)
  rowColors <- reactiveVal(NULL)
  numOfFeats <- reactiveVal(50)
  logScale <- reactiveVal(TRUE)
  updatePlot <- reactiveVal(TRUE)
  ## stores the R code needed to build the heatmap
  repCode <- reactiveVal(NULL)
  
  observe({
    req(reset())
    shinyjs::disable("optionbutton_state")
    hmColors(NULL)
    rowColors(NULL)
    numOfFeats(50)
    updatePlot(TRUE)
    repCode(NULL)
    shinyjs::reset("abhmplot")
  })
  
  
  ## update plot options based on GO click
  observeEvent(input$changeHeatmapSettings,{
    changed <- FALSE
    updatePlot(FALSE)
    if(!is.null(input$hmcol) && !input$hmcol %in% hmColors()){
      hmColors(input$hmcol)
      changed <- TRUE
    }
    if(!is.null(input$rowcol) && !input$rowcol %in% rowColors()){
      rowColors(input$rowcol)
      changed <- TRUE
    }
    if(!is.null(input$hmfeats) && !input$hmfeats %in% numOfFeats()){
      numOfFeats(input$hmfeats)
      changed <- TRUE
    }
    if(!is.null(input$hmlog) && !input$hmlog %in% logScale()){
      logScale(input$hmlog)
      changed <- TRUE
    }
    updatePlot(changed)
  })
  
  ## render the heatmap
  output$abhmplot <- plotly::renderPlotly({
    req(aggDat(), hmSort(), updatePlot())
    shinyjs::enable("optionbutton_state")
    plotTitle <- 'if'(is.null(hmFeatList()),
                      paste0("Top ", numOfFeats(),
                             " features sorted by ",
                             hmSort()),
                      paste0("Selected features"))
    plotTitle <- paste0(plotTitle,
                        " at ",
                        featLevel(),
                        " level")
    repCode(paste(
      paste0("#' ### Feature Heatmap"),
      paste0("", "#- fig.width = 7, fig.height = 9"),
      paste0("plotHeatmap(aggDat,"),
      paste0("\tfeatures = ", 
             'if'(is.null(hmFeatList()), "NULL", 
                  paste0("c(",
                         paste0("\"", hmFeatList(), "\"",
                                collapse = ", "), ")")),
             ","),
      paste0("\tlog = ", logScale(), ","),
      paste0("\tsort_by = \"", hmSort(), "\","),
      paste0("\tnfeat = ", numOfFeats(), ","),
      paste0("\tcol_by = ", 
             'if'(is.null(hmColors()), "NULL", 
                  paste0("c(",
                         paste0("\"", hmColors(), "\"", 
                                collapse = ", "),")")),
             ","),
      paste0("\trow_by = ", 
             'if'(is.null(rowColors()), "NULL", 
                  paste0("\"", rowColors(), "\"")), ","),
      paste0("\tplotTitle = \"",plotTitle,"\")\n\n"),
      sep = "\n"
    ))
    
    plotHeatmap(
      aggdat = aggDat(),
      features = hmFeatList(),
      log = logScale(),
      sort_by = hmSort(),
      nfeat = isolate(numOfFeats()),
      col_by = isolate(hmColors()),
      row_by = isolate(rowColors()),
      plotTitle = plotTitle
    )
  })
  
  return(repCode)
}



#' Beta Diversity module - UI
#'
#' @param id namespace identifier
#'
#' @author Janina Reeder
#'
#' @return box holding the ui code
betaDiversityUI <- function(id) {
  ns <- NS(id)
  
  box(
    title = "BETA DIVERSITY", solidHeader = TRUE, 
    collapsible = TRUE, width = 12,
    fluidRow(
      column(
        width = 11,
        shinycssloaders::withSpinner(
          plotly::plotlyOutput(ns("betaDiv"), 
                               width = "auto",
                               height = "auto"), 
          type = 3, color = "#424242", 
          color.background = "#fdfdfc"),
        shinyjs::disabled(
          shinyWidgets::dropdownButton(
            tags$h3("Plot Options"),
            selectInput(
              inputId = ns("Xbeta"),
              label = "X",
              choices = c("PC1", "PC2", "PC3", "PC4", 
                          "PC5", "PC6", "PC7", "PC8"),
              selected = "PC1"
            ),
            selectInput(
              inputId = ns("Ybeta"),
              label = "Y",
              choices = c("PC1", "PC2", "PC3", "PC4", 
                          "PC5", "PC6", "PC7", "PC8"),
              selected = "PC2"
            ),
            selectInput(
              inputId = ns("betacol"),
              label = "Color by",
              choices = ""
            ),
            checkboxInput(inputId = ns("confEllipse"),
                          label = "Add confidence ellipse",
                          value = FALSE),
            sliderInput(inputId = ns("confLevel"),
                        label = "Confidence Level",
                        min = 0.01, max = 0.99,
                        step = 0.01, value = 0.95),
            selectInput(inputId = ns("betashape"),
                        label = "Shape by",
                        choices = ""
            ),
            numericInput(inputId = ns("betasize"),
                         label = "Point size",
                         value = 8, min = 1
            ),
            sliderInput(
              inputId = ns("plotWidth"),
              label = "Adjust plot width",
              value = 650,
              min = 250,
              max = 1600,
              round = TRUE
            ),
            actionButton(ns("changeBetaSettings"), 
                         label = "GO", width = "50px"),
            circle = FALSE, status = "danger", 
            icon = icon("gear"), width = "300px",
            label = "Plot Options",
            inputId = ns("optionbutton")
          ))
      )
    ),
    fluidRow(
      column(width = 1),
      column(width = 10, class = "statsrow",
             DT::DTOutput(ns("statsdatatable"), 
                          width = "90%", height = "auto")
      )
    )
  )
}


#' Beta Diversity module - server 
#'
#' @param input shiny input
#' @param output shiny output
#' @param session shiny session
#' @param aggDat MRExperiment storing data
#' @param aggLevel aggregation level
#' @param colorOptions phenotype selection options for color 
#' @param shapeOptions phenotype selection options for shape
#' @param betadistance distance measured used for beta diversity analysis
#' @param betaSettings input choices for beta diversity
#' @param reset boolean reactive which resets the module if TRUE
#'
#' @author Janina Reeder
#'
#' @return R code needed to generate the beta diversity plot
betaDiversity <- function(input, output, session, 
                          aggDat, 
                          aggLevel,
                          colorOptions, 
                          shapeOptions,
                          betadistance,
                          betaSettings,
                          reset) {
  ns <- session$ns
  
  ## reatives storing plot options and whether plot update is required
  shapeChoice <- reactiveVal("None")
  colorChoice <- reactiveVal("No Color")
  sizeChoice <- reactiveVal(8)
  xbeta <- reactiveVal("PC1")
  ybeta <- reactiveVal("PC2")
  updatePlot <- reactiveVal(FALSE)
  ## stores the R code needed to build the plot
  repCode <- reactiveVal(NULL)
  ## reactive storing distance matrices (to be pulled from DB)
  distMat <- reactiveVal(NULL)
  pcaVals <- reactiveVal(NULL)
  
  adonisVar <- reactiveVal(NULL)
  adonisStrata <- reactiveVal(NULL)
  adonisText <- reactiveVal(NULL)
  adonisCode <- reactiveVal(NULL)
  plotComplete <- reactiveVal(FALSE)
  confInterval <- reactiveVal(NULL)
  
  observe({
    req(reset())
    shinyjs::disable("optionbutton_state")
    shapeChoice("None")
    colorChoice("No Color")
    sizeChoice(8)
    xbeta("PC1")
    ybeta("PC2")
    updatePlot(FALSE)
    repCode(NULL)
    distMat(NULL)
    pcaVals(NULL)
    adonisVar(NULL)
    adonisStrata(NULL)
    adonisText(NULL)
    adonisCode(NULL)
    plotComplete(FALSE)
  })
  
  ## update available choices for color and shape
  observe({
    updateSelectInput(session, "betacol", 
                      choices = c("No Color", colorOptions()))
    updateSelectInput(session, "betashape", 
                      choices = c("None",shapeOptions()))
  })
  
  observeEvent(input$betacol,{
    if(input$betacol == "No Color"){
      shinyjs::disable("confEllipse")
      shinyjs::disable("confLevel")
    } else {
      shinyjs::enable("confEllipse")
    }
  })
  
  observeEvent(input$optionbutton_state,{
    req(input$betacol == "No Color")
    shinyjs::disable("confEllipse")
    shinyjs::disable("confLevel")
  })
  
  observe({
    req(input$optionbutton_state)
    if(input$confEllipse){
      shinyjs::enable("confLevel")
    } else {
      shinyjs::disable("confLevel")
    }
  })
  
  observe({
    req(input$betacol != "No Color", input$confEllipse)
    confInterval(input$confLevel)
  })
  
  ## update adonisVar and adonisStrata if inputs change
  observeEvent(betaSettings()$adonisvar,{
    req(aggDat())
    if(!is.null(betaSettings()$adonisvar) && betaSettings()$adonisvar != ""){
      phenoData <- pData(aggDat())
      ## check if columns contains NA
      varNA <- (sum(is.na(
        phenoData[,betaSettings()$adonisvar])) |
          sum(phenoData[,betaSettings()$adonisvar] %in% "")) > 0
      strataNA <- FALSE
      if(!is.null(betaSettings()$adonisstrata)){
        strataNA <- (sum(is.na(
          phenoData[,betaSettings()$adonisstrata])) |
            sum(phenoData[,betaSettings()$adonisstrata] %in% "")) > 0
      }
      ## in case NAs are present, we need to check in with user
      if(varNA | strataNA){
        showModal(modalDialog(
          title = "Adonis Variables contains NA",
          paste0("NA values are present in the adonis variables:\n\n",
                 'if'(
                   varNA, 
                   paste0("\n\t",
                          betaSettings()$adonisvar),
                   ""),
                 'if'(
                   strataNA,
                   paste0("\n\t",
                          betaSettings()$adonisstrata),
                   ""),
                 "\nTo compute adonis statistics, all NAs need to be removed
                 or treated as \"NA\" strings. \n\n We advise to use adonis
                variables not containing NA values.\n\n
                Do you want to replace NA with \"NA\" or abort to review?"),
          footer = tagList(actionButton(ns("confirmReplace"), 
                                        "Replace"),
                           actionButton(ns("confirmReview"), 
                                        "Abort")
          ),
          easyClose = FALSE))
      } else {
        if(!varNA){
          adonisVar(phenoData[,betaSettings()$adonisvar]) 
        }
        if(!strataNA){
          if(!is.null(betaSettings()$adonisstrata)){
            adonisStrata(
              phenoData[,betaSettings()$adonisstrata])
          } else {
            adonisStrata(NULL)
          }
        }
      }
    } else {
      adonisText(NULL)
      adonisVar(NULL)
      adonisStrata(NULL)
      adonisCode(NULL)
    }
  })
  
  observeEvent(input$confirmReview,{
    removeModal()
    adonisVar(NULL)
    adonisStrata(NULL)
  })    
  
  
  observeEvent(input$confirmReplace,{
    removeModal()
    phenoData <- pData(aggDat())
    adonisVar('if'(
      is.na(phenoData[,betaSettings()$adonisvar]),
      "NA",
      phenoData[,betaSettings()$adonisvar]))
    
    if(!is.null(betaSettings()$adonisstrata)){
      adonisStrata('if'(
        is.na(
          phenoData[,betaSettings()$adonisstrata]),
        "NA",
        phenoData[,betaSettings()$adonisstrata]))
    }
  })
  
  observeEvent(c(aggDat(), betadistance()),{
    req(aggDat(), betadistance() != "")
    if(!(betadistance() %in% names(distMat()))){
      withProgress({
        newDist <- list(computeDistMat(aggDat(),betadistance()))
        incProgress(0.5)
        oldDists <- distMat()
        oldDists[[betadistance()]] <- newDist
        distMat(oldDists)
      }, message = "Calculating distance matrix and PCAs")
    }
    pcaVals(calculatePCAs(distMat()[[betadistance()]][[1]], 
                          c(xbeta(),ybeta())))
  }, priority = 20)
  
  output$statsdatatable <- DT::renderDT({
    req(aggDat(),distMat(), adonisVar())
    if(!is.null(betaSettings()$adonisvar) && betaSettings()$adonisvar != ""){ 
      beta_dis <- distMat()[[betadistance()]][[1]]
      req(labels(beta_dis) == rownames(pData(aggDat())))
      phenoData <- pData(aggDat())
      phenoData[,betaSettings()$adonisvar] <-
        adonisVar()
      
      formula <- stats::as.formula(paste("beta_dis ~ ",
                                         betaSettings()$adonisvar))
      x <- tryCatch(
        vegan::adonis(formula = formula, 
                      data = phenoData, 
                      strata = isolate(adonisStrata()),
                      permutations=500, 
                      parallel = 1),
        error = function(e){
          showModal(modalDialog(
            title = "Error running adonis",
            e,
            easyClose = TRUE
          ))
          return(NULL)
        }
      )
      
      tableCaption <- paste0("Adonis variance of ",
                             isolate(betaSettings()$adonisvar))
      if(!is.null(adonisStrata()))
        tableCaption <- paste0(tableCaption,
                               " with strata ",
                               isolate(betaSettings()$adonisstrata))
      
      adonisCode(paste(
        paste0("\nphenoData <- pData(aggDat)"),
        paste0("phenoData[,\"",betaSettings()$adonisvar,
               "\"] <- c(", paste0(
                 paste0("\"",adonisVar(),"\""),
                 collapse = ","),")"),
        paste0("formula <- as.formula(distMat ~ ",
               betaSettings()$adonisvar,")\n"),
        paste0("x <- vegan::adonis(formula = formula,"), 
        paste0("\tdata = phenoData,"), 
        paste0("\tstrata = ", 'if'(is.null(adonisStrata()), "NULL",
                                   paste0("c(",paste0(
                                     paste0("\"",adonisStrata(),"\""),
                                     collapse = ","),")")),","),
        paste0("\tpermutations=500,"), 
        paste0("\tparallel = 1)"),
        paste0("adonisData <- as.data.frame(x$aov.tab)"),
        paste0("adonisData[] <- sapply(adonisData, as.numeric)"),
        paste0("adonisData[] <- sapply(adonisData, round, digits = getOption(\"me.round_digits\"))"),
        paste0("\n\nif(doctype == \"html\"){"),
        paste0("\tDT::datatable(data = adonisData, 
                   class = \"stripe hover cell-border order-column\","),
        paste0("\t\tfilter = \"none\", style = \"bootstrap\","),
        paste0("\t\tcaption = \"",tableCaption,"\","),
        paste0("\t\toptions = list(scrollX = TRUE, paging = FALSE, 
                   digits = 4,dom = \"<tp>\"),"),
        paste0("\t\tescape = FALSE)"),
        paste0("} else {"),
        paste0("\tkable(adonisData)"),
        paste0("}\n\n"),
        sep = "\n"))
      
      
      adonisData <- as.data.frame(x$aov.tab)
      adonisData[] <- vapply(adonisData, as.numeric, numeric(3))
      adonisData[] <- vapply(adonisData, round, 
                             digits = getOption("me.round_digits"),
                             numeric(3))
      
      adonisText(paste0("R2: ", adonisData[1,"R2"],
                        "; Pr(>F): ", adonisData[1,"Pr(>F)"]))
      
      
      DT::datatable(
        data = adonisData, 
        class = "stripe hover cell-border order-column",
        filter = "none", style = "bootstrap",
        caption = tableCaption,
        options = list(
          scrollX = TRUE,
          paging = FALSE,
          digits = 4,
          dom = "<tp>"
        ),
        escape = FALSE
      )
    } else {
      adonisText(NULL)
      return (NULL)
    }
  })
  
  ## update plot options based on GO click
  observeEvent(input$changeBetaSettings,{
    changed <- FALSE
    updatePlot(FALSE)
    if(!input$Xbeta %in% xbeta()){
      xbeta(input$Xbeta)
      changed <- TRUE
    }
    if(!input$Ybeta %in% ybeta()){
      ybeta(input$Ybeta)
      changed <- TRUE
    }
    if(changed){
      req(distMat(), betadistance() != "")
      pcaVals(calculatePCAs(distMat()[[betadistance()]][[1]], 
                            c(xbeta(),ybeta())))
    }
    if(!input$betacol %in% colorChoice()){
      colorChoice(input$betacol)
      changed <- TRUE
    }
    if(!input$betashape %in% shapeChoice()){
      shapeChoice(input$betashape)
      changed <- TRUE
    }
    if(!input$betasize %in% sizeChoice()){
      sizeChoice(input$betasize)
      changed <- TRUE
    }
    if(input$confEllipse && !(input$confLevel != confInterval())){
      confInterval(input$confLevel)
      changed <- TRUE
    }
    if(!input$confEllipse && !is.null(confInterval())){
      confInterval(NULL)
      changed <- TRUE
    }
    updatePlot(changed)
  })
  
  observe({
    req(plotComplete())
    aT <- adonisText()
    if(is.null(aT))
      aT <- ""
    plotly::plotlyProxy("betaDiv", session) %>%
      plotly::plotlyProxyInvoke(
        "restyle",
        list(text = aT),
        0
      )
    shinyjs::js$resetAxes()
  })
  
  observeEvent(input$plotWidth, {
    req(pcaVals(), plotComplete())
    plotly::plotlyProxy("betaDiv", session) %>%
      plotly::plotlyProxyInvoke(
        "relayout",
        list(width = input$plotWidth)
      )
  })
  
  ## calls beta diversity plotting function
  output$betaDiv <- plotly::renderPlotly({
    req(pcaVals())
    shinyjs::enable("optionbutton_state")
    updatePlot()
    plotComplete(FALSE)
    plotTitle <- paste0(betadistance(),
                        " diversity at ", 
                        aggLevel()," level")
    color <- as.character(isolate(colorChoice()))
    if (stringr::str_detect(color, "No Color")) {
      color <- NULL
    }
    shape <- as.character(isolate(shapeChoice()))
    if (stringr::str_detect(shape, "None")) {
      shape <- NULL
    }
    
    
    pb <- plotBeta(
      aggdat = aggDat(),
      dist_method = tolower(betadistance()),
      pcas = isolate(pcaVals()),
      dim = isolate(c(xbeta(), ybeta())),
      col_by = color,
      shape_by = shape,
      plotTitle = plotTitle,
      pt_size = isolate(sizeChoice()),
      plotText = isolate(adonisText()),
      confInterval = isolate(confInterval())
    )
    shinyjs::delay(500, plotComplete(TRUE))
    pb
  })
  
  ## update the stored R code based on input choices
  observe({
    req(pcaVals())
    color <- as.character(colorChoice())
    if (stringr::str_detect(color, "No Color")) {
      color <- NULL
    }
    shape <- as.character(shapeChoice())
    if (stringr::str_detect(shape, "None")) {
      shape <- NULL
    }
    plotTitle <- paste0(betadistance(),
                        " diversity at ", 
                        aggLevel()," level")
    
    repCode(paste(
      paste0("#' ### Beta diversity"),
      paste0("", "#- fig.width = 7"),
      paste0("distMat <- computeDistMat(aggDat, \"", 
             tolower(betadistance()), "\")"),
      paste0("pcaVals <- calculatePCAs(distMat, 
                      c(\"", xbeta(), "\", \"", ybeta(), "\"))\n"),
      paste0("plotBeta(aggDat,"),
      paste0("\tdist_method = \"", tolower(betadistance()), "\","),
      paste0("\tpcas = pcaVals,"),
      paste0("\tdim = c(\"", xbeta(), "\", \"", ybeta(), "\"),"),
      paste0("\tcol_by = ", 
             'if'(is.null(color), "NULL", 
                  paste0("\"", color, "\"")), ","),
      paste0("\tshape_by = ", 
             'if'(is.null(shape), "NULL", 
                  paste0("\"", shape, "\"")), ","),
      paste0("\tplotTitle = \"",plotTitle, "\","),
      paste0("\tpt_size = \"",sizeChoice(), "\","),
      paste0("\tplotText = \"",adonisText(), "\","),
      paste0("\tconfInterval = ",confInterval(), ","),
      paste0("\tallowWebGL = FALSE)\n\n"),
      adonisCode(),"\n\n",
      sep = "\n"
    ))
  })
  
  return(repCode)
}
