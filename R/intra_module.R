
## Modules for intra sample analysis

#' Relative abundance plot module - UI
#'
#' @param id namespace identifier
#'
#' @author Janina Reeder
#'
#' @return box containing the ui code
relAbundanceUI <- function(id) {
  ns <- NS(id)
  
  box(
    title = "RELATIVE ABUNDANCE", solidHeader = TRUE, 
    collapsible = TRUE, width = 12,
    fluidRow(
      column(
        width = 11,
        shinycssloaders::withSpinner(
          plotly::plotlyOutput(
            ns("relabplot"), 
            width = "auto", 
            height = "auto"),
          type = 3, color = "#424242", color.background = "#fdfdfc"),
        br(),
        shinyjs::disabled(
          shinyWidgets::dropdownButton(
            tags$h3("Plot Options"),
            shinyWidgets::radioGroupButtons(
              inputId = ns("relAbChoice"),
              label = "Y axis",
              choices = c("Percentage", "Reads"),
              individual = TRUE,
              checkIcon = list(
                yes = tags$i(
                  class = "fa fa-circle",
                  style = "color: steelblue"
                ),
                no = tags$i(
                  class = "fa fa-circle-o",
                  style = "color: steelblue"
                )
              )
            ),
            div(id=ns("pspan"), class = "percentage",
                p("Normalization is required to show percentage")
            ),
            span(
              numericInput(
                inputId = ns("numofFeatures"),
                label = "Max number of Features to show",
                value = 10,
                min = 1
              ),
              sliderInput(
                inputId = ns("plotWidth"),
                label = "Adjust plot width",
                value = 650,
                min = 250,
                max = 1600,
                round = TRUE
              )
            ),
            circle = FALSE, status = "danger",
            icon = icon("gear"), width = "300px",
            label = "Plot Options",
            inputId = ns("optionbutton")
          ))
      )
    )
  )
}

#' Relative abundance plot module - server
#'
#' @param input shiny input
#' @param output shiny output
#' @param session shiny session
#' @param aggDat aggregated MRExperiment
#' @param featLevel chosen feature level (aggregation level)
#' @param intraSettings analysis input settings passed over to this module
#' @param normalizedData boolean indicating whether data has been normalized
#' @param reset boolean reactive which resets the module if TRUE
#' 
#'
#' @return list storing plot clicks and number of features displayed 
#' (passed to feature plot module) as well as the R code to make plot
relAbundance <- function(input, output, session, 
                         aggDat, 
                         featLevel, 
                         intraSettings,
                         normalizedData,
                         reset) {
  ns <- session$ns
  
  ## stores the plotly plot (needed to set source correctly)
  paPlot <- reactiveVal(NULL)
  ## stores the R code needed to build the relative abundance plot
  repCode <- reactiveVal(NULL)
  width <- reactiveVal(650)
  
  clickedFeature <- reactive(input$clickedFeature)
  
  observe({
    req(reset())
    shinyjs::disable("optionbutton_state")
    shinyjs::enable("relAbChoice")
    paPlot(NULL)
    repCode(NULL)
  })
  
  observe({
    width(input$plotWidth)
  })
  
  
  ## debouncing number of features to avoid multiple redraws
  nof <- reactive(input$numofFeatures)
  numOfFeats <- debounce(nof,1000)
  
  
  ## calls the plotting function
  observe({
    req(aggDat(), intraSettings(),numOfFeats())
    ylab <- "Reads"
    if(isFALSE(normalizedData())){
      shinyWidgets::updateRadioGroupButtons(session = session,
                                            inputId = "relAbChoice",
                                            selected = "Reads")
      shinyjs::removeClass(id = "pspan", class = "hide")
      shinyjs::disable("relAbChoice")
    } else {
      ylab <- input$relAbChoice
      shinyjs::addClass(id = "pspan", class = "hide")
      shinyjs::enable("relAbChoice")
    }
    facet1 <- intraSettings()$facetby
    facet2 <- intraSettings()$facetby2
    plotTitle <- paste0("Top ",
                        numOfFeats(),
                        " feature ", 
                        tolower(gsub("Reads","abundance",ylab)),
                        " at ",
                        featLevel()," level")
    if(!is.null(facet1)){
      plotTitle <- paste0(plotTitle,
                          " split by ",
                          facet1)
    } else if(!is.null(facet2)){
      plotTitle <- paste0(plotTitle,
                          " split by ",
                          facet2)
    }
    if(!is.null(facet1) & !is.null(facet2)){
      plotTitle <- paste0(plotTitle,
                          " and ",
                          facet2)
    }
    pa <- plotAbundance(aggDat(),
                        level = featLevel(),
                        x_var = intraSettings()$xvariable,
                        facet1 = facet1,
                        facet2 = facet2,
                        ind = seq_len(numOfFeats()), 
                        plotTitle = plotTitle,
                        ylab = ylab
    ) 
    pa$x$source <- "relabplot"
    paPlot(pa)
    shinyjs::enable("optionbutton_state")
  })
  
  output$relabplot <- plotly::renderPlotly({
    paPlot()
  })
  
  observeEvent(width(), {
    req(paPlot())
    plotly::plotlyProxy("relabplot", session) %>%
      plotly::plotlyProxyInvoke(
        "relayout",
        list(width = width())
      )
  })
  
  
  ## reactive stores plotly click events
  observe({
    pe <- plotly::event_data("plotly_click", source = "relabplot")
    ## get the name of the trace that was clicked (= feature name)
    shinyjs::js$getTraceName(isolate(paPlot()), pe$curveNumber)
  })
  
  ## update R code if inputs change
  observe({
    req(aggDat(), intraSettings())
    ylab <- input$relAbChoice
    facet1 <- intraSettings()$facetby
    facet2 <- intraSettings()$facetby2
    plotTitle <- paste0("Top ",
                        numOfFeats(),
                        " feature ", 
                        tolower(gsub("Reads","abundance",ylab)),
                        " at ",
                        featLevel()," level")
    if(!is.null(facet1)){
      plotTitle <- paste0(plotTitle,
                          " split by ",
                          facet1)
    } else if(!is.null(facet2)){
      plotTitle <- paste0(plotTitle,
                          " split by ",
                          facet2)
    }
    if(!is.null(facet1) & !is.null(facet2)){
      plotTitle <- paste0(plotTitle,
                          " and ",
                          facet2)
    }
    repCode(paste(
      paste0("#' ### Relative Abundance"),
      paste0("", "#- fig.width = 10"),  
      paste0("plotAbundance(aggDat,"),
      paste0("\tlevel = \"", featLevel(), "\","),
      paste0("\tx_var = \"", 
             intraSettings()$xvariable, "\","),
      paste0(
        "\tfacet1 = ", 
        'if'(is.null(facet1),"NULL",paste0("\"", facet1, "\"")),
        ","
      ),
      paste0(
        "\tfacet2 = ", 
        'if'(is.null(facet2),"NULL",paste0("\"",facet2, "\"")),
        ","
      ),
      paste0("\tind = 1:", numOfFeats(), ","),
      paste0("\tplotTitle = \"",plotTitle, "\","),
      paste0("\tylab = \"", ylab, "\")\n\n"),
      sep = "\n"
    ))
  })
  
  return(list(featName = clickedFeature, 
              ylabMode = reactive(input$relAbChoice),
              numOfFeats = numOfFeats, 
              repCode = repCode))
}


#' Feature plot module - UI
#'
#' @param id namespace identifier
#'
#' @return box holding the UI code
featAbundanceUI <- function(id) {
  ns <- NS(id)
  
  box(title = "FEATURE PLOT", solidHeader = TRUE, 
      collapsible = TRUE, width = 12,
      fluidRow(
        column(
          width = 11,
          shinycssloaders::withSpinner(
            plotly::plotlyOutput(ns("featplot"), 
                                 width = "auto", 
                                 height = "auto"),
            type = 3, color = "#424242", 
            color.background = "#fdfdfc"),
          br(),
          shinyjs::disabled(
            shinyWidgets::dropdownButton(
              tags$h3("Plot Options"),
              shinyWidgets::radioGroupButtons(
                inputId = ns("featChoice"),
                label = "Y axis",
                choices = c("Reads", "Percentage"),
                individual = TRUE,
                checkIcon = list(
                  yes = tags$i(
                    class = "fa fa-circle",
                    style = "color: steelblue"
                  ),
                  no = tags$i(
                    class = "fa fa-circle-o",
                    style = "color: steelblue"
                  )
                )
              ),          
              div(id=ns("pspan"), class = "percentage",
                  p("Normalization is required to show percentage")
              ),
              shinyWidgets::switchInput(
                inputId = ns("sP"),
                label = "Show Points",
                value = TRUE,
                size = "mini",
                labelWidth = "80px"
              ),
              shinyWidgets::switchInput(
                inputId = ns("logS"),
                label = "Log Scale",
                value = TRUE,
                size = "mini",
                labelWidth = "80px"
              ),
              sliderInput(
                inputId = ns("plotWidth"),
                label = "Adjust plot width",
                value = 650,
                min = 250,
                max = 1600,
                round = TRUE
              ),
              circle = FALSE, status = "danger", 
              icon = icon("gear"), width = "300px",
              label = "Plot Options",
              inputId = ns("optionbutton")
            ))
        )
      )
  )
}




#' Feature plot module - server
#'
#' @param input shiny input
#' @param output shiny output
#' @param session shiny session
#' @param aggDat aggregated MRExperiment
#' @param featLevel chosen feature level (aggregation level)
#' @param intraSettings analysis settings passed over from analysis input module
#' @param selectedFeat feature selected via drop down element of analysis input
#' @param featName plotly click event passed via relative abundance
#' @param numOfFeats number of features shown in relative abundance plot 
#' (affects plotly click data)
#' @param ylabMode character indication if raw \"Reads\" or \"Percentage\" 
#' should be shown
#' @param normalizedData boolean indicating whether data has been normalized
#' @param reset boolean reactive which resets the module if TRUE
#' 
#'
#' @author Janina Reeder
#' 
#' @return R code needed to build the feature plot
featAbundance <- function(input, output, session, 
                          aggDat, 
                          featLevel, 
                          intraSettings, 
                          selectedFeat, 
                          featName, 
                          numOfFeats,
                          ylabMode,
                          normalizedData,
                          reset) {
  ns <- session$ns
  
  ## selected feature for plot is updated both via dropdown (selectedFeat) 
  ## as well as plot click (featName)
  selectedFeature <- reactiveVal(NULL)
  ## stores R code needed to make the feature Plot
  repCode <- reactiveVal(NULL)
  width <- reactiveVal(650)
  
  observe({
    req(reset())
    shinyjs::disable("optionbutton_state")
    shinyjs::enable("featChoice")
    selectedFeature(NULL)
    repCode(NULL)
  })
  
  observe({
    width(input$plotWidth)
  })
  
  
  observe({
    req(ylabMode())
    shinyWidgets::updateRadioGroupButtons(
      session,"featChoice", selected = ylabMode())
  })
  
  observe({
    req(selectedFeat())
    selectedFeature(selectedFeat())
  })
  
  observe({
    selectedFeature(featName())
  })
  
  ## add/remove points using plotlyProxy to avoid redrawing plot
  observeEvent(input$sP, {
    req(aggDat(), intraSettings()$xvariable, 
        featLevel(), selectedFeature(), numOfFeats())
    show_points <- 'if'(input$sP, "all", FALSE)
    plotly::plotlyProxy("featplot", session) %>%
      plotly::plotlyProxyInvoke(
        "restyle",
        list(boxpoints = show_points)
      )
  })
  
  observeEvent(width(), {
    req(aggDat(), intraSettings()$xvariable, 
        featLevel(), selectedFeature(), numOfFeats())
    plotly::plotlyProxy("featplot", session) %>%
      plotly::plotlyProxyInvoke(
        "relayout",
        list(width = width())
      )
  })
  
  
  ## calls function to do a feature abundance plot
  output$featplot <- plotly::renderPlotly({
    req(aggDat(), intraSettings()$xvariable, 
        featLevel(), selectedFeature(), numOfFeats())
    shinyjs::enable("optionbutton_state")
    ylab <- "Reads"
    if(!normalizedData()){
      shinyWidgets::updateRadioGroupButtons(session = session,
                                            inputId = "featChoice",
                                            selected = "Reads")
      shinyjs::removeClass(id = "pspan", class = "hide")
      shinyjs::disable("featChoice")
    } else {
      ylab <- input$featChoice
      shinyjs::addClass(id = "pspan", class = "hide")
      shinyjs::enable("featChoice")
    }
    facet1 <- intraSettings()$facetby
    facet2 <- intraSettings()$facetby2
    plotTitle <- paste0(gsub("Reads","Abundance",ylab),
                        " of ",
                        selectedFeature())
    if(!is.null(facet1)){
      plotTitle <- paste0(plotTitle,
                          " split by ",
                          facet1)
    } else if(!is.null(facet2)){
      plotTitle <- paste0(plotTitle,
                          " split by ",
                          facet2)
    }
    if(!is.null(facet1) & !is.null(facet2)){
      plotTitle <- paste0(plotTitle,
                          " and ",
                          facet2)
    }
    plotSingleFeature(aggDat(),
                      x_var = intraSettings()$xvariable,
                      ind = seq_len(numOfFeats()),
                      plotTitle = plotTitle,
                      facet1 = intraSettings()$facetby,
                      facet2 = intraSettings()$facetby2,
                      feature = selectedFeature(),
                      ylab = ylab,
                      log = input$logS,
                      showPoints = isolate(input$sP)
    )
  })
  
  ## update R code needed to make plot
  observe({
    req(aggDat(), intraSettings()$xvariable, 
        featLevel(), selectedFeature(), numOfFeats())
    facet1 <- intraSettings()$facetby
    facet2 <- intraSettings()$facetby2
    plotTitle <- paste0(gsub("Reads","Abundance",input$featChoice),
                        " of ",
                        selectedFeature())
    if(!is.null(facet1)){
      plotTitle <- paste0(plotTitle,
                          " split by ",
                          facet1)
    } else if(!is.null(facet2)){
      plotTitle <- paste0(plotTitle,
                          " split by ",
                          facet2)
    }
    if(!is.null(facet1) & !is.null(facet2)){
      plotTitle <- paste0(plotTitle,
                          " and ",
                          facet2)
    }
    repCode(paste(
      paste0("#' ### Feature plot"),
      paste0("", "#- fig.width = 8"),
      paste0("plotSingleFeature(aggDat,"),
      paste0("\tx_var = \"", 
             intraSettings()$xvariable, "\","),
      paste0("\tind = 1:", numOfFeats(), ","),
      paste0("\tplotTitle = \"",plotTitle, "\","),
      paste0(
        "\tfacet1 = ", 
        'if'(is.null(intraSettings()$facetby),
             "NULL",
             paste0("\"", 
                    intraSettings()$facetby, "\"")
        ),
        ","
      ),
      paste0(
        "\tfacet2 = ", 
        'if'(is.null(intraSettings()$facetby2),
             "NULL",
             paste0("\"", 
                    intraSettings()$facetby2, "\"")
        ),
        ","
      ),
      paste0("\tfeature = \"", selectedFeature(), "\","),
      paste0("\tylab = \"", input$featChoice, "\","),
      paste0("\tlog = ", 
             'if'(is.null(input$log), "TRUE", input$log), ","),
      paste0("\tshowPoints = ", input$sP, ")\n\n"),
      sep = "\n"
    ))
  })
  
  return(repCode)
}

#' Alpha Diversity module - UI
#'
#' @param id namespace identifier
#'
#' @author Janina Reeder
#'
#' @return box holding the UI code
alphaDiversityUI <- function(id) {
  ns <- NS(id)
  
  box(title = "ALPHA DIVERSITY", solidHeader = TRUE, 
      collapsible = TRUE, width = 12,
      fluidRow(
        column(
          width = 11,
          shinycssloaders::withSpinner(
            plotly::plotlyOutput(
              ns("alphaDiv"), 
              width = "auto", 
              height = "auto"),
            type = 3, color = "#424242", 
            color.background = "#fdfdfc"),
          br(),
          shinyjs::disabled(
            shinyWidgets::dropdownButton(
              tags$h3("Plot Options"),
              selectInput(
                inputId = ns("alphaInd"),
                label = "Index",
                choices = c(
                  "Shannon", "Simpson", "Inverse Simpson",
                  "Richness"
                )
              ),
              selectInput(
                inputId = ns("alphacol"),
                label = "Color by",
                choices = ""
              ),
              sliderInput(
                inputId = ns("plotWidth"),
                label = "Adjust plot width",
                value = 650,
                min = 250,
                max = 1600,
                round = TRUE
              ),
              circle = FALSE, status = "danger", 
              icon = icon("gear"), width = "300px",
              label = "Plot Options",
              inputId = ns("optionbutton")
            ))
        )
      )
  )
}

#' Alpha Diversity module - server
#'
#' @param input shiny input
#' @param output shiny output
#' @param session shiny session
#' @param aggDat aggregated MRExperiment
#' @param featLevel chosen feature level (aggregation level)
#' @param intraSettings analysis settings as passed over from analysis 
#' input module
#' @param colorOptions phenotype selections: used for color choices
#' @param reset boolean reactive which resets the module if TRUE
#'
#' @author Janina Reeder
#' 
#' @return R code used to make the alpha diversity plot
alphaDiversity <- function(input, output, session, 
                           aggDat, 
                           featLevel, 
                           intraSettings, 
                           colorOptions,
                           reset) {
  ns <- session$ns
  
  width <- reactiveVal(650)
  ## stores index to pass on to alpha diversity calculation
  chosenIndex <- reactive({
    switch(input$alphaInd,
           "Shannon" = "shannon",
           "Simpson" = "simpson",
           "Inverse Simpson" = "invsimpson",
           "Richness" = "richness"
    )
  })
  
  ## handles color selection
  colorChoice <- reactiveVal("No Color")
  ## stores R code needed to build the plot
  repCode <- reactiveVal(NULL)
  
  observe({
    req(reset())
    shinyjs::disable("optionbutton_state")
    repCode(NULL)
    colorChoice("No Color")
  })
  
  observe({
    width(input$plotWidth)
  })
  
  ## update color options available
  observe({
    updateSelectInput(session, "alphacol", 
                      choices = c("No Color", colorOptions()))
  })
  
  observeEvent(input$alphacol, {
    colorChoice(input$alphacol)
  })
  
  observeEvent(width(), {
    req(aggDat(), featLevel(), chosenIndex(), 
        intraSettings(), colorChoice())
    plotly::plotlyProxy("alphaDiv", session) %>%
      plotly::plotlyProxyInvoke(
        "relayout",
        list(width = width())
      )
  })
  
  ## calls the alpha diversity plotting function
  output$alphaDiv <- plotly::renderPlotly({
    req(aggDat(), featLevel(), chosenIndex(), 
        intraSettings(), colorChoice())
    shinyjs::enable("optionbutton_state")
    facet1 <- intraSettings()$facetby
    facet2 <- intraSettings()$facetby2
    plotTitle <- paste0(stringr::str_to_sentence(chosenIndex()),
                        " diversity index at ", 
                        featLevel()," level")
    if(!is.null(facet1)){
      plotTitle <- paste0(plotTitle,
                          " split by ",
                          facet1)
    } else if(!is.null(facet2)){
      plotTitle <- paste0(plotTitle,
                          " split by ",
                          facet2)
    }
    if(!is.null(facet1) & !is.null(facet2)){
      plotTitle <- paste0(plotTitle,
                          " and ",
                          facet2)
    }
    color <- as.character(colorChoice())
    if (stringr::str_detect(color, "No Color")) {
      color <- NULL
    }
    plotAlpha(aggDat(),
              level = featLevel(),
              index = chosenIndex(),
              x_var = intraSettings()$xvariable,
              facet1 = intraSettings()$facetby,
              facet2 = intraSettings()$facetby2,
              col_by = color,
              plotTitle = plotTitle
    )
  })
  
  ## update R code based on input choices
  observe({
    req(aggDat(), featLevel(), chosenIndex(), 
        intraSettings(), colorChoice())
    color <- as.character(colorChoice())
    if (stringr::str_detect(color, "No Color")) {
      color <- NULL
    }   
    facet1 <- intraSettings()$facetby
    facet2 <- intraSettings()$facetby2
    plotTitle <- paste0(stringr::str_to_sentence(chosenIndex()),
                        " diversity index at ", 
                        featLevel()," level")
    if(!is.null(facet1)){
      plotTitle <- paste0(plotTitle,
                          " split by ",
                          facet1)
    } else if(!is.null(facet2)){
      plotTitle <- paste0(plotTitle,
                          " split by ",
                          facet2)
    }
    if(!is.null(facet1) & !is.null(facet2)){
      plotTitle <- paste0(plotTitle,
                          " and ",
                          facet2)
    }
    color <- as.character(colorChoice())
    if (stringr::str_detect(color, "No Color")) {
      color <- NULL
    }
    repCode(paste(
      paste0("#' ### Alpha diversity"),
      paste0("", "#- fig.width = 7"),
      paste0("plotAlpha(aggDat,"),
      paste0("\tlevel = \"", featLevel(), "\","),
      paste0("\tindex = \"", chosenIndex(), "\","),
      paste0("\tx_var = \"", 
             intraSettings()$xvariable, "\","),
      paste0(
        "\tfacet1 = ", 
        'if'(is.null(intraSettings()$facetby),
             "NULL",
             paste0("\"", 
                    intraSettings()$facetby, "\"")
        ),
        ","
      ),
      paste0(
        "\tfacet2 = ", 
        'if'(is.null(intraSettings()$facetby2),
             "NULL",
             paste0("\"", 
                    intraSettings()$facetby2, "\"")
        ),
        ","
      ),
      paste0("\tcol_by = ", 
             'if'(is.null(color), "NULL", 
                  paste0("\"", color, "\"")), ","),
      paste0("\tplotTitle = \"",plotTitle, "\")\n\n"),
      sep = "\n"
    ))
  })
  
  return(repCode)
}
