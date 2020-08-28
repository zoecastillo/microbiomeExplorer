
## Modules for feature sample analysis

#' Relative abundance plot module - UI
#'
#' @param id namespace identifier
#'
#' @author Janina Reeder
#'
#' @return box containing the ui code
avgAbundanceUI <- function(id) {
  ns <- NS(id)
  
  box(
    title = "AVERAGE RELATIVE ABUNDANCE", solidHeader = TRUE, 
    collapsible = TRUE, width = 12,
    fluidRow(
      column(
        width = 11,
        shinycssloaders::withSpinner(
          plotly::plotlyOutput(ns("avgabplot"), 
                               width = "auto", 
                               height = "auto"),
          type = 3, color = "#424242", color.background = "#fdfdfc"),
        br(),
        shinyjs::disabled(shinyWidgets::dropdownButton(
          tags$h3("Plot Options"),
          shinyWidgets::radioGroupButtons(
            inputId = ns("avgAbChoice"),
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
#' @param featureSettings analysis input settings passed over to this module
#' @param normalizedData boolean indicating whether data has been normalized
#' @param reset boolean reactive which resets the module if TRUE
#' 
#'
#' @return list storing plot clicks and number of features displayed 
#' (passed to feature plot module) as well as the R code to make plot
avgAbundance <- function(input, output, session, 
                         aggDat, 
                         featLevel,
                         featureSettings,
                         normalizedData,
                         reset) {
  ns <- session$ns
  
  ## stores the plotly plot (needed to set source correctly)
  paPlot <- reactiveVal(NULL)
  ## stores the R code needed to build the relative abundance plot
  repCode <- reactiveVal(NULL)
  width <- reactiveVal(650)
  
  
  observe({
    req(reset())
    shinyjs::disable("optionbutton_state")
    shinyjs::enable("avgAbChoice")
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
    req(aggDat(), featureSettings(),numOfFeats())
    ylab <- "Reads"
    if(isFALSE(normalizedData())){
      shinyWidgets::updateRadioGroupButtons(session = session,
                              inputId = "avgAbChoice",
                              selected = "Reads")
      shinyjs::removeClass(id = "pspan", class = "hide")
      shinyjs::disable("avgAbChoice")
    } else {
      ylab <- input$avgAbChoice
      shinyjs::addClass(id = "pspan", class = "hide")
      shinyjs::enable("avgAbChoice")
    }
    facet1 <- featureSettings()$facetby
    facet2 <- featureSettings()$facetby2
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
    pa <- plotAvgAbundance(aggDat(),
                        level = featLevel(),
                        facet1 = facet1,
                        facet2 = facet2,
                        ind = seq_len(numOfFeats()), 
                        plotTitle = plotTitle,
                        ylab = ylab
    ) 
    pa$x$source <- "avgabplot"
    paPlot(pa)
    shinyjs::enable("optionbutton_state")
  })
  
  output$avgabplot <- plotly::renderPlotly({
    paPlot()
  })
  
  observeEvent(width(), {
    req(paPlot())
    plotly::plotlyProxy("avgabplot", session) %>%
      plotly::plotlyProxyInvoke(
        "relayout",
        list(width = width())
      )
  })
  

  
  ## update R code if inputs change
  observe({
    req(aggDat(), featureSettings())
    ylab <- input$avgAbChoice
    facet1 <- featureSettings()$facetby
    facet2 <- featureSettings()$facetby2
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
      paste0("#' ### Average relative abundance"),
      paste0("", "#- fig.width = 10"),  
      paste0("plotAvgAbundance(aggDat,"),
      paste0("\tlevel = \"", featLevel(), "\","),
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
      paste0("\tind = seq_len(", numOfFeats(), "),"),
      paste0("\tplotTitle = \"",plotTitle, "\","),
      paste0("\tylab = \"", ylab, "\")\n\n"),
      sep = "\n"
    ))
  })
  
  return(repCode = repCode)
}

