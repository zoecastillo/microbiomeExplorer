

#' Longitudinal Analysis module UI
#'
#' @param id namespace identifier
#' 
#' @author Janina Reeder
#' 
#' @return row containing the UI elements
longResultsUI <- function(id) {
  ns <- NS(id)
  
  fluidRow(
    box(
      title = "LONGITUDINAL ANALYSIS", solidHeader = TRUE, 
      collapsible = TRUE, width = 12,
      fluidRow(
        column(
          width = 11, id = ns("longplot_column"),
          shinycssloaders::withSpinner(
            plotly::plotlyOutput(ns("longplot"), 
                                 width = "auto", 
                                 height = "auto"),
            type = 3, color = "#424242", 
            color.background = "#fdfdfc"),
          br(),
          shinyjs::disabled(
            shinyWidgets::dropdownButton(
              tags$h3("Plot Options"),
              shinyWidgets::radioGroupButtons(
                inputId = ns("longChoice"),
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
  )
}

#' Longitudinal analysis module server code
#'
#' @param input shiny input
#' @param output shiny output
#' @param session shiny session
#' @param aggDat aggregated MRExperiment
#' @param featLevel chosen feature level (aggregation level)
#' @param longSettings reactive storing values selected in analysis 
#' input interface
#' @param normalizedData reactive boolean indicating if data has been normalized
#' @param reset boolean reactive which resets the module if TRUE
#' 
#' @author Janina Reeder
#' 
#' @return list containing R code for analysis and for feature plots
longResults <- function(input, output, session, 
                        aggDat, 
                        featLevel, 
                        longSettings, 
                        normalizedData,
                        reset) {
  ns <- session$ns
  
  repCode <- reactiveVal(NULL)
  width <- reactiveVal(650)
  
  observe({
    req(reset())
    shinyjs::disable("optionbutton_state")
    shinyjs::enable("longChoice")
    shinyjs::js$removeInputElems()
    repCode(NULL)
  })
  
  observe({
    width(input$plotWidth)
  })
  
  
  observeEvent(width(), {
    req(aggDat(), longSettings()$featureselect)
    plotly::plotlyProxy("longplot", session) %>%
      plotly::plotlyProxyInvoke(
        "relayout",
        list(width = width())
      )
  })
  
  
  ## calls function to do a feature abundance plot
  output$longplot <- plotly::renderPlotly({
    req(aggDat(), longSettings()$featureselect, longSettings()$comparison,
        longSettings()$phenolevels)
    shinyjs::js$removeInputElems()
    shinyjs::enable("optionbutton_state")
    ylab <- "Reads"
    if(isFALSE(normalizedData())){
      shinyWidgets::updateRadioGroupButtons(session = session,
                                            inputId = "longChoice",
                                            selected = "Reads")
      shinyjs::removeClass(id = "pspan", class = "hide")
      shinyjs::disable("longChoice")
    } else {
      ylab <- input$longChoice
      shinyjs::addClass(id = "pspan", class = "hide")
      shinyjs::enable("longChoice")
    }
    plotTitle <- paste0(gsub("Reads","Abundance",ylab),
                        " of ",
                        longSettings()$featureselect)
    phenoId <- longSettings()$phenoid
    showLines <- TRUE
    if(phenoId == ""){
      phenoId <- "SAMPLE_ID"
      showLines <- FALSE
    }
    
    plotLongFeature(aggDat(),
                    x_var = longSettings()$comparison,
                    id_var = phenoId,
                    plotTitle = plotTitle,
                    feature = longSettings()$featureselect,
                    ylab = ylab,
                    log = input$logS,
                    showLines = showLines,
                    x_levels = longSettings()$phenolevels
    ) 
  })
  
  ## update R code needed to make plot
  observe({
    req(aggDat(), longSettings()$featureselect)
    plotTitle <- paste0(gsub("Reads","Abundance",input$longChoice),
                        " of ",
                        longSettings()$featureselect)
    
    repCode(paste(
      paste0("#' ### Feature plot"),
      paste0("", "#- fig.width = 8"),
      paste0("plotLongFeature(aggDat,"),
      paste0("\tx_var = \"", longSettings()$comparison, "\","),
      'if'(longSettings()$phenoid != "",
           paste0("\tid_var = \"", longSettings()$phenoid, "\","),
           paste0("\tshowLines = FALSE,")),
      paste0("\tplotTitle = \"",plotTitle, "\","),
      paste0("\tfeature = \"", longSettings()$featureselect, "\","),
      paste0("\tylab = \"", input$longChoice, "\","),
      paste0("\tlog = ", 
             'if'(is.null(input$logS), "TRUE", input$logS), ","),
      paste0("\tx_levels = c(", paste0("\"",
                                       longSettings()$phenolevels, 
                                       "\"",
                                       collapse = ","),
             "))\n\n"),
      sep = "\n"
    ))
  })
  
  return(repCode)
  
}
