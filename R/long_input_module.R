
#' Main diffanalysis input module.
#' Set up to handle diff analysis tabs in the app depending on given parameters
#'
#' @param id element identifier - namespace
#'
#' @author Janina Reeder
#'
#' @return box containing ui element
longInputUI <- function(id) {
  ns <- NS(id)
  
  box(
    width = 12, id = ns("analysisbox"),
    h4("ANALYSIS PARAMETERS"),
    selectizeInput(
      ns("featureselect"),
      label = "Select Feature", choices = "", 
      multiple = FALSE, 
      width = "250px"
    ),
    div(
      selectInput(
        ns("comparison"),
        label = "Longitudinal phenotype", choices = "", 
        multiple = FALSE, selectize = FALSE, width = "250px"
      ),
      selectizeInput(
        ns("phenolevels"),
        label = "Phenotype level order",
        choices = "", multiple = TRUE, width = "250px"
      ),
      selectInput(
        ns("phenoid"),
        label = "ID phenotype (optional)", choices = "", 
        multiple = FALSE, selectize = FALSE, width = "250px"
      )
    ),
    ## buttons are used to submit events. This ensures plots are not redrawn 
    ## while user still adjusts parameters
    div(
      id = ns("buttondiv"),
      fluidRow( 
        width = 12, id = "actionbuttonrow",
        ## update analysis outputs (plots/tables)
        actionButton(
          ns("updatebutton"), 
          icon = icon("far fa-chart-bar"),
          label = HTML("&nbsp;UPDATE"), width = "120px"),
        actionButton(
          ns("reportButton"), 
          label = HTML("<i class='far fa-bookmark'></i>&nbsp;&nbsp;REPORT"), 
          width = "120px")
      )
    )
  )
}


#' Server side for the analysis input module handling analysis control
#'
#' @param input shiny input
#' @param output shiny output
#' @param session shiny session
#' @param meData MRexperiment object storing all data
#' @param facetOptions named vector of available facet choices
#' @param reset reactive boolean determining if all inputs should be reset
#' @param aggDat aggregated MRexperiment
#'
#' @author Janina Reeder
#' 
#' @importFrom metagenomeSeq MRcounts
#' @importFrom Biobase pData fData
#'
#' @return list holding all chosen values and the selected feature
longInput <- function(input, output, session, 
                      meData, facetOptions = NULL,
                      reset, aggDat = reactive(NULL)) {
  ns <- session$ns
  
  ## reactive value storing choices made in the UI
  chosenValues <- reactiveVal(NULL)
  ## reactive storing all available phenotype levels
  phenoOpts <- reactiveVal(NULL)
  aggFeatures <- reactiveVal(NULL)
  
  observe({
    req(aggDat())
    aggFeatures(rownames(MRcounts(aggDat())))
  })
  
  ## reset all entries
  observe({
    req(reset())
    chosenValues(NULL)
    phenoOpts(NULL)
    aggFeatures(NULL)
  })
  
  ## update any inputs based on feature values
  observeEvent(aggFeatures(),{
    req(aggFeatures())
    featdata  <- data.frame(value = aggFeatures(), 
                            label = aggFeatures())
    selfeat <- ""
    if(!is.null(chosenValues()$featureselect))
      selfeat <- chosenValues()$featureselect
    updateSelectizeInput(session, "featureselect",
                         choices = featdata,
                         selected = selfeat,
                         options = list(placeholder = "Select Feature"),
                         server = TRUE)
  })
  
  ## initialize input elements
  observe({
    req(meData())
    updateSelectInput(session, "comparison", 
                      choices = c("", facetOptions()))
    updateSelectInput(session, "phenoid", 
                      choices = c("", names(pData(meData()))))
  })
  
  observe({
    req(input$comparison)
    req(input$comparison %in% names(pData(meData())))
    phenoOpts(levels(
      forcats::fct_explicit_na(
        factor(pData(meData())[, input$comparison]),
        na_level = "NA"))
    )
    updateSelectizeInput(session, "phenolevels",
                         choices = phenoOpts(),
                         options = list(placeholder = "Select levels"),
                         selected = "")
    remainingIds <- names(pData(meData()))[!(names(pData(meData())) 
                                             %in% input$comparison)]
    updateSelectInput(session, "phenoid", 
                      choices = c("", remainingIds))
  }) 
  
  observe({
    req(!is.null(input$comparison))
    req(input$comparison == "")
    updateSelectizeInput(session, "phenolevels",
                         choices = "",
                         options = list(placeholder = "Select levels"),
                         selected = "")
  })
  
  ## main control button: store input choices in chosenValues
  observeEvent(input$updatebutton, {
    comparison <- NULL
    if(input$comparison != "")
      comparison <- input$comparison
    cV <- list("featureselect" = input$featureselect,
               "comparison" = comparison,
               "phenolevels" = input$phenolevels,
               "phenoid" = input$phenoid)
    chosenValues(cV)
  })
  
  
  return(chosenValues)
}
