
#' Main intra analysis input module.
#' Set up to handle all analysis tabs in the app depending on given parameters
#'
#' @param id element identifier - namespace
#'
#' @author Janina Reeder
#'
#' @return box containing ui element
intraInputUI <- function(id) {
  ns <- NS(id)
  
  box(
    width = 12, id = ns("analysisbox"),
    h4("ANALYSIS PARAMETERS"),
    
    selectInput(
      ns("xvariable"),
      label = "Phenotype", choices = "", 
      multiple = FALSE, 
      selectize = FALSE, width = "250px"
    ),
    div(
      selectInput(
        ns("facetby"),
        label = "Facet columns by", choices = "", 
        multiple = FALSE, 
        selectize = FALSE, width = "250px"
      ),
      selectInput(
        ns("facetby2"),
        label = "Facet rows by", choices = "", 
        multiple = FALSE, 
        selectize = FALSE, width = "250px"
      )
    ),
    selectizeInput(
      ns("featureselect"),
      label = "Select Feature", choices = "", 
      multiple = FALSE, 
      width = "250px"
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
          label = HTML("&nbsp;UPDATE"), 
          width = "120px"),
        actionButton(
          ns("reportButton"), 
          label = HTML("<i class='far fa-bookmark'></i>&nbsp;&nbsp;REPORT"), 
          width = "120px")
      )
    )
  )
}


#' Server side for the intra analysis input module
#'
#' @param input shiny input
#' @param output shiny output
#' @param session shiny session
#' @param meData MRExperiment object storing all data
#' @param facetOptions named vector of available facet choices
#' @param reset reactive boolean determining if all inputs should be reset 
#' @param aggDat aggregated MRExperiment object (default is NULL)
#'
#' @author Janina Reeder
#' 
#' @importFrom metagenomeSeq MRcounts
#' @importFrom Biobase pData fData
#'
#' @return list holding all chosen values and the selected feature
intraInput <- function(input, output, session, 
                       meData, facetOptions = NULL, 
                       reset, aggDat = reactive(NULL)) {
  ns <- session$ns
  
  ## reactive value storing choices made in the UI
  chosenValues <- reactiveVal(NULL)
  ## specific reactive for just the feature selection
  ## this is required to only update the feature plot upon click/selection
  selectedFeat <- reactiveVal(NULL)
  ## reactive storing sorted features at given level
  aggFeatures <- reactiveVal(NULL)
  
  observe({
    req(aggDat())
    aggFeatures(rownames(MRcounts(aggDat())))
    shinyjs::enable("savebutton")
  })
  
  ## reset all entries
  observe({
    req(reset())
    chosenValues(NULL)
    selectedFeat(NULL)
    aggFeatures(NULL)
    shinyjs::disable("savebutton")
  })
  
  
  ## initialize input elements
  observe({
    req(meData())
    updateSelectInput(session, "xvariable", 
                      choices = names(pData(meData())),
                      selected = chosenValues()$xvariable)
    updateSelectInput(session, "facetby", 
                      choices = c("", facetOptions()), 
                      selected = chosenValues()$facetby)
    updateSelectInput(session, "facetby2", 
                      choices = c("", facetOptions()), 
                      selected = chosenValues()$facetby2)
  })
  
  ## react to main phenotype variable (facets should be distinct)
  observeEvent(input$xvariable, {
    req(facetOptions())
    facetopts <- facetOptions()
    facetopts <- facetopts[!facetopts %in% input$xvariable]
    updateSelectInput(session, "facetby", choices = c("", facetopts))
    updateSelectInput(session, "facetby2", choices = c("", facetopts))
  }, priority = 15, ignoreInit = TRUE)
  
  ## react to first facet variable (phenotype variables should be distinct)
  observeEvent(input$facetby, {
    req(facetOptions())
    facetopts <- facetOptions()
    facetopts <- facetopts[!facetopts %in% input$facetby]
    xopts <- names(pData(meData()))
    xopts <- xopts[!xopts %in% input$facetby]
    facetopts <- facetopts[!facetopts %in% input$xvariable]
    updateSelectInput(session, "xvariable", choices = xopts, 
                      selected = input$xvariable)
    updateSelectInput(session, "facetby2", choices = c("", facetopts), 
                      selected = input$facetby2)
  }, priority = 15, ignoreInit = TRUE)
  
  ## react to second facet variable (phenotype variables should be distinct)
  observeEvent(input$facetby2, {
    req(facetOptions())
    facetopts <- facetOptions()
    facetopts <- facetopts[!facetopts %in% input$facetby2]
    xopts <- names(pData(meData()))
    xopts <- xopts[!xopts %in% input$facetby2]
    facetopts <- facetopts[!facetopts %in% input$xvariable]
    updateSelectInput(session, "xvariable", choices = xopts, 
                      selected = input$xvariable)
    updateSelectInput(session, "facetby", choices = c("", facetopts), 
                      selected = input$facetby)
  }, priority = 15, ignoreInit = TRUE)
  
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
  
  ## main control button: store input choices in chosenValues
  observeEvent(input$updatebutton, {
    facet1 <- NULL
    if(!is.null(input$facetby)){
      if(input$facetby != "")
        facet1 <- input$facetby
    }
    facet2 <- NULL
    if(!is.null(input$facetby2)){
      if(input$facetby2 != "")
        facet2 <- input$facetby2
    }
    cV <- list("xvariable" = input$xvariable,
               "facetby" = facet1,
               "facetby2" = facet2)
    chosenValues(cV)
  })
  
  
  
  ## selected feature is stored specifically to update feature plot
  observeEvent(input$featureselect, {
    selectedFeat(input$featureselect)
  }, ignoreInit = TRUE)
  
  return(list(chosenValues = chosenValues, selectedFeat = selectedFeat))
}
