
#' Main correlation analysis input module.
#' Handles correlation analysis tab in the app
#'
#' @param id element identifier - namespace
#' @param type determines if 'feature' or 'pheno' correlation
#'
#' @author Janina Reeder
#'
#' @return box containing ui element
corrInputUI <- function(id, type) {
  ns <- NS(id)
  
  box(width = 12, id = ns("analysisbox"),
      h4("ANALYSIS PARAMETERS"),
      div(
        h5('if'(type == "feature", 
                "Feature vs Feature",
                "Feature vs Phenotype")),
        selectizeInput(
          ns("featurecorr1"),
          label = "Feature correlation (base)", choices = "", 
          multiple = FALSE, width = "250px"
        ),
        if(type == "feature") {
          selectizeInput(
            ns("featurecorr2"),
            label = "Feature correlation 2", choices = "", 
            multiple = FALSE, width = "250px"
          )
        } else {
          selectInput(
            ns("phenoselect"),
            label = "Phenotype correlation", choices = "", 
            multiple = FALSE, selectize = TRUE, width = "250px"
          )
        }
      ),
      div(
        selectInput(
          ns("facetby"),
          label = "Facet columns by", choices = "", multiple = FALSE, 
          selectize = FALSE, width = "250px"
        ),
        selectInput(
          ns("facetby2"),
          label = "Facet rows by", choices = "", multiple = FALSE, 
          selectize = FALSE, width = "250px"
        )
      ),
      selectInput(
        ns("method"),
        label = "Method", choices = getOption("me.corrmethods"), 
        multiple = FALSE, 
        selectize = FALSE, width = "250px"
      ),
      ## buttons are used to submit events. This ensures plots are not redrawn 
      ## while user still adjusts parameters
      div(id = ns("buttondiv"),
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
#' @param type of the correlation (feature vs phenotype)
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
corrInput <- function(input, output, session, type,
                      meData, facetOptions = NULL,
                      reset, aggDat = reactive(NULL)) {
  ns <- session$ns
  
  ## reactive value storing choices made in the UI
  chosenValues <- reactiveVal(NULL)
  ## reactive storing sorted features at given level
  aggFeatures <- reactiveVal(NULL)
  
  observe({
    req(aggDat())
    aggFeatures(rownames(MRcounts(aggDat())))
  })
  
  ## check which phenotypes are numeric
  numericPheno <- reactive({
    req(meData())
    names(pData(meData()))[vapply(pData(meData()), is.numeric, logical(1))]
  })
  
  ## reset all entries
  observe({
    req(reset())
    chosenValues(NULL)
    aggFeatures(NULL)
  })
  
  
  ## initialize input elements
  observe({
    req(meData())
    if(!is.null(chosenValues()$method))
      updateSelectInput(session, "method", 
                        selected = chosenValues()$method)
    updateSelectInput(session, "facetby", 
                      choices = c("", facetOptions()), 
                      selected = chosenValues()$facetby)
    updateSelectInput(session, "facetby2", 
                      choices = c("", facetOptions()), 
                      selected = chosenValues()$facetby2)
  })
  
  observe({
    req(meData(),numericPheno(), type() == "pheno")
    updateSelectInput(session, "phenoselect", 
                      choices = c("", numericPheno()), 
                      selected = chosenValues()$phenoselect)
  })
  
  # ## react to phenotype selection variable (facets should be distinct)
  observeEvent(input$phenoselect, {
    req(facetOptions())
    facetopts <- facetOptions()
    facetopts <- facetopts[!facetopts %in% input$phenoselect]
    updateSelectInput(session, "facetby", choices = c("", facetopts))
    updateSelectInput(session, "facetby2", choices = c("", facetopts))
  }, priority = 15, ignoreInit = TRUE)
  
  ## react to first facet variable (phenotype variables should be distinct)
  observeEvent(input$facetby, {
    req(facetOptions(), type())
    facetopts <- facetOptions()
    facetopts <- facetopts[!facetopts %in% input$facetby]
    if (type() == "pheno") {
      req(numericPheno())
      numPheno <- numericPheno()
      numPheno <- numPheno[!numPheno %in% input$facetby]
      facetopts <- facetopts[!facetopts %in% input$phenoselect]
      updateSelectInput(session, "phenoselect", choices = c("", numPheno), 
                        selected = input$phenoselect)
    }
    updateSelectInput(session, "facetby2", choices = c("", facetopts), 
                      selected = input$facetby2)
  }, priority = 15, ignoreInit = TRUE)
  
  ## react to second facet variable (phenotype variables should be distinct)
  observeEvent(input$facetby2, {
    req(facetOptions(), type())
    facetopts <- facetOptions()
    facetopts <- facetopts[!facetopts %in% input$facetby2]
    if (type() == "pheno") {
      req(numericPheno())
      numPheno <- numericPheno()
      numPheno <- numPheno[!numPheno %in% input$facetby2]
      facetopts <- facetopts[!facetopts %in% input$phenoselect]
      updateSelectInput(session, "phenoselect", choices = c("", numPheno),
                        selected = input$phenoselect)
    }
    updateSelectInput(session, "facetby", choices = c("", facetopts), 
                      selected = input$facetby)
  }, priority = 15, ignoreInit = TRUE)
  
  ## update any inputs based on feature values
  observeEvent(aggFeatures(),{
    req(aggFeatures(), type())
    corrdata  <- data.frame(value = aggFeatures(), 
                            label = aggFeatures())
    selcorr <- corrdata[1,]
    updateSelectizeInput(session, "featurecorr1",
                         choices = corrdata,
                         selected = selcorr,
                         server = TRUE)
    selcorr <- ""
    if(type() == "feature"){
      if(!is.null(chosenValues()$featurecorr2))
        selcorr <- chosenValues()$featurecorr2
      updateSelectizeInput(session, "featurecorr2",
                           choices = corrdata,
                           selected = selcorr,
                           server = TRUE,
                           options = list(placeholder = "Select feature"))
    }
  })
  
  
  ## main control button: store input choices in chosenValues
  observeEvent(input$updatebutton, {
    req(type())
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
    if(type() == "feature"){
      cV <- list("featurecorr1" = input$featurecorr1,
                 "featurecorr2" = input$featurecorr2,
                 "facetby" = facet1,
                 "facetby2" = facet2,
                 "method" = input$method)
    } else {
      cV <- list("featurecorr1" = input$featurecorr1,
                 "phenoselect" = input$phenoselect,
                 "facetby" = facet1,
                 "facetby2" = facet2,
                 "method" = input$method)
    }
    chosenValues(cV)
  })
  
  return(chosenValues)
}
