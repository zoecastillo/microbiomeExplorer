
#' Main diffanalysis input module.
#' Set up to handle diff analysis tabs in the app depending on given parameters
#'
#' @param id element identifier - namespace
#'
#' @author Janina Reeder
#'
#' @return box containing ui element
#' @export
diffInputUI <- function(id) {
  ns <- NS(id)
  
  box(width = 12, id = ns("analysisbox"),
      h4("ANALYSIS PARAMETERS"),
      selectInput(ns("method"),
                  label = "Method", choices = getOption("me.diffmethods"), multiple = FALSE, 
                  selectize = FALSE, width = "250px"
      ),
      div(
        selectInput(ns("comparison"),
                    label = "Comparison phenotype", choices = "", 
                    multiple = FALSE, selectize = FALSE, width = "250px"
        ),
        selectInput(ns("phenolevel1"),
                    label = "Comparison level 1",
                    choices = "", multiple = FALSE, width = "250px"
        ),
        selectInput(ns("phenolevel2"),
                    label = "Comparison level 2",
                    choices = "", multiple = FALSE, width = "250px"
        )
      ),
      ## buttons are used to submit events. This ensures plots are not redrawn 
      ## while user still adjusts parameters
      div(id = ns("buttondiv"),
          fluidRow( width = 12, id = "actionbuttonrow",
                    ## update analysis outputs (plots/tables)
                    actionButton(ns("updatebutton"), 
                                 icon = icon("far fa-chart-bar"),
                                 label = HTML("&nbsp;UPDATE"), width = "120px"),
                    actionButton(ns("reportButton"), 
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
#' @param meData MRExperiment object storing all data
#' @param facetOptions named vector of available facet choices
#' @param reset reactive boolean determining if all inputs should be reset
#'
#' @author Janina Reeder
#' 
#' @importFrom metagenomeSeq MRcounts
#' @importFrom Biobase pData fData
#'
#' @return list holding all chosen values and the selected feature
#' @export
diffInput <- function(input, output, session, 
                          meData, facetOptions = NULL,
                          reset) {
  ns <- session$ns
  
  ## reactive value storing choices made in the UI
  chosenValues <- reactiveVal(NULL)
  ## reactive storing all available phenotype levels
  phenoOpts <- reactiveVal(NULL)
  
  ## reset all entries
  observe({
    req(reset())
    chosenValues(NULL)
    phenoOpts(NULL)
  })
  
  ## initialize input elements
  observe({
    updateSelectInput(session, "comparison", 
                      choices = c("", facetOptions()))
  })
  
  observe({
    req(input$comparison)
    req(input$comparison %in% names(pData(meData())))
    phenoOpts(levels(
      forcats::fct_explicit_na(
        factor(pData(meData())[, input$comparison]),
        na_level = "NA"))
    )
    updateSelectInput(session, "phenolevel1",
                      choices = phenoOpts(),
                      selected = phenoOpts()[1])
    updateSelectInput(session, "phenolevel2",
                      choices = phenoOpts(),
                      selected = phenoOpts()[2])
  }) 
  
  observe({
    req(!is.null(input$comparison))
    req(input$comparison == "")
    updateSelectInput(session, "phenolevel1",
                      choices = "")
    updateSelectInput(session, "phenolevel2",
                      choices = "")
  })
  
  ## react to selecting first phenolevel (second must be different)
  observeEvent(input$phenolevel1, {
    req(input$phenolevel1 != "", input$phenolevel1 == input$phenolevel2)
    updateSelectInput(session, "phenolevel2",
                      selected = phenoOpts()[!(phenoOpts() %in% 
                                                 input$phenolevel1)][1])
  })
  
  ## react to selecting second phenolevel (first must be different)
  observeEvent(input$phenolevel2, {
    req(input$phenolevel2 != "", input$phenolevel1 == input$phenolevel2)
    updateSelectInput(session, "phenolevel1",
                      selected = phenoOpts()[!(phenoOpts() %in% 
                                                 input$phenolevel2)][1])
  })
  
  ## main control button: store input choices in chosenValues
  observeEvent(input$updatebutton, {
    comparison <- NULL
    if(input$comparison != "")
      comparison <- input$comparison
    cV <- list("method" = input$method,
               "comparison" = comparison,
               "phenolevel1" = input$phenolevel1,
               "phenolevel2" = input$phenolevel2)
    chosenValues(cV)
  })
  
  
  return(chosenValues)
}
