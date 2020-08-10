
## Main module for inter sample analysis

#'inter Analysis Module - UI
#'
#' @param id namespace identifier
#'
#' @author Janina Reeder
#'
#' @return fluidRow containing the ui code
#' 
#' @examples interAnalysisUI("interanalysis_id")
#' 
#' @export
interAnalysisUI <- function(id) {
  ns <- NS(id)
  
  fluidRow(
    width = 12,
    column(
      width = 3,
      betaInputUI(ns("betaInput")),
      heatmapInputUI(ns("heatmapInput"))
    ),
    column(
      width = 9,
      fluidRow(
        width = 11,
        box(
          width = 10,
          p("The inter sample page contains methods to analyze the microbial 
           composition and diversity between samples using PCA and interactive 
           heatmaps.")
        ),
        betaDiversityUI(ns("betaDiv")),
        abundanceHeatmapUI(ns("heatmap"))
      )
    )
  )
}


#' inter Analysis Module - server
#'
#' @param input shiny input
#' @param output shiny output
#' @param session shiny session
#' @param data the main data object returned from data_input_module
#' @param levelOpts available levels to aggregate on (depends on input data)
#' @param chosenLevel previously selected level (passed from different instance)
#' @param resetInput reactive boolean determining if reset is required
#' @param aggData the aggregated MRExperiment object
#'
#' @return reactive holding code to be used in reports
interAnalysis <- function(input, output, session, data, levelOpts,
                          chosenLevel, resetInput, aggData) {
  ## inter SAMPLE ANALYSIS
  ns <- session$ns
  
  betaSettings <- callModule(betaInput,
                             "betaInput",
                             meData = data$meData,
                             adonisOptions = data$adonisOptions,
                             reset = resetInput
  )
  
  ## reactive storing beta distance settings
  betaDist <- reactiveVal("")  
  ## reactives storing heatmap settings
  hmSort <- reactiveVal("")
  hmFeatList <- reactiveVal(NULL)
  
  ## update betaDist if a new measure was selected
  observe({
    req(betaDist() != betaSettings()$distance)
    betaDist(betaSettings()$distance)
  })
  
  observe({
    req(resetInput())
    betaDist("")
    hmSort("")
    hmFeatList(NULL)
  })
  
  betaDivRep <- callModule(betaDiversity,
                           "betaDiv",
                           aggDat = reactive(aggData$mrobj),
                           aggLevel = reactive(aggData$level),
                           colorOptions = data$facetOptions,
                           shapeOptions = data$shapeOptions,
                           betadistance = betaDist,
                           betaSettings = betaSettings,
                           reset = resetInput
  )
  
  
  heatmapSettings <- callModule(heatmapInput,
                                "heatmapInput",
                                meData = data$meData,
                                reset = resetInput,
                                aggDat = reactive(aggData$mrobj)
  )
  
  
  
  ## update sorting method (only) if a new value was selected
  observe({
    req(hmSort() != heatmapSettings()$sorting)
    hmSort(heatmapSettings()$sorting)
  }, priority = 20)
  
  ## stores selected features
  observe({
    hmFeatList(heatmapSettings()$featureselect)
  }, priority = 20)
  
  abHeatRep <- callModule(abundanceHeatmap,
                          "heatmap",
                          aggDat = reactive(aggData$mrobj),
                          featLevel = chosenLevel,
                          colorOptions = data$facetOptions,
                          levelOpts = levelOpts,
                          hmSort = hmSort,
                          hmFeatList = hmFeatList,
                          reset = resetInput
  )
  
  interRep <- reactive({
    list(
      "beta" = betaDivRep(),
      "abheat" = abHeatRep()
    )
  })
  
  
  return(interRep)
  
}