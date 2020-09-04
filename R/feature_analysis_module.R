
## Main module for feature sample analysis

#' feature Analysis Module - UI
#'
#' @param id namespace identifier
#'
#' @author Janina Reeder
#'
#' @return fluidRow containing the ui code
#' 
#' @examples featureAnalysisUI("featureanalysis_id")
#' 
#' @export
featureAnalysisUI <- function(id) {
  ns <- NS(id)
  
  fluidRow(
    width = 12,
    column(
      width = 3,
      featureInputUI(ns("featureInput"))
    ),
    column(
      width = 9,
      fluidRow(
        width = 11,
        box(width = 10,
            p("The feature sample page contains methods to analyze the microbial composition 
              and diversity within a sample")
        ),
        avgAbundanceUI(ns("avgFeaturePlot")),
      )
    )
  )
}


#' feature Analysis Module - server
#'
#' @param input shiny input
#' @param output shiny output
#' @param session shiny session
#' @param data the main data object returned from data_input_module
#' @param resetInput reactive boolean determining if reset is required
#' @param aggData the aggregated MRExperiment object
#' @param normalizedData boolean indicating if normalization was done
#' 
#' @author Janina Reeder
#' 
#'
#' @return reactive holding code to be used in reports
featureAnalysis <- function(input, output, session, data,
                          resetInput, aggData, normalizedData) {
  ## feature SAMPLE ANALYSIS
  ns <- session$ns
  
  
  featureSettings <- callModule(featureInput,
                              "featureInput",
                              meData = data$meData,
                              facetOptions = data$facetOptions,
                              reset = resetInput,
                              aggDat = reactive(aggData$mrobj)
  )
  
  
  featureRep <- callModule(avgAbundance,
                         "avgFeaturePlot",
                         aggDat = reactive(aggData$mrobj),
                         featLevel = reactive(aggData$level),
                         featureSettings = featureSettings$chosenValues,
                         normalizedData = normalizedData,
                         reset = resetInput
  )
  
  return(featureRep)
  
}