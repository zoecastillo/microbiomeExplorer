
## Main module for long sample analysis

#'Long Analysis Module - UI
#'
#' @param id namespace identifier
#'
#' @author Janina Reeder
#'
#' @return fluidRow containing the ui code
#' 
#' @examples longAnalysisUI("longanalysis_id")
#' 
#' @export
longAnalysisUI <- function(id) {
  ns <- NS(id)
  
  fluidRow(
    width = 12,
    column(
      width = 3,
      longInputUI(ns("longInput"))
    ),
    column(
      width = 9,
      fluidRow(
        width = 11,
        box(
          width = 10,
          p("Longitudinal analysis provides a module to compare the microbial 
            composition across time points or conditions (e.g. tissues).
            The order of the phenotype levels is preserved in the visualization.
            Use the optional phenotype ID parameter to add connections 
            between IDs over all selected time points/conditions.
            Individual lines can be highlighted via mouse clicks or
            the input area with option to choose colors. Select multiple
            connections by holding the \"Shift\" key. Double clicks
            near the edge of the plot remove selections.")
        ),
        longResultsUI(ns("longResults"))
      )
    )
  )
}


#' long Analysis Module - server
#'
#' @param input shiny input
#' @param output shiny output
#' @param session shiny session
#' @param data the main data object returned from data_input_module
#' @param levelOpts available levels to aggregate on (depends on input data)
#' @param chosenLevel previously selected level (passed from longerent instance)
#' @param resetInput reactive boolean determining if reset is required
#' @param aggData the aggregated MRExperiment object
#' @param normalizedData boolean indicating if normalization was done
#' 
#' @author Janina Reeder
#' 
#'
#' @return reactive holding code to be used in reports
longAnalysis <- function(input, output, session, data, levelOpts,
                         chosenLevel, resetInput, aggData, normalizedData) {
  ns <- session$ns
  
  longSettings <- callModule(longInput,
                             "longInput",
                             meData = data$meData,
                             facetOptions = data$facetOptions,
                             reset = resetInput,
                             aggDat = reactive(aggData$mrobj)
  )
  
  longResultsRep <- callModule(longResults,
                               "longResults",
                               aggDat = reactive(aggData$mrobj),
                               featLevel = chosenLevel,
                               longSettings = longSettings,
                               normalizedData = normalizedData,
                               reset = resetInput
  )
  
  return(longResultsRep)
  
}