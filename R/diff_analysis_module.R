
## Main module for diff sample analysis

#'Diff Analysis Module - UI
#'
#' @param id namespace identifier
#'
#' @author Janina Reeder
#'
#' @return fluidRow containing the ui code
#' 
#' @examples diffAnalysisUI("diffanalysis_id")
#' 
#' @export
diffAnalysisUI <- function(id) {
  ns <- NS(id)
  
  fluidRow(
        width = 12,
        column(
            width = 3,
            diffInputUI(ns("diffInput"))
        ),
        column(
            width = 9,
            fluidRow(
                width = 11,
                box(width = 10,
                    p("Differential abundance analysis allows to compare the microbiome composition
                    across certain conditions.")
                ),
                diffTableUI(ns("differentialTable"))
            )
        )
    )
}


#' diff Analysis Module - server
#'
#' @param input shiny input
#' @param output shiny output
#' @param session shiny session
#' @param data the main data object returned from data_input_module
#' @param levelOpts available levels to aggregate on (depends on input data)
#' @param chosenLevel previously selected level (passed from different instance)
#' @param resetInput reactive boolean determining if reset is required
#' @param aggData the aggregated MRExperiment object
#' @param normalizedData boolean indicating if normalization was done
#' 
#' @author Janina Reeder
#' 
#'
#' @return reactive holding code to be used in reports
diffAnalysis <- function(input, output, session, data, levelOpts,
                          chosenLevel, resetInput, aggData, normalizedData) {
  ## diff SAMPLE ANALYSIS
  ns <- session$ns
  
  diffSettings <- callModule(diffInput,
                             "diffInput",
                             meData = data$meData,
                             facetOptions = data$facetOptions,
                             reset = resetInput
  )
  
  diffTableRep <- callModule(diffTable,
                             "differentialTable",
                             aggDat = reactive(aggData$mrobj),
                             featLevel = chosenLevel,
                             diffSettings = diffSettings,
                             reset = resetInput,
                             normalized = normalizedData
  )
  
  ## store differential analysis code for reports
  diffRep <- reactive({
    list(
      "table" = diffTableRep$reptable(),
      "plot" = diffTableRep$repplot()
    )
  })
  
  return(diffRep)
  
}