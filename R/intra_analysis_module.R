
## Main module for intra sample analysis

#' Intra Analysis Module - UI
#'
#' @param id namespace identifier
#'
#' @author Janina Reeder
#'
#' @return fluidRow containing the ui code
#' @export
intraAnalysisUI <- function(id) {
  ns <- NS(id)
  
  fluidRow(
    width = 12,
    column(
      width = 3,
      intraInputUI(ns("intraInput"))
    ),
    column(
      width = 9,
      fluidRow(
        width = 11,
        box(width = 10,
            p("The intra sample page contains methods to analyze the microbial composition 
              and diversity within a sample")
        ),
        relAbundanceUI(ns("abundancePlot")),
        featAbundanceUI(ns("featurePlot")),
        alphaDiversityUI(ns("alphaDiv"))
      )
    )
  )
}


#' Intra Analysis Module - server
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
#'
#' @export
intraAnalysis <- function(input, output, session, data, levelOpts,
                          chosenLevel, resetInput, aggData, normalizedData) {
  ## INTRA SAMPLE ANALYSIS
  ns <- session$ns
  
  
  intraSettings <- callModule(intraInput,
                              "intraInput",
                              meData = data$meData,
                              facetOptions = data$facetOptions,
                              reset = resetInput,
                              aggDat = reactive(aggData$mrobj)
  )
  
  
  relAbReact <- callModule(relAbundance,
                           "abundancePlot",
                           aggDat = reactive(aggData$mrobj),
                           featLevel = reactive(aggData$level),
                           intraSettings = intraSettings$chosenValues,
                           normalizedData = normalizedData,
                           reset = resetInput
  )
  
  featAbRep <- callModule(featAbundance,
                          "featurePlot",
                          aggDat = reactive(aggData$mrobj),
                          featLevel = reactive(aggData$level),
                          intraSettings = intraSettings$chosenValues,
                          selectedFeat = intraSettings$selectedFeat,
                          featName = relAbReact$featName,
                          numOfFeats = relAbReact$numOfFeats,
                          ylabMode = relAbReact$ylabMode,
                          normalizedData = normalizedData,
                          reset = resetInput
  )
  
  alphDivRep <- callModule(alphaDiversity,
                           "alphaDiv",
                           colorOptions = data$facetOptions,
                           aggDat = reactive(aggData$mrobj),
                           featLevel = reactive(aggData$level),
                           intraSettings = intraSettings$chosenValues,
                           reset = resetInput
  )
  
  ## store intra analysis code for reports
  intraRep <- reactive({
    list(
      "relab" = relAbReact$repCode(),
      "featab" = featAbRep(),
      "alpha" = alphDivRep()
    )
  })
  
  return(intraRep)
  
}