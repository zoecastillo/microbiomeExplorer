
## Main module for corr sample analysis

#'corr Analysis Module - UI
#'
#' @param id namespace identifier
#'
#' @author Janina Reeder
#'
#' @return fluidRow containing the ui code
#' 
#' @examples corrAnalysisUI("coranalysis_id")
#' 
#' @export
corrAnalysisUI <- function(id) {
  ns <- NS(id)
  
  fluidRow(
    width = 12,
    column(
      width = 3,
      corrInputUI(ns("fCorrInput"), type = "feature"),
      corrInputUI(ns("pCorrInput"), type = "pheno")
    ),
    column(
      width = 9,
      fluidRow(
        width = 11,
        box(width = 10,
            p("The correlation page contains methods to investigate how specific
            taxa are correlated with each other or to numeric phenotypes.")
        ),
        featureCorrUI(ns("fCorrelation")),
        phenotypeCorrUI(ns("pCorrelation"))
      )
    )
  )
}


#' corr Analysis Module - server
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
corrAnalysis <- function(input, output, session, data, levelOpts,
                          chosenLevel, resetInput, aggData) {
  ## corr SAMPLE ANALYSIS
  ns <- session$ns
  
  fCorrSettings <- callModule(corrInput,
                              "fCorrInput",
                              type = reactive("feature"),
                              meData = data$meData,
                              facetOptions = data$facetOptions,
                              reset = resetInput,
                              aggDat = reactive(aggData$mrobj)
  )
  
  pCorrSettings <- callModule(corrInput,
                              "pCorrInput",
                              type = reactive("pheno"),
                              meData = data$meData,
                              facetOptions = data$facetOptions,
                              reset = resetInput,
                              aggDat = reactive(aggData$mrobj)
  )
  
  ## reactive values storing correlation input settings
  fCorrFeatBase <- reactiveVal("")
  pCorrFeatBase <- reactiveVal("")
  fCorrMethod <- reactiveVal("")
  pCorrMethod <- reactiveVal("")

  ## update base feature for comparisons only if changed
  observe({
      req(fCorrFeatBase() !=
              fCorrSettings()$featurecorr1)
      fCorrFeatBase(fCorrSettings()$featurecorr1)
  }, priority = 20)
  
  ## update base feature for comparisons only if changed
  observe({
    req(pCorrFeatBase() !=
          pCorrSettings()$featurecorr1)
    pCorrFeatBase(pCorrSettings()$featurecorr1)
  }, priority = 20)

  ## update correlation method only if changed
  observe({
      req(fCorrMethod() != fCorrSettings()$method)
      fCorrMethod(fCorrSettings()$method)
  }, priority = 20)

  ## update correlation method only if changed
  observe({
    req(pCorrMethod() != pCorrSettings()$method)
    pCorrMethod(pCorrSettings()$method)
  }, priority = 20)
  

  fCorrRep <- callModule(featureCorr, "fCorrelation",
      aggDat = reactive(aggData$mrobj),
      colorOptions = data$facetOptions,
      corFeatBase = fCorrFeatBase,
      corFeat2 = reactive(fCorrSettings()$featurecorr2),
      corFacet1 = reactive(fCorrSettings()$facetby),
      corFacet2 = reactive(fCorrSettings()$facetby2),
      corMethod = fCorrMethod,
      reset = resetInput
  )

  pCorrRep <- callModule(phenotypeCorr, "pCorrelation",
                        aggDat = reactive(aggData$mrobj),
                        colorOptions = data$facetOptions,
                        corFeatBase = pCorrFeatBase,
                        corPheno = reactive(pCorrSettings()$phenoselect),
                        corFacet1 = reactive(pCorrSettings()$facetby),
                        corFacet2 = reactive(pCorrSettings()$facetby2),
                        corMethod = pCorrMethod,
                        reset = resetInput
  )


  ## store correlation analysis code for reports
  corrRep <- reactive({
      list(
          "fcor" = fCorrRep(),
          "pcor" = pCorrRep()
      )
  })
  
  return(corrRep)
  
}