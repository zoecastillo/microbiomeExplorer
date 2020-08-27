

## MAIN SERVER

#' Set up the shiny server
#'
#' @param session 
#' @param input 
#' @param output 
#' 
#' @author Janina Reeder
#'
#' @import shiny
#'
#' @export
shinyServer(function(session, input, output) {
  
  
  ## tabs are initially disabled
  js$disableTab("PHENOTYPE")
  js$disableTab("FEATURES")
  js$disableTab("INTRA SAMPLE")
  js$disableTab("INTRA FEATURE")
  js$disableTab("INTER SAMPLE")
  js$disableTab("CORRELATION")
  js$disableTab("DIFFERENTIAL")
  js$disableTab("LONGITUDINAL")
  
  ## reactive steering reset of all analysis parts of the app
  resetInput <- reactiveVal(NULL)
  ## reactive steering reset of reports
  resetReports <- reactiveVal(FALSE)
  #pheno data modified/added
  addPheno <- reactiveVal(FALSE)
  
  ## reactives storing code to be passed to report generation tab
  dataSource <- reactiveVal(NULL)
  dataFilterRep <- reactiveVal(NULL)
  qcRep <- reactiveVal(NULL)
  phenoModRep <- reactiveVal(NULL)
  featureModRep <- reactiveVal(NULL)
  analysisRep <- reactiveVal(NULL)
  
  levelOpts <- reactiveVal(NULL)
  ## reactive storing aggregation level
  chosenLevel <- reactiveVal(NULL)
  
  ## reactive values storing aggregated data and aggregation level
  aggData <- reactiveValues(
    mrobj = NULL,
    level = NULL
  )
  normalizedData <- reactiveVal(FALSE)
  
  ## DATA FILTERING AND PREPROCESSING (SUBSETTING, NORMALIZATION)
  
  data <- callModule(
    microbiomeExplorer:::dataInput, "loadnfilter",
    dataSource = dataSource,
    dataFilterRep = dataFilterRep,
    qcRep = qcRep,
    addPheno = addPheno,
    resetReports = resetReports
  )
  
  ## enable tabs once data is available
  observeEvent(data$meData(),{
    resetInput(NULL)
    js$enableTab("PHENOTYPE")
    js$enableTab("FEATURES")
    js$enableTab("INTRA SAMPLE")
    js$enableTab("INTRA FEATURE")
    js$enableTab("INTER SAMPLE")
    js$enableTab("CORRELATION")
    js$enableTab("DIFFERENTIAL")
    js$enableTab("LONGITUDINAL")
    disable("intraAnalysis-intraInput-analysisbox")
    disable("featureAnalysis-featureInput-analysisbox")
    disable("interAnalysis-betaInput-analysisbox")
    disable("interAnalysis-heatmapInput-analysisbox")
    disable("corrAnalysis-fCorrInput-analysisbox")
    disable("corrAnalysis-pCorrInput-analysisbox")
    disable("diffAnalysis-diffInput-analysisbox")
    disable("longAnalysis-longInput-analysisbox")
    featName <- colnames(fData(data$meData()))
    if(length(colnames(fData(data$meData()))) == 0){
      featName <- "unavailable"
    } else {
      featName <- featName[!(featName %in% getOption("me.featurenames"))]
      safeIndex <- sapply(featName, function(f){
        length(levels(as.factor(fData(data$meData())[[f]]))) > 1
      })
      featName <- featName[safeIndex]
    }
    if(!identical(levelOpts(),featName)){
      genusindex <- grepl("enus",featName)
      if(sum(genusindex) == 1){
        chosenLevel(featName[genusindex])
      } else if(length(featName) > 1){
        chosenLevel(featName[length(featName)-1])
      } else {
        chosenLevel(featName[1])
      }
      levelOpts(featName)
    }
    resetInput(TRUE)
  }, ignoreNULL = TRUE, ignoreInit = TRUE)
  
  
  ## PHENOTABLE MODULE
  callModule(
    microbiomeExplorer:::phenotypeTable, 
    "phenotab", data$meData, phenoModRep, addPheno)
  ## FEATURE TABLE MODULE
  callModule(
    microbiomeExplorer:::featureTable, 
    "featuretab", data$meData, featureModRep)
  
  ## store all data loading, filtering and modifications for report
  preprocessRep <- reactive({
    req(data$meData())
    c(
      phenoModRep(), featureModRep(), dataFilterRep(),
      paste(paste0("", "#+ mrdata, echo=TRUE"),
            paste0("\nmeData\n"),
            sep = "\n")
    )
  })
  
  
  ## AGGREGATION
  
  aggResults <- callModule(
    microbiomeExplorer:::aggregationTab, 
    "aggregateTab",
    resetInput, levelOpts, chosenLevel, data$meData)
  
  observe({
    aggData$mrobj <- aggResults$mrobj()
    aggData$level <- chosenLevel()
    normalizedData(aggResults$normalizedData())
  })
  
  observeEvent(aggResults$aggCode(),{
    analysisRep(unlist(c(analysisRep(),aggResults$aggCode())))
  })
  
  observe({
    req(aggData$mrobj)
    # enable analysis input modules when data is available
    enable("intraAnalysis-intraInput-analysisbox")
    enable("featureAnalysis-featureInput-analysisbox")
    enable("interAnalysis-betaInput-analysisbox")
    enable("interAnalysis-heatmapInput-analysisbox")
    enable("corrAnalysis-fCorrInput-analysisbox")
    enable("corrAnalysis-pCorrInput-analysisbox")
    enable("diffAnalysis-diffInput-analysisbox")
    enable("longAnalysis-longInput-analysisbox")
    js$enableAll()
  })
  
  ## INTRA SAMPLE ANALYSIS
  intraRep <- callModule(
    microbiomeExplorer:::intraAnalysis, "intraAnalysis", data, levelOpts,
    chosenLevel, resetInput, aggData, normalizedData)
  
  ## INTRA FEATURE ANALYSIS
  featureRep <- callModule(
    microbiomeExplorer:::featureAnalysis, "featureAnalysis", data,
    resetInput, aggData, normalizedData)
  
  ## INTER SAMPLE ANALYSIS
  interRep <- callModule(
    microbiomeExplorer:::interAnalysis, "interAnalysis", data, levelOpts,
    chosenLevel, resetInput, aggData)
  
  ## CORRELATION ANALYSIS
  corrRep <- callModule(
    microbiomeExplorer:::corrAnalysis, 
    "corrAnalysis", data, levelOpts,
    chosenLevel, resetInput, aggData)
  
  ## DIFFERENTIAL ANALYSIS
  diffRep <- callModule(
    microbiomeExplorer:::diffAnalysis, "diffAnalysis", data, levelOpts,
    chosenLevel, resetInput, aggData, normalizedData)
  
  ## DIFFERENTIAL ANALYSIS
  longRep <- callModule(
    microbiomeExplorer:::longAnalysis, "longAnalysis", data, levelOpts,
    chosenLevel, resetInput, aggData, normalizedData)
  
  ## ANALYSIS REPORTS
  
  ## capture report button clicks via reactive value (multiple buttons)
  reportClick <- reactiveVal(NULL)
  onclick("intraAnalysis-intraInput-reportButton", reportClick("intra"))
  onclick("featureAnalysis-featureInput-reportButton", reportClick("feature"))
  onclick("interAnalysis-betaInput-reportButton", reportClick("beta"))
  onclick("interAnalysis-heatmapInput-reportButton", reportClick("heatmap"))
  onclick("corrAnalysis-fCorrInput-reportButton", reportClick("fcorr"))
  onclick("corrAnalysis-pCorrInput-reportButton", reportClick("pcorr"))
  onclick("diffAnalysis-diffInput-reportButton", reportClick("diff"))
  onclick("longAnalysis-longInput-reportButton", reportClick("long"))
  
  ## reactive storing which code chunks contain aggregation code
  aggIndex <- reactiveVal(NULL)
  
  ## handle report button clicks: add new code to analysisRep
  observeEvent(reportClick(), {
    analysisRep(unlist(
      c(
        analysisRep(),
        switch(reportClick(),
               "intra" = intraRep(),
               "feature" = featureRep(),
               "beta" = interRep()$beta,
               "heatmap" = interRep()$abheat,
               "fcorr" = corrRep()$fcor,
               "pcorr" = corrRep()$pcor,
               "diff" = diffRep(),
               "long" = longRep()
        )
      )
    ))
    analysisRep(analysisRep()[!is.null(analysisRep())])
    ## capture aggregation chunks (cannot be deselected)
    aggIndex(stringr::str_starts(analysisRep(), "aggDat|meData"))
    ## omit duplicated code
    duplicatedEntries <- duplicated(analysisRep())
    analysisRep(analysisRep()[aggIndex() | !duplicatedEntries])
    reportClick(NULL)
  })
  
  observe({
    req(resetReports())
    dataFilterRep(NULL)
    qcRep(NULL)
    phenoModRep(NULL)
    featureModRep(NULL)
    analysisRep(NULL)
    resetReports(FALSE)
  })
  
  
  ## REPORT TAB
  callModule(
    microbiomeExplorer:::reportList, "reportlist",
    dataSource = dataSource,
    preprocessRep = preprocessRep,
    qcRep = qcRep,
    analysisRep = analysisRep,
    aggIndex = aggIndex,
    reset = resetInput
  )
})
