## Module handling data import/upload and filtering


#' Module handling file upload for the application: UI
#' In a deployed version this module should be replaced with database access
#'
#' @param id module identifier
#' 
#' @return div holding ui elements
#'
#' @author Janina Reeder
fileUploadUI <- function(id) {
  ns <- NS(id)
  
  div(
    ## inputs restricted to certain filetypes defined via global.R
    fileInput(
      ns("datafile"),
      label = "Upload Feature Count Data", 
      accept = getOption("me.filetypes")),
    shinyjs::disabled(
      div(
        id = ns("phenodiv"),
        fileInput(
          ns("phenofile"),
          label = "Link Phenotype Info", 
          accept = getOption("me.delim"))
      )),
    shinyjs::disabled(
      div(
        id = ns("featdiv"),
        fileInput(ns("featfile"),
                  label = "Add Taxonomy Info", 
                  accept = getOption("me.delim"))
      ))
  )
}

#' Module handling file upload for the application: server
#'
#' @param input module input
#' @param output module output
#' @param session app session
#' @param meData main reactive storing the MRexperiment data
#' @param meName main reactive storing the filename uploaded
#' @param initializeData reactiveVal keeping track of new uploads to reset data
#' @param addPheno reactiveVal keeping track of phenodata changes
#' @param dataSource reactive Value storing commands for loading data
#' @param resetFile indicating if module should be reset
#' 
#' @return boolean denoting successful upload of a file
#'
#' @author Janina Reeder
fileUpload <- function(input, output, session, 
                       meData, meName, initializeData, 
                       addPheno, dataSource, 
                       resetFile = reactive(NULL)) {
  ns <- session$ns
  
  fileUploaded <- reactiveVal(NULL)
  
  ## reset all input elements
  observe({
    req(resetFile())
    shinyjs::reset("datafile")
    shinyjs::reset("phenofile")
    shinyjs::reset("featfile")
    fileUploaded(NULL)
  })
  
  observe({
    req(meData())
    shinyjs::enable("phenodiv")
  })
  
  ## loading a data file: uses getFileType & readData from functions.R
  observeEvent(input$datafile, {
    withProgress(message = "Data Uploading", {
      fileext <- tools::file_ext(input$datafile$name)
      filetype <- getFileType(fileext)
      shinyjs::enable("phenodiv")
      meName(input$datafile$name)
      meData(filterMEData(readData(input$datafile$datapath, 
                                   type = filetype), 
                          minpresence = 1, minfeats = 2, minreads = 2))
      dataSource(
        paste0("meData <- filterMEData(readData(\"", 
               input$datafile$datapath, "\", type = \"", 
               filetype,
               "\"),minpresence = 1, minfeats = 2, minreads = 2)", "\n",
               "##Data file name: ",input$datafile$name))
      if (filetype %in% c("TAB", "CSV")) {
        shinyjs::enable("featdiv")
      } else if (filetype == "BIOM" & ncol(fData(meData())) < 2) {
        shinyjs::enable("featdiv")
      } else { 
        shinyjs::disable("featdiv")
      }
      ## UPLOAD COMPLETED
      initializeData(TRUE)
      shinyjs::reset("phenofile")
      shinyjs::reset("featfile")
      fileUploaded(TRUE)
    })
  })
  
  ## handling change in pheno file;
  observeEvent(input$phenofile, {
    if(ncol(pData(meData())) > 1){
      showModal(modalDialog(
        title = "Overwrite phenotype data",
        "Phenotype data is already available for this project.\n
                Do you want to replace or extend the existing phenotype
                data?",
        footer = tagList(actionButton(ns("confirmReplace"), "Replace"),
                         actionButton(ns("confirmExtend"), "Extend"),
                         modalButton("Cancel")
        ),
        easyClose = FALSE))
    } else {
      bufMEDate <- addPhenoData(meData(), input$phenofile$datapath)
      if (is.null(bufMEDate)) {
        showModal(modalDialog(
          title = "Error merging phenotype information",
          "An error occurred when attempting to merge the phenotype 
                information with the counts data.\n
                The phenotable provided must have a row for every sample column 
                in the counts data.",
          easyClose = TRUE))
      } else {
        dataSource(paste0(dataSource(), "\n", 
                          "meData <- addPhenoData(meData, \"", 
                          input$phenofile$datapath, "\")", "\n",
                          "##Pheno file name: ",input$phenofile$name))
        meData(bufMEDate)
        addPheno(TRUE)
      }
    }
  })
  
  observeEvent(input$confirmReplace,{
    removeModal()
    bufMEDate <- addPhenoData(meData(), input$phenofile$datapath)
    if (is.null(bufMEDate)) {
      showModal(modalDialog(
        title = "Error merging phenotype information",
        "An error occurred when attempting to merge the phenotype 
                information with the counts data.\n
                The phenotable provided must have a row for every sample column 
                in the counts data.",
        easyClose = TRUE))
    } else {
      dataSource(paste0(dataSource(), "\n", 
                        "meData <- addPhenoData(meData, \"", 
                        input$phenofile$datapath, "\")", "\n",
                        "##Pheno file name: ",input$phenofile$name))
      meData(bufMEDate)
      addPheno(TRUE)
    }
  })
  
  observeEvent(input$confirmExtend,{
    removeModal()
    bufMEDate <- extendPhenoData(meData(), input$phenofile$datapath)
    if (is.null(bufMEDate)) {
      showModal(modalDialog(
        title = "Error merging phenotype information",
        "An error occurred when attempting to merge the phenotype 
                information with the counts data.\n
                The phenotable provided must have a row for every sample column 
                in the counts data.",
        easyClose = TRUE))
    } else {
      dataSource(paste0(dataSource(), "\n", 
                        "meData <- extendPhenoData(meData, \"", 
                        input$phenofile$datapath, "\")", "\n",
                        "##Pheno file name: ",input$phenofile$name))
      meData(bufMEDate)
      addPheno(TRUE)
    }
  })
  
  
  ## handling change in feature file (taxonomy);
  observeEvent(input$featfile, {
    bufMEDate <- addFeatData(meData(), input$featfile$datapath)
    if (is.null(bufMEDate)) {
      showModal(modalDialog(
        title = "Error merging feature information",
        "An error occurred when attempting to merge the feature 
                information with the counts data.\n
                The feature information provided must have a corresponding row 
                for every row in the counts data.",
        easyClose = TRUE))
    } else {
      dataSource(paste0(dataSource(), "\n", 
                        "meData <- addFeatData(meData, \"", 
                        input$featfile$datapath, "\")", "\n",
                        "##Pheno file name: ",input$featfile$name))
      meData(bufMEDate)
    }
  })
  
  return(fileUploaded)
}

## END UPLOAD DATA FROM FILE


## DATA QC AND FILTERING

#' Main Data input UI where the user selects files to upload to the app or 
#' connects to database
#'
#' @param id module identifier
#' 
#' @return fluidRow holding UI interface
#'
#' @author Janina Reeder
#' 
#' @examples dataInputUI("datainput_id")
#' 
#' @export
dataInputUI <- function(id) {
  ns <- NS(id)
  
  fluidRow(
    column(
      width = 4, id = ns("datauploadcol"),
      h4("UPLOAD DATA", style = "padding-top: 45px;"),
      fileUploadUI(ns("fileupload")),
      hr(),
      shinyjs::disabled(div(
        id = ns("filterdiv"),
        h4("FILTER DATA", style = "padding-top: 20px;"),
        h5("Features", style = "padding-top:0;"),
        sliderInput(
          ns("featuresams"),
          label = "Minimum Sample Presence",
          value = 1, min = 1, max = 1000, 
          round = TRUE, width = "250px"),
        h5("Samples"),
        sliderInput(
          ns("minfeats"),
          label = "Minimum Number of Features",
          value = 2, min = 2, max = 1000, 
          round = TRUE, width = "250px"),
        sliderInput(
          ns("minreads"),
          label = "Minimum Number of Reads",
          value = 2, min = 2, max = 1000, 
          round = TRUE, width = "250px"),
        box(
          collapsible = TRUE, collapsed = TRUE, width = 10,
          style = "padding:0; margin-top: -40px;",
          numericInput(
            ns("readnum"),
            label = "Set min #reads",
            value = 2, min = 2,
            width = "100px")
        )
      )),
      shinyjs::disabled(
        div(
          id = ns("phenodiv"), style = "margin-top: 100px",
          h5("Phenotypes", 
             style = "padding-top:0; margin-top: -11px;"),
          selectizeInput(
            ns("removepheno"),
            label = "Remove Phenotypes",
            choices = "", 
            multiple = FALSE, 
            options = list(placeholder = 
                             "Select Phenotype"), 
            width = "250px")
        )),
      shinyjs::disabled(
        div(
          id = ns("buttondiv"),
          fluidRow(
            width = 12, id = "actionbuttonrow2",
            ## filter updates meData for further analysis
            actionButton(
              ns("filterbutton"), 
              icon = icon("fas fa-filter"),
              label = HTML("&nbsp;FILTER"), 
              width = "120px"),
            ## resets back to original data 
            actionButton(
              ns("resetbutton"), 
              icon = icon("fas fa-redo-alt"),
              label = HTML("&nbsp;RESET"), 
              width = "120px")
          ),
          hr(),
          h4("OBTAIN DATA", style = "padding-top: 50px;"),
          fluidRow(
            width = 12, id = "actionbuttonrow2",
            ## downloads data to user's computer
            shinyjs::disabled(
              downloadButton(
                ns("savebutton"), 
                label = "GET DATA")),
            ## add QC plots to report
            actionButton(
              ns("reportButton"), 
              label = HTML("<i class='far fa-bookmark'></i>&nbsp;&nbsp;REPORT"), 
              width = "120px")
          )
        ))
    ),
    column(
      width = 8,
      uiOutput(ns("projecthead")),
      fluidRow(
        width = 10, id = "histogramrow",
        column(
          width = 5,
          plotly::plotlyOutput(
            ns("featpres"), 
            width = "200px", 
            height = "200px")
        ),
        column(
          width = 5,
          plotly::plotlyOutput(
            ns("libsize"), 
            width = "200px", 
            height = "200px")
        )
      ),
      shinycssloaders::withSpinner(
        plotly::plotlyOutput(
          ns("dataqcplot"), 
          width = "675px", 
          height = "550px"),
        type = 3, color = "#424242", 
        color.background = "#fdfdfc"),
      uiOutput(ns("plotoptions")),
      fluidRow(
        width = 12, id = "sampleoturow",
        plotly::plotlyOutput(
          ns("sampleotus"), 
          width = "650px", 
          height = "400px")
      ),
      uiOutput(ns("samplebaroptions"))
    )
  )
}


#' Helper function to filter phenodata for interesting phenotypes to be used 
#' for filtering or subsetting
#'
#' @param MRobj the MRexperiment storing the data
#'
#' @author Janina Reeder
#' 
#' @importFrom Biobase pData
#'
#' @return list of named vectors with names being pData column headers and 
#' values being unique entries; columns with only one entry or those with 
#' different values for each samples are omitted
getFilterChoices <- function(MRobj) {
  maxvalues <- nrow(pData(MRobj))
  ## Keeping only those that have more than one unique value,
  ## but less than the number of samples
  allchoices <- lapply(names(pData(MRobj)), function(s) {
    uniqvalues <- unique(pData(MRobj)[[s]])
    if (length(uniqvalues) <= 1 | length(uniqvalues) == maxvalues) {
      return(NULL)
    }
    return(uniqvalues)
  })
  names(allchoices) <- names(pData(MRobj))
  if("SAMPLE_ID" %in% names(allchoices)){
    allchoices[["SAMPLE_ID"]] <- pData(MRobj)[["SAMPLE_ID"]]
  }
  allchoices <- allchoices[vapply(allchoices, function(a) {
    !is.null(a)
  }, logical(1))]
  return(allchoices)
}


#' Main Data input server where the user selects files to upload to the app or 
#' connects to database
#'
#' @param input module input
#' @param output module output
#' @param session app session
#' @param dataSource reactive Value storing commands for loading data
#' @param dataFilterRep reactive Value storing commands for filtering data
#' @param qcRep reactive Value storing commands for producing qc plots
#' @param addPheno reactive boolean keeping track of phenodata changes
#' @param resetReports reactive boolean indicating whether reports need to be reset
#'
#' @author Janina Reeder
#'
#' @return list of reactives containing the uploaded and filtered data as well 
#' as the filterChoices on phenotypes
#'
#' @importFrom metagenomeSeq MRcounts 
#' @importFrom Biobase pData
#' @import shiny
#' @import shinydashboard
dataInput <- function(input, output, session, 
                      dataSource, dataFilterRep, qcRep,
                      addPheno, resetReports) {
  ns <- session$ns
  
  ## main data reactive to be passed to analysis tabs
  meData <- reactiveVal(NULL)
  meName <- reactiveVal("") 
  ## backup of original data: used when resetting
  origData <- reactiveVal(NULL)
  #complete data update
  initializeData <- reactiveVal(FALSE)
  filterChoices <- reactiveVal(NULL)
  remainingFilters <- reactiveVal(NULL)
  facetOptions <- reactiveVal(NULL)
  shapeOptions <- reactiveVal(NULL)
  colorChoice <- reactiveVal("No color")
  logChoice <- reactiveVal("none")
  featureSams <- reactiveVal(1)
  minReads <- reactiveVal(2)
  minFeats <- reactiveVal(2)
  qcCode <- reactiveVal(NULL)
  sortPheno <- reactiveVal(NULL)
  phenoValues <- reactiveVal(NULL)
  
  ## essentially meData, but only updated if plot changes are required
  plotData <- reactiveVal(NULL)
  phenoFilters <- reactiveVal(NULL)
  ## phenotypes selected during filtering steps
  selectedChoices <- reactiveVal(NULL)
  
  callModule(fileUpload, "fileupload", 
             meData = meData, 
             meName = meName,  
             initializeData = initializeData, 
             addPheno = addPheno,
             dataSource = dataSource)
  
  
  observeEvent(filterChoices(),{
    filterNames <- names(filterChoices())
    filterNames <- filterNames[!(filterNames %in% "SAMPLE_ID")]
    if(length(filterNames) == 0){
      shapeOptions(NULL)
      facetOptions(NULL)
    } else {
      shapeIndices <- vapply(filterNames, function(n){
        length(unique(pData(meData())[,n])) <= getOption("me.shapenum")
      }, logical(1))
      if(length(shapeIndices) > 0){
        sO <- filterNames[shapeIndices]
        if(any(!(sO %in% shapeOptions())) |
           any(!(shapeOptions() %in% sO))){
          shapeOptions(sO)
        }
      }
      facetIndices <- vapply(filterNames, function(n){
        length(unique(pData(meData())[,n])) <= getOption("me.facetnum")
      }, logical(1))
      if(length(facetIndices) > 0){
        fO <- filterNames[facetIndices]
        if(any(!(fO %in% facetOptions())) | 
           any(!(facetOptions() %in% fO))){
          facetOptions(fO)
        }
      }
    }
    remainingFilters(filterChoices()[!names(filterChoices()) %in% 
                                       selectedChoices()])
  })
  
  observe({
    updateSelectInput(session, "qc_color",
                      choices = c("No color", facetOptions()),
                      selected = isolate(colorChoice()))
    updateSelectInput(session, "bar_color",
                      choices = c("No color", facetOptions()),
                      selected = isolate(colorChoice()))
    updateSelectInput(session, "sort_pheno",
                      choices = facetOptions())
  })
  
  ## update input elements based on meData
  observeEvent(initializeData(),{
    req(initializeData(), meData())
    shinyjs::enable("filterdiv")
    shinyjs::enable("buttondiv")
    shinyjs::enable("savebutton")
    origData(meData())
    colorChoice("No color")  
    logChoice("none")
    featureSams(1)
    minFeats(2)
    minReads(2)    
    removeUI(".phenoremoverow", multiple = TRUE, immediate = TRUE)
    dataFilterRep(NULL)
    phenoFilters(NULL)  
    filterChoices(NULL)
    remainingFilters(NULL)    
    facetOptions(NULL)
    shapeOptions(NULL)
    plotData(meData())
    resetReports(TRUE)
    shinyjs::js$parsePhenoFilters()
    updateSliderInput(session, "featuresams", 
                      value = 1, min = 1,
                      max = max(rowSums(MRcounts(meData()) > 0)))
    updateSliderInput(session, "minfeats", 
                      value = 2, min = 2,
                      max = max(colSums(MRcounts(meData()) > 0)))
    updateSliderInput(session, "minreads", 
                      value = 2, min = 2,
                      max = max(colSums(MRcounts(meData()))))
    shinyjs::enable("phenodiv")
    filterChoices(getFilterChoices(meData()))
    updateSelectizeInput(session, "removepheno",
                         choices = names(filterChoices()), 
                         selected = "",
                         options = list(placeholder = "Select Phenotype"))
    initializeData(FALSE)
  })
  
  observeEvent(addPheno(),{
    req(addPheno(), meData())
    if (ncol(pData(meData())) > 1) {
      shinyjs::enable("phenodiv")
      filterChoices(getFilterChoices(meData()))
      updateSelectizeInput(session, "removepheno",
                           choices = names(filterChoices()), 
                           selected = "",
                           options = list(placeholder = "Select Phenotype")) 
      isolate({
        if(!colorChoice() %in% names(filterChoices()))
          colorChoice("No color")
      })
      plotData(meData())
    } else {
      shinyjs::disable("phenodiv")
    }
    addPheno(FALSE)
  })
  
  ## workaround to allow the user to select a specific value
  observeEvent(input$readnum,{
    updateSliderInput(session, "minreads", value = input$readnum)
  })
  
  ## show the qc data in log scale
  observeEvent(input$logoptions, {
    prevLog <- logChoice()
    minf <- input$minfeats
    minr <- input$minreads
    
    if(input$logoptions == "none"){
      if(prevLog == "x axis"){
        minr <- exp(minr)
      } else if(prevLog == "y axis"){
        minf <- exp(minf)
      } else if(prevLog == "both"){
        minr <- exp(minr)
        minf <- exp(minf)
      }
      updateSliderInput(session, "minfeats", min = minFeats(),
                        max = max(colSums(MRcounts(meData()) > 0)),
                        value = minf)
      updateSliderInput(session, "minreads", min = minReads(),
                        max = max(colSums(MRcounts(meData()))),
                        value = minr)
    } else if(input$logoptions == "x axis"){
      if(prevLog == "none"){
        minr <- log(minr)
      } else if(prevLog == "y axis"){
        minf <- exp(minf)
        minr <- log(minr)
      } else if(prevLog == "both"){
        minf <- exp(minf)
      }
      updateSliderInput(session, "minfeats", min = minFeats(),
                        max = max(colSums(MRcounts(meData()) > 0)),
                        value = minf)
      updateSliderInput(session, "minreads", min = log(minReads()),
                        max = log(max(colSums(MRcounts(meData())))),
                        value = minr)
      
    } else if(input$logoptions == "y axis"){     
      if(prevLog == "none"){
        minf <- log(minf)
      } else if(prevLog == "x axis"){
        minr <- exp(minr)
        minf <- log(minf)
      } else if(prevLog == "both"){
        minr <- exp(minr)
      }
      updateSliderInput(session, "minfeats", min = log(minFeats()),
                        max = log(max(colSums(MRcounts(meData()) > 0))),
                        value = minf)
      updateSliderInput(session, "minreads", min = minReads(),
                        max = max(colSums(MRcounts(meData()))),
                        value = minr)
    } else if(input$logoptions == "both"){
      if(prevLog == "none"){
        minf <- log(minf)
        minr <- log(minr)
      }
      if(prevLog == "x axis"){
        minf <- log(minf)
      }
      if(prevLog == "y axis"){
        minr <- log(minr)
      }
      updateSliderInput(session, "minfeats", min = log(minFeats()),
                        max = log(max(colSums(MRcounts(meData()) > 0))),
                        value = minf)
      updateSliderInput(session, "minreads", min = log(minReads()),
                        max = log(max(colSums(MRcounts(meData())))),
                        value = minr)
    }
    logChoice(input$logoptions)
  }, ignoreInit = TRUE)
  
  
  ## add pheno filter UI element
  observeEvent(input$removepheno, {
    req(input$removepheno,remainingFilters())
    insertUI(
      selector = paste0("#", ns("removepheno")),
      where = "afterEnd",
      ui = fluidRow(
        class = "phenoremoverow", 
        id = ns(paste0("ui_", input$removepheno)),
        column(width = 1,
               actionLink(ns(paste0("rmbut_", input$removepheno)),
                          label = "", icon = icon("times-circle-o"),
                          style = paste0("color: darkred"), class = "rmpclass")
        ),
        column(
          width = 10, id = ns("checkboxcol"),
          selectInput(
            inputId = ns(paste0("rp_", input$removepheno)),
            label = paste0(input$removepheno, 
                           ": Select values to remove"),
            multiple = TRUE,
            choices = as.character(
              remainingFilters()[[input$removepheno]]),
          )
        )
      ),
      immediate = TRUE
    )
    selectedChoices(c(selectedChoices(), input$removepheno))
    shinyjs::js$parsePhenoFilters()
  })
  
  observeEvent(selectedChoices(),{
    remainingFilters(filterChoices()[!names(filterChoices()) %in% 
                                       selectedChoices()])
    updateSelectizeInput(session, "removepheno", 
                         choices = c(names(remainingFilters())), 
                         selected = "",
                         options = list(placeholder = "Select Phenotype"))
  })
  
  ## observes remove button clicks handled via javascript in UI
  observe({
    ## last button is changed via Shiny.setInputValue within 'on' click 
    ## function (see UI)
    req(input$last_btn)
    isolate({
      sc <- selectedChoices()
      selectedChoices(sc[!sc %in% gsub(".*rmbut_", "", input$last_btn)])
    })
    ## jquery (used by removeUI) can't handle dots in names, 
    ##so we need to explicitely include control characters
    removeUI(paste0("#", gsub("rmbut", "ui", 
                              gsub("\\.", "\\\\.", input$last_btn))), 
             immediate = TRUE)
    shinyjs::js$parsePhenoFilters()
  })
  
  
  ## process filtering operation on data: if origData hasn't been saved yet, 
  ##store current meData as origData
  observeEvent(input$filterbutton, {
    mustabort <- FALSE
    if (is.na(input$featuresams)) {
      showModal(modalDialog(
        title = "Input Error",
        "No value entered for minimum sample presence.\n 
                Please revise.",
        easyClose = TRUE
      ))
    } else {
      ## parse input elements for max allowed values: shiny doesn't 
      ## control text inputs.
      if (input$featuresams > max(rowSums(MRcounts(meData()) > 0))) {
        mustabort <- TRUE
      }
      
      if (input$minfeats > max(colSums(MRcounts(meData()) > 0))) {
        mustabort <- TRUE
      }
      
      if (input$minreads > max(colSums(MRcounts(meData())))) {
        mustabort <- TRUE
      }
      if (mustabort == TRUE) {
        showModal(modalDialog(
          title = "Filtering Error",
          "One of the filter values chosen is larger than the overall 
                    maximum. All data will be removed.\n Please revise.",
          easyClose = TRUE
        ))
      } else {
        ## PERFORM FILTERING BASED ON FEATURES AND READS
        prevfeatsams <- featureSams()
        prevminfeats <- minFeats()
        prevminreads <- minReads()
        featureSams(input$featuresams)
        minFeats(input$minfeats)
        minReads(input$minreads)
        minf <- minFeats()
        minr <- minReads()
        if(logChoice() %in% c("both","x axis"))
          minr <- exp(minr)
        if(logChoice() %in% c("both","y axis"))
          minf <- exp(minf)
        if(prevfeatsams != featureSams() | 
           prevminfeats != minFeats() | 
           prevminreads != minReads()){
          bufData <- filterMEData(meData(), 
                                  minpresence = featureSams(), 
                                  minfeats = minf,
                                  minreads = minr)
          if (nrow(pData(bufData)) == 0) {
            showModal(modalDialog(
              title = "Error filtering data",
              "No data remains when applying size filters. Please revise ...",
              easyClose = TRUE
            ))
          } else {
            meData(bufData) 
            plotData(meData())
            updateSliderInput(session, "featuresams", 
                              min = input$featuresams)
            updateSliderInput(session, "minfeats", 
                              min = input$minfeats)
            updateSliderInput(session, "minreads", 
                              min = input$minreads)
            dataFilterRep(
              paste0(dataFilterRep(),
                     "\nmeData <- filterMEData(meData,minpresence = ", 
                     featureSams(),
                     ", minfeats = ", minf, ", minreads = ", 
                     minr, ")"))
          }
        }
        ## FILTER ON PHENOTABLE
        ## go through all open filters and adjust meData
        if (!is.null(input$parsedFilters)) {
          phenovalstoremove <- lapply(
            input$parsedFilters, 
            function(f) {
              selectedFilterVals <- input[[gsub("loadnfilter-ui", "rp", f)]]
              selectedFilterVals[selectedFilterVals == ""] <- "NA"
              selectedFilterVals
            }
          )
          names(phenovalstoremove) <- gsub("loadnfilter-ui_", "", 
                                           input$parsedFilters)
          bufData <- filterByPheno(meData(), phenovalstoremove)
          if (nrow(pData(bufData)) == 0) {
            showModal(
              modalDialog(
                title = "Error filtering data",
                "No data remains when applying phenotype filters. Please revise ...",
                easyClose = TRUE
              ))
          } else {
            removeUI(".phenoremoverow", multiple = TRUE, 
                     immediate = TRUE)
            selectedChoices(NULL)
            meData(bufData)
            filterChoices(getFilterChoices(meData()))
            remainingFilters(filterChoices())
            updateSelectizeInput(session, "removepheno",
                                 choices = names(filterChoices()), 
                                 selected = "",
                                 options = list(placeholder = "Select Phenotype"))
            isolate({
              if(!colorChoice() %in% names(pData(meData())))
                colorChoice("No color")
            })
            plotData(meData())
            #Adjust code for reports
            phenoFilters(c(phenoFilters(),phenovalstoremove))
            duplPF <- duplicated(phenoFilters())
            phenoFilters(phenoFilters()[!duplPF])
            
            dataFilterRep(
              paste0(dataFilterRep(),
                     "\nfiltList <- ",
                     paste0("list(",
                            paste0('if'(lengths(phenoFilters()) == 1,"\"",""),
                                   phenoFilters(), 
                                   'if'(lengths(phenoFilters()) == 1,"\"",""),
                                   collapse= ","),")"),
                     "\nnames(filtList) <- ",
                     paste0("c(",
                            paste0("\"",
                                   names(phenoFilters()),
                                   "\"",
                                   collapse=", "),
                            ")", collapse = ", "),
                     "\nmeData <- filterByPheno(meData, filtList)"
              )
            )
          }
        }
      }
    }
  })
  
  ## reset any filtering operations: check to proceed
  observeEvent(input$resetbutton, {
    showModal(modalDialog(
      tagList(
        p("Reset will set data back to last saved version. All 
                  modifications will be lost.\nPlease confirm reset!")
      ),
      title = "Resetting Data",
      footer = tagList(
        actionButton(ns("confirmReset"), "Reset"),
        modalButton("Cancel")
      )
    ))
  })
  
  ## proceed approved: perform reset
  observeEvent(input$confirmReset, {
    colorChoice("No color")
    logChoice("none")
    minFeats(2)
    minReads(2)
    featureSams(1)
    qcCode(NULL)
    updateSliderInput(session, "featuresams", value = 1, min = 1)
    updateSliderInput(session, "minfeats", value = 2, min = 2)
    updateSliderInput(session, "minreads", value = 2, min = 2)
    removeUI(".phenoremoverow", multiple = TRUE, immediate = TRUE)
    dataFilterRep(NULL)
    phenoFilters(NULL)
    meData(origData())
    filterChoices(getFilterChoices(meData()))
    remainingFilters(filterChoices())
    updateSelectizeInput(session, "removepheno",
                         choices = names(filterChoices()), 
                         selected = "",
                         options = list(placeholder = "Select Phenotype"))
    updateSelectInput(session, "qc_color", selected = "No color")
    plotData(meData())
    removeModal()
  })
  
  ## Render header
  output$projecthead <- renderUI({
    if(!is.null(meData())){
      div(
        h2(id = "projecttitle", meName()),
        br(),
        div(
          p("Performing quality control (QC) is an important step prior to downstream analyses. 
          The following QC plots highlight the distribution of present features and library sizes 
          as well as the relationship between these two. 
          Colors can be added based on phenotype information, e.g. sequencing or library preparation batches.
          Data can be subset based on phenotype information or filtered to feature/library constraints.
          It is important for the user to consider the appropriate filtering approaches. Downstream analyses 
          require normalization to account for sequencing variability. Currently, 
          proportions and cumulative sum scaling (CSS) are implemented (Paulson, 2013).
          To enable inputs in the analysis sections, aggregate to a specific taxonomy level. If no feature data
          is available, choose \"unavailable\" which will allow analysis at the raw counts level.")
        ),
        h4("QC OVERVIEW")
      )
    } else{
      div(id = "datahelp",
          h2("Data input"),
          div(
            h5("Microbiome Explorer accepts several different data upload formats:"),
            tags$li("MRexperiment-class objects stored as RDATA or RDS files"),
            tags$li("The Biological Observation Matrix (BIOM) formatted files"),
            tags$li("Raw counts files"),
            p(""),
            p("A ",strong("counts"), " file is required and can be uploaded in delimited form (csv, tsv). The required format
             is such that each ",strong("sample is a column"), " and each unique ",strong("feature is a row"), " in the data set.
            The counts matrix should have no row names and store the information on the features/OTUs in its first columns. 
            All other columns names should correspond to sample names."),
            p("A tabular ", strong("phenotype"), "file can be linked to the counts data. Each ", strong("sample is a row"),
              "with the names of the rows in the phenotype data corresponding to the names of the columns in the count
            data. Appending several phenotype files by subsequent uploads is possible.
            If no phenotype file is given, the names of the columns of the count data are
            used as the phenotype data."),
            p("A ", strong("feature"), " data file must be provided if aggregation to a particular phylogenetic level is desired.
            Each unique ",strong("feature is a row"), " which must correspond to the ones in the counts data.
            Each column is a taxonomy level.
            If the feature data file is omitted, analysis can only be done at the raw counts level."),
            p("An MRexperiment-class object from metagenomeSeq already combines these different data files into one data structure based on 
           an extended eSet-class.") 
          )
      )
    }
  })
  
  ## render plot options
  output$plotoptions <- renderUI({
    req(plotData())
    selection <- isolate('if'(colorChoice() != "", 
                              colorChoice(), "No Color"))
    shinyWidgets::dropdownButton(
      tags$h3("Plot Options"),
      selectInput(
        ns("qc_color"),
        label = "Color By",
        choices = c("No color", isolate(facetOptions())),
        selected = selection,
        multiple = FALSE,
        selectize = FALSE
      ),      
      sliderInput(
        inputId = ns("plotWidth"),
        label = "Adjust plot width",
        value = 550,
        min = 250,
        max = 1600,
        round = TRUE
      ),
      shinyWidgets::radioGroupButtons(
        inputId = ns("logoptions"),
        label = "Log scale on",
        choices = c("none", 
                    "x axis", 
                    "y axis",
                    "both"),
        status = "primary",
        checkIcon = list(
          yes = icon("ok", 
                     lib = "glyphicon"),
          no = icon("remove",
                    lib = "glyphicon"))
      ),
      circle = FALSE, status = "danger", icon = icon("gear"), 
      width = "315px",
      label = "Plot Options"
    )
  })
  
  ## adjust color of plot
  observeEvent(input$qc_color, {
    colorChoice(input$qc_color)
    updateSliderInput(session, "minreads",
                      value = input$minreads - 1)
  }, ignoreInit = TRUE)
  
  ## store qc code in reactive Value
  observeEvent(input$reportButton, {
    qcRep(qcCode())
  })
  
  ## update plot based on min features or min reads
  observeEvent(c(input$minfeats, input$minreads), {
    req(plotData(), input$minfeats, input$minreads)
    if (nrow(plotData()) > 0) {
      minf <- input$minfeats
      minr <- input$minreads
      featthresh <- 2
      readthresh <- 2
      if(logChoice() %in% c("both","x axis"))
        readthresh <- log(2)
      if(logChoice() %in% c("both","y axis"))
        featthresh <- log(2)
      featalpha <- 0
      featalphafill <- 0
      if (minf > featthresh) {
        featalpha <- 1
        featalphafill <- 0.2
      }
      readalpha <- 0
      readalphafill <- 0
      if (minr > readthresh) {
        readalpha <- 1
        readalphafill <- 0.2
      }
      
      plotly::plotlyProxy("dataqcplot", session) %>%
        plotly::plotlyProxyInvoke("deleteTraces", 0)
      plotly::plotlyProxy("dataqcplot", session) %>%
        plotly::plotlyProxyInvoke(
          "addTraces",
          list(
            x = c(0, max(colSums(MRcounts(plotData())))),
            y = c(minf, minf),
            type = "scatter",
            mode = "lines",
            line = list(
              color = paste0("rgba(119, 0, 0,",
                             featalpha, ")"),
              width = 1, dash = "dash"
            ),
            fill = "tozeroy",
            fillcolor = paste0("rgba(219,219,219,",
                               featalphafill, ")")
          ),
          0
        )
      
      plotly::plotlyProxy("dataqcplot", session) %>%
        plotly::plotlyProxyInvoke("deleteTraces", 1)
      plotly::plotlyProxy("dataqcplot", session) %>%
        plotly::plotlyProxyInvoke(
          "addTraces",
          list(
            y = c(0, max(colSums(MRcounts(plotData()) > 0))),
            x = c(minr, minr),
            type = "scatter",
            mode = "lines",
            line = list(
              color = paste0("rgba(119, 0, 0,",
                             readalpha, ")"),
              width = 1, dash = "dash"
            ),
            fill = "tozerox",
            fillcolor = paste0("rgba(219,219,219,",
                               readalphafill, ")")
          ),
          1
        )
    }
  }, ignoreInit = TRUE)
  
  observeEvent(input$plotWidth, {
    req(plotData())
    req(nrow(plotData()) > 0)
    plotly::plotlyProxy("dataqcplot", session) %>%
      plotly::plotlyProxyInvoke(
        "relayout",
        list(width = input$plotWidth)
      )
  })
  
  observeEvent(input$barWidth, {
    req(plotData())
    req(nrow(plotData()) > 0)
    plotly::plotlyProxy("sampleotus", session) %>%
      plotly::plotlyProxyInvoke(
        "relayout",
        list(width = input$barWidth)
      )
  })
  
  
  ## render the QC plot
  output$dataqcplot <- plotly::renderPlotly({
    req(plotData())
    req(nrow(plotData()) > 0)
    minf <- minFeats()
    minr <- minReads()
    featthresh <- 2
    readthresh <- 2
    if(logChoice() %in% c("both","x axis"))
      readthresh <- log(2)
    if(logChoice() %in% c("both","y axis"))
      featthresh <- log(2)
    if (minf <= featthresh) {
      minf <- 0
    }
    if (minr <= readthresh) {
      minr <- 0
    }
    qcCode(paste(
      paste0("#' ## QC plots"),
      paste0("", "#+ qcplot, fig.width = 7"),
      paste0("makeQCPlot(meData, col_by = \"", colorChoice(), "\","),
      paste0("\tlog = \"", logChoice(),"\","),
      paste0("\tfilter_feat = ", minf, ","),
      paste0("\tfilter_read = ", minr, ","),
      paste0("\tallowWebGL = FALSE)"),
      paste0("\n"),            
      paste0("", "#+ qcsample, fig.width = 7"),
      paste0("plotlySampleBarplot(meData,"),
      paste0("\tcol_by = \"",input$bar_color, "\","),
      paste0("\tpwidth = 600,"),
      paste0("\tsortbyfreq = ",input$sortbyfreq, ","),
      paste0("\tpheno_sort = ",
             'if'(is.null(sortPheno()),"NULL",
                  paste0("\"",sortPheno(),"\"")), ","),
      paste0("\tx_levels = ",
             'if'(is.null(phenoValues()),"NULL",
                  paste0("c(", paste0("\"",
                                      phenoValues(), 
                                      "\"",
                                      collapse = ","),")")),
             ")\n"),
      paste0("#'", ""),
      paste0(""),
      paste0("", "#+ histogram1, fig.width = 3, fig.height = 3"),
      paste0("plotlyHistogram(histvalue = colSums(MRcounts(meData)>0),"),
      paste0("\tplotTitle = \"No. Present Features\","),
      paste0("\txaxisTitle = \"features\","),
      paste0("\tyaxisTitle = \"frequency\")"),
      paste0("\n"),
      paste0("", "#+ histogram2, fig.width = 3, fig.height = 3"),
      paste0("plotlyHistogram(histvalue = colSums(MRcounts(meData)),"),
      paste0("\tplotTitle = \"Library Size\","),
      paste0("\txaxisTitle = \"Library Size\","),
      paste0("\tyaxisTitle = \"frequency\")"),
      paste0("#'", ""),
      paste0(""),
      sep = "\n"))
    
    
    microbiomeExplorer::makeQCPlot(plotData(),
                                   col_by = colorChoice(),
                                   log = logChoice(),
                                   filter_feat = minf,
                                   filter_read = minr,
                                   pwidth = isolate(input$plotWidth))
    
  })
  
  
  ## renders present features histogram; updates with meData
  output$featpres <- plotly::renderPlotly({
    req(plotData())
    if (nrow(plotData()) > 0) {
      microbiomeExplorer::plotlyHistogram(
        histvalue = colSums(MRcounts(plotData()) > 0),
        plotTitle = "Feature distribution",
        xaxisTitle = "features",
        yaxisTitle = "frequency"
      )
    }
  })
  
  ## renders library size histogram; updates with meData
  output$libsize <- plotly::renderPlotly({
    req(plotData())
    if (nrow(plotData()) > 0) {
      microbiomeExplorer::plotlyHistogram(
        histvalue = colSums(MRcounts(plotData())),
        plotTitle = "Read distribution",
        xaxisTitle = "library size",
        yaxisTitle = "frequency"
      )
    }
  })
  
  ## render plot options
  output$samplebaroptions <- renderUI({
    req(plotData())
    selection <- isolate('if'(colorChoice() != "", 
                              colorChoice(), "No Color"))
    shinyWidgets::dropdownButton(
      tags$h3("Plot Options"),
      selectInput(
        ns("bar_color"),
        label = "Color By",
        choices = c("No color", isolate(facetOptions())),
        selected = selection,
        multiple = FALSE,
        selectize = FALSE
      ),      
      sliderInput(
        inputId = ns("barWidth"),
        label = "Adjust plot width",
        value = 550,
        min = 250,
        max = 1600,
        round = TRUE
      ),
      shinyWidgets::switchInput(
        inputId = ns("sortbyfreq"),
        label = "Sort by frequency",
        size = "mini",
        value = FALSE,
        labelWidth = "100px"
      ),
      shinyWidgets::switchInput(
        inputId = ns("sortbypheno"),
        label = "Sort by phenotype",
        size = "mini",
        value = FALSE,
        labelWidth = "100px"
      ),
      shinyjs::hidden(selectInput(
        ns("sort_pheno"),
        label = "Sort By",
        choices = isolate(facetOptions()),
        multiple = FALSE,
        selectize = FALSE
      )),
      shinyjs::hidden(selectizeInput(
        ns("pheno_values"),
        label = "in order of",
        choices = "",
        options = list(placeholder = "Select levels"),
        selected = "",
        multiple = TRUE
      )),
      shinyjs::hidden(actionButton(
        ns("changeFactorOrder"), 
        label = "reorder", 
        width = "100px")),
      circle = FALSE, status = "danger", icon = icon("gear"), 
      width = "315px",
      label = "Plot Options"
    )
  })
  
  output$sampleotus <- plotly::renderPlotly({
    req(plotData())
    
    if(nrow(plotData()) > 0){
      microbiomeExplorer::plotlySampleBarplot(
        plotData(),
        col_by = input$bar_color,
        pwidth = 600,
        sortbyfreq = input$sortbyfreq,
        pheno_sort = sortPheno(),
        x_levels = phenoValues())
    }
  })
  
  observe({
    if(isTRUE(input$sortbyfreq)){
      shinyjs::hide("sort_pheno")
      shinyjs::hide("pheno_values")
      shinyjs::hide("changeFactorOrder")
      shinyWidgets::updateSwitchInput(session, "sortbypheno", value = FALSE)
    }
    else{
      if(isTRUE(input$sortbypheno)){
        shinyjs::show("sort_pheno")
        shinyjs::show("pheno_values")
        shinyjs::show("changeFactorOrder")
      } else {
        shinyjs::hide("sort_pheno")
        shinyjs::hide("pheno_values")
        shinyjs::hide("changeFactorOrder")
      }
    }
  })
  
  observeEvent(input$sort_pheno,{
    updateSelectizeInput(session,"pheno_values", 
                         choices = unique(plotData()[[input$sort_pheno]]))
  }, ignoreNULL = TRUE)
  
  observeEvent(input$changeFactorOrder,{
    sortPheno('if'(length(input$sort_pheno) > 0, input$sort_pheno,NULL))
    phenoValues('if'(length(input$pheno_values) > 0, input$pheno_values,NULL))
  })
  
  observeEvent(input$sortbypheno,{
    if(isFALSE(input$sortbypheno)){
      sortPheno(NULL)
      phenoValues(NULL)
    } else {
      shinyWidgets::updateSwitchInput(session, "sortbyfreq", value = FALSE)
    }
  })
  
  ## download modified data as RDS
  output$savebutton <- downloadHandler(
    filename = function() {
      fext <- tools::file_ext(meName())
      if(fext == ""){
        paste0(meName(),".rds")
      } else {
        gsub(fext, "rds", meName())
      }
    },
    content = function(file) {
      saveRDS(meData(), file)
    }
  )
  
  adonisOptions <- reactive({
    fC <- filterChoices()
    fC[!(names(fC) %in% "SAMPLE_ID")]
  })
  
  # current data after filtering is returned and passed to other tabs
  return(list(meData = meData, 
              facetOptions = facetOptions, 
              shapeOptions = shapeOptions,
              adonisOptions = adonisOptions,
              qcCode = qcCode, 
              dataName = meName))
}


