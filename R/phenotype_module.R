

#' Phenotype table UI module
#'
#' @param id namespace identifier
#'
#' @author Janina Reeder
#'
#' @return fluidRow holding the ui code
#' 
#' @examples phenotypeTableUI("phenotype_id")
#' 
#' @export
phenotypeTableUI <- function(id) {
  ns <- NS(id)
  
  fluidRow(
    column(width = 3, id = ns("phenointeractioncol"),
           h4("COMBINE PHENOTYPE COLUMNS"),
           selectInput(
             ns("pheno1"),
             label = "Phenotype 1",
             choices = "", multiple = FALSE, selectize = FALSE, 
             width = "250px"),
           selectInput(
             ns("pheno2"),
             label = "Phenotype 2",
             choices = "", multiple = FALSE, selectize = FALSE, 
             width = "250px"),
           textInput(
             ns("interactionname"),
             label = "Column Name",
             value = "", width = "250px"),
           fluidRow(
             width = 12, id = "actionbuttonrow",
             actionButton(
               ns("createbutton"), 
               icon = icon("fas fa-plus"),
               label = HTML("&nbsp;ADD"), width = "100px"),
             actionButton(
               ns("savebutton"), 
               icon = icon("far fa-save"),
               label = HTML("&nbsp;SAVE"), width = "100px")
           ),
           br(),
           box(
             width = 12,
             title = "ADJUST DATATYPES (optional)",
             selectInput(
               ns("phenotype"),
               label = "Adjust Phenotype",
               choices = "", multiple = FALSE, selectize = FALSE, 
               width = "250px"),
             selectInput(
               ns("phenodatatype"),
               label = "Change data type",
               choices = "", multiple = FALSE, selectize = FALSE, 
               width = "250px"),
             fluidRow(
               width = 12, id = "actionbuttonrow2",
               actionButton(
                 ns("changebutton"), 
                 icon("fas fa-exchange-alt"),
                 label = HTML("&nbsp;ADJUST"), width = "90px"),
               actionButton(
                 ns("savebutton2"),
                 icon = icon("far fa-save"),
                 label = HTML("&nbsp;SAVE"), width = "90px")
             ),
             collapsed = TRUE,
             collapsible = TRUE
           )
    ),
    column(
      width = 9,
      box(
        width = 10,
        h2("PHENOTYPE INFORMATION"),
        p("Create new interaction phenotypes by combining two columns and adjust  
         column data types if needed. New phenotype information can be added  
         via \"Load & Filter\" tab. Table settings allow paging through sections
          of the data, choosing how many entries to display or searching for 
         specific entries.Modifications must be saved in order to be available 
         in the analysis sections.")
      ),
      box(
        width = 12, 
        div(
          id = ns("downloaddiv"),
          downloadButton(ns("download_button"),"Download"),
          DT::DTOutput(ns("phenodatatable"), width = "100%")
        )
      )
    )
  )
}

#' Helper function returning the code used to modify the phenotable as a string
#'
#' @param name interaction name
#' @param pheno1 first interaction phenotype
#' @param pheno2 second interaction phenotype
#' 
#' @author Janina Reeder
#'
#' @return String storing code to perform modification
getPhenoModCode <- function(name, pheno1, pheno2) {
  paste(paste0("\nnew_pheno <- interaction(pData(meData)[,c(\"",
               pheno1, "\",\"", pheno2, "\")])"),
        paste0("mutatedRows <- row.names(pData(meData))"),
        paste0("mutatedData <- dplyr::mutate(pData(meData), \"", 
               name, "\" = new_pheno)"),
        paste0("row.names(mutatedData) <- mutatedRows"),
        paste0("meData <- addPhenoData(meData,mutatedData)"),
        sep = "\n"
  )
}

adjustPhenoColCode <- function(selectedCols){
  paste(paste0("\nmutatedRows <- row.names(pData(meData))"),
        paste0("mutatedData <- pData(meData)[, ",
               paste0("c(",paste0(selectedCols, collapse=", "),
                      ")", collapse = ", "), "]"),
        paste0("rownames(mutatedData) <- mutatedRows"),
        paste0("meData <- addPhenoData(meData,mutatedData)"),
        sep = "\n")
}

#' Helper function returning the code used to modify the data types of the
#' pheno table
#'
#' @param phenotype name of the phenotype column header
#' @param datatype variable type to assign to column
#' 
#' @author Janina Reeder
#'
#' @return String storing code to perform modification
getPhenoChanges <- function(phenotype, datatype){
  paste(paste0("\nnew_pheno <- pData(meData)"),
        paste0("pfdata <- new_pheno[,\"",phenotype,"\"]"),
        paste0("pfdata <- switch(\"",datatype,"\","),
        paste0("\t\"factor\" = as.factor(pfdata),"),
        paste0("\t\"integer\" = as.integer(pfdata),"),
        paste0("\t\"numeric\" = as.numeric(pfdata),"),
        paste0("\t\"logical\" = as.logical(pfdata))"),
        paste0("new_pheno[,\"",phenotype,"\"] <- pfdata"),
        paste0("meData <- addPhenoData(meData,new_pheno)\n\n"),
        sep = "\n"
  )
}



#' Phenotype table server module
#'
#' @param input shiny input
#' @param output shiny output
#' @param session shiny session
#' @param meData MRExperiment storing the data
#' @param phenoModRep reactive Value storing any phenotable modifications made
#' @param addPheno reactive boolean keeping track of pheno data modifications
#'
#' @return phenotype table server fragment - no return value
#'
#' @author Janina Reeder
#' 
#' @importFrom Biobase pData
#' @importFrom rlang :=
phenotypeTable <- function(input, output, session, meData, 
                           phenoModRep, addPheno) {
  ns <- session$ns
  variable <- NULL
  ## stores the currently shown phenotype table
  phenoFrame <- reactiveVal(NULL)
  colTypes <- reactiveVal(NULL)
  ## captures whether any changes were made
  mutated <- reactiveVal(FALSE)
  phenoChangeCode <- reactiveVal(NULL)
  phenoInteractCode <- reactiveVal(NULL)
  phenoColChangeCode <- reactiveVal(NULL)
  phenoChanged <- reactiveVal(FALSE)
  
  ## initializes the local pheno frame and input fields
  observe({
    req(meData())
    isolate({
      phenoFrame(pData(meData()))
      colTypes(tibble::as_tibble(phenoFrame()) %>%
                 dplyr::mutate_if(lubridate::is.POSIXt, as.Date) %>%
                 dplyr::summarise_all(class) %>% 
                 tidyr::gather(variable, class))
    })
    updateSelectInput(session, "pheno1", choices = names(pData(meData())))
    updateSelectInput(session, "pheno2", choices = names(pData(meData())))
    updateSelectInput(session, "phenotype", 
                      choices = names(pData(meData())))
  })
  
  observeEvent(input$phenotype,{
    req(colTypes())
    df <- colTypes()
    updateSelectInput(session, "phenodatatype",
                      choices = c("factor","integer","numeric","logical", 
                                  "Date"),
                      selected = df[df$variable %in% input$phenotype,
                                    "class"])
  })
  
  observeEvent(input$changebutton,{
    shinyjs::disable("savebutton2")
    df <- phenoFrame()
    pfdata <- phenoFrame()[,input$phenotype]
    pfdata <- switch(input$phenodatatype,
                     "factor" = as.factor(pfdata),
                     "integer" = as.integer(pfdata),
                     "numeric" = as.numeric(pfdata),
                     "logical" = as.logical(pfdata),
                     "Date" = tryCatch(as.Date(pfdata), 
                                       error = function(e) pfdata))
    df[,input$phenotype] <- pfdata
    phenoFrame(df)
    phenoChangeCode(c(phenoChangeCode(),
                      getPhenoChanges(phenotype = input$phenotype,
                                      datatype = input$phenodatatype)))
    shinyjs::enable("savebutton2")
  })
  
  ## add new interaction column to phenotable
  observeEvent(input$createbutton, {
    req(phenoFrame())
    shinyjs::disable("savebutton")
    mustabort <- FALSE
    
    if (is.null(input$pheno1) | is.null(input$pheno2) | 
        input$pheno1 == input$pheno2) {
      showModal(modalDialog(
        title = "Error Creating New Phenotype",
        "Two different phenotypes must be selected in order to create 
                a new one. \nPlease revise.",
        easyClose = TRUE))
      mustabort <- TRUE
    }
    interactname <- parseInteractionName(input$interactionname)
    if (startsWith(interactname, "ERROR")) {
      showModal(modalDialog(
        title = "Error Creating New Phenotype",
        paste0("The interaction name is not valid: ", 
               gsub("ERROR: ", "", interactname), "\nPlease revise."),
        easyClose = TRUE))
      mustabort <- TRUE
    }
    if (startsWith(interactname, "WARNING")) {
      showModal(modalDialog(
        title = "Warning",
        paste0("The interaction name contains invalid characters which 
                       have been adjusted: ", 
               gsub("WARNING: ", "", interactname)),
        easyClose = TRUE))
    }
    if (!mustabort) {
      new_pheno <- interaction(phenoFrame()[, c(input$pheno1, 
                                                input$pheno2)])
      ## store old row names to copy over to modified frame
      mutatedRows <- rownames(phenoFrame())
      mutatedData <- dplyr::mutate(phenoFrame(), 
                                   !!interactname := new_pheno)
      rownames(mutatedData) <- mutatedRows
      phenoInteractCode(getPhenoModCode(interactname,input$pheno1,input$pheno2))
      ## mark that table has been modified
      mutated(TRUE)
      ## display modified phenotype table
      phenoFrame(mutatedData)
    }
    shinyjs::enable("savebutton")
  })
  
  ## make changes permanent
  observeEvent(c(input$savebutton, input$savebutton2), {
    req(meData())
    selectedCols <- vapply(
      seq_along(input$phenodatatable_state$columns),
      function(i) input$phenodatatable_state$columns[[i]]$visible,
      logical(1))
    selectedCols <- selectedCols[-1]
    if (sum(selectedCols) != ncol(phenoFrame())) {
      ## store rownames to pass over to modified frame
      mutatedRows <- rownames(phenoFrame())
      mutatedData <- phenoFrame()[, selectedCols]
      rownames(mutatedData) <- mutatedRows
      phenoColChangeCode(adjustPhenoColCode(selectedCols))
      ## mark as modified
      mutated(TRUE)
      ## store modification
      phenoFrame(mutatedData)
    }
    if (mutated()) {
      phenoChanged(TRUE)
      phenoModRep(paste0(c(phenoModRep(),
                           phenoInteractCode()),
                         collapse = "\n"))
      phenoInteractCode(NULL)
      mutated(FALSE)
    }
    if (!is.null(phenoChangeCode())) {
      phenoChanged(TRUE)
      phenoModRep(paste0(c(phenoModRep(),
                           phenoChangeCode()),
                         collapse = "\n"))
      phenoChangeCode(NULL)
    }
    if (!is.null(phenoColChangeCode())) {
      phenoChanged(TRUE)
      phenoModRep(paste0(c(phenoModRep(),
                           phenoColChangeCode()),
                         collapse = "\n"))
      phenoColChangeCode(NULL)
    }
    if(phenoChanged()){
      ## pass modifications over to base data structure
      meData(addPhenoData(meData(), phenoFrame()))
      addPheno(TRUE)
      phenoChanged(FALSE)
    }
  })
  
  ## render pheno data as a table
  output$phenodatatable <- DT::renderDT({
    req(phenoFrame())
    shinyjs::js$moveButton(ns("downloaddiv"),ns("download_button"))
    DT::datatable(
      data = phenoFrame(), 
      class = "stripe hover cell-border order-column",
      filter = "none", style = "bootstrap",
      callback = DT::JS("$('div.dwnld_pheno').append($('#phenotab-download_button'));"),
      extensions = c("FixedColumns", "ColReorder", "Buttons"),
      plugins = "ellipsis",
      options = list(
        scrollX = TRUE,
        fixedColumns = list(leftColumns = 1, rightColumns = 0),
        paging = TRUE,
        colReorder = TRUE,
        stateSave = TRUE,
        stateLoadParams = DT::JS("function (settings, data) 
                                         {return false;}"),
        columnDefs = list(list(
          targets = "_all",
          render = DT::JS("$.fn.dataTable.render.ellipsis( 50 , false)"))),
        buttons = list(
          list(
            extend = "colvis",
            text = "Select columns"
          )), 
        dom = '<"dwnld_pheno"B>ftlp'
      )
    )
  })
  
  
  
  
  output$download_button <- downloadHandler(
    filename = function() {
      paste0("phenodata-", Sys.Date(), ".csv")
    },
    content = function(file) {
      readr::write_csv(phenoFrame(), file)
    }
  )
}
