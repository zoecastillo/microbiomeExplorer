

#' report tab ui
#'
#' @param id namespace identifier
#' 
#' @return fluidRow holding ui elements
#'
#' @author Janina Reeder
#' 
#' @examples reportListUI("reportlist_id")
#' 
#' @export
reportListUI <- function(id) {
    ns <- NS(id)
    
    fluidRow(
        column(
            width = 4,
            h4("REPORT SETTINGS"),
            textInput(
                ns("repfile"),
                label = "File name",
                value = "me_report"
            ),
            textInput(
                ns("reptitle"),
                label = "Report title",
                value = "MicrobiomeExplorer Report"
            ),
            textInput(
                ns("repauthor"),
                label = "Author",
                value = ""
            ),
            textAreaInput(
                ns("intro_text"), "Introductory text (optional)",
                placeholder = "Any text entered here will be presented at the 
                top of the generated report"
            ),
            shinyWidgets::prettyCheckbox(
                inputId = ns("repcontents"),
                label = "Include table of contents",
                value = TRUE,
                icon = icon("check"),
                status = "success",
                animation = "rotate",
                outline = TRUE
            ),
            h4("REPORT FORMAT", style = "padding-bottom: 0;"),
            shinyWidgets::prettyCheckboxGroup(
                inputId = ns("repformat"),
                label = "",
                shape = "square",
                icon = icon("check"),
                choices = getOption("me.reportformatshort"),
                selected = "HTML",
                status = "success",
                animation = "rotate",
                outline = TRUE
            ),
            fluidRow(
                width = 12,
                actionButton(
                    ns("generatebutton"), 
                    icon = icon("fas fa-file-alt"),
                    label = HTML("&nbsp;GENERATE"), 
                    width = "110px"),
                shinyjs::disabled(
                    actionButton(ns("exportbutton"), 
                                 icon = icon("fas fa-download"),
                                 label = HTML("&nbsp;EXPORT"), 
                                 width = "90px"))
            ),
            fluidRow(
                width = 12,
                downloadButton(ns("downloadbutton"), 
                               label = "Download", 
                               style = "visibility: hidden;")
            )
        ),
        column(
            width = 8,
            h2("REPORT CONTENTS"),
            uiOutput(ns("dataRowBox"), width = "90%"),
            uiOutput(ns("qcRowBox"), width = "90%"),
            uiOutput(ns("reportListBox"), width = "90%")
        )
    )
}

#' Report tab module server
#'
#' @param input module input
#' @param output module output
#' @param session app session
#' @param dataSource R code to obtain data for rendering
#' @param preprocessRep R code containing preprocessing steps of data
#' @param qcRep R Code to generate QC plots
#' @param analysisRep R Code to generate all analyses saved to reports
#' @param aggIndex boolean value representing aggregation steps in analysisRep
#' @param reset boolean reactive which resets the module if TRUE
#' 
#' @return report list server fragment - no return value
#' 
#' @author Janina Reeder
reportList <- function(input, output, session, 
                       dataSource,
                       preprocessRep, 
                       qcRep, 
                       analysisRep, 
                       aggIndex,
                       reset) {
    ns <- session$ns
    
    ## stores info on whether qc plots should be included in final report
    qcSelected <- reactiveVal(FALSE)
    ## stores info on which analysis rows should be included in final report
    selectedRows <- reactiveVal(NULL)
    
    
    dataCode <- reactive({
        req(preprocessRep(), dataSource())
        paste0(dataSource(),
               paste0(preprocessRep(), collapse = "\n"),
               sep = "\n\n")
    })
    
    ## server code for data row module
    observe({
        req(dataCode())
        callModule(reportRow, "datarow", type = TRUE, content = dataCode)
    })
    
    ## server code for qc plot module
    observe({
        req(qcRep())
        qcSelected(callModule(reportRow, "qcrow", 
                              type = FALSE, content = qcRep))
    })
    
    ## ui call for data row
    output$dataRowBox <- renderUI({
        req(preprocessRep())
        div(
            box(
                width = 11,
                reportRowUI(ns("datarow"), type = TRUE)
            )
        )
    })
    
    ## ui call for qc plot row
    output$qcRowBox <- renderUI({
        req(qcRep())
        div(
            box(
                width = 11,
                reportRowUI(ns("qcrow"), type = FALSE)
            )
        )
    })
    
    ## calling server rows for all analysis code parts
    observe({
        req(aggIndex(), analysisRep())
        selectedRows(lapply(
            seq_len(length(analysisRep())),
            function(i) callModule(reportRow, paste0("rr", i), 
                                   type = aggIndex()[i], 
                                   content = reactive(analysisRep()[i]))
        ))
    })
    
    ## calling ui code for all analysis code parts
    output$reportListBox <- renderUI({
        req(aggIndex(), analysisRep())
        div(box(
            width = 11,
            lapply(seq_len(length(analysisRep())), 
                   function(i) reportRowUI(ns(paste0("rr", i)), 
                                           type = aggIndex()[i]))
        ))
    })
    
    ## get a temp directory to use for rendering reports
    report_dir <- tempdir()
    ## reactive keeping track of whether the report is ready
    report_ready <- reactiveVal(FALSE)
    fileName <- reactive(gsub(" ","",input$repfile))
    
    ## prepare R code for rendering
    myreport <- reactive({
        mr <- list(unname(unlist(preprocessRep())))
        if (qcSelected()) {
            mr <- append(mr, qcRep())
        }
        if (!is.null(selectedRows())) {
            mr <- append(
                mr, 
                lapply(
                    seq_len(length(selectedRows())), 
                    function(si) {
                        if (is.null(selectedRows()[[si]])) {
                            unlist(analysisRep()[si])
                        } else if (isTRUE(selectedRows()[[si]])) {
                            unlist(analysisRep()[si])
                        } else {
                            NULL
                        }
                    }))
        }
        mr <- mr[vapply(mr, function(m) !is.null(m), logical(1))]
        labelchunks <- vapply(mr, function(m) {
            sum(stringr::str_count(m, pattern = "#-"))
        }, numeric(1))
        sumofchunks <- sum(labelchunks > 0)
        chunklocs <- mr[labelchunks == 1]
        ## each chunk should have a unique id
        mr[labelchunks == 1] <- vapply(seq_len(sumofchunks), function(i) {
            stringr::str_replace(chunklocs[i], "#-", 
                                 paste0("#+ codechunk", i, ","))
        }, character(1))
        mr
    })
    
    chosenFormat <- reactive({
        vapply(input$repformat,
               switch,
               "HTML" = "html_document",
               "PDF" = "pdf_document",
               "DOC" = "word_document",
               "PPT" = "powerpoint_presentation",
               FUN.VALUE = character(1))
    })
    
    ## start building the reports
    observeEvent(input$generatebutton, {
        report_ready(FALSE)
        req(chosenFormat())
        shinyjs::disable("exportbutton")
        
        showModal(modalDialog(
            title = "Creating Report",
            p("Report is rendering, this may take a moment. To download 
              generated report, click Export when report is ready")
        ))
        
        cat(file = stderr(), paste("-- Sending report to", report_dir, "--\n"))
        
        file.remove(list.files(report_dir, paste0(fileName(), "[.]"),
                               full.names = TRUE
        ))
        try(
            file.remove(
                file.path(report_dir, paste0(fileName(), c(".Rmd", ".html")))
            ),
            silent = TRUE
        )
        cat(file = stderr(), paste("-- Removing", 
                                   file.path(report_dir, 
                                             paste0(fileName(), 
                                                    ".Rmd", "--\n"))))
        
        ## main function generating report
        x <- try(generateReport(myreport(),
                                filename = fileName(),
                                dir = report_dir,
                                title = input$reptitle,
                                author = input$repauthor,
                                date = format(Sys.time(), "%d %B, %Y"),
                                data.source = dataSource(),
                                output = chosenFormat(),
                                intro_text = input$intro_text,
                                toc = input$repcontents),
                 silent = TRUE
        )
        
        cat(file = stderr(), "gen report try value:", class(x))
        removeModal()
        if (is(x, "try-error")) {
            showNotification(duration = 3, type = "error", 
                             ui = "Report could not be generated")
            validate(
                need(!is(x, "try-error"), 
                     message = "Report could not be generated")
            )
        } else {
            report_ready(TRUE)
            shinyjs::enable("exportbutton")
        }
    })
    
    ## downloding generated report (this button is invisible)
    output$downloadbutton <- downloadHandler(
        filename = function() {
            paste0(fileName(), ".zip")
        },
        content = function(file) {
            paths <- c(
                file.path(report_dir, paste0(fileName(), ".Rmd")),
                vapply(
                    chosenFormat(), switch,
                    "html_document" = file.path(
                        report_dir, 
                        paste0(fileName(), ".html")),
                    "pdf_document" = file.path(
                        report_dir, 
                        paste0(fileName(), ".pdf")),
                    "word_document" = file.path(
                        report_dir, 
                        paste0(fileName(), ".docx")),
                    "powerpoint_presentation" = file.path(
                        report_dir, 
                        paste0(fileName(), ".pptx")),
                    FUN.VALUE = character(1)
                )
            )
            paths <- c(paths, dir(file.path(report_dir, 
                                            paste0(fileName(), "_files"), 
                                            "figure-markdown_strict"),
                                  pattern = ".png",
                                  full.names = TRUE
            ))
            utils::zip(file, paths)
        }
    )
    
    ## call invisible download button
    observeEvent(input$exportbutton, {
        path <- file.path(report_dir, paste0(fileName(), ".Rmd"))
        if (!file.exists(path)) {
            showNotification(ui = "Could not download report", 
                             type = "error", duration = 3)
            stop()
        }
        req(file.exists(path))
        shinyjs::runjs("document.getElementById('reportlist-downloadbutton').
                       click();")
    })
}
