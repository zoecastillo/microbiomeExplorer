

#' Differential Analysis module UI
#'
#' @param id namespace identifier
#' 
#' @author Janina Reeder
#' 
#'
#' @return row containing the UI elements
#'
#' @export
diffTableUI <- function(id) {
  ns <- NS(id)
  
  fluidRow(
    box(title = "DIFFERENTIAL ANALYSIS", solidHeader = TRUE, 
        collapsible = TRUE, width = 12,
        fluidRow(
          column(width = 1),
          column(width = 10,
                 div(id = ns("downloaddiv"),
                     shinyjs::hidden(downloadButton(ns("download_button"),
                                                    "Download")),
                     DT::DTOutput(ns("diffdatatable"), width = "100%")
                 )
          )
        )
    ),
    box(width = 12,
        fluidRow(
          column(width = 1),
          column(width = 10,
                 plotly::plotlyOutput(ns("clickedFeature"), 
                                      width = "100%", height = "100%"),
                 br(),
                 shinyjs::disabled(shinyWidgets::dropdownButton(tags$h3("Plot Options"),
                                                                shinyWidgets::switchInput(
                                                    inputId = ns("sP"),
                                                    label = "Show Points",
                                                    size = "mini",
                                                    labelWidth = "80px"
                                                  ),
                                                  shinyWidgets::switchInput(
                                                    inputId = ns("logS"),
                                                    label = "Log Scale",
                                                    value = TRUE,
                                                    size = "mini",
                                                    labelWidth = "80px"
                                                  ),
                                                  circle = FALSE, 
                                                  status = "danger", 
                                                  icon = icon("gear"), 
                                                  width = "300px",
                                                  label = "Plot Options",
                                                  inputId = ns("optionbutton")
                 ))
          )
        )
    )
  )
}

#' Differential analysis module server code
#'
#' @param input shiny input
#' @param output shiny output
#' @param session shiny session
#' @param aggDat aggregated MRExperiment
#' @param featLevel chosen feature level (aggregation level)
#' @param diffSettings reactive storing values selected in analysis 
#' input interface
#' @param reset boolean reactive which resets the module if TRUE
#' @param normalized boolean reactive indicating if data has been normalized
#' 
#' @author Janina Reeder
#' 
#' @import metagenomeSeq
#' 
#' @return list containing R code for analysis and for feature plots
#'
#' @export
diffTable <- function(input, output, session, 
                      aggDat, featLevel, diffSettings, 
                      reset, normalized) {
  ns <- session$ns
  
  ## reactive storing analysis code for reports
  repCode <- reactiveVal(NULL)
  ## reactive storing feature plot code for reports
  featListRep <- reactiveVal(NULL)
  normalizeCode <- reactiveVal(NULL)
  diffResults <- reactiveVal(NULL)
  plotHeight <- reactiveVal(350)
  
  observe({
    req(reset())
    shinyjs::disable("optionbutton_state")
    diffResults(NULL)
  })
  
  diffAggData <- reactive({
    req(aggDat(), featLevel())
    # FORCE NORMALIZATION IF NOT PERFORMED YET
    if (isFALSE(normalized())) {
      diffAggDat <- aggDat()
      ## need trycatch for only one feature level remaining
      p <- tryCatch({
        cumNormStatFast(diffAggDat)
      },
      error = function(e){
        return(0.5)
      })
      if(is.na(sum(normFactors(aggDat()))))
        diffAggDat <- cumNorm(diffAggDat, p)
      normalizeCode(paste0(c("diffAggDat <- aggDat",
                             "p <- tryCatch({",
                             "\tcumNormStatFast(diffAggDat)",
                             "},",
                             "error = function(e){",
                             "\treturn(0.5)",
                             "})",
                             "if(is.na(sum(normFactors(diffAggDat))))",
                             "\tdiffAggDat <- cumNorm(diffAggDat, p)",
                             "diffResults <- runDiffTest(diffAggDat,"),
                           collapse = "\n"))
      diffAggDat
    } else {
      normalizeCode("diffResults <- runDiffTest(aggDat,")
      aggDat()
    }
  })
  
  # performs the differential analysis based on given data, level and settings
  observe({
    req(diffAggData(), featLevel(), 
        diffSettings()$comparison, 
        diffSettings()$method,
        diffSettings()$phenolevel1,
        diffSettings()$phenolevel2)
    
    if(diffSettings()$method == "DESeq2"){
      if(make.names(diffSettings()$phenolevel1) ==
         make.names(diffSettings()$phenolevel2)){
        showModal(modalDialog(
          title = "DeSeq2 Issue",
          "DeSeq2 uses make.names() on the comparison phenotypes.\n 
                The chosen values lead to non-unique values. \n
                Please update the phenotable for this project in order
                to run DeSeq2 differential analysis using the selected comparison levels.",
          easyClose = TRUE
        ))
        return(NULL)
      }
    }
    
    tableCaption <- paste0(diffSettings()$method,
                           " comparison of ",
                           diffSettings()$comparison,
                           ": ",
                           diffSettings()$phenolevel1,
                           " vs ",
                           diffSettings()$phenolevel2)
    
    repCode(paste(
      paste0("#' ### Differential analysis"),
      paste0("#' #### ", diffSettings()$method,": ",
             diffSettings()$phenolevel1, " vs ",
             diffSettings()$phenolevel2),
      paste0("", "#- fig.width = 8"),
      paste0(normalizeCode()),
      paste0("\tlevel = \"", featLevel(), "\","),
      paste0("\tphenotype = \"", 
             diffSettings()$comparison, "\","),
      paste0(
        "\tphenolevels = c(\"",
        diffSettings()$phenolevel1,
        "\", \"",
        diffSettings()$phenolevel2, "\"),"
      ),
      paste0("\tmethod = \"", 
             diffSettings()$method, "\")"),
      paste0("\n"),
      paste0("if(doctype == \"html\"){"),
      paste0("\tDT::datatable(data = diffResults, 
                   class = \"stripe hover cell-border order-column\","),
      paste0("\t\tfilter = \"none\", style = \"bootstrap\","),
      paste0("\t\tcaption = \"",tableCaption,"\","),
      paste0("\t\toptions = list(scrollX = TRUE,paging = TRUE, 
                   digits = 4,dom = \"<Bftp>\"),"),
      paste0("\t\trownames = FALSE)"),
      paste0("} else {"),
      paste0("\tprint(\"data too large to display in pdf\")"),
      paste0("}"),
      sep = "\n"
    ))
    
    withProgress({
      diffResults(runDiffTest(
        aggdat = diffAggData(),
        level = featLevel(),
        phenotype = diffSettings()$comparison,
        phenolevels = c(diffSettings()$phenolevel1, 
                        diffSettings()$phenolevel2),
        method = diffSettings()$method
      ))
    }, message = "Calculating ...")
  })
  
  ## render the table showing the results; this needs to be client side, 
  ## so sorting doesn't remove selections
  output$diffdatatable <- DT::renderDT({
    req(diffResults())
    shinyjs::js$moveButton(ns("downloaddiv"),ns("download_button"))
    tableCaption <- isolate(paste0(diffSettings()$method,
                                   " comparison of ",
                                   diffSettings()$comparison,
                                   ": ",
                                   diffSettings()$phenolevel1,
                                   " vs ",
                                   diffSettings()$phenolevel2))
    shinyjs::show("download_button")
    DT::datatable(
      data = diffResults(), 
      class = "stripe hover cell-border order-column",
      filter = "none", style = "bootstrap",
      callback = DT::JS("$('div.dwnld_diff').append($('#diffAnalysis-differentialTable-download_button'));"),
      caption = tableCaption,
      extensions = c("ColReorder", "Buttons", "Select"),
      options = list(
        scrollX = TRUE,
        paging = TRUE,
        colReorder = TRUE,
        stateSave = TRUE,
        stateLoadParams = DT::JS("function (settings, data) 
                                         {return false;}"),
        digits = 4,
        columnDefs = list(list(className = "select-checkbox", 
                               targets = 0)),
        select = list(style = "multi", selector = "tr"),
        buttons = list(
          list(
            extend = "colvis",
            text = "Select columns"
          )),
        dom = '<"dwnld_diff"B>ftlp'
      ),
      rownames = "",
      selection = "none")
  }, server = FALSE)
  
  output$download_button <- downloadHandler(
    filename = function() {
      gsub(" ","",paste0(diffSettings()$phenolevel1,
            "_vs_",
            diffSettings()$phenolevel2, "_", Sys.Date(), ".csv"))
    },
    content = function(file) {
      readr::write_csv(diffResults(), file)
    }
  )
  
  ## add/remove points in feature plots
  observeEvent(input$sP, {
    req(input$diffdatatable_rows_selected, diffResults(), aggDat(), 
        featLevel(), diffSettings())
    featListRep(sapply(featListRep(), 
                       function(f){
                         gsub("showPoints = [A-Z].+",
                              paste0("showPoints = ", input$sP, ")"),f)
                         })
                
                )
    showPoints <- 'if'(input$sP, "all", FALSE)
    plotly::plotlyProxy("clickedFeature", session) %>%
      plotly::plotlyProxyInvoke(
        "restyle",
        list(boxpoints = showPoints)
      )
  })


  observeEvent(plotHeight(), {
    req(input$diffdatatable_rows_selected, diffResults(), aggDat(), 
        featLevel(), diffSettings())
    plotly::plotlyProxy("clickedFeature", session) %>%
      plotly::plotlyProxyInvoke(
        "relayout",
        list(height = plotHeight())
      )
  },ignoreInit = TRUE)


  ## renders a plot if datatable rows is clicked and stores results for report
  output$clickedFeature <- plotly::renderPlotly({
    req(input$diffdatatable_rows_selected, diffResults(), 
        aggDat(), featLevel(), diffSettings())
    isolate(featListRep(NULL))
    shinyjs::enable("optionbutton_state")
    ## determine how many rows the plots should be distributed upon
    numofrows <- ceiling(length(input$diffdatatable_rows_selected) / 2)
    
    plotList <- lapply(input$diffdatatable_rows_selected, function(i) {
      isolate({
        featListRep(c(
          featListRep(),
          paste(
            paste0("#' ### Feature plot"),
            paste0("", "#- fig.width = 7"),
            paste0("plotSingleFeature(aggDat,"),
            paste0("\tfeature = \"", 
                   diffResults()[i, featLevel()], "\","),
            paste0("\tx_var = \"", 
                   diffSettings()$comparison, "\","),
            paste0("\tylab = \"Reads\","), 
            paste0("\tlog = ", input$logS, ","),
            paste0("\tshowPoints = ", input$sP, ")\n\n"),
            sep = "\n"
          )
        ))
      })
      
      plotSingleFeature(
        aggdat = aggDat(),
        feature = diffResults()[i, featLevel()],
        x_var = diffSettings()$comparison,
        ylab = "Reads",
        xlab = diffResults()[i, featLevel()],
        log = input$logS,
        showPoints = isolate(input$sP),
        fixedHeight = numofrows * 350,
        x_levels = c(diffSettings()$phenolevel1, 
                     diffSettings()$phenolevel2)
      )
    })

    plotly::subplot(plotList, nrows = numofrows, 
                    margin = 0.1, titleY = TRUE, titleX = TRUE) %>%
      plotly::layout(showlegend = FALSE)
  })
  
  return(list(reptable = repCode, repplot = featListRep))
}
