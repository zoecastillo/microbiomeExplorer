## Module handling correlation analysis in the Shiny application

#' Feature correlation analysis module UI
#'
#' @param id namespace identifier
#' @author Janina Reeder
#'
#' @return box containing the UI elements
featureCorrUI <- function(id) {
    ns <- NS(id)
    
    box(title = "FEATURE CORRELATION", solidHeader = TRUE, 
        collapsible = TRUE, width = 12,
        fluidRow(
            column(width = 11,
                   shinycssloaders::withSpinner(
                       plotly::plotlyOutput(
                           ns("featcorrplot"), 
                           width = "auto", 
                           height = "auto"),
                       type = 3, color = "#424242", 
                       color.background = "#fdfdfc"),
                   br(),
                   shinyjs::disabled(
                       shinyWidgets::dropdownButton(
                           tags$h3("Plot Options"),
                           selectInput(
                               inputId = ns("fcol"),
                               label = "Color by",
                               choices = ""
                           ),
                           shinyWidgets::switchInput(
                               inputId = ns("logS"),
                               label = "Log Scale",
                               value = TRUE,
                               size = "mini",
                               labelWidth = "120px"
                           ),
                           shinyWidgets::switchInput(
                               inputId = ns("regL"),
                               label = "Regression Line",
                               value = TRUE,
                               size = "mini",
                               labelWidth = "120px"
                           ),
                           sliderInput(
                               inputId = ns("plotWidth"),
                               label = "Adjust plot width",
                               value = 550,
                               min = 250,
                               max = 1600,
                               round = TRUE
                           ),
                           circle = FALSE, 
                           status = "danger", 
                           icon = icon("gear"), 
                           width = "300px",
                           label = "Plot Options",
                           inputId = ns("optionbutton")
                       ))
            )
        ),
        fluidRow(
            column(width = 1),
            column(width = 10, class = "statsrow",
                   DT::DTOutput(ns("statsdatatable"), 
                                width = "90%", height = "auto")
            )
        )
    )
}

#' Feature correlation analysis server module
#'
#' @param input module input
#' @param output module output
#' @param session app session
#' @param aggDat aggregated MRExperiment
#' @param colorOptions reactive storing filters available via data input
#' @param corFeatBase first correlation feature
#' @param corFeat2 second correlation feature
#' @param corFacet1 first correlation facet
#' @param corFacet2 second correlation facet
#' @param corMethod correlation method to use
#' @param reset boolean reactive which resets the module if TRUE
#' 
#' @author Janina Reeder
#' @return R code used to do the correlation analysis (character)
featureCorr <- function(input, output, session,
                        aggDat,
                        colorOptions,
                        corFeatBase,
                        corFeat2,
                        corFacet1,
                        corFacet2,
                        corMethod,
                        reset) {
    ns <- session$ns
    
    observe({
        updateSelectInput(session, "fcol", 
                          choices = c("No Color", colorOptions()))
    })
    
    ## stores color selection
    colorChoice <- reactiveVal("No Color")
    ## stores featcorrplot
    corPlot <- reactiveVal(NULL)
    ## stores statsdatatable
    corStats <- reactiveVal(NULL)
    ## stores R code for reports
    repCode <- reactiveVal(NULL)
    width <- reactiveVal(550)
    
    observe({
        req(reset())
        shinyjs::disable("optionbutton_state")
        colorChoice("No Color")
        corPlot(NULL)
        corStats(NULL)
        repCode(NULL)
    })
    
    
    observeEvent(input$fcol, {
        colorChoice(input$fcol)
    })
    
    observe({
        width(input$plotWidth)
    })
    
    observeEvent(width(), {
        req(corPlot())
        plotly::plotlyProxy("featcorrplot", session) %>%
            plotly::plotlyProxyInvoke(
                "relayout",
                list(width = width())
            )
    })
    
    
    observe({
        req(aggDat(), corFeatBase(), corFeat2(),colorChoice())
        color <- as.character(colorChoice())
        if (stringr::str_detect(color, "No Color")) {
            color <- NULL
        }
        plotTitle <- paste0(stringr::str_to_sentence(corMethod()),
                            " correlation of ",
                            corFeatBase(), " vs ",
                            corFeat2())
        if(!is.null(corFacet1())){
            plotTitle <- paste0(plotTitle,
                                " split by ",
                                corFacet1())
        } else if(!is.null(corFacet2())){
            plotTitle <- paste0(plotTitle,
                                " split by ",
                                corFacet2())
        }
        if(!is.null(corFacet1()) & !is.null(corFacet2())){
            plotTitle <- paste0(plotTitle,
                                " and ",
                                corFacet2())
        }
        ## main function to perform feature-feature correlation
        withProgress({
            cf <- corrFeature(aggDat(),
                              feat1 = corFeatBase(),
                              feat2 = corFeat2(),
                              log = input$logS,
                              facet1 = corFacet1(),
                              facet2 = corFacet2(),
                              method = corMethod(),
                              addRegression = input$regL,
                              plotTitle = plotTitle,
                              col_by = color
            )
            corPlot(cf$plot)
            setProgress(0.9)
            
            statsframe <- dplyr::bind_rows(
                lapply(
                    names(cf$stats), 
                    function(n) {
                        s <- cf$stats[[n]]
                        list(
                            "facet" = n,
                            "method" = corMethod(),
                            "estimate" = round(s$estimate, 
                                               getOption("me.round_digits")),
                            "p" = round(s$p.value, 
                                        getOption("me.round_digits")),
                            "lower CI" = round(s$lower, 
                                               getOption("me.round_digits")),
                            "upper CI" = round(s$upper, 
                                               getOption("me.round_digits"))
                        )
                    }))
        }, message = "Analyzing correlation ...")
        
        ## code needed to do the analysis
        repCode(paste(
            paste0("#' ### Feature correlation"),
            paste0("", "#- fig.width = 10"),
            paste0("cf <- corrFeature(aggDat,"),
            paste0("\tfeat1 = \"", corFeatBase(), "\","),
            paste0("\tfeat2 = \"", corFeat2(), "\","),
            paste0("\tlog = ", input$logS, ","),
            paste0(
                "\tfacet1 = ", 'if'(is.null(corFacet1()),
                                    "NULL",
                                    paste0("\"", corFacet1(), "\"")
                ),
                ","
            ),
            paste0(
                "\tfacet2 = ", 'if'(is.null(corFacet2()),
                                    "NULL",
                                    paste0("\"", corFacet2(), "\"")
                ),
                ","
            ),
            paste0("\tmethod = \"", corMethod(), "\","),
            paste0("\taddRegression = ", input$regL, ","),
            paste0("\tplotTitle = \"", plotTitle, "\","),
            paste0("\tcol_by = ", 'if'(is.null(color), 
                                       "NULL", 
                                       paste0("\"", color, "\"")),","),
            paste0("\tallowWebGL = FALSE)\n\n"),
            paste0("cf$plot"),
            paste0("\n\n"),
            paste0("statsframe <- dplyr::bind_rows(lapply(names(cf$stats), 
                   function(n){"),
            paste0("\ts <- cf$stats[[n]]"),
            paste0("\tlist(\"facet\" = n,"),
            paste0("\t\t\"method\" = \"",corMethod(),"\","),
            paste0("\t\t\"estimate\" = round(s$estimate, getOption(\"me.round_digits\")),"),
            paste0("\t\t\"p\" = round(s$p.value, getOption(\"me.round_digits\")),"),
            paste0("\t\t\"lower CI\" = round(s$lower, getOption(\"me.round_digits\")),"),
            paste0("\t\t\"upper CI\" = round(s$upper, getOption(\"me.round_digits\")))"),
            paste0("}))\n\n"),
            paste0("if(doctype == \"html\"){"),
            paste0("\tDT::datatable(data = statsframe, 
                   class = \"stripe hover cell-border order-column\","),
            paste0("\t\tfilter = \"none\", style = \"bootstrap\","),
            paste0("\t\toptions = list(scrollX = TRUE, paging = TRUE, 
                   digits = 4,dom = \"<Bftp>\"),"),
            paste0("\t\trownames = FALSE, escape = FALSE)"),
            paste0("} else {"),
            paste0("\tkable(statsframe)"),
            paste0("}\n\n"),
            sep = "\n"
        ))
        corStats(statsframe)
    })
    
    
    output$featcorrplot <- plotly::renderPlotly({
        req(corPlot())
        shinyjs::enable("optionbutton_state")
        corPlot()
    })
    
    output$statsdatatable <- DT::renderDT({
        req(corStats())
        
        DT::datatable(
            data = corStats(), class = "stripe hover cell-border order-column",
            filter = "none", style = "bootstrap",
            options = list(
                scrollX = TRUE,
                paging = FALSE,
                digits = 4,
                dom = "<tp>"
            ),
            rownames = NULL,
            escape = FALSE
        )
    })
    
    return(repCode)
}



#' Phenotype correlation analysis module
#'
#' @param id namespace identifier
#' @author Janina Reeder
#'
#' @return box containing the UI element
phenotypeCorrUI <- function(id) {
    ns <- NS(id)
    
    box(title = "PHENOTYPE CORRELATION", solidHeader = TRUE, 
        collapsible = TRUE, width = 12,
        fluidRow(
            column(
                width = 11,
                shinycssloaders::withSpinner(
                    plotly::plotlyOutput(ns("phenocorrplot"), 
                                         width = "auto", 
                                         height = "auto"),
                    type = 3, color = "#424242", color.background = "#fdfdfc"),
                br(),
                shinyjs::disabled(
                    shinyWidgets::dropdownButton(
                        tags$h3("Plot Options"),
                        selectInput(
                            inputId = ns("pcol"),
                            label = "Color by",
                            choices = ""
                        ),
                        shinyWidgets::switchInput(
                            inputId = ns("logS"),
                            label = "Log Scale",
                            value = TRUE,
                            size = "mini",
                            labelWidth = "120px"
                        ),
                        shinyWidgets::switchInput(
                            inputId = ns("regL"),
                            label = "Regression line",
                            value = TRUE,
                            size = "mini",
                            labelWidth = "120px"
                        ),
                        sliderInput(
                            inputId = ns("plotWidth"),
                            label = "Adjust plot width",
                            value = 550,
                            min = 250,
                            max = 1600,
                            round = TRUE
                        ),
                        circle = FALSE, status = "danger", 
                        icon = icon("gear"), width = "300px",
                        label = "Plot Options",
                        inputId = ns("optionbutton")
                    ))
            )
        ),
        fluidRow(
            column(width = 1),
            column(width = 10, class = "statsrow",
                   column(
                       width = 10, class = "statsrow",
                       DT::DTOutput(ns("statsdatatable"), 
                                    width = "90%", height = "auto")
                   )
            )
        )
    )
}

#' Phenotype correlation analysis server module
#'
#' @param input shiny input
#' @param output shiny output
#' @param session shiny session
#' @param aggDat aggregated MRExperiment
#' @param colorOptions reactive storing filters available via data input
#' @param corFeatBase first correlation feature
#' @param corPheno correlation phenotype
#' @param corFacet1 first correlation facet
#' @param corFacet2 second correlation facet
#' @param corMethod correlation method to use
#' @param reset boolean reactive which resets the module if TRUE
#' 
#' @author Janina Reeder
#' @return R code used to do the correlation analysis (character)
phenotypeCorr <- function(input, output, session,
                          aggDat,
                          colorOptions,
                          corFeatBase,
                          corPheno,
                          corFacet1,
                          corFacet2,
                          corMethod,
                          reset) {
    ns <- session$ns
    
    observe({
        updateSelectInput(session, "pcol", 
                          choices = c("No Color", colorOptions()))
    })
    
    ## stores color selection
    colorChoice <- reactiveVal("No Color")
    ## stores featcorrplot
    corPlot <- reactiveVal(NULL)
    ## stores statsdatatable
    corStats <- reactiveVal(NULL)
    ## stores R code for reports
    repCode <- reactiveVal(NULL)
    width <- reactiveVal(550)
    
    observe({
        req(reset())
        shinyjs::disable("optionbutton_state")
        colorChoice("No Color")
        corPlot(NULL)
        corStats(NULL)
        repCode(NULL)
    })
    
    observeEvent(input$pcol, {
        colorChoice(input$pcol)
    })
    
    observe({
        width(input$plotWidth)
    })
    
    observeEvent(width(), {
        req(corPlot())
        plotly::plotlyProxy("phenocorrplot", session) %>%
            plotly::plotlyProxyInvoke(
                "relayout",
                list(width = width())
            )
    })
    
    
    observe({
        req(aggDat(), corFeatBase(), corPheno(), colorChoice())
        color <- as.character(colorChoice())
        if (stringr::str_detect(color, "No Color")) {
            color <- NULL
        }
        plotTitle <- paste0(stringr::str_to_sentence(corMethod()),
                            " correlation of ",
                            corFeatBase(), " vs ",
                            corPheno())
        if(!is.null(corFacet1())){
            plotTitle <- paste0(plotTitle,
                                " split by ",
                                corFacet1())
        } else if(!is.null(corFacet2())){
            plotTitle <- paste0(plotTitle,
                                " split by ",
                                corFacet2())
        }
        if(!is.null(corFacet1()) & !is.null(corFacet2())){
            plotTitle <- paste0(plotTitle,
                                " and ",
                                corFacet2())
        }
        ## main function to perform feature-phenotype correlation
        withProgress({
            cp <- corrPhenotype(aggDat(),
                                feature = corFeatBase(),
                                phenotype = corPheno(),
                                log = input$logS,
                                facet1 = corFacet1(),
                                facet2 = corFacet2(),
                                method = corMethod(),
                                addRegression = input$regL,
                                plotTitle = plotTitle,
                                col_by = color
            )
            corPlot(cp$plot)
            setProgress(0.9)
            
            statsframe <- dplyr::bind_rows(
                lapply(
                    names(cp$stats), 
                    function(n) {
                        s <- cp$stats[[n]]
                        list(
                            "facet" = n,
                            "method" = corMethod(),
                            "estimate" = round(s$estimate, 
                                               getOption("me.round_digits")),
                            "p" = round(s$p.value, 
                                        getOption("me.round_digits")),
                            "lower CI" = round(s$lower, 
                                               getOption("me.round_digits")),
                            "upper CI" = round(s$upper, 
                                               getOption("me.round_digits"))
                        )
                    }))
            corStats(statsframe)
        }, message = "Analyzing correlation ...")
        
        ## code needed to do the analysis
        repCode(paste(
            paste0("#' ### Phenotype correlation"),
            paste0("", "#- fig.width = 10"),
            paste0("cp <- corrPhenotype(aggDat,"),
            paste0("\tfeature = \"", corFeatBase(), "\","),
            paste0("\tphenotype = \"", corPheno(), "\","),
            paste0("\tlog = ", input$logS, ","),
            paste0(
                "\tfacet1 = ", 'if'(is.null(corFacet1()),
                                    "NULL",
                                    paste0("\"", corFacet1(), "\"")
                ),
                ","
            ),
            paste0(
                "\tfacet2 = ", 'if'(is.null(corFacet2()),
                                    "NULL",
                                    paste0("\"", corFacet2(), "\"")
                ),
                ","
            ),
            paste0("\tmethod = \"", corMethod(), "\","),
            paste0("\taddRegression = ", input$regL, ","),
            paste0("\tplotTitle = \"", plotTitle, "\","),
            paste0("\tcol_by = ", 'if'(is.null(color), "NULL", 
                                       paste0("\"", color, "\"")),","),
            paste0("\tallowWebGL = FALSE)\n\n"),
            paste0("cp$plot"),
            paste0("\n\n"),
            paste0("statsframe <- dplyr::bind_rows(lapply(names(cp$stats), 
                   function(n){"),
            paste0("\ts <- cp$stats[[n]]"),
            paste0("\tlist(\"facet\" = n,"),
            paste0("\t\t\"method\" = \"",corMethod(),"\","),
            paste0("\t\t\"estimate\" = round(s$estimate, getOption(\"me.round_digits\")),"),
            paste0("\t\t\"p\" = round(s$p.value, getOption(\"me.round_digits\")),"),
            paste0("\t\t\"lower CI\" = round(s$lower, getOption(\"me.round_digits\")),"),
            paste0("\t\t\"upper CI\" = round(s$upper, getOption(\"me.round_digits\")))"),
            paste0("}))\n\n"),
            paste0("if(doctype == \"html\"){"),
            paste0("\tDT::datatable(data = statsframe, 
                   class = \"stripe hover cell-border order-column\","),
            paste0("\t\tfilter = \"none\", style = \"bootstrap\","),
            paste0("\t\toptions = list(scrollX = TRUE,paging = TRUE, 
                   digits = 4,dom = \"<Bftp>\"),"),
            paste0("\t\trownames = FALSE, escape = FALSE)"),
            paste0("} else {"),
            paste0("\tkable(statsframe, escape = FALSE)"),
            paste0("}"),
            sep = "\n"
        ))
    })
    
    output$phenocorrplot <- plotly::renderPlotly({
        req(corPlot())
        shinyjs::enable("optionbutton_state")
        corPlot()
    })
    
    output$statsdatatable <- DT::renderDT({
        req(corStats())
        
        DT::datatable(
            data = corStats(), class = "stripe hover cell-border order-column",
            filter = "none", style = "bootstrap",
            options = list(
                scrollX = TRUE,
                paging = FALSE,
                digits = 4,
                dom = "<tp>"
            ),
            rownames = NULL,
            escape = FALSE
        )
    })
    
    return(repCode)
}
