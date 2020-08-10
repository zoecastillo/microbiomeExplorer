
#' Report row module consisting of a checkbox, image and description/R code area
#'
#' @param id namespace identifier
#' @param type boolean indicating if a selector checkbox should be added
#' 
#' @author Janina Reeder
#' @return div holding the UI code
reportRowUI <- function(id, type) {
    ns <- NS(id)
    
    div(
        fluidRow(
            width = 11, id = "reportrow",
            column(
                width = 1,
                if (!type) {
                    shinyWidgets::prettyCheckbox(
                        inputId = ns("includeCheck"),
                        label = "",
                        value = TRUE,
                        icon = icon("check"),
                        status = "success",
                        animation = "rotate",
                        outline = TRUE
                    )
                }
            ),
            column(
                width = 4,
                plotOutput(
                    ns("reportImage"), height = "175px", width = "175px")
            ),
            column(width = 1),
            column(
                width = 6,
                uiOutput(ns("reportText"), height = "200px")
            )
        ),
        hr()
    )
}


#' Report Row
#'
#' @param input module input
#' @param output module output
#' @param session app session
#' @param type boolean indicating whether checkbox should be included
#' @param content R code to show
#' 
#' @author Janina Reeder
#' @return reactive boolean indicating whether row is selected
reportRow <- function(input, output, session, type, content) {
    ns <- session$ns
    
    ## preprocess code to only display relevant parts
    codecontent <- reactive({
        req(content())
        cc <- lapply(stringr::str_split(content(), "\n")[[1]], function(c) {
            'if'(stringr::str_starts(c, "#"), "", c)
        })
        cc[cc != ""]
    })
    
    ## determine which icon to show
    imagename <- reactive({
        req(codecontent())
        if (stringr::str_starts(codecontent()[1], 
                                "makeQCPlot")) {
            "qcpic"
        } else if (stringr::str_starts(codecontent()[1], 
                                       "plotSingleFeature")) {
            "featplot"
        } else if (stringr::str_starts(codecontent()[1], 
                                       "plotAbundance")) {
            "relAb"
        } else if (stringr::str_starts(codecontent()[1], 
                                       "plotAlpha")) {
            "alphadiv"
        } else if (stringr::str_starts(codecontent()[1], 
                                       "plotHeatmap")) {
            "heatmap"
        } else if (stringr::str_starts(codecontent()[1], 
                                       "distMat")) {
            "betadiv"
        } else if (stringr::str_starts(codecontent()[1], 
                                       "c[fp]")) {
            "correlation"
        } else if (stringr::str_starts(codecontent()[1], 
                                       "diff")) {
            "difftest"
        } else if (stringr::str_starts(codecontent()[1], 
                                       "plotLongFeature")) {
            "longitudinal"
        } else {
            "blank"
        }
    })
    
    ## render the icon
    output$reportImage <- renderImage({
        req(imagename())
        filename <- normalizePath(file.path(
            "www/",
            paste(imagename(), ".png", sep = "")
        ))
        
        # Return a list containing the filename
        list(
            src = filename,
            width = "175px",
            height = "175px"
        )
    }, deleteFile = FALSE)
    
    ## render the code sections
    output$reportText <- renderUI({
        req(codecontent())
        div(pre(id = "rcode", codecontent()),
            class = "reportDivBox"
        )
    })
    
    return(input$includeCheck)
}
