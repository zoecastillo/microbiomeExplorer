## Module handling aggregating the MRExperiment data

#' Aggregation module ui function
#'
#' @param id namespace identifier
#' @author Janina Reeder
#' 
#' @return box holding aggregation input elements
#' 
#' @export
aggregationTabUI <- function(id) {
    ns <- NS(id)

    box(
        id = ns("aggregation_box"),
        style = "padding-top: 0; padding-left: 20px;",
        width = 12,
        title = textOutput(ns("agg_title")),
        collapsible = TRUE,
        collapsed = FALSE,
        
        
        fluidRow(width = 12, 
                 shinyjs::disabled(
                     div(id = ns("normaggdiv"),
                         column(
                             width = 3,
                             selectInput(ns("normalizedata"),
                                         label = "Normalize Data",
                                         choices = c("Proportion", "CSS", "Do not normalize"),
                                         multiple = FALSE, 
                                         selectize = FALSE,
                                         width = "150px")
                         ),
                         column(
                             width = 3,
                             selectInput(ns("featurelevel"),
                                         label = "Feature level", 
                                         choices = "",
                                         multiple = FALSE, 
                                         selectize = FALSE, 
                                         width = "250px"
                             )
                         ),
                         column(
                             width = 2,
                             div(style = "margin-top: 24px;",
                                 actionButton(ns("aggregatebutton"), 
                                              icon = icon("fas fa-compress"),
                                              label = HTML("&nbsp;AGGREGATE"),
                                              width = "120px"))
                         ),
                         column(
                             width = 2,
                             div(style = "margin-top: 24px;",
                                 shinyjs::disabled(downloadButton(ns("savebutton"),
                                                                  label = "GET DATA",
                                                                  width = "120px"))
                             )
                         )
                     )
                 )
        )    
    )
}

#' Aggregation module server function
#'
#' @param input shiny input
#' @param output shiny output
#' @param session shiny session
#' @param resetInput boolean updated to TRUE if new data is available
#' @param levelOpts available levels to aggregate on (depends on input data)
#' @param chosenLevel previously selected level (passed from different instance)
#' @param meData the main MRexperiment object
#' 
#' @author Janina Reeder
#' 
#' @return reactive list holding aggregated object, aggregation code and boolean on normalization
#' 
#' @export
aggregationTab <- function(input, output, session, 
                           resetInput, levelOpts, chosenLevel, meData) {
    ns <- session$ns
    
    aggCode <- reactiveVal(NULL)
    aggMRobj <- reactiveVal(NULL)
    normalizedData <- reactiveVal(FALSE)
    normUpdate <- reactiveVal(FALSE)
    
    observeEvent(meData(),{
        if(isFALSE(normUpdate())){
            if(is.null(fData(meData()))){
                shinyjs::disable("normaggdiv")
            } else {
                shinyjs::enable("normaggdiv")
            }
            aggCode(NULL)
            aggMRobj(NULL)
            normalizedData(FALSE)
        } else {
            normUpdate(FALSE)
        }
    }, ignoreNULL = TRUE)
    
    observeEvent(levelOpts(),{
        updateSelectInput(session, "featurelevel", 
                          choices = levelOpts(), 
                          selected = chosenLevel())
    })

    
    observeEvent(input$aggregatebutton, {
        req(meData())
        aggC <- ""
        if(input$normalizedata %in% c("Proportion", "CSS")){
            normUpdate(TRUE)
            meData(normalizeData(meData(), norm_method = input$normalizedata))
            aggC <- paste0("\nmeData <- normalizeData(meData, norm_method = \"",
                                 input$normalizedata, "\")")
            normalizedData(TRUE)
        } else if(all(is.na(normFactors(meData())))){
            normalizedData(FALSE)
        } else {
            normalizedData(TRUE)
        }
        aggMRobj(aggFeatures(meData(), level = input$featurelevel))
        chosenLevel(input$featurelevel)
        aggCode(paste0(aggC,"\naggDat <- aggFeatures(meData, level = \"",
                       input$featurelevel, "\")\n"))
        shinyjs::js$collapse(ns("aggregation_box"))
    })
    
    output$agg_title <- renderText({
        if (is.null(aggMRobj())) {
            return("Aggregation")
        } else {
            return(paste0("Aggregation level: ",chosenLevel()))
        }
    })
    
    observe({
        if(!is.null(aggMRobj())){
            shinyjs::enable("savebutton")
        } else {
            shinyjs::disable("savebutton")
        }
    })
    
    ## download aggregated counts as csv
    output$savebutton <- downloadHandler(
        filename = function() {
            paste0("featureCounts.csv")
        },
        content = function(file) {
            rawcounts <- as.data.frame(MRcounts(aggMRobj()))
            rawcounts <- rawcounts %>%
                tibble::rownames_to_column(var = "Features")
            readr::write_csv(rawcounts, file)
        }
    )
 
    return(list(aggCode = aggCode,
                mrobj = aggMRobj,
                normalizedData = normalizedData))

}
