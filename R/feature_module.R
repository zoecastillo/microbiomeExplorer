

#' Helper function to replace any un-annotated features with the term unknown
#'
#' @param featcol vector of entried to be replaced where needed (fData column)
#'
#' @author Janina Reeder
#'
#' @return modified featcol
#' 
#' @examples 
#' data("mouseData", package = "metagenomeSeq")
#' featcol <- fData(mouseData)[["genus"]]
#' featcol[featcol == "NA"] <- NA
#' replaceWithUnknown(featcol)
#' 
#' @export
replaceWithUnknown <- function(featcol) {
    featcol <- stringr::str_replace_na(featcol, replacement = "unknown")
    featcol <- stringr::str_replace_all(featcol, "^[kpcofgs]{1}__$", "unknown")
    featcol[stringr::str_length(featcol) == 0] <- "unkown"
    return(featcol)
}

#' Helper function which rolls down annotated from closest higher order with 
#' annotation
#'
#' @param featrow vector of entries to be replaced where needed (fData row)
#'
#' @author Janina Reeder
#'
#' @return modified featurerow
#' 
#' @examples 
#' data("mouseData", package = "metagenomeSeq")
#' featrow <- fData(mouseData)[5,]
#' rollDownFeatures(featrow)
#' 
#' @export
rollDownFeatures <- function(featrow) {
    featrow <- stringr::str_replace_na(t(featrow), replacement = "NA")
    featrow[stringr::str_length(featrow) == 0] <- "NA"
    rdlocs <- stringr::str_detect(featrow, "^[kpcofgs]{1}__$|^NA$")
    if (sum(rdlocs) > 0) {
        ri <- min(which(rdlocs))
        featrow[rdlocs] <- paste0("unknown_", featrow[ri - 1])
    }
    return(featrow)
}

splitColumns <- function(featData){
    if(is.null(featData)){
        return(NULL)
    }
    if(length(featData) == 0){
        return(NULL)
    }

    featList <- lapply(featData, function(f){
        stringr::str_trim(
            stringr::str_split(f,";")[[1]]
        )
    })
    numofcols <- max(lengths(featList))
    featList %>% purrr::map_dfr(~ as.data.frame(t(.)))
}


#' Feature table UI module
#'
#' @param id namespace identifier
#'
#' @author Janina Reeder
#'
#' @return fluidRow containing the UI code for feature tables
#' 
#' @examples featureTableUI("feature_id")
#' 
#' @export
featureTableUI <- function(id) {
    ns <- NS(id)
    
    fluidRow(
        column(
            width = 3, id = ns("featannocol"),
            shinyjs::hidden(
                div(
                    id = ns("splitdiv"),
                    h4("SPLIT COLUMNS"),
                    selectInput(
                        ns("splittaxonomy"),
                        label = "Choose taxonomy column",
                        choices = "", 
                        multiple = FALSE, selectize = FALSE, width = "250px"),
                    actionButton(
                        ns("splitbutton"), 
                        icon = icon("fas fa-angle-double-right"),
                        label = HTML("&nbspSPLIT"), 
                        width = "100px")
                )
            ),
            h4("ANNOTATE BLANK VALUES"),
            selectInput(
                ns("featureanno"),
                label = "Method",
                choices = c("Roll down taxonomy", "Mark as unknown"), 
                multiple = FALSE, selectize = FALSE, width = "250px"),
            div(
                id = ns("buttondiv"),
                fluidRow(
                    width = 12, id = "actionbuttonrow",
                    actionButton(
                        ns("annobutton"), 
                        icon = icon("fas fa-angle-double-right"),
                        label = HTML("&nbsp;ASSIGN"), 
                        width = "100px"),
                    actionButton(
                        ns("resetbutton"),
                        icon = icon("fas fa-redo-alt"),
                        label = HTML("&nbsp;RESET"), 
                        width = "100px")
                ),
                fluidRow(
                    width = 12, id = "actionbuttonrow2",
                    shinyjs::disabled(
                        actionButton(
                            ns("savebutton"), 
                            icon = icon("far fa-save"),
                            label = HTML("&nbsp;SAVE"),
                            width = "100px"))
                )
            )
        ),
        column(
            width = 9,
            box(
                width = 10,
                h2("FEATURE OVERVIEW"),
                p("Available feature taxonomy for the counts data.  
                Table settings allow paging through sections of the data, 
                choosing how many entries to display
                 or searching for specific entries.
                Unannotated features can be marked as \"Unknown\" or obtain
                annotation via the next available higher taxonomy level in a 
                roll down mechanism. Modifications must be saved in order
                  to be available in the analysis sections.")
            ),
            box(
                width = 12,
                div(
                    id = ns("downloaddiv"),
                    downloadButton(ns("download_button"),"Download"),
                    DT::DTOutput(ns("featuredatatable"), width = "100%")
                )
            )
        )
    )
}

#' Helper function returning the fData modifications as strings for 
#' report generation
#'
#' @param featureanno type of feature annotation; values are "Mark unknown" or 
#' "Roll down"
#'
#' @return String containing R code performing the modification
getFeatModCode <- function(featureanno) {
    if (featureanno == "Mark as unknown") {
        paste(paste0("bufrownames <- row.names(fData(meData))"),
              paste0("df <- as.data.frame(apply(fData(meData),2, 
                   replaceWithUnknown))"),
              paste0("rownames(df) <- bufrownames"),
              paste0("meData <- addFeatData(meData,df)"),
              sep = "\n"
        )
    } else if (featureanno == "Roll down taxonomy") {
        paste(paste0("bufcolnames <- names(fData(meData))"),
              paste0("df <- as.data.frame(t(apply(fData(meData),1, 
                   rollDownFeatures)))"),
              paste0("names(df) <- bufcolnames"),
              paste0("meData <- addFeatData(meData,df)"),
              sep = "\n"
        )
    }
}

#' Helper function returning the fData modifications as strings for 
#' report generation
#'
#' @param splittaxonomy name of column to split on
#'
#' @return String containing R code performing the modification
getFeatSplitCode <- function(splittaxonomy) {
    paste(paste0("bufrownames <- row.names(fData(meData))"),
          paste0("df <- splitColumns(fData(meData)[[\"",
                 splittaxonomy,"\"]])"),
          paste0("rownames(df) <- bufrownames"),
          paste0("meData <- addFeatData(meData,df)"),
          sep = "\n"
    )
}

#' Feature table module server code
#'
#' @param input shiny input
#' @param output shiny output
#' @param session shiny session
#' @param meData MRExperiment storing the data
#' @param featureModRep reactiveValue storing modifications performed on fData
#'
#' @return feature table server fragment - no return value
#'
#' @author Janina Reeder
#' 
#' @importFrom Biobase fData
featureTable <- function(input, output, session, meData, featureModRep) {
    ns <- session$ns
    
    ## stores the fData of the given MRExperiment
    featFrame <- reactiveVal(NULL)
    ## keeps track of whether annotation was performed
    annotated <- reactiveVal(FALSE)
    taxsplit <- reactiveVal(FALSE)
    
    ## initialize featFrame when meData becomes available
    observe({
        req(meData())
        isolate(featFrame(fData(meData())))
    })
    
    observe({
        req(featFrame())
        if(ncol(featFrame()) <= 2){
            shinyjs::disable("annobutton")
        } else {
            shinyjs::enable("annobutton")
        }
        if(nrow(featFrame()) > 0){
            semicolon_index <- grepl(";",featFrame())
            if(any(semicolon_index)){
                shinyjs::show("splitdiv")
                updateSelectInput(session,"splittaxonomy",
                                  choices = names(featFrame())[semicolon_index])
            } else {
                shinyjs::hide("splitdiv")
            }
        } else {
            shinyjs::hide("splitdiv")
        }
    })
    
    observeEvent(input$splitbutton, {
        req(input$splittaxonomy %in% names(featFrame()))
        bufrownames <- row.names(featFrame())
        shinyjs::addClass("featuredatatable", "transparent")
        df <- splitColumns(featFrame()[[input$splittaxonomy]])
        rownames(df) <- bufrownames
        shinyjs::removeClass("featuredatatable", "transparent")
        ## update stored feature frame
        featFrame(df)
        taxsplit(TRUE)
        shinyjs::enable("savebutton")
    })
    
    ## perform annotation op
    observeEvent(input$annobutton, {
        if (input$featureanno == "Mark as unknown") {
            ## everything will just be unknown
            bufrownames <- row.names(featFrame())
            df <- as.data.frame(apply(featFrame(), 2, replaceWithUnknown))
            rownames(df) <- bufrownames
            ## update stored feature frame
            featFrame(df)
            ## mark as annotated
            annotated(TRUE)
        } else if (input$featureanno == "Roll down taxonomy") {
            ## use annotation of parent
            bufcolnames <- names(featFrame())
            ## Mark datatable as tranparent to show the user the roll down is 
            ## in progress (css switch)
            shinyjs::addClass("featuredatatable", "transparent")
            ## roll down features row by row; used stringr, so fairly fast. 
            df <- as.data.frame(t(apply(featFrame(), 1, rollDownFeatures)))
            ## Return to normal
            shinyjs::removeClass("featuredatatable", "transparent")
            names(df) <- bufcolnames
            ## update stored feature frame
            featFrame(df)
            ## mark as annotated
            annotated(TRUE)
        }
        shinyjs::enable("savebutton")
    })
    
    ## revert annotation changes
    observeEvent(input$resetbutton, {
        req(meData())
        if (annotated() || taxsplit()) {
            ## go back to original dataset
            featFrame(fData(meData()))
            annotated(FALSE)
            taxsplit(FALSE)
            featureModRep(NULL)
            shinyjs::disable("savebutton")
        }
    })
    
    ## make changes permanent
    observeEvent(input$savebutton, {
        req(meData())
        if (annotated() || taxsplit()) {
            ## adjust original dataset
            meData(addFeatData(meData(), featFrame()))
            if(taxsplit()){
                featureModRep(getFeatSplitCode(input$splittaxonomy))
                taxsplit(FALSE)
            }
            if(annotated()){
                featureModRep(paste(featureModRep(),
                                    getFeatModCode(input$featureanno),
                                    sep = "\n"))
                annotated(FALSE)
            }
            shinyjs::disable("savebutton")
        }
    })
    
    ## render table showing feature data
    output$featuredatatable <- DT::renderDT({
        req(meData())
        shinyjs::js$moveButton(ns("downloaddiv"),ns("download_button"))
        DT::datatable(
            data = featFrame(), class = "stripe hover cell-border order-column",
            filter = "none", style = "bootstrap",
            callback = DT::JS("$('div.dwnld_feat').append($('#featuretab-download_button'));"),
            extensions = c("FixedColumns", "ColReorder", "Buttons"),
            options = list(
                scrollX = TRUE,
                paging = TRUE,
                colReorder = TRUE,
                stateSave = TRUE,
                stateLoadParams = DT::JS("function (settings, data) {
                                         return false;}"),
                buttons = list(
                    list(
                        extend = "colvis",
                        text = "Select columns"
                    )),
                dom = '<"dwnld_feat"B>ftlp'
            )
        )
    })
    
    output$download_button <- downloadHandler(
        filename = function() {
            paste0("featuredata-", Sys.Date(), ".csv")
        },
        content = function(file) {
            readr::write_csv(featFrame(), file)
        }
    )
}
