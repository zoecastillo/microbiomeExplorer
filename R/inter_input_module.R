
#' Main beta analysis input module.
#' Set up to handle all analysis tabs in the app depending on given parameters
#'
#' @param id element identifier - namespace
#'
#' @author Janina Reeder
#'
#' @return box containing ui element
#' @export
betaInputUI <- function(id) {
    ns <- NS(id)

    box(width = 12, id = ns("analysisbox"),
        h4("ANALYSIS PARAMETERS"),
        
        div(
          h5("Beta Diversity"),
          selectInput(ns("distance"),
                      label = "Distance (beta)",
                      choices = c(
                        "Bray", "Canberra", "Jaccard", "Euclidean",
                        "Manhattan", "Clark", "Kulczynski", "Gower",
                        "altGower", "Morisita", "Horn", "Mountford",
                        "Raup", "Binomial", "Chao", "Cao", "Mahalanobis"
                      ),
                      selected = "Bray", multiple = FALSE,
                      selectize = FALSE, width = "250px"
          ),
          box(collapsible = TRUE, collapsed = TRUE, width = 11,
              style = "padding:0;",
              title = "*",
              p("Adonis refers to \"permutational multivariate analysis of variance using distance matrices\"
                 from the vegan package. The adonis variable specifies the column of the pheno data holding
                the independant variable whereas strata (optional) defines the groups within which to constrain 
                permutations. For more details and descriptions of the specific dissimilarity matrices, 
                please refer to the vegan package.")
          ),
          selectInput(ns("adonisvar"),
                      label = "Adonis variable (Optional)*",
                      choices = "",
                      selected = "", multiple = FALSE,
                      selectize = FALSE, width = "250px"
          ),
          selectInput(ns("adonisstrata"),
                      label = "Adonis strata (Optional)*",
                      choices = "",
                      selected = "", multiple = FALSE,
                      selectize = FALSE, width = "250px"
          )
        ),
        ## buttons are used to submit events. This ensures plots are not redrawn 
        ## while user still adjusts parameters
        div(id = ns("buttondiv"),
            fluidRow( width = 12, id = "actionbuttonrow",
                ## update analysis outputs (plots/tables)
                actionButton(ns("updatebutton"), 
                             icon = icon("far fa-chart-bar"),
                             label = HTML("&nbsp;UPDATE"), 
                             width = "120px"),
                actionButton(ns("reportButton"), 
                             label = HTML("<i class='far fa-bookmark'></i>&nbsp;&nbsp;REPORT"), 
                             width = "120px")
            )
        )
    )
}


#' Server side for the analysis input module handling analysis control
#'
#' @param input shiny input
#' @param output shiny output
#' @param session shiny session
#' @param meData MRExperiment object storing all data
#' @param adonisOptions phenodata colums ready for adonis analysis
#' @param reset reactive boolean determining if all inputs should be reset 
#' @author Janina Reeder
#' 
#' @importFrom metagenomeSeq MRcounts
#' @importFrom Biobase pData fData
#'
#' @return list holding all chosen values and the selected feature
#' @export
betaInput <- function(input, output, session, 
                          meData, adonisOptions, reset) {
    ns <- session$ns

    ## reactive value storing choices made in the UI
    chosenValues <- reactiveVal(NULL)

    ## reset all entries
    observe({
        req(reset())
        chosenValues(NULL)
    })

    ## initialize input elements
    observe({
      req(meData(), adonisOptions())
      isolate({
        updateSelectInput(session, "adonisvar", 
                          choices = c("", names(adonisOptions())), 
                          selected = chosenValues()$adonisvar)
        updateSelectInput(session, "adonisstrata", 
                          choices = c("", names(pData(meData()))), 
                          selected = chosenValues()$adonisstrata)
      })
    })
 
    ## main control button: store input choices in chosenValues
    observeEvent(input$updatebutton, {
      adonisvar <- NULL
      if(input$adonisvar != "")
        adonisvar <- input$adonisvar
      adonisstrata <- NULL
      if(input$adonisstrata != "")
        adonisstrata <- input$adonisstrata
      cV <- list("distance" = input$distance,
                 "adonisvar" = adonisvar,
                 "adonisstrata" = adonisstrata)
      chosenValues(cV)
    })

    return(chosenValues)
}


#' Heatmap analysis input module.
#' Set up to handle all analysis tabs in the app depending on given parameters
#'
#' @param id element identifier - namespace
#'
#' @author Janina Reeder
#'
#' @return box containing ui element
#' @export
heatmapInputUI <- function(id) {
  ns <- NS(id)
  
  box(width = 12, id = ns("analysisbox"), style = "padding-top: 0;",
      div(
        h5("Heatmap"),
        selectInput(ns("sorting"),
                    label = "Top features sorted by",
                    choices = c("Variance", "Fano", "MAD"),
                    selected = "Variance", multiple = FALSE,
                    selectize = FALSE, width = "250px"
        )
      ),
      selectizeInput(ns("featureselect"),
                     label = "Feature selection (optional)", choices = "", 
                     multiple = TRUE, 
                     width = "250px"
      ),
      ## buttons are used to submit events. This ensures plots are not redrawn 
      ## while user still adjusts parameters
      div(id = ns("buttondiv"),
          fluidRow( width = 12, id = "actionbuttonrow",
                    ## update analysis outputs (plots/tables)
                    actionButton(ns("updatebutton"), 
                                 icon = icon("far fa-chart-bar"),
                                 label = HTML("&nbsp;UPDATE"), 
                                 width = "120px"),
                    actionButton(ns("reportButton"), 
                                 label = HTML("<i class='far fa-bookmark'></i>&nbsp;&nbsp;REPORT"), 
                                 width = "120px")
          )
      )
  )
}


#' Server side for the analysis input module handling analysis control
#'
#' @param input shiny input
#' @param output shiny output
#' @param session shiny session
#' @param meData MRExperiment object storing all data
#' @param reset reactive boolean determining if all inputs should be reset 
#' @param aggDat aggregated MRExperiment object (default is NULL)
#'
#' @author Janina Reeder
#' 
#' @importFrom metagenomeSeq MRcounts
#' @importFrom Biobase pData fData
#'
#' @return list holding all chosen values and the selected feature
#' @export
heatmapInput <- function(input, output, session, 
                          meData, reset, aggDat = reactive(NULL)) {
  ns <- session$ns
  
  ## reactive value storing choices made in the UI
  chosenValues <- reactiveVal(NULL)
  ## reactive storing sorted features at given level
  aggFeatures <- reactiveVal(NULL)

  observe({
    req(aggDat())
    aggFeatures(rownames(MRcounts(aggDat())))
  })
  
  ## reset all entries
  observe({
    req(reset())
    chosenValues(NULL)
    aggFeatures(NULL)
  })
  
  ## update any inputs based on feature values
  observeEvent(aggFeatures(),{
    req(aggFeatures())
    selfeat <- ""
    if(!is.null(chosenValues()$featureselect))
      selfeat <- chosenValues()$featureselect
    corrdata  <- data.frame(value = aggFeatures(), 
                            label = aggFeatures())
    updateSelectizeInput(session, "featureselect",
                         choices = corrdata,
                         selected = selfeat,
                         options = list(placeholder = "Select Feature"),
                         server = TRUE)
  })
  
  
  ## main control button: store input choices in chosenValues
  observeEvent(input$updatebutton, {
    cV <- list("sorting" = input$sorting,
               "featureselect" = input$featureselect)
    chosenValues(cV)
  })
  
  return(chosenValues)
}

