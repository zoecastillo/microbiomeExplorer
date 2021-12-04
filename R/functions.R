### NON-PLOTTING, DATA-PROCESSING AND HELPER FUNCTIONS


#' Main function to start the Microbiome Explorer Shiny app via a command line
#' call
#'
#' @return the shiny application
#' 
#'
#' @export
runMicrobiomeExplorer <- function(){
  shiny::runApp(appDir = system.file("shiny", package = "microbiomeExplorer"))
}


## DATA PROCESSING ############################################################

#' Add phenotype data to object.
#'
#' This function adds phenotype data to the phenoData slot in an MRexperiment
#' object.
#'
#' @param MRobj An MRexperiment object.
#' @param phenodata Phenotype data frame or file path.
#'
#' @return An updated MRexperiment object.
#'
#' @importFrom metagenomeSeq loadPhenoData MRcounts
#' @importFrom Biobase phenoData<- pData<-
#'
#' @export
addPhenoData <- function(MRobj, phenodata = NULL) {
  if (is.null(phenodata)) {
    return(MRobj)
  }
  if (is.character(phenodata)) {
    phenodata <- utils::read.delim(phenodata, 
                               sep = 'if'(
                                 tools::file_ext(phenodata) == "csv", 
                                 ",", "\t"),
                               row.names = 1)
  }
  ord <- match(colnames(MRcounts(MRobj)), rownames(phenodata))
  if (any(is.na(ord))) {
    return(NULL)
  }
  phenodata <- phenodata[ord, ,drop = FALSE]
  phenoData(MRobj) <- Biobase::AnnotatedDataFrame(phenodata)
  if (!"SAMPLE_ID" %in% names(pData(MRobj))) {
    pData(MRobj) <- data.frame(SAMPLE_ID = rownames(pData(MRobj)),
                               pData(MRobj))
  }
  return(MRobj)
}

#' Extends existing phenodata for an object
#'
#' This function adds phenotype data to the phenoData slot in an MRexperiment
#' object.
#'
#' @param MRobj An MRexperiment object.
#' @param phenodata Phenotype data frame or file path.
#'
#' @return An updated MRexperiment object.
#'
#' @importFrom metagenomeSeq loadPhenoData MRcounts
#' @importFrom Biobase phenoData<-
#'
#' @export
extendPhenoData <- function(MRobj, phenodata = NULL) {
  if (is.null(phenodata)) {
    return(MRobj)
  }
  if (is.character(phenodata)) {
    phenodata <- as.data.frame(
      utils::read.delim(phenodata, 
                    sep = 'if'(
                      tools::file_ext(phenodata) == "csv", 
                      ",", "\t"),
                    row.names = 1))
  }
  if (!"SAMPLE_ID" %in% names(phenodata)) {
    phenodata <- data.frame(SAMPLE_ID = rownames(phenodata),
                            phenodata)
  }
  phenodata <- dplyr::left_join(pData(MRobj),phenodata)
  rownames(phenodata) <- phenodata$SAMPLE_ID
  phenoData(MRobj) <- Biobase::AnnotatedDataFrame(phenodata)
  
  return(MRobj)
}


#' Add feature data to MRobj.
#'
#' This function adds feature data to the featureData slot in an MRexperiment
#' object.
#'
#' @param MRobj An MRexperiment object.
#' @param featdata Feature data frame or file path.
#'
#' @return An updated MRexperiment object.
#'
#' @import metagenomeSeq 
#' @importFrom Biobase featureData<-
#'
#' @export
addFeatData <- function(MRobj, featdata = NULL) {
  if (is.null(featdata)) {
    return(MRobj)
  }
  if (is.character(featdata)) {
    featdata <- utils::read.delim(featdata,
                           sep = 'if'(
                             tools::file_ext(featdata) == "csv", ",", "\t"), 
                           row.names = 1,
                           stringsAsFactors = FALSE
    )
  }
  ord <- match(rownames(MRcounts(MRobj)), rownames(featdata))
  if (sum(is.na(ord)) > 0) {
    return(NULL)
  }
  featdata <- featdata[ord, ]
  featureData(MRobj) <- Biobase::AnnotatedDataFrame(featdata)
  return(MRobj)
}


#' Reads in data
#'
#' This function reads in an MRexperiment object saved as an RDS file, a Biom
#' file, or a tab - delimited count matrix with features as rows
#' and samples as columns.
#'
#' @param filepath Relative or absolute file path of data object.
#' @param type The type of file to be read; default is "RDS", other options are
#'  "RDATA", "BIOM", "TAB", "CSV"
#'
#' @return An MRexperiment object.
#'
#' @importFrom metagenomeSeq biom2MRexperiment loadMeta newMRexperiment
#' @importFrom Biobase pData<- fData<-
#' @importFrom methods is
#'
#' @export
readData <- function(filepath, type = "RDS") {
  if (type == "RDS") {
    bufobj <- readr::read_rds(filepath)
    MRobj <- newMRexperiment(counts = MRcounts(bufobj))
    pData(MRobj) <- pData(bufobj)
    fData(MRobj) <- fData(bufobj)
    if (!is.null(shiny::getDefaultReactiveDomain())) {
      incProgress(0.9)
    }
  } else if (type == "RDATA") {
    bufobj <- get(load(filepath))
    MRobj <- newMRexperiment(counts = MRcounts(bufobj))
    pData(MRobj) <- pData(bufobj)
    fData(MRobj) <- fData(bufobj)
    if (!is.null(shiny::getDefaultReactiveDomain())) {
      incProgress(0.9)
    }
  } else if (type == "BIOM") {
    if (!is.null(shiny::getDefaultReactiveDomain())) {
      incProgress(message = "Reading biom file", 
                  detail = "\nConversion in progress ...")
    }
    Bobj <- biomformat::read_biom(filepath)
    if (!is.null(shiny::getDefaultReactiveDomain())) {
      incProgress(0.3)
    }
    MRobj <- biom2MRexperiment(Bobj)
    if (!is.null(shiny::getDefaultReactiveDomain())) {
      incProgress(0.3)
    }
    featdata <- biomformat::observation_metadata(Bobj)
    if (!is.null(shiny::getDefaultReactiveDomain())) {
      incProgress(0.3)
    }
  } else if (type == "TAB") {
    countdata <- loadMeta(filepath)
    MRobj <- newMRexperiment(countdata$counts)
  } else if (type == "CSV") {
    countdata <- loadMeta(filepath, sep = ",")
    MRobj <- newMRexperiment(countdata$counts)
  } else {
    return(NULL)
  }
  if(!is(MRobj, "MRexperiment")){
    return(NULL)
  } 
  ## add required columns if they aren't available yet
  if (!"SAMPLE_ID" %in% names(pData(MRobj))) {
    pData(MRobj) <- data.frame(SAMPLE_ID = rownames(pData(MRobj)),
                               pData(MRobj))
  }
  
  if(!is.null(pData(MRobj))){
    pData(MRobj) <- pData(MRobj) %>%
      tibble::rownames_to_column("samname") %>%
      dplyr::mutate_if(is.character, as.factor) %>%
      tibble::column_to_rownames("samname")
  }
  
  if(!is.null(fData(MRobj))){
    fData(MRobj) <- fData(MRobj) %>%
      tibble::rownames_to_column("samname") %>%
      dplyr::mutate_if(is.character, as.factor) %>%
      tibble::column_to_rownames("samname")
  }
  
  return(MRobj)
}

#' Helper function assigning different file extensions to specific short texts 
#' identifying the types
#'
#' @param fileext the file extension found after '.'
#'
#' @author Janina Reeder
#'
#' @return character string for the filetype
getFileType <- function(fileext) {
  switch(tolower(fileext),
         "rds" = "RDS",
         "rda" = "RDATA",
         "rdata" = "RDATA",
         "tsv" = "TAB",
         "txt" = "TAB",
         "csv" = "CSV",
         "csv2" = "CSV",
         "biom" = "BIOM",
         "biom2" = "BIOM"
  )
}



#' Function to filter the MRexperiment data by numerical parameters
#'
#' @param MRobj MRExperiment object to filter
#' @param minpresence minimum sample presence per feature
#' @param minfeats minimum number of features per sample
#' @param minreads minimum number of reads per sample
#'
#' @author Janina Reeder
#'
#' @return the filtered MRobj
#' 
#' @examples 
#' data("mouseData", package = "metagenomeSeq")
#' filterMEData(MRobj = mouseData, minpresence = 4, minfeats = 300)
#' 
#' @export
filterMEData <- function(MRobj, minpresence = 1, 
                         minfeats = 2, minreads = 2) {
  filteredObj <- MRobj
  while(TRUE){
    if(nrow(filteredObj) == 0 || ncol(filteredObj) == 0) 
      break
    rowsToRemove <- rowSums(MRcounts(filteredObj) > 0) < minpresence
    colsToRemove <- colSums(MRcounts(filteredObj) > 0) < minfeats |
      colSums(MRcounts(filteredObj)) < minreads
    if(any(rowsToRemove) || any(colsToRemove)){
        filteredObj <- filteredObj[!rowsToRemove,!colsToRemove]
    } else {
      break
    }
  }

  return(filteredObj)
}



#' Function to filter the MRexperiment by certain phenotype values
#'
#' @param MRobj the MRexperiment to subset
#' @param rm_phenovalues list of named vectors with names corresponding to 
#' column names in pData and values representing phenotypes within the column
#'
#' @author Janina Reeder
#'
#' @return the filtered MRobj
#' 
#' @examples 
#' data("mouseData", package = "metagenomeSeq")
#' filterByPheno(MRobj = mouseData, 
#'   rm_phenovalues = list("diet" = c("BK"),"mouseID" = c("PM1","PM10")))
#' 
#' @export
filterByPheno <- function(MRobj, rm_phenovalues) {
  samstoremove <- rowSums(
    dplyr::bind_cols(lapply(names(rm_phenovalues), 
                            function(r) {
                              naind <- which(rm_phenovalues[[r]] == "NA")
                              rm_phenovalues[[r]][naind] <- NA
                              pData(MRobj)[, r] %in% rm_phenovalues[[r]]
                            })
    )
  )
  return(MRobj[, !samstoremove])
}


## PHENOTABLE 


#' Helper function used to build a correct interactionName based on the chosen 
#' columns
#'
#' @param interactionName as chosen by user. This may not be good to store 
#' internally
#'
#' @return updated interactionName or warning/error string
parseInteractionName <- function(interactionName) {
  if (is.null(interactionName) | interactionName == "") {
    return("ERROR: No interaction name specified")
  }
  if (interactionName != make.names(interactionName)) {
    return(paste0("WARNING: ", 
                  gsub("\\.+", "_", make.names(interactionName))))
  }
  return(interactionName)
}


#' Makes header for R script
#'
#' This function makes the header for the report R script to be rendered by
#' knitr into Rmarkdown and rendered into HTML, PDF, or Word.
#'
#' This was adapted from https://yihui.name/knitr/demo/stitch/
#'
#' @param title Title of the report.
#' @param author Author of the report.
#' @param date Date of the report.
#' @param data.source R code used to obtain the dataset
#' @param output Output of Rmarkdown file. 
#' @param toc Table of contents. Default is TRUE.
#'
#' @return A character vector where each element is a line in the R script.
createHeader <- function(title = "MicrobiomeExplorer Report", 
                         author = "",
                         date = "",
                         data.source = "",
                         output = getOption("me.reportformat"),
                         toc = TRUE) {
  cat(file = stderr(), "-- creating header --\n")
  
  if (is.null(output)) {
    output <- "html_document"
  }
  # output <- match.arg(output)
  output <- c(output, "md_document")
  
  cat(file = stderr(), "-- setting toc ---")
  if (toc) {
    output2 <- unlist(lapply(output, function(o){
      paramlist <- c(paste0(o, ":"),
                     "  toc: true")
      if(o %in% "powerpoint_presentation"){
        paramlist <- c(paramlist,"  slide_level: 3")
      }
      return (paramlist)
    }))
  } else {
    output2 <- paste0(output, ": default")
  }
  
  cat(file = stderr(), "-- data processing ---")
  
  ep_project <- 'if'(stringr::str_starts(data.source, "proj"), TRUE, FALSE)
  
  header <- c(
    "---",
    paste0("title: ", title),
    paste0("author: ", author),
    paste0("date: ", date),
    paste0("output:"),
    paste0("  ", output2),
    paste0("always_allow_html: true"),
    paste0("font-size: 11pt"),
    "---"
  )
  
  header <- paste0("#' ", header)
  header <- c(
    header, "", "#+ setup, include=FALSE",
    "knitr::opts_chunk$set(echo=FALSE, message=FALSE,
              warning=FALSE, cache=FALSE)"
  )
  header <- c(
    header, "#'", "",
    "#+ load_package, cache=FALSE",
    "library(knitr)",
    "library(microbiomeExplorer)",
    "library(magrittr)",
    "library(metagenomeSeq)",
    "library(shiny)",
    'if'(ep_project, "library(DataSetDB)", ""),
    "#'", "",
    "",
    "",
    "#' This report has been generated from the Microbiome Explorer app",
    "#' Data is stored externally and must be available to render", 
    "#' the Rmarkdown document\n\n",
    "#'", "",
    "doctype <- opts_knit$get(\"rmarkdown.pandoc.to\")",
    "options(me.round_digits = 3)",
    "options(me.modebar = list(list(\"toImage\", \"zoom2d\", \"pan2d\",", 
    "\t\"select2d\", \"lasso2d\", \"resetScale2d\")))"
  )
  header
}

#' Generates report
#'
#' This function generates the pieces of the report, which includes the R
#' script, Rmarkdown file, and any Rmarkdown outputs.
#'
#' Adapted from https://yihui.name/knitr/demo/stitch/
#'
#' @param rcode A named list where each element corresponds to a different
#' analysis (Alpha diversity, Beta diversity). The name of the list is used to
#' denote the first part of the code chunks in each analysis section
#' (alpha, beta). Each element is itself a list of R commands corresponding to
#' a code chunk. 
#' @param filename Name of output files. Default is "report".
#' @param dir Directory of output. Default is "out".
#' @param title Title of the report.
#' @param author Author of the report.
#' @param date Date of the report.
#' @param data.source R code used to obtain the dataset
#' @param output Output of Rmarkdown file. Options defined in global.R
#' @param toc Table of contents. Default is TRUE.
#' @param intro_text Introductory text to include with the report (optional)
#'
#' @return A character vector where each element is a line in the R script.
generateReport <- function(rcode,
                           filename = "report",
                           dir = "out",
                           title = "MicrobiomeExplorer Report",
                           author = "",
                           date = "`r format(Sys.time(), '%d %B, %Y')`",
                           data.source = "",
                           output = c("html_document"),
                           toc = TRUE,
                           intro_text = NULL) {
  cat(file = stderr(), "-- entering my report generator --\n")
  
  ## Generate body code and combine with header
  header <- createHeader(title, author, date, data.source, 
                         output, toc)
  cat(file = stderr(), "-- header created! --\n")
  
  if (!is.null(intro_text) && length(intro_text) > 0) {
    intro_text <- paste("#' ## Introduction\n", intro_text, "\n\n  ")
    intro_text <- vapply(intro_text, stringr::str_replace_all,
                         pattern = "\n", replacement = "\n#' ", 
                         FUN.VALUE = character(1))
  } else {
    intro_text <- NULL
  }
  
  cat(file = stderr(), "-- build intro text --\n")
  
  ep_project <- 'if'(stringr::str_starts(data.source, "proj"), TRUE, FALSE)
  dataloading <- c("#' ## Data loading and preparation",
                   data.source)
  
  report <- c(header, intro_text, dataloading, 
              unlist(rcode, use.names = FALSE))
  
  cat(file = stderr(), "-- merged code --\n")
  
  ## handle files, processing and output
  rscript <- file.path(dir, paste0(filename, ".R"))
  cat(report, sep = "\n", file = rscript)
  knitr::spin(rscript, knit = FALSE)
  rmarkdown::render(rscript, output_format = output)
  out.ext <- c("Rmd", 
               gsub("powerpoint", "pptx", 
                    gsub("word", "docx", 
                         vapply(strsplit(output, split = "_"), `[[`, 1,
                                FUN.VALUE = character(1)))))
  print(c(rscript, paste(paste(dir, filename, sep = "/"), 
                         out.ext, sep = ".")))
}


## HELPER FUNCTIONS

#' Sets up a dataframe used by several plotting functions by joining
#' the required data with relevant phenotype data
#'
#' @param df dataframe storing plotting data values
#' @param phenoTable pData of the MRexperiment; all following parameters must be 
#' a column of the phenoTable
#' @param x_var main plotting variable
#' @param facet1 column-based faceting (can be NULL)
#' @param facet2 row-based faceting (can be NULL)
#' @param col_by coloring factor (can be NULL)
#' @param col_name character to be used as name for col_by
#' @param id_var variable used to connect samples longitudinally (can be NULL)
#'
#' @return dataframe obtained by joining df and relevant columns of phenoTable
buildPlottingDF <- function(df, phenoTable, x_var = NULL, 
                            facet1 = NULL, facet2 = NULL, 
                            col_by = NULL, col_name = col_by,
                            id_var = NULL){
  sample_data <- data.frame(
    samname = rownames(phenoTable),
    X1 = phenoTable[, x_var],
    X2 = phenoTable[, facet1],
    X3 = phenoTable[, facet2],
    X4 = phenoTable[, col_by],
    X5 = phenoTable[, id_var]
  )
  colnames(sample_data) <- c("samname", x_var, facet1, facet2, col_name, id_var)
  sample_data[] <- lapply(sample_data, function(s) gsub("^$","NA",s))
  sample_data[] <- lapply(sample_data, factor, exclude = NULL)
  sample_data[] <- lapply(sample_data, 
                          forcats::fct_explicit_na, na_level = "NA")
  suppressMessages(dplyr::left_join(df, sample_data))
}




#' Helper function to account for issues plotly has with very small widths
#' (these end up being 1 and cover the entire plotting area)
#'
#' @param df2 dataframe storing plotting data
#' @param facets column facets
#' @param x_var x variable
#' @param drop passed on as .drop to dplyr::group_by
#'
#' @return widths for each facet
getWidths <- function(df2, facets, x_var, drop = TRUE){
  numinfacet <- df2 %>%
    dplyr::group_by_at(c(facets, x_var), .drop = drop)  %>%
    dplyr::tally() %>%
    dplyr::group_by_at(facets) %>%
    dplyr::summarise(maxnum = length(x_var))
  totalwidth <- sum(numinfacet$maxnum)
  widths <- numinfacet$maxnum/totalwidth
  widths <- ifelse(widths < 0.025, 0.025, widths)
  if(sum(widths) > 1){
    adj_ind <- which(widths > 0.025)
    diff <- (sum(widths) - 1) / length(adj_ind)
    widths[adj_ind] <- widths[adj_ind] - diff
  }
  widths
}

#' Function to find a non-empty facet in the last row. This will be the
#' one to be connected to the plot legend to avoid duplicates within
#'
#' @param df2 plotting data frame
#' @param facets column facets
#' @param facet2s row facets
#'
#' @return the name of the column-based facet which can be used as legend
getLegendLevel <- function(df2, facets, facet2s){
  legendlevel <- "nofacets"
  grouptally <- df2 %>%
    dplyr::group_by_at(c(facets, facet2s)) %>% 
    dplyr::tally() 
  if(nrow(grouptally) > 0)
    legendlevel <- grouptally[nrow(grouptally),"facets"]
  legendlevel
}


#' Creates an empty plotly plot using the given labels on the x and y axis
#'
#' @param xaxis_text x axis label
#' @param ylab y axis label
#'
#' @return call to plotly_empty
buildEmptyPlotlyPlot <- function(xaxis_text, ylab){
  plotly::plotly_empty(type = "scatter", 
                       mode = "markers") %>%
    plotly::layout(
      xaxis = list(
        title = list(
          text = xaxis_text,
          font = list(
            size = 12
          )
        ),
        showgrid = FALSE,
        visible = TRUE,
        tickfont = list(
          size = 9,
          color = "#5b5b5b"
        )
      ),
      yaxis = list(
        title = list(
          text = ylab,
          font = list(
            size = 12
          )
        ),
        visible = TRUE,
        showgrid = TRUE,
        showticklabels = TRUE,
        tickfont = list(
          size = 9,
          color = "#5b5b5b"
        )
      )
    )
}

#' Adds a layout call based on plotly::layout
#'
#' @param plotTitle plot title to use
#' @param xaxis_text x axis label to use
#' @param ylab y axis label to use
#' @param .data plotly data object to apply the layout call to
#'
#' @return plotly::layout call
add_plotly_layout <- function(.data, plotTitle, xaxis_text, ylab){
  plotly::layout(.data,
    #barmode = "stack",
    title = list(
      text = plotTitle,
      font = list(
        size = 14
      )
    ),
    xaxis = list(
      title = list(
        text = xaxis_text,
        font = list(
          size = 12
        )
      ),
      showgrid = FALSE,
      visible = TRUE,
      tickangle = 45,
      tickfont = list(
        size = 9,
        color = "#5b5b5b"
      )
    ),
    yaxis = list(
      title = list(
        text = ylab,
        font = list(
          size = 12
        )
      ),
      visible = TRUE,
      showgrid = TRUE,
      showticklabels = TRUE,
      tickfont = list(
        size = 9,
        color = "#5b5b5b"
      )
    ),
    margin = list(
      b = 150,r = 50, t = 25, l = 50
    ),
    legend = list(
      font = list(
        family = "sans-serif",
        size = 11,
        color = "#000"
      ),
      bordercolor = "#E2E2E2",
      borderwidth = 2
    )
  )
}

#' Adds a config call based on plotly::config
#'
#' @param .data plotly data object to apply the config call to
#'
#' @return plotly::config call
add_plotly_config <- function(.data){
  plotly::config(.data,
    toImageButtonOptions = list(
      format = "svg"
    ),
    displaylogo = FALSE,
    modeBarButtons = list(list("toImage", "zoom2d", "pan2d", 
                               "select2d", "lasso2d", "resetScale2d"))
  )
}