.onLoad <- function(libname, pkgname){
  
  op <- options()
  
  op.microbiomeExplorer <- list(
    ## ALLOWED FILE EXTENSIONS (letter case is ignored)
    me.filetypes = c(
      ".biom", ".biom2", ".txt", ".tsv", ".csv",
      ".RDS", ".rds", ".rda", ".Rdata"
    ),
    me.delim = c(".tsv", ".csv", ".txt"),
    ## used to filter out columns in the feature annotation section
    ## that cannot be used for aggregation
    ## adjust as needed
    me.featurenames = c("Taxonomy","taxonomy","TAXONOMY",
                          "confidence","Confidence","CONFIDENCE",
                          "clusterCenter","clustercenter","CLUSTERCENTER"),
    me.corrmethods = c("spearman", "pearson", "kendall"),
    #removing ZILN as an option to avoid error in lmFit asking 
    #for identical ncol and nrow instead of equal
    me.diffmethods = c("DESeq2", "Kruskal-Wallis", "limma"),
    me.facetnum = 15,
    me.shapenum = 6,
    me.minwebgl = 500,
    ## COLOR PALETTE: Adjust as needed for different color options
    ## first should be the neutral color used if no color is specified
    me.colorvalues = c("#0C1707", "#9B110E", "#899DA4", "#FDD262", "#046C9A",
                    "#C7B19C", "#9986A5", "#4E2A1E", "#DC863B", "#273046", 
                    "#F4B5BD", "#446455", "#EBCC2A",  "#0B775E", "#C93312", 
                    "#9C964A"),
    me.round_digits = 3,
    me.modebar = list(list("toImage", "zoom2d", "pan2d", 
                                "select2d", "lasso2d", "resetScale2d")),
    ## AVAILABLE REPORT FORMATS (depends on pandoc version and pdflatex)
    me.reportformat = c("html_document", 
                       "pdf_document", 
                       "word_document",
                       "powerpoint_presentation"),
    me.reportformatshort = c("HTML", 
                             "PDF", 
                             "DOC", 
                             "PPT")
  )
  
  toset <- !(names(op.microbiomeExplorer) %in% names(op))
  if (any(toset)) options(op.microbiomeExplorer[toset])
  
  invisible()
}