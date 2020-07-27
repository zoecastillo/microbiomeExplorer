library(shinyjs)
library(shinydashboard)
library(shinycssloaders)
library(shinyWidgets)
library(magrittr)
library(microbiomeExplorer)
library(metagenomeSeq)
library(rmarkdown)
library(DESeq2)

#MAX ALLOWED FILE SIZE for file upload
#modify as needed, but note that the app will become slow with larger data
options(shiny.maxRequestSize = 10 * 1024^2)
options(warn = -1)




