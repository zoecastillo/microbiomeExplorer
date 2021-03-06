% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_input_module.R
\name{getFilterChoices}
\alias{getFilterChoices}
\title{Helper function to filter phenodata for interesting phenotypes to be used 
for filtering or subsetting}
\usage{
getFilterChoices(MRobj)
}
\arguments{
\item{MRobj}{the MRexperiment storing the data}
}
\value{
list of named vectors with names being pData column headers and 
values being unique entries; columns with only one entry or those with 
different values for each samples are omitted
}
\description{
Helper function to filter phenodata for interesting phenotypes to be used 
for filtering or subsetting
}
\author{
Janina Reeder
}
