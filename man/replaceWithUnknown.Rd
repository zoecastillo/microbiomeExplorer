% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/feature_module.R
\name{replaceWithUnknown}
\alias{replaceWithUnknown}
\title{Helper function to replace any un-annotated features with the term unknown}
\usage{
replaceWithUnknown(featcol)
}
\arguments{
\item{featcol}{vector of entried to be replaced where needed (fData column)}
}
\value{
modified featcol
}
\description{
Helper function to replace any un-annotated features with the term unknown
}
\examples{
data("mouseData", package = "metagenomeSeq")
featcol <- fData(mouseData)[["genus"]]
featcol[featcol == "NA"] <- NA
replaceWithUnknown(featcol)

}
\author{
Janina Reeder
}
