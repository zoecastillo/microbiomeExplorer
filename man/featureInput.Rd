% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/feature_input_module.R
\name{featureInput}
\alias{featureInput}
\title{Server side for the feature analysis input module}
\usage{
featureInput(
  input,
  output,
  session,
  meData,
  facetOptions = NULL,
  reset,
  aggDat = reactive(NULL)
)
}
\arguments{
\item{input}{shiny input}

\item{output}{shiny output}

\item{session}{shiny session}

\item{meData}{MRExperiment object storing all data}

\item{facetOptions}{named vector of available facet choices}

\item{reset}{reactive boolean determining if all inputs should be reset}

\item{aggDat}{aggregated MRExperiment object (default is NULL)}
}
\value{
list holding all chosen values and the selected feature
}
\description{
Server side for the feature analysis input module
}
\author{
Janina Reeder
}
