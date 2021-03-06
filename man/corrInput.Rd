% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/corr_input_module.R
\name{corrInput}
\alias{corrInput}
\title{Server side for the analysis input module handling analysis control}
\usage{
corrInput(
  input,
  output,
  session,
  type,
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

\item{type}{of the correlation (feature vs phenotype)}

\item{meData}{MRExperiment object storing all data}

\item{facetOptions}{named vector of available facet choices}

\item{reset}{reactive boolean determining if all inputs should be reset}

\item{aggDat}{aggregated MRExperiment object (default is NULL)}
}
\value{
list holding all chosen values and the selected feature
}
\description{
Server side for the analysis input module handling analysis control
}
\author{
Janina Reeder
}
