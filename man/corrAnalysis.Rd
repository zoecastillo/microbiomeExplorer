% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/corr_analysis_module.R
\name{corrAnalysis}
\alias{corrAnalysis}
\title{corr Analysis Module - server}
\usage{
corrAnalysis(
  input,
  output,
  session,
  data,
  levelOpts,
  chosenLevel,
  resetInput,
  aggData
)
}
\arguments{
\item{input}{shiny input}

\item{output}{shiny output}

\item{session}{shiny session}

\item{data}{the main data object returned from data_input_module}

\item{levelOpts}{available levels to aggregate on (depends on input data)}

\item{chosenLevel}{previously selected level (passed from different instance)}

\item{resetInput}{reactive boolean determining if reset is required}

\item{aggData}{the aggregated MRExperiment object}
}
\value{
reactive holding code to be used in reports
}
\description{
corr Analysis Module - server
}
