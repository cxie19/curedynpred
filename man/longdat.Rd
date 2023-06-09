% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/longdat.R
\docType{data}
\name{longdat}
\alias{longdat}
\title{A toy data set used for examples}
\format{
A data frame with 500 observations on the following 6 variables.
\describe{
\item{patient.id}{patient id}
\item{mes.time}{measurement time points for the biomarker}
\item{measure}{repeatedly measure biomarker values}
\item{trt}{a binary variable with 0 for treatment A and 1 for treatment B}
\item{event.time}{observed survival time}
\item{event}{an event indicator with 1 for dead and 0 for censored}
}
}
\source{
{Generated from JMFHC to serve as an example.}
}
\usage{
longdat
}
\description{
This toy data set with a cure fraction contains 500 patients' time-to-event outcomes, baseline information,
and a biomarker's repeated measured values.
}
\examples{
data(longdat)
}
\keyword{datasets}
