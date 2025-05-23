% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/train_dat.R
\docType{data}
\name{train_dat}
\alias{train_dat}
\title{A sample training data set used for examples}
\format{
A data frame with 500 observations on the following 6 variables.
\describe{
\item{id}{patient id}
\item{measure.time}{measurement time points for the biomarker}
\item{measure}{repeatedly measure biomarker values}
\item{x1}{a binary variable with values 0 and 1}
\item{x2}{a continous variable following a standard normal distribution}
\item{event.time}{observed survival time}
\item{event}{an event indicator with 1 for event of interest and 0 for censored}
}
}
\source{
{Generated from JMFHC to serve as an example.}
}
\usage{
train_dat
}
\description{
This training data set with a cure fraction contains 500 patients' time-to-event outcomes, baseline information (X1 and X2),
and a biomarker's repeated measured values (measure).
}
\examples{
data(train_dat)
}
\keyword{datasets}
