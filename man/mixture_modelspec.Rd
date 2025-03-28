% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulator.R
\name{mixture_modelspec}
\alias{mixture_modelspec}
\title{Ensemble Setup}
\usage{
mixture_modelspec(...)
}
\arguments{
\item{...}{either a list of valid simulation objects or individual objects
passed to the function}
}
\value{
A object of class tssim.mixture ready for ensembling,
}
\description{
Ensemble Setup
}
\details{
The function performs certain checks on the inputs to ensure they
conform to the simulation models in the package and are of the same length.
}
