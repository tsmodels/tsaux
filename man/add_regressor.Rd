% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulator.R
\name{add_regressor}
\alias{add_regressor}
\alias{add_regressor.issm.component}
\title{Regressor Component}
\usage{
add_regressor(x, ...)

\method{add_regressor}{issm.component}(x, xreg = NULL, pars = NULL, ...)
}
\arguments{
\item{x}{an object of class issm.component or other supported class.}

\item{...}{additional parameters.}

\item{xreg}{a matrix of regressors.}

\item{pars}{regressors coefficients.}
}
\value{
An object of class issm.component updated with the regressor
components.
}
\description{
Regressor Component
}
