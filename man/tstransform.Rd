% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transform.R
\name{tstransform}
\alias{tstransform}
\title{General transformation function}
\usage{
tstransform(method = "box-cox", lambda = NULL, lower = 0, upper = 1, ...)
}
\arguments{
\item{method}{valid methods are currently \dQuote{box-cox},
\dQuote{logit}, \dQuote{soft-logit} and \dQuote{sigmoid}.}

\item{lambda}{parameter in the Box Cox transformation.}

\item{lower}{lower bound for the transformations.}

\item{upper}{upper bound for the transformations.}

\item{\dots}{additional arguments taken by the transformations.}
}
\value{
A list with the transform and inverse functions.
}
\description{
Includes the Box Cox, logit, softplus-logit and sigmoif transforms.
Returns a list of functions for the transform and its inverse.
}
\author{
Alexios Galanos
}
