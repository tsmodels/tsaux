% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transform.R
\name{softlogit}
\alias{softlogit}
\title{The softplus logit transformation}
\usage{
softlogit(lower = 0, upper = 1, ...)
}
\arguments{
\item{lower}{lower bound of the variable.}

\item{upper}{upper bound of the variable.}

\item{...}{not currently used.}
}
\value{
A list with the transform and inverse functions.
}
\description{
The softplus logit transformation is an alternative to the logit transform for bounded
outcomes with positive output.
}
\examples{

y = cumprod(c(1, 1 + rnorm(100,0.01, 0.005)))
B = softlogit(lower = 0,  upper = 15)
yt = B$transform(y)
ye = B$inverse(yt)
}
\author{
Alexios Galanos
}
