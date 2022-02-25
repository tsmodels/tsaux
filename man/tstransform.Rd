\name{tstransform}
\alias{tstransform}
\title{
General transformation function
}
\description{
Includes the Box Cox and logit transforms. Returns a function for transform
and inverse.
}
\usage{
tstransform(method = "box-cox", lambda = NULL, lower = 0, upper = 1, ...)
}
\arguments{
  \item{method}{
Currently, valid methods are \dQuote{box-cox} and \dQuote{logit}.
}
\item{lambda}{
Parameter in the Box Cox transformation.
}
\item{lower}{
Lower bound for the transformations.
}
\item{upper}{
Upper bound for the transformations.
}
\item{\dots}{
Additional arguments taken by the transformations.
}
}
\details{
Returns a list with functions for transform and its inverse.
}
\value{
A list of functions
}
\author{
Alexios Galanos
}
