\name{seasonality_test}
\alias{seasonality_test}
\title{
Simple Seasonality Test
}
\description{
Checks for the presence of seasonality based on the QS test of Gomez  and  Maravall (1996).
}
\usage{
seasonality_test(x, frequency = NULL)
}
\arguments{
\item{x}{
An (xts) vector (usually of a stationary series).
}
\item{frequency}{
Ovverides any frequency automatically identified in the index of x.
}
}
\details{
Given the identified frequency of the xts vector (using the \code{\link{sampling_frequency}}),
the function checks for seasonality at that frequency. The frequency can be ovveriden by
directly supplying a frequency argument, in which case y does not need to be a xts vector.
}
\value{
Logical.
}
\references{
Gomez, Victor and Maravall, Agustin, \emph{Programs TRAMO and SEATS}, 1986.
}
\author{
Alexios Galanos
}
