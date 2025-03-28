% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transform.R
\name{box_cox}
\alias{box_cox}
\title{Box-Cox transform specification}
\usage{
box_cox(lambda = NA, lower = 0, upper = 1.5, multivariate = FALSE, ...)
}
\arguments{
\item{lambda}{the power parameters. If NA then it will automatically
calculate the optimal parameter using the method of Guerrero (for univariate
case) else for the multivariate case, the method of Velilla (1993) which
implemented in the \code{car} package of John Fox. This targets a
transformation to multivariate normality. If any of the inputs has a
frequency other than 1, then an stl decomposition is first applied and the
seasonal component removed prior to the estimation in order to avoid
confounding the estimation by seasonality. It is also possible to pass a
vector equal to the number of columns of the dataset (with numeric values
mixed with NAs which will calculate the univariate optimal lambda).}

\item{lower}{optional parameter lower bound for cases when it is calculated.}

\item{upper}{optional parameter upper bound for cases when it is calculated.}

\item{multivariate}{flag for the multivariate case. If lambda is a single
parameter, then that is applied to all series (including NA which results in
the multivariate transformation described above).}

\item{...}{not currently used.}
}
\value{
A list with the transform and inverse functions.
}
\description{
Creates a specification for the Box Cox transformation.
}
\details{
The function returns a list of 2 functions called \dQuote{transform} and
\dQuote{inverse} which can be called with a data object and a frequency to
calculate the transformed values. It is meant to be used in the transform
argument of the model specifications in the ts universe of models. The
auto_lambda function uses the method of Guerrero(1993).
}
\note{
The returned transform function will take additional argument \dQuote{frequency}
which determines whether a series is seasonal or not. When estimating lambda (when
setting this to NA), a series with frequency > 1 will first be de-seasonalized
using an STL decomposition.
}
\examples{

y = cumprod(c(1, 1 + rnorm(100,0.01, 0.005)))
B = box_cox(lambda = NA)
yt = B$transform(y, frequency = 1)
lambda = attr(yt,"lambda")
ye = B$inverse(yt, lambda)

}
\references{
\insertRef{Box1964}{tsaux}\cr
\insertRef{Velilla1993}{tsaux}\cr
\insertRef{Guerrero1993}{tsaux}\cr
}
\author{
Alexios Galanos for the BoxCox function.\cr John Fox for the
powerTransform function used in the multivariate case.
}
