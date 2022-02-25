\name{box_cox}
\alias{box_cox}
\alias{auto_lambda}
\title{
Box-Cox transform specification
}
\description{
Creates a specification for the Box Cox transformtion.
}
\usage{
box_cox(lambda = NA, lower = 0, upper = 1.5, multivariate = FALSE, ...)
auto_lambda(y, lower = 0, upper = 1, nonseasonal_length = 2, frequency = 1, ...)
}
\arguments{
\item{lambda}{
The power parameters. If NA then it will automatically calculate the
optimal parameter using the method of Guerrero (for univariate case)
else for the multivariate case, the method of Velilla (1993)
which implemented in the \code{car} package of John Fox. This targets a
transformation to multivariate normality. If any of the inputs has a 
frequency other than 1, then an stl decomposition is first applied and
the seasonal component removed prior to the estimation in order to
avoid confounding the estimation by seasonality. It is also possible
to pass a vector equal to the number of columns of the dataset (with
numeric values mixed with NAs which will calculate the univariate
optimal lambda).
}
\item{lower}{
Optional parameter lower bound for cases when it is calculated.
}
\item{upper}{
Optional paraemter upper bound for cases when it is calculated.
}
\item{y}{
An xts vector of strictly positive values.
}
\item{frequency}{
The seasonal period of the data.
}
\item{nonseasonal_length}{
The measurement periodicity of the data. The maximum of this or the
frequency argument is taken when calculating the optimal lambda.
}
\item{multivariate}{
Flag for the multivariate case. If lambda is a single parameter,
then that is applied to all series (including NA which results in
the multivariate transformation described above).
}
\item{\ldots}{
Not currently used.
}
}
\details{
The function returns a list of 2 functions called \dQuote{transform}
and \dQuote{inverse} which can be called with a data object and
a frequency to calculate the transformed values. It is meant to be used
in the transform argument of the model specifications in the ts universe
of models.
The auto_lambda function uses the method of Guerrero(1993).
}
\value{
A list with the transform and invtransform functions. For the auto_lambda
function the optimal lambda is returned.
}
\references{
Box, G. E. P. and Cox, D. R. (1964),\emph{An analysis of transformations}. JRSS B \bold{26} 211--246.\cr
Guerrero, V.M. (1993), \emph{Time-series analysis supported by power transformations}. Journal of Forecasting, \bold{12}, 37--48.\cr
Velilla, S. (1993), \emph{A note on the multivariate Box-Cox transformation to normality}. Statistics and Probability Letters, \bold{17}, 259-263.
}
\author{
Alexios Galanos for the BoxCox function.\cr
John Fox for the powerTransform function used in the multivariate case.
}
\examples{
y = cumprod(c(1, 1 + rnorm(100,0.01, 0.005)))
B = box_cox(lambda = NA)
yt = B$transform(y, frequency = 1)
lambda = attr(yt,"lambda")
ye = B$inverse(yt, lambda)
}