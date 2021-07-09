\name{auto_regressors}
\alias{auto_regressors}
\title{
Automatic Detection of Outliers, Trends Breaks and Temporary Changes
}
\description{
A wrapper function for the \emph{tso} function from the tsoutliers package which
takes as input a univariate xts object and returns an xts object with any 
identified outliers, trend breaks and/or temporary changes to be used as regressors
during estimation. Providing a forecast horizon in the function will also return
the values of those variables for the forecast horizon (see details).
}
\usage{
auto_regressors(y, frequency = 1, xreg = NULL, h = NULL, forc_dates = NULL, 
args.tsmethod = list(), ...)
}
\arguments{
\item{y}{
A univariate xts object.
}
\item{frequency}{
The frequency of the time series. If the frequency is 1 then seasonal estimation
will be turned off in the auto.arima function.
}
\item{xreg}{
An optional xts object of known regressors to include in the auto.arima estimation.
}
\item{h}{
An optional value for the forecast horizon (if planning to also use for 
prediction).
}
\item{forc_dates}{
An optional vector of Date to be used for indexing the series when h is not NULL.
If this is not provided then the sampling frequency of the series will be estimated
in order to generate this.
}
\item{args.tsmethod}{
Additional arguments to the auto.arima function.
}
\item{\dots}{
Any additional arguments passed to the tso function (refer to the tso 
function documentation).
}
}
\details{
For generating future values of the identified outliers, the filter function is
used with additive outliers having a filter value of 0, trend changes a value of
1, and temporary changes have value between 0 and 1.
}
\value{
Either an xts object if any outliers where identified else a NULL object.
}
\author{
Alexios Galanos for this wrapper function.\cr
Javier López-de-Lacalle for the tsoutliers package.
}
\examples{
library(xts)
set.seed(200)
y = cumprod(c(100,(1+rnorm(100,0.01, 0.04))))
y = xts(y, as.Date(1:101, origin = as.Date("2000-01-01")))
outlier1 = rep(0, 101)
outlier1[20] = 0.35
outlier2 = rep(0, 101)
outlier2[40] = 0.25
outlier2 = as.numeric(filter(outlier2, filter = 0.6, method = "recursive"))
y = y + y*xts(outlier1, index(y))
y = y + y*xts(outlier2, index(y))
# may need some tweaking of the tso options.
xreg = auto_regressors(y, frequency = 1, h = 20, check.rank = TRUE, discard.cval = 4)
head(xreg)
tail(xreg)
}
