\name{auto_regressors}
\alias{auto_regressors}
\title{
Automatic Detection of Outliers, Trends Breaks and Temporary Changes
}
\description{
A wrapper function for either the \emph{UC} or \emph{tso} functions from 
the UComp and tsoutliers packages respectively. Takes as input a univariate xts 
object and returns a list with an xts object with any 
identified outliers, trend breaks and/or temporary changes to be used as 
regressors during estimation as well as any ar and/or ma orders. Providing a 
forecast horizon in the function will also return
the values of those variables for the forecast horizon (see details).
}
\usage{
auto_regressors(y, frequency = 1, lambda = NULL, xreg = NULL, h = NULL, 
forc_dates = NULL, use = c("UC","tso"), args.tsmethod = list(), ...)
}
\arguments{
\item{y}{
A univariate xts object.
}
\item{frequency}{
The frequency of the time series. If the frequency is 1 then seasonal estimation
will be turned off in the auto.arima or UC function.
}
\item{lambda}{
An optional Box Cox transformation parameter. The routines are then run on the
transformed dataset.
}
\item{xreg}{
An optional xts object of known regressors to include in the estimation.
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
\item{use}{
The function to use. The tso function makes use of the appoach described in 
Chen & Liu (1993), whilst the UC function makes use of the approach described
in de Jong, P. & Penzer, J. (1998) based on smoothations from the disturbance
smoother.
}
\item{args.tsmethod}{
Additional arguments to the auto.arima function.
}
\item{\dots}{
Any additional arguments passed to the tso or UC functions (refer to the 
documentation of those functions). For the UC function, it is suggested that
the use pass additional information, particularly if dealing with series with
multiple seasonality (which require the period argument with harmonics for each
frequency).
}
}
\details{
For generating future values of the identified outliers, the filter function is
used with additive outliers having a filter value of 0, trend changes a value of
1, and temporary changes have value between 0 and 1.
}
\value{
A list with an xts outlier matrix (if any where identified) as well as the 
order of the ar and ma terms identified.
}
\author{
Alexios Galanos for this wrapper function.\cr
Diego J. Pedregal for the UComp package.\cr
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
x = auto_regressors(y, frequency = 1, h = 20, check.rank = TRUE, discard.cval = 4)
head(x$xreg)
tail(x$xreg)
}
