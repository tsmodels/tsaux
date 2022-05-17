\name{auto_regressors}
\alias{auto_regressors}
\title{
Automatic Detection of Outliers, Trends Breaks and Temporary Changes
}
\description{
A wrapper function for \emph{tso} from  the tsoutliers packages. Takes as input 
a univariate xts object and returns a list with an xts object with any 
identified outliers, trend breaks and/or temporary changes to be used as 
regressors during estimation as well initial coefficients (see details).
}
\usage{
auto_regressors(y, frequency = 1, lambda = NULL, forc_dates = NULL, 
sampling = NULL, h = 0, stlm_opts = list(etsmodel = "AAN"), 
auto_arima_opts = list(max.p = 1, max.q = 1, d = 1, allowdrift = FALSE), 
return_table = FALSE, method = c("sequential", "full"), ...)
}
\arguments{
\item{y}{
A univariate xts object.
}
\item{frequency}{
The frequency of the time series. If the frequency is 1 then seasonal estimation
will be turned off. Will also accept multiple seasonal frequencies.
}
\item{sampling}{
The sampling frequency the series. If h>0 and forc_dates is not provided, then
this is required in order to generate future time indices (valid values are
days, months, hours, mins, secs etc).
}
\item{lambda}{
An optional Box Cox transformation parameter. The routines are then run on the
transformed dataset.
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
\item{stlm_opts}{
Additional arguments to the stlm function.
}
\item{auto_arima_opts}{
Additional arguments to the auto.arima function in the tso routine.
}
\item{return_table}{
Whether to return a data.table instead with the anomalies detected rather
than an xts matrix with the pre-processed and ready to use anomalies.
}
\item{method}{
Whether to apply a sequential identification of anomalies using STL decomposition
in order to only pass the stationary residuals to the tso function, else 
to pass the series directly to the tso package.
}
\item{\dots}{
Any additional arguments passed to the tso functions (refer to the 
documentation of the tsoutliers package).
}
}
\details{
For generating future values of the identified outliers, the filter function is
used with additive outliers having a filter value of 0, trend changes a value of
1, and temporary changes have value between 0 and 1.
For the sequential method, the routine first interpolates any missing values, 
followed by an optional Box Cox transformation, and then elimination 
(and identification) of any outliers during the first pass. The cleaned series 
is then run through an stl filter (if any frequency is greater than 1) in 
order to deseasonalize the data (with multiple seasonality supported), after 
which the deseasonlized series is passed to the tso function where any additive 
outliers (AO), temporary shifts (TC) or level shift (LS) are identified. 
Additive outliers from this stage are added to any identified outliers from 
the initial stage. For each regressor, initial parameter values are returned 
together with the regressor matrix which should be passed to the estimation 
routine. This is critically important since in the absence of good parameter 
scaling, initial values are key to good convergence. Care should be taken with 
regards to any automatic Box Cox parameter estimation. In the presence of 
large outliers or level shifts, this is likely to be badly estimated which is 
why we do not allow automatic calculation of this, but instead place the burden 
on the user to decide what is a reasonable value (if any).
If a Box Cox transformation is used in the estimation routine, then it is important
to use the same lambda parameter in this function in order to get sensible results.
Again, avoid automatic Box Cox calculations throughout when you suspect significant
contamination of the series by outliers and breaks.
For the full method, the series is directly passed to the tso function of the
tsoutliers package.
Finally, it should be noted that this function is still experimental, and may
change in the future.
}
\value{
A list with an xts outlier matrix (if any where identified) as well as a vector
of initial parameter for use in the initialization of the optimizer.
}
\author{
Alexios Galanos for this wrapper function.\cr
Rob Hyndman for the forecast package.\cr
Javier López-de-Lacalle for the tsoutliers package.
}
\examples{
library(xts)
set.seed(200)
y = cumprod(c(100,(1+rnorm(100,0.01, 0.04))))
y = xts(y, as.Date(1:101, origin = as.Date("2000-01-01")))
yclean = y
outlier1 = rep(0, 101)
outlier1[20] = 0.35
outlier2 = rep(0, 101)
outlier2[40] = 0.25
outlier2 = as.numeric(filter(outlier2, filter = 0.6, method = "recursive"))
y = y + y*xts(outlier1, index(y))
y = y + y*xts(outlier2, index(y))
# may need some tweaking of the tso options.
x = auto_regressors(y, frequency = 1, sampling = "days", h = 20, 
check.rank = TRUE, discard.cval = 4)
head(x$xreg)
tail(x$xreg)
min(which(x$xreg[,1]==1))
min(which(x$xreg[,2]==1))
plot(as.numeric(y), type = "l", ylab = "")
lines(as.numeric(yclean) + (x$xreg \%*\% x$init)[1:101], col = 2)
}
