\name{tslinear}
\alias{tslinear}
\title{
Linear Time Series Filter
}
\description{
Provides a quick linear time series filter with trend and seasonality.
}
\usage{
tslinear(y, trend = FALSE, seasonal = FALSE, xreg = NULL, frequency = 1, ...)
}
\arguments{
\item{y}{
A vector.
}
\item{trend}{
Whether to include a linear trend.
}
\item{seasonal}{
Whether to include seasonal dummies.
}
\item{xreg}{
An optional matrix of regressors.
}
\item{frequency}{
The frequency of the series (required if seasonal is TRUE).
}
\item{\dots}{
Not currently used.
}
}
\details{
Estimates a simple linear time series model with trend, seasonal and regressors.
}
\value{
An object of class \dQuote{tslinear} which also inherits \dQuote{lm}.
}
\author{
Alexios Galanos
}
