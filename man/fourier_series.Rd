\name{fourier_series}
\alias{fourier_series}
\title{
Fourier terms for modelling seasonality
}
\description{
Returns a matrix containing terms from a Fourier series, up to order K
}
\usage{
fourier_series(dates, period = NULL, K = NULL)
}
\arguments{
\item{dates}{
A Date vector representing the length of the series for which
the fourier terms are required.
}
\item{period}{
The frequency of the underlying series, if NULL will try to infer
it from the difference in the Date vector.
}
\item{K}{
Maximum order of the Fourier terms.
}
}
\value{
A matrix of size N (length of dates) by 2*K.
}
