\name{anomalies}
\alias{additive_outlier}
\alias{temporary_change}
\alias{level_shift}
\alias{anomalies}
\title{
Anomaly Creation
}
\description{
Creates specific types of anomalies given a series.
}
\usage{
additive_outlier(y, time = 1, parameter = 0.5, add = TRUE)
temporary_change(y, time = 1, parameter = 0.5, alpha = 0.7, add = TRUE) 
level_shift(y, time = 1, parameter = 0.5, add = TRUE)
}
\arguments{
\item{y}{
a univariate xts object or numeric series.
}
\item{time}{
the time index at which the anomaly takes place.
}
\item{parameter}{
the coefficient on the anomaly (the percent of the value of y at the specified 
time index representing the anomaly).
}
\item{alpha}{
the AR(1) coefficient for the temporary change which determines how quickly
the effect decays.
}
\item{add}{
whether to contaminate the series (add the anomaly to the series) else will
return a matrix with the anomaly (without the effect of the parameter).
}
}
\details{
These functions allow the generation of anomalies and may be chained together.
}
\value{
Either the contaminated series else a matrix of the anomaly.
}
\author{
Alexios Galanos for this wrapper function.\cr
}