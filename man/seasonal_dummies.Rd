\name{seasonal_dummies}
\alias{seasonal_dummies}
\title{
Seasonal Dummies
}
\description{
Creates a matrix of seasonal dummies.
}
\usage{
seasonal_dummies(y = NULL, n = nrow(y), seasons = 12)
}
\arguments{
\item{y}{
An optional data series.
}
\item{n}{
If y is missing, then the length of the series is required.
}
\item{seasons}{
The number of seasons in a cycle.
}
}
\details{
Generates seasons-1 dummy variables.
}
\value{
Either a matrix (if y is missing or y is not an xts vector) or
an xts matrix (when y is an xts vector).
}
\author{
Alexios Galanos
}
\examples{
head(seasonal_dummies(n=100, seasons=12))
}
