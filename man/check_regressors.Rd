\name{check_xreg}
\alias{check_xreg}
\alias{check_newxreg}
\title{
Checks on regressor matrix.
}
\description{
Used internally by other packages, these functions provides some commonly
used validation checks on regressor matrices in both in and out of sample.
}
\usage{
check_xreg(xreg, valid.index)
check_newxreg(newdata, xnames, h = 1, forc_dates = NULL)

}
\arguments{
\item{xreg}{
An xts matrix of named regressors.
}
\item{valid.index}{
A vector of dates against which the xreg matrix index is compared for validity.
}
\item{newdata}{
An xts matrix of out of named sample regressors.
}
\item{xnames}{
Names of regressors used insample.
}
\item{h}{
The forecast horizon
}
\item{forc_dates}{
An optional vector of forecast dates. This is used if newdata is not an 
xts matrix in which case it formats the data into such using the forc_dates
vector.
}
}
\value{
Returns the xts input matrix if checks are passed else raises an error. 
}
\author{
Alexios Galanos
}