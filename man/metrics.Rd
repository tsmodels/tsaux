% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metrics.R
\name{mape}
\alias{mape}
\alias{rmape}
\alias{smape}
\alias{mase}
\alias{mslre}
\alias{mis}
\alias{msis}
\alias{bias}
\alias{wape}
\alias{wslre}
\alias{wse}
\alias{pinball}
\alias{crps}
\title{Forecast Performance Metrics}
\usage{
mape(actual, predicted)

bias(actual, predicted)

mslre(actual, predicted)

mase(actual, predicted, original_series = NULL, frequency = 1)

mis(actual, lower, upper, alpha)

wape(actual, predicted, weights)

wslre(actual, predicted, weights)

wse(actual, predicted, weights)

pinball(actual, distribution, alpha = 0.1)

crps(actual, distribution)

smape(actual, predicted)

msis(actual, lower, upper, original_series, frequency = 1, alpha)
}
\arguments{
\item{actual}{the actual values corresponding to the forecast period.}

\item{predicted}{the predicted values corresponding to the forecast period.}

\item{original_series}{the actual values corresponding to the training
period.}

\item{frequency}{the seasonal frequency of the series used in the model.}

\item{lower}{the lower distributional forecast for the quantile
corresponsing to the coverage ratio alpha (i.e. alpha/2).}

\item{upper}{the upper distributional forecast for the quantile
corresponding to the coverage ratio alpha (i.e. 1 - alpha/2).}

\item{alpha}{the distributional coverage.}

\item{weights}{a vector of weights for generating weighted metrics. If the
actual and predicted inputs are univariate, this should be equal to the
length of the actual series and calculates a time weighted average, else the
weights should be of length equal to the number of series in a multivariate
case in which case a cross-sectional average is calculated.}

\item{distribution}{the forecast distribution (returned in the distribution
slot of the prediction object). This is used in the continuous ranked
probability score (crps) of Gneiting et al (2005), and calculated using the
function from the scoringRules package.}
}
\value{
A numeric value.
}
\description{
Functions to calculate a number of performance metrics.
}
\note{
The bias metric returns the percent bias. The rmape is the rescaled
measure for mape based on the paper by Swanson et al.
}
\references{
Tofallis (2015). A better measure of relative prediction
accuracy for model selection and model estimation, Journal of the
Operational Research Society \emph{66(8)}, 1352-1362.\cr Hyndman, Rob J and
Koehler, Anne B (2006). Another look at measures of forecast accuracy,
International Journal of Forecasting\emph{22(4)}, 679-688.\cr Gneiting,
Tilmann, Raftery, Adrian E, Westveld, Anton H and Goldman, Tom (2005).
Calibrated Probabilistic Forecasting Using Ensemble Model Output Statistics
and Minimum CRPS Estimation, Monthly Weather Review, \emph{133(5)},
1098-1118.\cr Gneiting, Tilmann and Raftery, Adrian E (2007). Strictly
proper scoring rules, prediction, and estimation, Journal of the American
statistical Association \emph{102(477)}, 359-378.\cr Swanson, D.A., Tayman,
J. and Bryan, T.M. (2011) MAPE-R: a rescaled measure of accuracy for
cross-sectional subnational population forecasts. J Pop Research \emph{28},
225–243.
}
\author{
Alexios Galanos
}
