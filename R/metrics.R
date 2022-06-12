#' Forecast Performance Metrics
#' 
#' Functions to calculate a number of performance metrics.
#' 
#' 
#' @aliases mape rmape smape mase mslre mis msis bias wape wslre wse
#' pinball crps
#' @param actual the actual values corresponding to the forecast period.
#' @param predicted the predicted values corresponding to the forecast period.
#' @param original_series the actual values corresponding to the training
#' period.
#' @param frequency the seasonal frequency of the series used in the model.
#' @param lower the lower distributional forecast for the quantile
#' corresponsing to the coverage ratio alpha (i.e. alpha/2).
#' @param upper the upper distributional forecast for the quantile
#' corresponding to the coverage ratio alpha (i.e. 1 - alpha/2).
#' @param alpha the distributional coverage.
#' @param distribution the forecast distribution (returned in the distribution
#' slot of the prediction object). This is used in the continuous ranked
#' probability score (crps) of Gneiting et al (2005), and calculated using the
#' function from the scoringRules package.
#' @param weights a vector of weights for generating weighted metrics. If the
#' actual and predicted inputs are univariate, this should be equal to the
#' length of the actual series and calculates a time weighted average, else the
#' weights should be of length equal to the number of series in a multivariate
#' case in which case a cross-sectional average is calculated.
#' @return A numeric value.
#' @export
#' @rdname metrics
#' @note The bias metric returns the percent bias. The rmape is the rescaled
#' measure for mape based on the paper by Swanson et al.
#' @author Alexios Galanos
#' @references Tofallis (2015). A better measure of relative prediction
#' accuracy for model selection and model estimation, Journal of the
#' Operational Research Society \emph{66(8)}, 1352-1362.\cr Hyndman, Rob J and
#' Koehler, Anne B (2006). Another look at measures of forecast accuracy,
#' International Journal of Forecasting\emph{22(4)}, 679-688.\cr Gneiting,
#' Tilmann, Raftery, Adrian E, Westveld, Anton H and Goldman, Tom (2005).
#' Calibrated Probabilistic Forecasting Using Ensemble Model Output Statistics
#' and Minimum CRPS Estimation, Monthly Weather Review, \emph{133(5)},
#' 1098-1118.\cr Gneiting, Tilmann and Raftery, Adrian E (2007). Strictly
#' proper scoring rules, prediction, and estimation, Journal of the American
#' statistical Association \emph{102(477)}, 359-378.\cr Swanson, D.A., Tayman,
#' J. and Bryan, T.M. (2011) MAPE-R: a rescaled measure of accuracy for
#' cross-sectional subnational population forecasts. J Pop Research \emph{28},
#' 225–243.
mape <- function(actual, predicted)
{
    actual <- as.numeric(actual)
    predicted <- as.numeric(predicted)
    n <- length(actual)
    if (length(predicted) != n) stop("\nactual and predicted must be of the same length.")
    if (any(is.na(actual))) {
        exc <- which(is.na(actual))
        actual <- actual[-exc]
        predicted <- predicted[-exc]
    }
    n <- length(actual)
    if (n == 0) return(NA)
    metric <- sum(abs((actual - predicted)/actual))/n
    return(metric)
}

#' @export
#' @rdname metrics
bias <- function(actual, predicted)
{
    actual <- as.numeric(actual)
    predicted <- as.numeric(predicted)
    n <- length(actual)
    if (length(predicted) != n) stop("\nactual and predicted must be of the same length.")
    if (any(is.na(actual))) {
        exc <- which(is.na(actual))
        actual <- actual[-exc]
        predicted <- predicted[-exc]
    }
    n <- length(actual)
    if (n == 0) return(NA)
    metric <- sum((predicted - actual)/abs(actual))/n
    return(metric)
}

#' @export
#' @rdname metrics
mslre <- function(actual, predicted)
{
    actual <- as.numeric(actual)
    predicted <- as.numeric(predicted)
    n <- length(actual)
    if (length(predicted) != n) stop("\nactual and predicted must be of the same length.")
    if (any(is.na(actual))) {
        exc <- which(is.na(actual))
        actual <- actual[-exc]
        predicted <- predicted[-exc]
    }
    n <- length(actual)
    if (n == 0) return(NA)
    
    ap <- actual/predicted
    if (any((ap) <= 0)) {
        if (all(ap <= 0)) {
            return(NA)
        } else {
            ap <- ap[which(ap > 0)]
        }
        n <- length(ap)
    }
    metric <- sum(log(ap)^2)/n
    return(metric)
}

#' @export
#' @rdname metrics
mase <- function(actual, predicted, original_series = NULL, frequency = 1)
{
    actual <- as.numeric(actual)
    predicted <- as.numeric(predicted)
    original_series <- as.numeric(original_series)
    if (any(is.na(actual))) {
        exc <- which(is.na(actual))
        actual <- actual[-exc]
        predicted <- predicted[-exc]
    }
    if (length(actual) == 0) return(NA)
    if (any(is.na(original_series))) {
        original_series <- original_series[-which(is.na(original_series))]
    }
    if (length(original_series) == 0) return(NA)
    n <- length(actual)
    if (length(predicted) != n) stop("\nactual and predicted must be of the same length.")
    if (is.null(original_series)) stop("\noriginal_series cannot be NULL for mase calculation.")
    N <- length(original_series)
    if (N < frequency) stop("\nlength of original series < frequency.")
    scaling_factor <- mean(abs(original_series[(frequency + 1):N] - original_series[1:(N - frequency)]))
    error <- mean(abs(predicted - actual))
    metric <- error/scaling_factor
    return(metric)
}

#' @export
#' @rdname metrics
mis <- function(actual, lower, upper, alpha)
{
    actual <- as.numeric(actual)
    lower <- as.numeric(lower)
    upper <- as.numeric(upper)
    alpha <- as.numeric(alpha)[1]
    if (any(is.na(actual))) {
        exc <- which(is.na(actual))
        actual <- actual[-exc]
        lower <- lower[-exc]
        upper <- upper[-exc]
    }
    n <- length(actual)
    if (n == 0) {
        metric <- as.numeric(NA)        
    } else {
        if (length(lower) != n | length(upper) != n ) stop("\nactual, lower and upper must be of the same length.")
        metric <- (sum(upper - lower) + 2/alpha * (sum((lower - actual) * (actual < lower)) + sum((actual - upper) * (actual > upper))))/n
    }
    return(metric)
}

#' @export
#' @rdname metrics
wape <- function(actual, predicted, weights)
{
    actual <- as.matrix(actual)
    predicted <- as.matrix(predicted)
    weights <- matrix(as.numeric(weights), ncol = 1)
    n <- NROW(actual)
    m <- NCOL(actual)
    if (m == 1) {
        if (NROW(weights) != n) stop("\nwhen actual is a univariate series, the weights must be of length equal to NROW(actual).")
    } else {
        if (NROW(weights) != m) stop("\nwhen actual is a multivariate series, the weights must be of length equal to NCOL(actuals).")
    }
    if (NROW(actual) != NROW(predicted)) stop("\nactual and predicted length not equal.")
    if (any(weights <= 0)) stop("\nweights must be strictly positive")
    # Normalizing weights to sum to 1
    weights <- weights/sum(weights)
    abs_error <- abs(predicted - actual)/actual
    weighted_metric <- mean(abs_error %*% weights, na.rm = T)
    return(weighted_metric)
}

#' @export
#' @rdname metrics
wslre <- function(actual, predicted, weights)
{
    actual <- as.matrix(actual)
    predicted <- as.matrix(predicted)
    weights <- matrix(as.numeric(weights), ncol = 1)
    n <- NROW(actual)
    m <- NCOL(actual)
    if (m == 1) {
        if (NROW(weights) != n) stop("\nwhen actual is a univariate series, the weights must be of length equal to NROW(actual).")
    } else {
        if (NROW(weights) != m) stop("\nwhen actual is a multivariate series, the weights must be of length equal to NCOL(actuals).")
    }
    if (NROW(actual) != NROW(predicted)) stop("\nactual and predicted length not equal.")
    if (any(weights <= 0)) stop("\nweights must be strictly positive")
    # Normalizing weights to sum to 1
    weights <- weights/sum(weights)
    weighted_metric <- mean(log(predicted/actual)^2 %*% weights, na.rm = T)
    return(weighted_metric)
}

#' @export
#' @rdname metrics
wse <- function(actual, predicted, weights)
{
    actual <- as.matrix(actual)
    predicted <- as.matrix(predicted)
    weights <- matrix(as.numeric(weights), ncol = 1)
    n <- NROW(actual)
    m <- NCOL(actual)
    if (m == 1) {
        if (NROW(weights) != n) stop("\nwhen actual is a univariate series, the weights must be of length equal to NROW(actual).")
    } else {
        if (NROW(weights) != m) stop("\nwhen actual is a multivariate series, the weights must be of length equal to NCOL(actuals).")
    }
    if (NROW(actual) != NROW(predicted)) stop("\nactual and predicted length not equal.")
    if (any(weights <= 0)) stop("\nweights must be strictly positive")
    # Normalizing weights to sum to 1
    weights <- weights/sum(weights)
    weighted_metric <- mean((predicted - actual)^2 %*% weights, na.rm = T)
    return(weighted_metric)
}

#' @export
#' @rdname metrics
pinball <- function(actual, distribution, alpha = 0.1){
    qtile <- apply(distribution, 2, quantile, q)
    loss <- rep(0, length(actual))
    loss[actual >= qtile] <- (1 - alpha/2) * (actual - qtile)
    loss[actual < qtile] <- (alpha/2) * (qtile - actual)
    metric <- 2 * mean(loss, na.rm = TRUE)
    return(metric)
}

#' @export
#' @rdname metrics
crps <- function(actual, distribution) 
{
    metric <- crps_sample(as.numeric(actual), t(distribution))
    metric <- mean(metric, na.rm = TRUE)
    return(metric)
}

rmape <- function(actual, predicted)
{
    actual <- as.numeric(actual)
    predicted <- as.numeric(predicted)
    abs_pe <- abs( (actual - predicted)/actual )
    lambda <- auto_lambda(abs_pe, lower = 1e-12, upper = 1.5)
    ape_t <- ((1 + abs_pe)^lambda - lambda)/lambda
    mape_t <- sum(abs(ape_t), na.rm = T)/length(ape_t)
    metric <- (lambda * (mape_t + 1))^(1/lambda) - 1
    return(metric)
}

#' @export
#' @rdname metrics
smape <- function(actual, predicted)
{
    actual <- as.numeric(actual)
    predicted <- as.numeric(predicted)
    metric <- mean( (abs(actual - predicted))/(abs(actual) + abs(predicted)), na.rm = T)
    return(metric)
}

#' @export
#' @rdname metrics
msis <- function(actual, lower, upper, original_series, frequency = 1, alpha)
{
    h <- length(upper)
    n <- length(original_series)
    a <- sum( (upper - lower) + (2/alpha) * (lower - actual) * as.integer(actual < lower) + (2/alpha) * (actual - upper) * as.integer(actual > upper), na.rm = T)
    b <- (1/(n - frequency)) * sum(abs(original_series[(frequency + 1):n] - original_series[1:(n - frequency)]), na.rm = T)
    metric <- (1/h) * (a/b)
    return(metric)
}