#' Forecast Performance Metrics
#'
#' Functions to calculate a number of performance metrics.
#'
#' @details
#' The following performance metrics are implemented:
#'
#' \describe{
#'   \item{\emph{Mean Absolute Percentage Error (MAPE)}}{Measures the average percentage deviation of predictions from actual values.
#'     \deqn{ MAPE = \frac{1}{n} \sum_{t=1}^{n} \left| \frac{A_t - P_t}{A_t} \right| }
#'     where \eqn{A_t} is the actual value and \eqn{P_t} is the predicted value.}
#'
#'   \item{\emph{Rescaled Mean Absolute Percentage Error (RMAPE)}}{A transformation of MAPE using a Box-Cox transformation for scale invariance (Swanson et al.).}
#'
#'   \item{\emph{Symmetric Mean Absolute Percentage Error (SMAPE)}}{An alternative to MAPE that symmetrizes the denominator.
#'     \deqn{ SMAPE = \frac{2}{n} \sum_{t=1}^{n} \frac{|A_t - P_t|}{|A_t| + |P_t|} }}
#'
#'   \item{\emph{Mean Absolute Scaled Error (MASE)}}{Compares the absolute error to the mean absolute error of a naive seasonal forecast.
#'     \deqn{ MASE = \frac{\frac{1}{n} \sum_{t=1}^{n} |P_t - A_t|}{\frac{1}{N-s} \sum_{t=s+1}^{N} |A_t - A_{t-s}|} }
#'     where \eqn{s} is the seasonal period.}
#'
#'   \item{\emph{Mean Squared Logarithmic Relative Error (MSLRE)}}{Measures squared log relative errors to penalize large deviations.
#'     \deqn{ MSLRE = \frac{1}{n} \sum_{t=1}^{n} \left( \log(1 + A_t) - \log(1 + P_t) \right)^2 }}
#'
#'   \item{\emph{Mean Interval Score (MIS)}}{Evaluates the accuracy of prediction intervals.
#'     \deqn{ MIS = \frac{1}{n} \sum_{t=1}^{n} (U_t - L_t) + \frac{2}{\alpha} [(L_t - A_t) I(A_t < L_t) + (A_t - U_t) I(A_t > U_t)] }
#'     where \eqn{L_t} and \eqn{U_t} are the lower and upper bounds of the interval.}
#'
#'   \item{\emph{Mean Scaled Interval Score (MSIS)}}{A scaled version of MIS, dividing by the mean absolute seasonal error.
#'     \deqn{ MSIS = \frac{1}{h} \sum_{t=1}^{h} \frac{(U_t - L_t) + \frac{2}{\alpha} [(L_t - A_t) I(A_t < L_t) + (A_t - U_t) I(A_t > U_t)]}{\frac{1}{N-s} \sum_{t=s+1}^{N} |A_t - A_{t-s}|} }}
#'
#'   \item{\emph{Bias}}{Measures systematic overestimation or underestimation.
#'     \deqn{ Bias = \frac{1}{n} \sum_{t=1}^{n} (P_t - A_t) }}
#'
#'   \item{\emph{Weighted Absolute Percentage Error (WAPE)}}{A weighted version of MAPE.
#'     \deqn{ WAPE = \sum_{t=1}^{n} \mathbf{w} \frac{|P_t - A_t|}{A_t} }
#'     where \eqn{\mathbf{w}} is the weight vector.}
#'
#'   \item{\emph{Weighted Squared Logarithmic Relative Error (WSLRE)}}{A weighted version of squared log relative errors.
#'     \deqn{ WSLRE = \sum_{t=1}^{n} \mathbf{w} (\log(P_t / A_t))^2 }}
#'
#'   \item{\emph{Weighted Squared Error (WSE)}}{A weighted version of squared errors.
#'     \deqn{ WSE = \sum_{t=1}^{n} \mathbf{w} \left( \frac{P_t}{A_t} \right)^2 }}
#'
#'   \item{\emph{Pinball Loss}}{A scoring rule used for quantile forecasts.
#'     \deqn{ \text{Pinball} = \frac{1}{n} \sum_{t=1}^{n} \left[ \tau (A_t - Q^\tau_t) I(A_t \geq Q^\tau_t) + (1 - \tau) (Q^\tau_t - A_t) I(A_t < Q^\tau_t) \right] }
#'     where \deqn{Q^\tau_t} is the predicted quantile at level \deqn{\tau}.}
#'
#'   \item{\emph{Continuous Ranked Probability Score (CRPS)}}{A measure of probabilistic forecast accuracy.
#'     \deqn{ CRPS = \frac{1}{n} \sum_{t=1}^{n} \int_{-\infty}^{\infty} (F_t(y) - I(y \geq A_t))^2 dy }
#'     where \eqn{F_t(y)} is the cumulative forecast distribution.}
#' }
#' @aliases mape rmape smape mase mslre mis msis bias wape wslre wse pinball crps
#' @param actual the actual values corresponding to the forecast period.
#' @param predicted the predicted values corresponding to the forecast period.
#' @param original_series the actual values corresponding to the training period.
#' @param frequency the seasonal frequency of the series used in the model.
#' @param lower the lower distributional forecast for the quantile corresponding to the coverage ratio alpha (i.e. alpha/2).
#' @param upper the upper distributional forecast for the quantile corresponding to the coverage ratio alpha (i.e. 1 - alpha/2).
#' @param alpha the distributional coverage.
#' @param distribution the forecast distribution (returned in the distribution slot of the prediction object). This is used in the continuous ranked probability score (crps) of Gneiting et al. (2005), and calculated using the function from the `scoringRules` package.
#' @param weights a vector of weights for generating weighted metrics. If the actual and predicted inputs are univariate, this should be equal to the length of the actual series and calculates a time-weighted average; otherwise, the weights should be of length equal to the number of series in a multivariate case, in which case a cross-sectional average is calculated.
#' @returns A numeric value.
#' @export
#' @rdname metrics
#' @note The RMAPE is the rescaled measure for MAPE based on the paper by Swanson et al.
#' @author Alexios Galanos
#' @references
#' \insertRef{Tofallis2015}{tsaux}\cr
#' \insertRef{Hyndman2006}{tsaux}\cr
#' \insertRef{Gneiting2005}{tsaux}\cr
#' \insertRef{Gneiting2007}{tsaux}\cr
#' \insertRef{Swanson2011}{tsaux}\cr
mape <- function(actual, predicted) {
    actual <- as.numeric(actual)
    predicted <- as.numeric(predicted)
    if (length(predicted) != length(actual)) stop("actual and predicted must be the same length.")
    valid <- actual != 0  # Avoid division by zero
    if (sum(valid) == 0) return(NA)  # If all actual values are zero, return NA
    metric <- mean(abs((actual[valid] - predicted[valid]) / actual[valid]), na.rm = TRUE)
    return(metric)
}

#' @export
#' @rdname metrics
bias <- function(actual, predicted) {
    actual <- as.numeric(actual)
    predicted <- as.numeric(predicted)
    if (length(predicted) != length(actual)) stop("actual and predicted must be of the same length.")
    valid <- !is.na(actual)
    if (sum(valid) == 0) return(NA)
    metric <- mean(predicted[valid] - actual[valid], na.rm = TRUE)
    return(metric)
}

#' @export
#' @rdname metrics
mslre <- function(actual, predicted) {
    actual <- as.numeric(actual)
    predicted <- as.numeric(predicted)
    valid <- actual > 0 & predicted > 0
    if (sum(valid) == 0) return(NA)  # If no valid data, return NA
    ap <- (log(1 + actual[valid]) - log(1 + predicted[valid]))^2
    metric <- mean(ap, na.rm = TRUE)
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
wape <- function(actual, predicted, weights) {
    actual <- as.matrix(actual)
    predicted <- as.matrix(predicted)
    weights <- matrix(as.numeric(weights), ncol = 1)
    n <- NROW(actual)
    m <- NCOL(actual)
    # Check weight dimensions
    if (m == 1) {
        if (NROW(weights) != n) stop("For univariate series, weights must have length equal to NROW(actual).")
    } else {
        if (NROW(weights) != m) stop("For multivariate series, weights must have length equal to NCOL(actual).")
    }
    if (NROW(actual) != NROW(predicted)) stop("actual and predicted must have the same length.")
    # Ensure all weights are strictly positive and non-zero sum
    if (any(weights <= 0)) stop("Weights must be strictly positive.")
    if (sum(weights) == 0) stop("Weights must not sum to zero.")
    # Normalize weights to sum to 1
    weights <- weights / sum(weights)
    # Handle zero values in `actual`
    valid <- which(actual != 0, arr.ind = TRUE)
    if (NROW(valid) > 0) {
        valid <- sort(unique(valid[,1]))
        abs_error <- abs(predicted[valid,] - actual[valid,]) / actual[valid,]
        weighted_metric <- sum(abs_error %*% weights, na.rm = TRUE)
        return(weighted_metric)
    } else {
        return(NA)
    }
}

#' @export
#' @rdname metrics
wslre <- function(actual, predicted, weights) {
    actual <- as.matrix(actual)
    predicted <- as.matrix(predicted)
    weights <- matrix(as.numeric(weights), ncol = 1)

    n <- NROW(actual)
    m <- NCOL(actual)

    # Check weight dimensions
    if (m == 1) {
        if (NROW(weights) != n) stop("For univariate series, weights must have length equal to NROW(actual).")
    } else {
        if (NROW(weights) != m) stop("For multivariate series, weights must have length equal to NCOL(actual).")
    }
    if (NROW(actual) != NROW(predicted)) stop("actual and predicted must have the same length.")

    # Ensure all weights are strictly positive and non-zero sum
    if (any(weights <= 0)) stop("Weights must be strictly positive.")
    if (sum(weights) == 0) stop("Weights must not sum to zero.")

    # Normalize weights to sum to 1
    weights <- weights / sum(weights)

    # Handle zero values in `actual` and `predicted` to prevent log errors
    valid <- actual > 0 & predicted > 0
    if (sum(valid) == 0) return(NA)  # If no valid data, return NA

    valid <- which(actual != 0, arr.ind = TRUE)
    if (NROW(valid) > 0) {
        valid <- sort(unique(valid[,1]))
        log_errors <- log(predicted[valid,] / actual[valid,])^2
        weighted_metric <- sum(log_errors %*% weights, na.rm = TRUE)
        return(weighted_metric)
    } else {
        return(NA)
    }
}

#' @export
#' @rdname metrics
wse <- function(actual, predicted, weights) {
    actual <- as.matrix(actual)
    predicted <- as.matrix(predicted)
    weights <- matrix(as.numeric(weights), ncol = 1)

    n <- NROW(actual)
    m <- NCOL(actual)

    # Check weight dimensions
    if (m == 1) {
        if (NROW(weights) != n) stop("For univariate series, weights must have length equal to NROW(actual).")
    } else {
        if (NROW(weights) != m) stop("For multivariate series, weights must have length equal to NCOL(actual).")
    }
    if (NROW(actual) != NROW(predicted)) stop("actual and predicted must have the same length.")

    # Ensure all weights are strictly positive and non-zero sum
    if (any(weights <= 0)) stop("Weights must be strictly positive.")
    if (sum(weights) == 0) stop("Weights must not sum to zero.")

    # Normalize weights to sum to 1
    weights <- weights / sum(weights)

    valid <- which(actual != 0, arr.ind = TRUE)
    if (NROW(valid) > 0) {
        valid <- sort(unique(valid[,1]))
        squared_errors <- (predicted[valid,] / actual[valid,])^2
        weighted_metric <- sum(squared_errors %*% weights, na.rm = TRUE)
        return(weighted_metric)
    } else {
        return(NA)
    }
}

#' @export
#' @rdname metrics
pinball <- function(actual, distribution, alpha = 0.1) {
    actual <- as.numeric(actual)
    distribution <- as.matrix(distribution)
    if (NCOL(distribution) != length(actual)) stop("actual and ncol(distribution) must be of the same length.")
    valid <- !is.na(actual)
    if (sum(valid) == 0) return(NA)
    qtile <- apply(distribution, 2, quantile, probs = alpha)
    qtile <- qtile[valid]
    actual <- actual[valid]
    loss <- ifelse(actual >= qtile, alpha * (actual - qtile),
                   (1 - alpha) * (qtile - actual))
    metric <- mean(loss, na.rm = TRUE)
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

#' @export
#' @rdname metrics
rmape <- function(actual, predicted)
{
    actual <- as.numeric(actual)
    predicted <- as.numeric(predicted)
    abs_pe <- abs( (actual - predicted)/actual)
    lambda <- box_cox_auto(abs_pe, lower = 1e-12, upper = 1.5)
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
    metric <- 2 * mean((abs(actual - predicted))/(abs(actual) + abs(predicted)), na.rm = T)
    return(metric)
}

#' @export
#' @rdname metrics
msis <- function(actual, lower, upper, original_series, frequency = 1, alpha) {
    h <- length(upper)
    n <- length(original_series)
    if (n <= frequency) return(NA)  # Prevent division by zero
    a <- sum((upper - lower) + (2 / alpha) * (lower - actual) * (actual < lower) +
                 (2 / alpha) * (actual - upper) * (actual > upper), na.rm = TRUE)
    b <- mean(abs(original_series[(frequency + 1):n] - original_series[1:(n - frequency)]), na.rm = TRUE)
    metric <- (1 / h) * (a / b)
    return(metric)
}
