mape <- function(actual, predicted)
{
    actual <- as.numeric(actual)
    predicted <- as.numeric(predicted)
    n <- length(actual)
    if (length(predicted) != n) stop("\nactual and predicted must be of the same length.")
    metric <- sum(abs((actual - predicted)/actual))/n
    return(metric)
}

bias <- function(actual, predicted)
{
    actual <- as.numeric(actual)
    predicted <- as.numeric(predicted)
    n <- length(actual)
    if (length(predicted) != n) stop("\nactual and predicted must be of the same length.")
    metric <- sum((predicted - actual)/abs(actual))/n
    return(metric)
}

mslre <- function(actual, predicted)
{
    actual <- as.numeric(actual)
    predicted <- as.numeric(predicted)
    n <- length(actual)
    if (length(predicted) != n) stop("\nactual and predicted must be of the same length.")
    metric <- sum(log(predicted/actual)^2)/n
    return(metric)
}

mase <- function(actual, predicted, original_series = NULL, frequency = 1)
{
    actual <- as.numeric(actual)
    predicted <- as.numeric(predicted)
    original_series <- as.numeric(original_series)
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

mis <- function(actual, lower, upper, alpha)
{
    actual <- as.numeric(actual)
    lower <- as.numeric(lower)
    upper <- as.numeric(upper)
    alpha <- as.numeric(alpha)[1]
    n <- length(actual)
    if (length(lower) != n | length(upper) != n ) stop("\nactual, lower and upper must be of the same length.")
    metric <- (sum(upper - lower) + 2/alpha * (sum((lower - actual) * (actual < lower)) + sum((actual - upper) * (actual > upper))))/n
    return(metric)
}

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
    weighted_metric <- mean(abs_error %*% weights)
    return(weighted_metric)
}

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
    weighted_metric <- mean(log(predicted/actual)^2 %*% weights)
    return(weighted_metric)
}

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
    weighted_metric <- mean((predicted - actual)^2 %*% weights)
    return(weighted_metric)
}

# pinball <- function(actual, lower, upper, q){
#
#     x <- ifelse(actual - q < 0, 1, 0)
#     score <- (actual - forecast) * (tau - indicator)
#     return(score)
# }
