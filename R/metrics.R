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

pinball <- function(actual, distribution, q = 0.95){
    qtile <- apply(distribution, 2, quantile, q)
    loss <- rep(0, length(actual))
    loss[actual >= qtile] <- q * (actual - qtile)
    loss[actual < qtile] <- (1 - q) * (qtile - actual)
    metric <- 2 * mean(loss, na.rm = TRUE)
    return(metric)
}

crps <- function(actual, distribution) 
{
    metric <- crps_sample(as.numeric(actual), t(distribution))
    metric <- mean(metric, na.rm = TRUE)
    return(metric)
}
