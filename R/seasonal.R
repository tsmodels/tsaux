seasonal_frequencies <- function(x, taper = 0, detrend = FALSE, demean = TRUE, plot = FALSE, ...)
{
    tmp <- ts(as.numeric(x), frequency = 1)
    sp <- spec.pgram(tmp, taper = taper, detrend = detrend, demean = demean, plot = FALSE)
    sp$spec <- 2 * sp$spec
    sp$spec[sp$freq == 0.5] <- sp$spec[sp$freq == 0.5]/2
    sol <- data.frame(freq = sp$freq, spec = sp$spec)
    sol <- sol[order(-sol$spec),]
    sol$time <- 1/sol$f
    if (plot) plot(y = sp$spec, x = sp$freq, type = "h")
    return(sol)
}

#' Simple Seasonality Test
#'
#' Checks for the presence of seasonality based on the QS test of Gomez and
#' Maravall (1996).
#'
#' Given the identified frequency of the xts vector (using the
#' \code{\link{sampling_frequency}}), the function checks for seasonality at
#' that frequency. The frequency can be overridden by directly supplying a
#' frequency argument, in which case y does not need to be a xts vector.
#'
#' @param x an (xts) vector (usually of a stationary series).
#' @param frequency overrides any frequency automatically identified in the
#' index of x.
#' @return Logical.
#' @export
#' @rdname seasonality_test
#' @author Alexios Galanos
#' @references
#' \insertRef{Gomez1995}{tsaux}
seasonality_test <- function(x, frequency = NULL){
    # Based on the QS test (Gomez  and  Maravall  1996)
    if (is.null(frequency)) {
        if (!is.xts(x)) stop("\nx must be an xts or zoo vector when frequency is NULL.")
        period <- na.omit(sampling_sequence(sampling_frequency(x)))
        period <- period[which(period > 1)][1]
    } else{
        period <- frequency
    }
    N <- NROW(x)
    # Used for determining whether the time series is seasonal
    tcrit <- 1.645
    if (N < (3 * period)) {
        test_seasonal <- FALSE
    } else {
        xacf <- acf(as.numeric(x), plot = FALSE, lag.max = 2 * period)$acf[-1, 1, 1]
        clim <- tcrit/sqrt(N) * sqrt(cumsum(c(1, 2 * xacf^2)))
        test_seasonal <- abs(xacf[period]) > clim[period]
        if (is.na(test_seasonal)) {
            test_seasonal <- FALSE
        }
    }
    return(test_seasonal)
}



#' Fourier terms for modelling seasonality
#'
#' Returns a matrix containing terms from a Fourier series, up to order K
#'
#'
#' @param dates a Date vector representing the length of the series for which
#' the fourier terms are required.
#' @param period frequency of the underlying series, if NULL will try to
#' infer it from the difference in the Date vector.
#' @param K maximum order of the Fourier terms.
#' @return A matrix of size N (length of dates) by 2*K.
#' @export
#' @rdname fourier_series
fourier_series <- function(dates, period = NULL, K = NULL) {
    N <- length(dates)
    if (is.character(period) | is.null(period)) {
        pd <- sampling_frequency(dates)
        period <- na.omit(sampling_sequence(pd))
        period <- period[which(period > 1)]
        if (length(period) == 0)
            return(NULL)
    }
    if (is.null(K)) {
        K <- pmax(1, period/2 - 1)
    }

    if (length(period) != length(K)) {
        stop("Need to provide K for each period")
    }
    if (any(K > period/2)) {
        stop("K is too big (greater than period/2)")
    }

    # Create arguments of trig functions and labels
    kseq <- numeric(0)
    name_seq <- character(0)
    for (k in seq_along(K)) {
        if (K[k] > 0) {
            kseq <- c(kseq, 1:K[k]/period[k])
            name_seq <- c(name_seq, paste(1:K[k], period[k], sep = "-"))
        }
    }
    # Remove duplicates for multicolinearity
    no_dupe <- !duplicated(kseq)
    kseq <- kseq[no_dupe]
    name_seq <- name_seq[no_dupe]
    # Grab number of columns of design mat
    bigK <- length(kseq)
    # Create design matrix
    design <- matrix(NA_real_, nrow = N, ncol = 2 * bigK)
    labels <- character(length = 2 * bigK)
    for (j in seq_along(kseq)) {
        design[, 2 * j - 1] <- sinpi(2 * kseq[j] * (1:N))
        design[, 2 * j] <- cospi(2 * kseq[j] * (1:N))
        labels[2 * j - 1] <- paste0("Sin", name_seq[j])
        labels[2 * j] <- paste0("Cos", name_seq[j])
    }
    colnames(design) <- labels
    design <- design[, !is.na(colSums(design)), drop = FALSE]
    attr(design, "K") <- K
    attr(design, "period") <- period
    return(design)
}



#' Seasonal Dummies
#'
#' Creates a matrix of seasonal dummies.
#'
#' Generates seasons-1 dummy variables.
#'
#' @param y optional data series.
#' @param n if y is missing, then the length of the series is required.
#' @param seasons number of seasons in a cycle.
#' @returns Either a matrix (if y is missing or y is not an xts vector) or an
#' xts matrix (when y is an xts vector).
#' @export
#' @rdname seasonal_dummies
#' @author Alexios Galanos
#' @examples
#'
#' head(seasonal_dummies(n=100, seasons=12))
#'
seasonal_dummies <- function(y = NULL, n = nrow(y), seasons = 12)
{
    if (is.null(y)) {
        if (is.null(n)) stop("\neither provide n or y")
    } else {
        n <- NROW(y)
    }
    dummy_matrix <- matrix(0, ncol = seasons - 1, nrow = n)
    dummy <- 1:(seasons - 1)
    for (i in 1:(seasons - 1)) {
        dummy_matrix[dummy == i, i] <- 1
    }
    colnames(dummy_matrix) <- paste0("s-",1:(seasons - 1))
    if (!is.null(y) & is.xts(y)) {
        dummy_matrix <- xts(dummy_matrix, index(y))
    }
    return(dummy_matrix)
}
