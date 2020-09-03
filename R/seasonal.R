seasonality_test <- function(x, frequency = NULL){
    # Based on the QS test (Gomez  and  Maravall  1996)
    if (is.null(frequency)) {
        if (!is.xts(x)) stop("\nx must be an xts or zoo vector when frequency is NULL.")
        period <- na.omit(sampling_sequence(sampling_frequency(x)))
        period <- period[which(period > 1)]
    } else{
        period <- frequency
    }
    N <- NROW(x)
    # Used for determining whether the time series is seasonal
    tcrit <- 1.645
    if (N < (3 * period)) {
        test_seasonal <- FALSE
    } else {
        xacf <- acf(as.numeric(x), plot = FALSE)$acf[-1, 1, 1]
        clim <- tcrit/sqrt(N) * sqrt(cumsum(c(1, 2 * xacf^2)))
        test_seasonal <- abs(xacf[period]) > clim[period]
        if (is.na(test_seasonal)) {
            test_seasonal <- FALSE
        }
    }
    return(test_seasonal)
}

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
