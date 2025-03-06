logit_transform <- function(x, lower = 0, upper = 1.0, ...) {
    if (!is.numeric(lower) || !is.numeric(upper)) {
        stop("\nlower and upper must be numeric values.")
    }
    if (length(lower) != 1 || length(upper) != 1) {
        stop("\nlower and upper must be scalars (single values).")
    }
    if (lower >= upper) {
        stop("\nlower must be strictly less than upper.")
    }
    if (any(x > upper, na.rm = TRUE) | any(x < lower, na.rm = TRUE))
        stop("\nrange of x is outside of lower and upper")
    return(log((x - lower) / (upper - x)))
}

logit_inverse <- function(x, lower = 0, upper = 1.0, ...) {
    return((upper - lower)/(1 + exp(-x)) + lower )
}

#' The logit transformation
#'
#' @description The logit transformation as an alternative to the Box Cox for bounded
#' outcomes.
#' @name logit
#' @aliases logit
#' @param lower lower bound of the variable.
#' @param upper upper bound of the variable.
#' @param ... not currently used.
#' @returns A list with the transform and inverse functions.
#' @export
#' @rdname logit
#' @author Alexios Galanos
logit <- function(lower = 0, upper = 1.0, ...)
{
    f <- function(y, ...){
        logit_transform(y, lower, upper)
    }
    fi <- function(y, ...){
        logit_inverse(y, lower, upper)
    }
    return(list(transform = f, inverse = fi))
}




#' Box-Cox transform specification
#'
#' Creates a specification for the Box Cox transformation.
#'
#' The function returns a list of 2 functions called \dQuote{transform} and
#' \dQuote{inverse} which can be called with a data object and a frequency to
#' calculate the transformed values. It is meant to be used in the transform
#' argument of the model specifications in the ts universe of models. The
#' auto_lambda function uses the method of Guerrero(1993).
#' @name box_cox
#' @aliases box_cox
#' @param lambda the power parameters. If NA then it will automatically
#' calculate the optimal parameter using the method of Guerrero (for univariate
#' case) else for the multivariate case, the method of Velilla (1993) which
#' implemented in the \code{car} package of John Fox. This targets a
#' transformation to multivariate normality. If any of the inputs has a
#' frequency other than 1, then an stl decomposition is first applied and the
#' seasonal component removed prior to the estimation in order to avoid
#' confounding the estimation by seasonality. It is also possible to pass a
#' vector equal to the number of columns of the dataset (with numeric values
#' mixed with NAs which will calculate the univariate optimal lambda).
#' @param lower optional parameter lower bound for cases when it is calculated.
#' @param upper optional parameter upper bound for cases when it is calculated.
#' @param multivariate flag for the multivariate case. If lambda is a single
#' parameter, then that is applied to all series (including NA which results in
#' the multivariate transformation described above).
#' @param ... not currently used.
#' @returns A list with the transform and inverse functions.
#' @note
#' The returned transform function will take additional argument \dQuote{frequency}
#' which determines whether a series is seasonal or not. When estimating lambda (when
#' setting this to NA), a series with frequency > 1 will first be de-seasonalized
#' using an STL decomposition.
#' @export
#' @rdname box_cox
#' @author Alexios Galanos for the BoxCox function.\cr John Fox for the
#' powerTransform function used in the multivariate case.
#' @references
#' \insertRef{Box1964}{tsaux}\cr
#' \insertRef{Velilla1993}{tsaux}\cr
#' \insertRef{Guerrero1993}{tsaux}\cr
#' @examples
#'
#' y = cumprod(c(1, 1 + rnorm(100,0.01, 0.005)))
#' B = box_cox(lambda = NA)
#' yt = B$transform(y, frequency = 1)
#' lambda = attr(yt,"lambda")
#' ye = B$inverse(yt, lambda)
#'
box_cox <- function(lambda = NA, lower = 0, upper = 1.5, multivariate = FALSE, ...)
{
    .lambda <- lambda
    .lower <- lower
    .upper <- upper
    if (multivariate) {
        f <- function(y, lambda = .lambda, frequency = 1, ...){
            n <- NCOL(y)
            if (length(frequency) == 1) {
                frequency <- rep(frequency, n)
            } else if (length(frequency) != n ) {
                stop("\nlength of frequency no equal ncol y")
            }
            ixs <- which(frequency > 1)
            # seasonally adjust series
            ynew <- y
            if (length(ixs) > 0) {
                for (i in ixs) {
                    tmp <- ts(as.numeric(y[,i]), frequency = frequency[i])
                    tmp <- stlplus(tmp, s.window = "periodic", n.p = frequency[i])$data$trend
                    ynew[,i] <- as.numeric(tmp)
                }
            }
            if (length(lambda) ==  1) {
                if (is.na(lambda)) {
                    lambda <- powerTransform(coredata(ynew))$lambda
                    # constrain to lower lambda
                    if (any(lambda < .lower)) {
                        lambda[which(lambda < .lower)] <- .lower
                    }
                    if (any(lambda > .upper)) {
                        lambda[which(lambda > .upper)] <- .upper
                    }
                    out <- do.call(cbind, lapply(1:ncol(y), function(i){ box_cox_transform(y[,i], lambda = lambda[i]) }))
                    colnames(out) <- colnames(y)
                    attr(out, "lambda") <- lambda
                    return(out)
                } else {
                    lambda <- rep(lambda, n)
                    out <- do.call(cbind, lapply(1:ncol(y), function(i){ box_cox_transform(y[,i], lambda = lambda[i]) }))
                    colnames(out) <- colnames(y)
                    attr(out, "lambda") <- lambda
                    return(out)
                }
            } else {
                if (length(lambda) == n) {
                    if (any(is.na(lambda))) {
                        use <- which(is.na(lambda))
                        for (i in use) lambda[i] <- box_cox_auto(y[,use], lower = .lower[1], upper = .upper[1], frequency = frequency[1])
                        out <- do.call(cbind, lapply(1:ncol(y), function(i){ box_cox_transform(y[,i], lambda = lambda[i]) }))
                        colnames(out) <- colnames(y)
                        attr(out, "lambda") <- lambda
                        return(out)
                    } else {
                        out <- do.call(cbind, lapply(1:ncol(y), function(i){ box_cox_transform(y[,i], lambda = lambda[i]) }))
                        colnames(out) <- colnames(y)
                        attr(out, "lambda") <- lambda
                        return(out)
                    }
                } else {
                    stop("\nlength of lambda should either be 1 or equal to number of columns of data")
                }
            }
        }
        fi <- function(y, lambda, ...){
            out <- do.call(cbind, lapply(1:ncol(y), function(i){ box_cox_inverse(y[,i], lambda = lambda[i]) }))
            colnames(out) <- colnames(y)
            return(out)
        }
    } else {
        f <- function(y, lambda = .lambda, frequency = 1, ...){
            if (NCOL(y) > 1) stop("\ny is multivariate. Call the function with multivariate = TRUE argument.")
            if (frequency > 1) {
                tmp <- ts(as.numeric(y), frequency = frequency)
                tmp <- stlplus(tmp, s.window = "periodic", n.p = frequency)$data$trend
            } else {
                tmp <- y
            }
            if (is.na(lambda)) lambda <- box_cox_auto(tmp, lower = .lower, upper = .upper, frequency = frequency)
            out <- box_cox_transform(y, lambda = lambda)
            attr(out, "lambda") <- lambda
            return(out)
        }
        fi <- function(y, lambda, ...){
            box_cox_inverse(y, lambda)
        }
    }
    return(list(transform = f, inverse = fi))
}

box_cox_inverse <- function(y, lambda = 1)
{
    if (lambda < 0) {
        stop("\nlambda must be strictly positive for the implemented Box-Cox transformation.")
    }
    if (lambda == 0) {
        y <- exp(y)
    } else {
        y <- (y * lambda + 1)^(1/lambda)
    }
    return(y)
}

box_cox_transform <- function(y, lambda = 1) {
    # Ensure y is strictly positive
    if (any(y <= 0, na.rm = TRUE)) {
        stop("\nBox-Cox transformation requires strictly positive values (y > 0).")
    }
    if (!is.numeric(lambda) || length(lambda) != 1) {
        stop("\nlambda must be a single numeric value.")
    }
    if (lambda < 0) {
        stop("\nlambda must be strictly positive for the implemented Box-Cox transformation.")
    }
    if (lambda == 0) {
        y <- log(y)
    } else {
        y <- (y^lambda - 1) / lambda
    }

    return(y)
}

box_cox_auto <- function(y, lower = 0, upper = 1, nonseasonal_length = 2, frequency = 1)
{
    if (any(y <= 0, na.rm = TRUE)) {
        lower <- max(lower, 0)
    }
    if (length(y) <= 2 * frequency) {
        stop("\nfrequency x 2 >= length(y). Not possible to calculate.")
    }
    return(guerrero(y, lower, upper, nonseasonal_length, frequency))
}


guerrero <- function(x, lower = 0, upper = 1, nonseasonal_length = 2, frequency = 1)
{
    return(optimize(guer_cv, c(lower, upper), x = x, nonseasonal_length = nonseasonal_length, frequency = frequency)$minimum)
}

guer_cv <- function(lambda, x, nonseasonal_length = 2, frequency)
{
    period <- round(max(nonseasonal_length, frequency))
    nobsf <- NROW(x)
    nyr <- floor(nobsf/period)
    nobst <- floor(nyr * period)
    x_mat <- matrix(x[(nobsf - nobst + 1):nobsf], period, nyr)
    x_mean <- apply(x_mat, 2, mean, na.rm = TRUE)
    x_sd <- apply(x_mat, 2, sd, na.rm = TRUE)
    x_rat <- x_sd/x_mean^(1 - lambda)
    return(sd(x_rat, na.rm = TRUE)/mean(x_rat, na.rm = TRUE))
}

#' The softplus logit transformation
#'
#' The softplus logit transformation is an alternative to the logit transform for bounded
#' outcomes with positive output.
#' @name softlogit
#' @aliases softlogit
#' @param lower lower bound of the variable.
#' @param upper upper bound of the variable.
#' @param ... not currently used.
#' @returns A list with the transform and inverse functions.
#' @export
#' @rdname softlogit
#' @author Alexios Galanos
#' @examples
#'
#' y = cumprod(c(1, 1 + rnorm(100,0.01, 0.005)))
#' B = softlogit(lower = 0,  upper = 15)
#' yt = B$transform(y)
#' ye = B$inverse(yt)
softlogit <- function(lower = 0, upper = 1.0, ...)
{
    f <- function(y, ...){
        softlogit_transform(y, lower, upper)
    }
    fi <- function(y, ...){
        softlogit_inverse(y, lower, upper)
    }
    return(list(transform = f, inverse = fi))
}

softlogit_transform <- function(x, lower = 0, upper = 1) {
    if (!is.numeric(lower) || !is.numeric(upper)) {
        stop("\nlower and upper must be numeric values.")
    }
    if (length(lower) != 1 || length(upper) != 1) {
        stop("\nlower and upper must be scalars (single values).")
    }
    if (lower >= upper) {
        stop("\nlower must be strictly less than upper.")
    }
    if (any(x <= lower | x >= upper)) {
        stop("x must be strictly between lower and upper bounds.")
    }
    logit_x <- log((x - lower) / (upper - x))
    softplus_x <- log(1 + exp(logit_x))
    return(softplus_x)
}

softlogit_inverse <- function(y, lower = 0, upper = 1) {
    x_recovered <- (lower + upper * (exp(y) - 1)) / exp(y)
    return(x_recovered)
}



#' The sigmoid transformation
#'
#' The sigmoid function is a smooth, S-shaped function that maps any real-valued input
#' into a bounded interval, typically  (0,1) . It is widely used in probability modeling,
#' logistic regression, and neural networks as an activation function.
#' @name sigmoid
#' @aliases sigmoid
#' @param lower lower bound of the variable.
#' @param upper upper bound of the variable.
#' @param ... not currently used.
#' @returns A list with the transform and inverse functions.
#' @export
#' @rdname sigmoid
#' @author Alexios Galanos
#' @examples
#'
#' y = cumprod(c(1, 1 + rnorm(100,0.01, 0.005)))
#' B = sigmoid()
#' yt = B$transform(y)
#' ye = B$inverse(yt)
sigmoid <- function(lower = 0.0, upper = 1.0, ...)
{
    f <- function(y, ...){
        sigmoid_transform(y, lower, upper)
    }
    fi <- function(y, ...){
        sigmoid_inverse(y, lower, upper)
    }
    return(list(transform = f, inverse = fi))
}

sigmoid_transform <- function(x, lower = 0, upper = 1) {
    if (!is.numeric(lower) || !is.numeric(upper)) {
        stop("\nlower and upper must be numeric values.")
    }
    if (length(lower) != 1 || length(upper) != 1) {
        stop("\nlower and upper must be scalars (single values).")
    }
    if (lower >= upper) {
        stop("\nlower must be strictly less than upper.")
    }
    lower + (upper - lower) / (1 + exp(-x))
}

sigmoid_inverse <- function(x, lower = 0, upper = 1) {
    if (any(x <= lower | x >= upper)) stop("Input must be strictly between lower and upper bounds")
    log((x - lower) / (upper - x))
}

#' General transformation function
#'
#' Includes the Box Cox, logit, softplus-logit and sigmoif transforms.
#' Returns a list of functions for the transform and its inverse.
#'
#' @param method valid methods are currently \dQuote{box-cox},
#' \dQuote{logit}, \dQuote{soft-logit} and \dQuote{sigmoid}.
#' @param lambda parameter in the Box Cox transformation.
#' @param lower lower bound for the transformations.
#' @param upper upper bound for the transformations.
#' @param \dots additional arguments taken by the transformations.
#' @returns A list with the transform and inverse functions.
#' @export
#' @rdname tstransform
#' @author Alexios Galanos
tstransform <- function(method = "box-cox", lambda = NULL, lower = 0, upper = 1,
                        ...)
{
    method = match.arg(method[1], c("box-cox","logit","softplus-logit","sigmoid"))
    f <- switch(method[1],
                "box-cox" = box_cox(lambda = lambda, lower = lower, upper = upper, ...),
                "logit" = logit(lower = lower, upper = upper, ...),
                "softplus-logit" = softlogit(lower = lower, upper = upper, ...),
                "sigmoid" = sigmoid(lower = lower, upper = upper, ...))
    return(f)
}
