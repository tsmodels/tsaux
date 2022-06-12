#' Box-Cox transform specification
#' 
#' Creates a specification for the Box Cox transformtion.
#' 
#' The function returns a list of 2 functions called \dQuote{transform} and
#' \dQuote{inverse} which can be called with a data object and a frequency to
#' calculate the transformed values. It is meant to be used in the transform
#' argument of the model specifications in the ts universe of models. The
#' auto_lambda function uses the method of Guerrero(1993).
#' 
#' @aliases box_cox auto_lambda
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
#' @param y an xts vector of strictly positive values.
#' @param frequency seasonal period of the data.
#' @param nonseasonal_length measurement periodicity of the data. The
#' maximum of this or the frequency argument is taken when calculating the
#' optimal lambda.
#' @param multivariate flag for the multivariate case. If lambda is a single
#' parameter, then that is applied to all series (including NA which results in
#' the multivariate transformation described above).
#' @param ... not currently used.
#' @return A list with the transform and invtransform functions. For the
#' auto_lambda function the optimal lambda is returned.
#' @export
#' @rdname box_cox
#' @author Alexios Galanos for the BoxCox function.\cr John Fox for the
#' powerTransform function used in the multivariate case.
#' @references Box, G. E. P. and Cox, D. R. (1964),\emph{An analysis of
#' transformations}. JRSS B \bold{26} 211--246.\cr Guerrero, V.M. (1993),
#' \emph{Time-series analysis supported by power transformations}. Journal of
#' Forecasting, \bold{12}, 37--48.\cr Velilla, S. (1993), \emph{A note on the
#' multivariate Box-Cox transformation to normality}. Statistics and
#' Probability Letters, \bold{17}, 259-263.
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
                    tmp <- tmp - stlplus(tmp, s.window = "periodic", n.p = frequency[i])$data$seasonal
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
            if (is.na(lambda)) lambda <- box_cox_auto(y, lower = .lower, upper = .upper, frequency = frequency)
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
        y[y > -1/lambda] <- NA
    }
    if (lambda == 0) {
        y <- exp(y)
    } else {
        xx <- y * lambda + 1
        y <- sign(xx) * abs(xx)^(1/lambda)
    }
    return(y)
}

box_cox_transform <- function(y, lambda = 1)
{
    if (lambda < 0) {
        y[y < 0] <- NA
    }
    if (lambda == 0) {
        y <- log(y)
    } else {
        y <- (sign(y) * abs(y)^lambda - 1)/lambda
    }
    return(y)
}

box_cox_auto <- function(y, lower = 0, upper = 1, nonseasonal_length = 2, frequency = 1)
{
    if (any(y <= 0, na.rm = TRUE)) {
        lower <- max(lower, 0)
    }
    if (length(y) <= 2 * frequency) {
        return(1)
    }
    return(guerrero(y, lower, upper, nonseasonal_length, frequency))
}

#' @export
#' @rdname box_cox
auto_lambda <- function(y, lower = 0, upper = 1, nonseasonal_length = 2, frequency = 1, ...)
{
    return( box_cox_auto(y = y, lower = lower, upper = upper, nonseasonal_length = nonseasonal_length, 
                 frequency = frequency) )
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



#' The logit transformation
#' 
#' The logit transformation as an alternative to the Box Cox for bounded
#' outcomes.
#' 
#' The logit function without a data argument returns a list of the 2 functions
#' (transform and inverse) with hard coded min and max bounds to be used in
#' other functions.
#' 
#' @aliases logit_transform logit logit_inverse
#' @param x a numeric vector (or object which can be coerced to such).
#' @param lower lower bound of the variable.
#' @param upper upper bound of the variable.
#' @param ... not currently used.
#' @return The transformed value or list of transformation functions.
#' @export
#' @rdname logit_transform
#' @author Alexios Galanos
logit_transform <- function(x, lower = 0, upper = 1.0, ...) {
    if (any(x > upper, na.rm = TRUE) | any(x < lower, na.rm = TRUE)) stop("\nrange of x is outside of lower and upper")
    return( -1.0 * log( ((upper - lower)/(x - lower)) - 1.0) )
}

logit_inverse <- function(x, lower = 0, upper = 1.0, ...) {
    return((upper - lower)/(1 + exp(-x)) + lower )
}

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



#' General transformation function
#' 
#' Includes the Box Cox and logit transforms. Returns a function for transform
#' and inverse.
#' 
#' Returns a list with functions for transform and its inverse.
#' 
#' @param method valid methods are currently \dQuote{box-cox} and
#' \dQuote{logit}.
#' @param lambda parameter in the Box Cox transformation.
#' @param lower lower bound for the transformations.
#' @param upper upper bound for the transformations.
#' @param \dots additional arguments taken by the transformations.
#' @return A list of functions
#' @export
#' @rdname tstransform
#' @author Alexios Galanos
tstransform <- function(method = "box-cox", lambda = NULL, lower = 0, upper = 1, 
                        ...)
{
    method = match.arg(method[1], c("box-cox","logit"))
    f <- switch(method[1],
                "box-cox" = box_cox(lambda = lambda, lower = lower, upper = upper, ...),
                "logit" = logit(lower = lower, upper = upper, ...))
    return(f)
}
