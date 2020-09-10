box_cox <- function(lambda = NA, lower = 0, upper = 1.5, multivariate = FALSE)
{
    .lambda <- lambda
    .lower <- lower
    .upper <- upper
    if (multivariate) {
        f <- function(y, lambda = .lambda, frequency = 1){
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
                    tmp <- tmp - stl(tmp, "periodic")$time.series[,"seasonal"]
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
        fi <- function(y, lambda){
            out <- do.call(cbind, lapply(1:ncol(y), function(i){ box_cox_inverse(y[,i], lambda = lambda[i]) }))
            colnames(out) <- colnames(y)
            return(out)
        }
    } else {
        f <- function(y, lambda = .lambda, frequency = 1){
            if (NCOL(y) > 1) stop("\ny is multivariate. Call the function with multivariate = TRUE argument.")
            if (is.na(lambda)) lambda <- box_cox_auto(y, lower = .lower, upper = .upper, frequency = frequency)
            out <- box_cox_transform(y, lambda = lambda)
            attr(out, "lambda") <- lambda
            return(out)
        }
        fi <- function(y, lambda){
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
    
    if (lambda == 1) {
        attr(y, "lambda") <- 1
        return(y)
    }
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
