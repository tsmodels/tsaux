####################################################
tslinear = function(y, trend = FALSE, seasonal = FALSE, xreg = NULL, frequency = 1, ...)
{
    if (NCOL(y) != 1) stop("\nonly univariate series allowed for y")
    good <- rep(1, NROW(y))
    has_missing <- FALSE
    if (any(is.na(y))) {
        good[which(is.na(y))] <- 0
        has_missing <- TRUE
    }
    y <- matrix(coredata(y), ncol = 1)
    colnames(y) <- "y"
    if (trend) {
        trend <- matrix(1:NROW(y), ncol = 1)
        colnames(trend) <- "trend"
    } else{
        trend <- NULL
    }
    if (seasonal) {
        if (frequency == 1) {
            stop("Non-seasonal time series")
        }
        seasonal <- matrix(0, ncol = frequency - 1, nrow = NROW(y))
        N = seq(1, nrow(y)) %% frequency
        N[which(N == 0)] <- frequency
        for (i in 1:(frequency - 1)) {
            seasonal[N == paste(i), i] <- 1
        }
        colnames(seasonal) <- paste0("s",1:(frequency - 1))
    } else {
        seasonal <- NULL
    }
    if (!is.null(xreg)) {
        if (NCOL(xreg) == 1) xreg <- matrix(xreg, ncol = 1)
        xreg <- check_xreg(xreg, index(y))
        if (is.null(colnames(xreg))) colnames(xreg) <- paste0("x",1:ncol(xreg))
    }
    data <- cbind(y, trend, seasonal, xreg)
    colnames(data) <- make.names(colnames(data))
    form <- as.formula(paste0("y~",paste0(colnames(data)[-1],collapse = "+")))
    fit <- lm(form, data = as.data.frame(data), na.action = na.exclude)
    fitted_value <- fitted(fit)
    if (has_missing) {
        p <- predict(fit, newdata = as.data.frame(data[which(good == 0),-1,drop = FALSE]))
        fitted_value[which(good == 0)] <- p
    }
    fit$data <- data
    responsevar <- deparse(form[[2]])
    fit$residuals <- residuals(fit)
    fit$x <- fit$residuals
    fit$x[!is.na(fit$x)] <- model.frame(fit)[, responsevar]
    fit$fitted.values <- fitted_value
    fit$method <- "Linear regression model"
    class(fit) <- c("tslinear", class(fit))
    return(fit)
}

check_xreg = function(xreg, valid.index)
{
    if (is.null(xreg)) return(xreg)
    n <- length(valid.index)
    if (NROW(xreg) != n) {
        stop("\nxreg does not have the same number of rows as y")
    }
    if (!all.equal(index(xreg),valid.index)) {
        stop("\nxreg time index does not match that of y")
    }
    if (any(is.na(xreg))) {
        stop("\nNAs found in xreg object")
    }
    if (is.null(colnames(xreg))) {
        colnames(xreg) <- paste0("x",1:ncol(xreg))
    }
    return(xreg)
}
