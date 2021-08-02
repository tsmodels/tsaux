auto_regressors <- function(y, frequency = 1, lambda = NULL, xreg = NULL, h = NULL, forc_dates = NULL, use = c("UC","tso"), args.tsmethod = list(), ...)
{
    if (!is.xts(y)) stop("\ny must be an object of class xts")
    if (NCOL(y) > 1) stop("\nonly univariate time series handled.")
    frequency <- as.integer(abs(frequency[1]))
    if (!is.null(h)) {
        h <- as.integer(abs(h))
        if (is.null(forc_dates)) {
            s <- sampling_frequency(index(y))
            forc_dates <- future_dates(tail(index(y),1), s, h)
        } else {
            if (length(forc_dates) != h) stop("\nforc_dates lenth must equal h")
        }
    }
    if (!is.null(lambda)) {
        y <- box_cox_transform(y, lambda = lambda)
    }
    use <- match.arg(use[1], c("UC","tso"))
    oxreg <- ma <- ar <- NULL
    if (use == "UC") {
        arglist <- list(...)
        if (is.null(arglist$model)) {
            if (frequency > 1) {
                arglist$model <- "?/equal/?"
            } else {
                arglist$model <- "?/none/?"
            }
        }
        if (frequency == 1) {
            arglist$periods <- NULL
        } else {
            if (is.null(arglist$periods)) {
                arglist$periods <- frequency
            }
        }
        if (is.null(arglist$outlier)) {
            arglist$outlier <- 3.5
        }
        arglist$y <- ts(as.numeric(y), frequency = frequency)
        arglist$u <- xreg
        arglist$h <- 0
        model <- do.call(UC, args = arglist, quote = TRUE)
        out <- model$comp
        cnames <- colnames(out)
        oxreg <- NULL
        if (any(grepl("AO",cnames))) {
            aoi <- which(grepl("AO",cnames))
            tmp <- xts(matrix(0, ncol = length(aoi), nrow = nrow(y)), index(y))
            colnames(tmp) <- cnames[aoi]
            for (i in 1:length(aoi)) tmp[which(out[,aoi[i]] != 0), i] <- 1
            if (!is.null(h)) {
                newtmp <- xts(matrix(0, ncol = ncol(tmp), nrow = h), forc_dates)
                colnames(newtmp) <- colnames(tmp)
                tmp <- rbind(tmp, newtmp)
            }
            oxreg <- cbind(oxreg, tmp)
        }
        if (any(grepl("LS",cnames))) {
            aoi <- which(grepl("LS",cnames))
            tmp <- xts(matrix(0, ncol = length(aoi), nrow = nrow(y)), index(y))
            colnames(tmp) <- cnames[aoi]
            for (i in 1:length(aoi)) tmp[which(out[,aoi[i]] != 0), i] <- 1
            if (!is.null(h)) {
                newtmp <- xts(matrix(1, ncol = ncol(tmp), nrow = h), forc_dates)
                colnames(newtmp) <- colnames(tmp)
                tmp <- rbind(tmp, newtmp)
            }
            oxreg <- cbind(oxreg, tmp)
        }
        coefnames <-  as.vector(rownames(coef(model)))
        if (any(grepl("AR\\([0-9]\\)",coefnames))) {
            ar <- length(which(grepl("AR\\([0-9]\\)",coefnames)))
        }
        if (any(grepl("MA\\([0-9]\\)",coefnames))) {
            ma <- length(which(grepl("MA\\([0-9]\\)",coefnames)))
        }
    } else {
        if (frequency == 1) seasonal <- FALSE else seasonal <- TRUE
        yts <- ts(as.numeric(y), frequency = frequency)
        args.tsmethod$seasonal <- seasonal
        model <- tso(yts, args.tsmethod = args.tsmethod, xreg = coredata(xreg), ...)
        if (!is.null(model$fit$xreg)) {
            oxreg <- xts(model$fit$xreg, index(y))
            cnames <- colnames(model$fit$xreg)
            colnames(oxreg) <- cnames
            if (any(cnames == "drift")) {
                oxreg <- oxreg[,-which(cnames == "drift")]
                if (NCOL(oxreg) == 0) return(NULL)
                cnames <- colnames(oxreg)
            }
            if (!is.null(h)) {
                oxreg <- rbind(oxreg, xts(matrix(as.numeric(NA), ncol = ncol(oxreg), nrow = h), forc_dates))
                if (any(grepl("AO",cnames))) {
                    ix <- which(grepl("AO",cnames))
                    oxreg[,ix] <- na.fill(oxreg[,ix], fill = 0)
                }
                if (any(grepl("LS",cnames))) {
                    ix <- which(grepl("LS",cnames))
                    oxreg[,ix] <- na.fill(oxreg[,ix], fill = 1)
                }
                if (any(grepl("TC",cnames))) {
                    ix <- which(grepl("TC",cnames))
                    for (j in 1:length(ix)) {
                        tmp <- oxreg[,ix[j]]
                        m <- which(tmp == 1)
                        delta <- as.numeric(tmp[m + 1])
                        z <- rep(0, nrow(oxreg))
                        z[m] <- 1
                        z <- filter(z, filter = delta, method = "recursive")
                        z <- as.numeric(z)
                        oxreg[,ix[j]] <- z
                    }
                }
            }
        } else {
            oxreg <- NULL
        }
        coefnames <- names(coef(model$fit))
        if (any(grepl("^ar[0-9]",coefnames))) {
            ar <- length(which(grepl("^ar[0-9]",coefnames)))
        }
        if (any(grepl("^ma[0-9]",coefnames))) {
            ma <- length(which(grepl("^ma[0-9]",coefnames)))
        }
    }
    return(list(xreg = oxreg, ar = ar, ma = ma))
}
