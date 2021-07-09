auto_regressors <- function(y, frequency = 1, xreg = NULL, h = NULL, forc_dates = NULL, args.tsmethod = list(), ...)
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
    return(oxreg)
}
