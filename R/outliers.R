#' Automatic Detection of Outliers, Trends Breaks and Temporary Changes
#'
#' @description
#' A wrapper function for \emph{tso} from the \code{\link[tsoutliers]{tso}} package.
#' Takes as input a univariate xts object and returns a list with an xts object with any
#' identified outliers, trend breaks and/or temporary changes to be used as
#' regressors during estimation as well initial coefficients (see details).
#' @details
#' For generating future values of the identified outliers, the filter function
#' is used with additive outliers having a filter value of 0, trend changes a
#' value of 1, and temporary changes have value between 0 and 1. For the
#' sequential method, the routine first interpolates any missing values,
#' followed by an optional Box Cox transformation, and then elimination (and
#' identification) of any outliers during the first pass. The cleaned series is
#' then run through an stl filter (if any frequency is greater than 1) in order
#' to deseasonalize the data (with multiple seasonality supported), after which
#' the deseasonalized series is passed to the tso function where any additive
#' outliers (AO), temporary shifts (TC) or level shift (LS) are identified.
#' Additive outliers from this stage are added to any identified outliers from
#' the initial stage. For each regressor, initial parameter values are returned
#' together with the regressor matrix which should be passed to the estimation
#' routine. This is critically important since in the absence of good parameter
#' scaling, initial values are key to good convergence. Care should be taken
#' with regards to any automatic Box Cox parameter estimation. In the presence
#' of large outliers or level shifts, this is likely to be badly estimated
#' which is why we do not allow automatic calculation of this, but instead
#' place the burden on the user to decide what is a reasonable value (if any).
#' If a Box Cox transformation is used in the estimation routine, then it is
#' important to use the same lambda parameter in this function in order to get
#' sensible results. Again, avoid automatic Box Cox calculations throughout
#' when you suspect significant contamination of the series by outliers and
#' breaks. For the full method, the series is directly passed to the tso
#' function of the tsoutliers package. Finally, it should be noted that this
#' function is still experimental, and may change in the future.
#'
#' @param y a univariate xts object.
#' @param frequency the frequency of the time series. If the frequency is 1
#' then seasonal estimation will be turned off. Will also accept multiple
#' seasonal frequencies.
#' @param sampling the sampling frequency the series. If h>0 and forc_dates is
#' not provided, then this is required in order to generate future time indices
#' (valid values are days, months, hours, mins, secs etc).
#' @param lambda an optional Box Cox transformation parameter. The routines are
#' then run on the transformed dataset.
#' @param h an optional value for the forecast horizon (if planning to also use
#' for prediction).
#' @param forc_dates an optional vector of Date to be used for indexing the
#' series when h is not NULL. If this is not provided then the sampling
#' frequency of the series will be estimated in order to generate this.
#' @param stlm_opts additional arguments to the stlm function.
#' @param auto_arima_opts additional arguments to the auto.arima function in
#' the tso routine.
#' @param return_table whether to return a data.table instead with the
#' anomalies detected rather than an xts matrix with the pre-processed and
#' ready to use anomalies.
#' @param method whether to apply a sequential identification of anomalies
#' using STL decomposition in order to only pass the stationary residuals to
#' the tso function, else to pass the series directly to the tso package.
#' @param \dots any additional arguments passed to the tso functions (refer to
#' the documentation of the tsoutliers package).
#' @returns A list with an xts outlier matrix (if any where identified) as well
#' as a vector of initial parameter for use in the initialization of the
#' optimizer.
#' @export
#' @rdname auto_regressors
#' @author Alexios Galanos for this wrapper function.\cr Rob Hyndman for the
#' forecast package.\cr Javier Lopez-de-Lacalle for the tsoutliers package.
#' @examples
#'
#' library(xts)
#' set.seed(200)
#' y = cumprod(c(100,(1+rnorm(100,0.01, 0.04))))
#' y = xts(y, as.Date(1:101, origin = as.Date("2000-01-01")))
#' yclean = y
#' outlier1 = rep(0, 101)
#' outlier1[20] = 0.35
#' outlier2 = rep(0, 101)
#' outlier2[40] = 0.25
#' outlier2 = as.numeric(filter(outlier2, filter = 0.6, method = "recursive"))
#' y = y + y*xts(outlier1, index(y))
#' y = y + y*xts(outlier2, index(y))
#' # may need some tweaking of the tso options.
#' x = auto_regressors(y, frequency = 1, sampling = "days", h = 20,
#' check.rank = TRUE, discard.cval = 4)
#' head(x$xreg)
#' tail(x$xreg)
#' min(which(x$xreg[,1]==1))
#' min(which(x$xreg[,2]==1))
#' #plot(as.numeric(y), type = "l", ylab = "")
#' #lines(as.numeric(yclean) + (x$xreg %*% x$init)[1:101], col = 2)
#'
auto_regressors <- function(y, frequency = 1, lambda = NULL, forc_dates = NULL, sampling = NULL, h = 0,
                            stlm_opts = list(etsmodel = "AAN"),
                            auto_arima_opts = list(max.p = 1, max.q = 1, d = 1, allowdrift = FALSE),
                            return_table = FALSE, method = c("sequential","full"), ...)
{
    # missing values initial
    method <- match.arg(method[1], c("sequential","full"))
    type <- NULL
    if (length(y) <= (2 * min(frequency))) {
        frequency <- 1
        warning("\nlength of y less than 2xfrequency. Setting frequency to 1.")
    }
    if (h > 0) {
        if (is.null(sampling) & is.null(forc_dates)) {
            stop("\nh>0 but neither forc_dates or the sampling frequency of the data has been provided.")
        }
        if (is.null(forc_dates)) {
            forc_dates <- future_dates(tail(index(y),1), frequency = sampling, n = h)
        }
        newindex <- c(index(y), forc_dates)
    } else {
        newindex <- index(y)
    }
    if (any(is.na(y))) {
        y <- xts(na.interp(ts(y, frequency = frequency[1])), index(y))
    }
    # box cox
    if (!is.null(lambda)) {
        y <- box_cox_transform(y, lambda = lambda)
    }
    # outliers
    if (method == "full") {
        ynew <- y
        init_pars <- NULL
        tso_opts <- list()
        scaler <- max(ynew)
        tso_opts$y <- ts(as.numeric(ynew)/scaler, frequency = frequency)
        tso_opts$args.tsmethod <- auto_arima_opts
        tso_opts <- c(tso_opts, list(...))
        mod <- do.call(tso, args = tso_opts, quote = TRUE)
        xreg <- mod$fit$xreg
        aox <- NULL
        tcx <- NULL
        lsx <- NULL
        oxreg <- NULL
        if (!is.null(xreg)) {
            cnames <- colnames(xreg)
            if (any(grepl("^AO", cnames))) {
                init_pars <- c(init_pars, scaler * coef(mod$fit)[cnames[grepl("^AO",cnames)]])
                aox <- c(aox, as.numeric(gsub("^AO","",cnames[grepl("^AO", cnames)])))
            }
            if (!is.null(aox)) {
                aoxtmp <- xts(matrix(0, ncol = length(aox), nrow = length(newindex)), newindex)
                colnames(aoxtmp) <- paste0("AO",aox)
                for (i in 1:length(aox)) {
                    aoxtmp[aox[i],i] <- 1
                }
                oxreg <- cbind(oxreg, aoxtmp)
            }
            if (any(grepl("^TC", cnames))) {
                init_pars <- c(init_pars, scaler * coef(mod$fit)[cnames[grepl("^TC",cnames)]])
                tcx <- c(tcx, as.numeric(gsub("^TC","",cnames[grepl("^TC", cnames)])))
                xtmp <- xreg[,grepl("^TC", cnames), drop = FALSE]
                tcxtmp <- xts(matrix(0, ncol = length(tcx), nrow = length(newindex)), newindex)
                colnames(tcxtmp) <- paste0("TC",tcx)
                delta <- rep(0, length(tcx))
                for (i in 1:length(tcx)) {
                    delta[i] <- xtmp[tcx[i] + 1, i]
                    z <- rep(0, length(newindex))
                    z[tcx[i]] <- 1
                    z <- filter(z, filter = delta[i], method = "recursive")
                    z <- as.numeric(z)
                    tcxtmp[,i] <- as.numeric(z)
                }
                oxreg <- cbind(oxreg, tcxtmp)
            }
            if (any(grepl("^LS", cnames))) {
                init_pars <- c(init_pars, scaler * coef(mod$fit)[cnames[grepl("^LS",cnames)]])
                lsx <- c(lsx, as.numeric(gsub("^LS","",cnames[grepl("^LS", cnames)])))
                lsxtmp <- xts(matrix(0, ncol = length(lsx), nrow = length(newindex)), newindex)
                colnames(lsxtmp) <- paste0("LS",lsx)
                for (i in 1:length(lsx)) {
                    z <- rep(0, length(newindex))
                    z[lsx[i]:length(z)] <- 1
                    lsxtmp[,i] <- as.numeric(z)
                }
                oxreg <- cbind(oxreg, lsxtmp)
            }
            if (return_table) {
                dt <- as.numeric(sapply(1:ncol(oxreg), function(i) substr(colnames(oxreg)[i],3,nchar(colnames(oxreg)[i]))))
                dt <- index(y)[dt]
                rtable <- data.table(type = substr(colnames(oxreg),1,2), date = dt, filter = 0, coef = init_pars)
                rtable[type == "TC", filter := delta]
                rtable[type == "LS", filter := 1]
                return(rtable)
            } else {
                return(list(xreg = oxreg, init = init_pars))
            }
        } else {
            if (return_table) {
                return(NULL)
            } else {
                return(list(xreg = NULL, init = NULL))
            }
        }
    } else {
        out_proposal <- tsoutliers(ts(y, frequency = frequency[1]))
        ynew <- y
        init_pars <- NULL

        if (length(out_proposal$index) > 0) {
            init_pars <- c(init_pars, ynew[out_proposal$index] - out_proposal$replacements)
            names(init_pars) <- paste0("AO", out_proposal$index)
            ynew[out_proposal$index] <- out_proposal$replacements
            # identify the additive coefficient
        }
        # stlm
        if (any(frequency > 1)) {
            stlm_opts$y <- msts(as.numeric(ynew), ts.frequency = frequency[1], seasonal.periods = frequency)
            smod <- do.call(stlm, args = stlm_opts, quote = TRUE)
            ynew <- xts(seasadj(smod$stl), index(ynew))
        }
        tso_opts <- list()
        scaler <- max(ynew)
        tso_opts$y <- ts(as.numeric(ynew)/scaler, frequency = 1)
        tso_opts$args.tsmethod <- auto_arima_opts
        tso_opts <- c(tso_opts, list(...))
        mod <- do.call(tso, args = tso_opts, quote = TRUE)
        xreg <- mod$fit$xreg
        aox <- NULL
        tcx <- NULL
        lsx <- NULL
        oxreg <- NULL
        if (length(out_proposal$index) > 0) {
            aox <- c(aox, out_proposal$index)
            names(aox) <- paste0("AO",out_proposal$index)
        }
        if (!is.null(xreg)) {
            cnames <- colnames(xreg)
            if (any(grepl("^AO", cnames))) {
                init_pars <- c(init_pars, scaler * coef(mod$fit)[cnames[grepl("^AO",cnames)]])
                aox <- c(aox, as.numeric(gsub("^AO","",cnames[grepl("^AO", cnames)])))
            }
            if (!is.null(aox)) {
                aoxtmp <- xts(matrix(0, ncol = length(aox), nrow = length(newindex)), newindex)
                colnames(aoxtmp) <- paste0("AO",aox)
                for (i in 1:length(aox)) {
                    aoxtmp[aox[i],i] <- 1
                }
                oxreg <- cbind(oxreg, aoxtmp)
            }
            if (any(grepl("^TC", cnames))) {
                init_pars <- c(init_pars, scaler * coef(mod$fit)[cnames[grepl("^TC",cnames)]])
                tcx <- c(tcx, as.numeric(gsub("^TC","",cnames[grepl("^TC", cnames)])))
                xtmp <- xreg[,grepl("^TC", cnames), drop = FALSE]
                tcxtmp <- xts(matrix(0, ncol = length(tcx), nrow = length(newindex)), newindex)
                colnames(tcxtmp) <- paste0("TC",tcx)
                delta <- rep(0, length(tcx))
                for (i in 1:length(tcx)) {
                    delta[i] <- xtmp[tcx[i] + 1, i]
                    z <- rep(0, length(newindex))
                    z[tcx[i]] <- 1
                    z <- filter(z, filter = delta[i], method = "recursive")
                    z <- as.numeric(z)
                    tcxtmp[,i] <- as.numeric(z)
                }
                oxreg <- cbind(oxreg, tcxtmp)
            }
            if (any(grepl("^LS", cnames))) {
                init_pars <- c(init_pars, scaler * coef(mod$fit)[cnames[grepl("^LS",cnames)]])
                lsx <- c(lsx, as.numeric(gsub("^LS","",cnames[grepl("^LS", cnames)])))
                lsxtmp <- xts(matrix(0, ncol = length(lsx), nrow = length(newindex)), newindex)
                colnames(lsxtmp) <- paste0("LS",lsx)
                for (i in 1:length(lsx)) {
                    z <- rep(0, length(newindex))
                    z[lsx[i]:length(z)] <- 1
                    lsxtmp[,i] <- as.numeric(z)
                }
                oxreg <- cbind(oxreg, lsxtmp)
            }
        } else {
            if (!is.null(aox)) {
                aoxtmp <- xts(matrix(0, ncol = length(aox), nrow = length(newindex)), newindex)
                colnames(aoxtmp) <- paste0("AO",aox)
                for (i in 1:length(aox)) {
                    aoxtmp[aox[i],i] <- 1
                }
                oxreg <- cbind(aoxtmp, oxreg)
            }
        }
        # one final check for the AO
        if (!is.null(oxreg)) {
            cnames <- colnames(oxreg)
            if (any(duplicated(cnames))) {
                exc <- which(duplicated(cnames))
                oxreg <- oxreg[,-exc]
                init_pars <- init_pars[-exc]
            }
        }
        if (return_table) {
            dt <- as.numeric(sapply(1:ncol(oxreg), function(i) substr(colnames(oxreg)[i],3,nchar(colnames(oxreg)[i]))))
            dt <- index(y)[dt]
            rtable <- data.table(type = substr(colnames(oxreg),1,2), date = dt, filter = 0, coef = init_pars)
            rtable[type == "TC", filter := delta]
            rtable[type == "LS", filter := 1]
            return(rtable)
        } else {
            return(list(xreg = oxreg, init = init_pars))
        }
    }

}


#' Automatic Cleaning of Outliers and Temporary Changes
#'
#' A wrapper function for \emph{tso} from the tsoutliers packages. Takes as
#' input a univariate xts object and returns a series decomtaminated from
#' outliers and temporary changes.
#'
#' Calls the \code{\link[tsaux]{auto_regressors}} function to obtain the matrix of
#' regressors and coefficients which are then used to decomtaminate the series.
#' If lambda is not NULL, the series is first transformed to perform the
#' decontamination and then back transformed afterwards.
#'
#' @param y a univariate xts object.
#' @param frequency the frequency of the time series. If the frequency is 1
#' then seasonal estimation will be turned off. Will also accept multiple
#' seasonal frequencies.
#' @param lambda an optional Box Cox transformation parameter. The routines are
#' then run on the transformed dataset.
#' @param stlm_opts additional arguments to the stlm function.
#' @param auto_arima_opts additional arguments to the auto.arima function in
#' the tso routine.
#' @param types the types of anomalies to search and decontaminate series from.
#' Defaults to Additive outliers and temporary changes. Can be enhanced with
#' trend breaks but not suggested for the purpose of forecasting.
#' @param method whether to apply a sequential identification of anomalies
#' using STL decomposition in order to only pass the stationary residuals to
#' the tso function, else to pass the series directly to the tso package.
#' @param \dots any additional arguments passed to the tso functions (refer to
#' the documentation of the tsoutliers package).
#' @return A xts vector.
#' @export
#' @rdname auto_clean
#' @author Alexios Galanos for this wrapper function.\cr Rob Hyndman for the
#' forecast package.\cr Javier López-de-Lacalle for the tsoutliers package.
auto_clean <- function(y, frequency = 1, lambda = NULL, types = c("AO","TC"),
                       stlm_opts = list(etsmodel = "AAN"),
                       auto_arima_opts = list(max.p = 1, max.q = 1, d = 1, allowdrift = FALSE),
                       method = c("sequential", "full"), ...)
{
    mod <- auto_regressors(y, frequency = frequency, lambda = lambda, forc_dates = NULL, sampling = NULL, h = 0,
                                       stlm_opts = stlm_opts, auto_arima_opts = auto_arima_opts,
                                       return_table = FALSE, types = types, method = method, ...)
    if (!is.null(mod$xreg)) {
        x <- coredata(mod$xreg) %*% mod$init
        if (!is.null(lambda)) {
            ynew <- box_cox_transform(y, lambda)
            ynew <- coredata(ynew) - x
            ynew <- box_cox_inverse(ynew, lambda)
        } else {
            ynew <- coredata(y) - x
        }
        ynew <- xts(ynew, index(y))
        return(ynew)
    } else {
        return(y)
    }
}


#' Anomaly Creation
#'
#' Creates specific types of anomalies given a series.
#'
#' These functions allow the generation of anomalies and may be chained
#' together.
#'
#' @aliases additive_outlier temporary_change level_shift
#' @param y a univariate xts object or numeric series.
#' @param time the time index at which the anomaly takes place.
#' @param parameter the coefficient on the anomaly (the percent of the value of
#' y at the specified time index representing the anomaly).
#' @param alpha the AR(1) coefficient for the temporary change which determines
#' how quickly the effect decays.
#' @param add whether to contaminate the series (add the anomaly to the series)
#' else will return a matrix with the anomaly (without the effect of the
#' parameter).
#' @returns Either the contaminated series else a matrix of the anomaly.
#' @export
#' @rdname anomalies
#' @author Alexios Galanos for this wrapper function.\cr
additive_outlier <- function(y, time = 1, parameter = 0.5, add = TRUE)
{
    time <- as.integer(time[1])
    parameter <- as.numeric(parameter[1])
    n <- NROW(y)
    if (!(time %in% (1:n))) stop("\ntime not inside length of y")
    v <- rep(0, n)
    v[time] <- 1
    out <- filter(v, filter = 0, method = "recursive", init = 0)
    if (add) {
        coredata(y) <- coredata(y) + parameter * (coredata(y) * out)
        return(y)
    } else {
        out <- matrix(as.numeric(out), ncol = 1)
        colnames(out) <- paste0("AO_",time)
        if (is.xts(y)) out <- xts(out, index(y))
        return(out)
    }
}

#' @export
#' @rdname anomalies
temporary_change <- function(y, time = 1, parameter = 0.5, alpha = 0.7, add = TRUE)
{
    time <- as.integer(time[1])
    parameter <- as.numeric(parameter[1])
    alpha <- as.numeric(alpha[1])
    n <- NROW(y)
    if (!(time %in% (1:n))) stop("\ntime not inside length of y")
    v <- rep(0, n)
    v[time] <- 1
    out <- filter(v, filter = alpha, method = "recursive", init = 0)
    if (add) {
        coredata(y) <- coredata(y) + parameter * (coredata(y) * out)
        return(y)
    } else {
        out <- matrix(as.numeric(out), ncol = 1)
        colnames(out) <- paste0("TC_",time)
        if (is.xts(y)) out <- xts(out, index(y))
        return(out)
    }
}

#' @export
#' @rdname anomalies
level_shift <- function(y, time = 1, parameter = 0.5, add = TRUE)
{
    time <- as.integer(time[1])
    parameter <- as.numeric(parameter[1])
    n <- NROW(y)
    if (!(time %in% (1:n))) stop("\ntime not inside length of y")
    v <- rep(0, n)
    v[time] <- 1
    out <- filter(v, filter = 1, method = "recursive", init = 0)
    if (add) {
        coredata(y) <- coredata(y) + parameter * (coredata(y) * out)
        return(y)
    } else {
        out <- matrix(as.numeric(out), ncol = 1)
        colnames(out) <- paste0("AO_",time)
        if (is.xts(y)) out <- xts(out, index(y))
        return(out)
    }
}
