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
