find_frequency <- function(x)
{
    n <- length(x)
    x <- as.ts(x)
    # x should be the transformed variable else if the data is not variance stabilized
    # the frequency will be badly determined
    x <- residuals(tslinear(x, trend = TRUE))
    n.freq <- 500
    spec <- spec.ar(c(na.contiguous(x)), plot = FALSE, n.freq = n.freq)
    if (max(spec$spec) > 10) {
        period <- floor(1/spec$freq[which.max(spec$spec)] + 0.5)
        if (period == Inf) {
            j <- which(diff(spec$spec) > 0)
            if (length(j) > 0) {
                nextmax <- j[1] + which.max(spec$spec[(j[1] + 1):n.freq])
                if (nextmax < length(spec$freq)) {
                    period <- floor(1/spec$freq[nextmax] + 0.5)
                }
                else {
                    period <- 1L
                }
            }
            else {
                period <- 1L
            }
        }
    }
    else {
        period <- 1L
    }
    return(as.integer(period))
}


#' @title Infers the sampling frequency of a time series
#' @description
#' Given either a vector of time indices or an xts object will infer the
#' sampling frequency.
#' @param x either an xts object (or one which has an index attribute) else a
#' vector of class Date or POSIX based time index
#' @returns the sampling period (character).
#' @export
#' @rdname sampling_frequency
#' @examples
#' w <- sampling_frequency(seq(as.Date("2010-01-01"), as.Date("2011-01-01"), by="weeks"))
#' m <- sampling_frequency(seq(as.POSIXct("2010-01-01 12:00:00"),
#' as.POSIXct("2010-01-02 12:00:00"), by="15 mins"))
#'
sampling_frequency <- function(x)
{
    if (is(x, "Date") || length(grep("POSIX", class(x))) > 0) {
        dates <- x
    } else {
        dates <- index(x)
    }
    u <- min(diff(dates))
    count <- attr(u, 'units')
    if (count == 'days') {
        u <- round(u)
        daily   <- c(1, 2, 3)
        weekly  <- c(4, 5, 6, 7)
        monthly <- c(27, 28, 29, 30, 31, 32)
        yearly  <-  355:370
        if (u %in% daily) {
            period <- "days"
            attr(period,"date_class") <- "Date"
        } else if (u %in% weekly) {
            period <- "weeks"
            attr(period,"date_class") <- "Date"
        } else if (u %in% monthly) {
            period <- "months"
            attr(period,"date_class") <- "Date"
        } else if (u %in% yearly) {
            period <- "years"
            attr(period,"date_class") <- "Date"
        } else {
            period <- "unknown"
            attr(period,"date_class") <- "POSIXct"
        }
    } else if (count == "hours") {
        period <- paste0(u, " hours")
        attr(period,"date_class") <- "POSIXct"
    } else if (count == "mins") {
        period <- paste0(u, " mins")
        attr(period,"date_class") <- "POSIXct"
    } else if (count == "secs") {
        period <- paste0(u," secs")
        attr(period,"date_class") <- "POSIXct"
    } else {
        period <- "unknown"
        attr(period,"date_class") <- "POSIXct"
    }
    if (period == "unknown") warning("\ncould not determine sampling frequency")
    return(period)
}

#' @title Sampling frequency sequence
#' @description
#' Given a sampling period, the function will return the proportion of units of
#' that period in secs, mins, hours, days, weeks, months and years, but will
#' return NA for periods of higher frequency i.e. for a period of days it will
#' return NA for secs, mins and hours. The function serves as a helper for
#' seasonal periodicity calculations.
#' @param period the period returned by a call to the function
#' \code{\link[tsaux]{sampling_frequency}}.
#' @returns A named numeric vector.
#' @export
#' @rdname sampling_sequence
#' @author Alexios Galanos
#' @examples
#' w <- sampling_sequence(sampling_frequency(seq(as.Date("2010-01-01"),
#' as.Date("2011-01-01"), by="weeks")))
#' m <- sampling_sequence(sampling_frequency(seq(as.POSIXct("2010-01-01 12:00:00"),
#' as.POSIXct("2010-01-02 12:00:00"), by="15 mins")))
#'
sampling_sequence <- function(period)
{
    # [secs, mins, hrs, days, weeks, months, year]
    if (period == "days") {
        out <- c(NA, NA, NA, 1, 7, 30.4167, 365.25)
    }
    if (period == "weeks") {
        out <- c(NA, NA, NA, NA, 1, 4.34524, 52.1429)
    }
    if (period == "months") {
        out <- c(NA, NA, NA, NA, NA, 1, 12)
    }
    if (period == "years") {
        out <- c(NA, NA, NA, NA, NA, NA, 1)
    }
    if (grepl("hours", period)) {
        split <- strsplit(period," ")[[1]]
        if (length(split) > 1) {
            f <- as.numeric(split[1])
        } else{
            f <- 1
        }
        out <- c(NA, NA, 1, 24, 168, 730.001, 8760)/f
    }
    if (grepl("mins", period)) {
        split <- strsplit(period," ")[[1]]
        if (length(split) > 1) {
            f <- as.numeric(split[1])
        } else {
            f <- 1
        }
        out <- c(NA, 1, 60, 1440, 10080, 43800, 525600)/f
    }
    if (grepl("secs", period)) {
        split <- strsplit(period," ")[[1]]
        if (length(split) > 1) {
            f <- as.numeric(split[1])
        } else {
            f <- 1
        }
        out <- c(1, 60, 3600, 86400, 604800, 2.628e+6, 3.154e+7)/f
    }
    names(out) <- c("secs", "mins", "hours", "days", "weeks", "months", "years")
    return(out)
}
