#' End of Week Date
#' 
#' Returns the last day of the week from a Date given a choice of week days.
#' 
#' Given a Date (such as 2019-01-02) and a day of 7, will return the Date for
#' the Sunday at or immediately after that. The week starting day is Monday
#' (1). A simple use case is when one wants to aggregate daily data to a
#' regular weekly sequence.
#' 
#' @param date a Date vector
#' @param day a value between 1 (Monday) and 7 (Sunday).
#' @return Date object
#' @export
#' @rdname calendar_eow
#' @author Alexios Galanos
calendar_eow <- function(date, day = 7)
{
    if (!is(date, "Date")) date <- as.Date(date)
    if (!day %in% 1:7) stop("\ninvalid day (must be 1-7)")
    date <- as.Date(date)
    w <- 7 - wday(date, week_start = day) + 1
    if (any(w == 7)) w[which(w == 7)] <- 0
    newdate <- date + w
    return(newdate)
}



#' End of Month Date
#' 
#' Returns the last day of the month from a Date within the month.
#' 
#' Given a Date (such as 2019-01-02), will return the last Date within that
#' year month.
#' 
#' @param date a Date vector
#' @return Date object
#' @export
#' @rdname calendar_eom
#' @author Alexios Galanos
calendar_eom <- function(date)
{
    if (!is(date, "Date")) date <- as.Date(date)
    # Add a month, then subtract a day:
    date.lt <- as.POSIXlt(date, format = "%Y-%m-%d", tz = tz(date))
    mon <- date.lt$mon + 2
    year <- date.lt$year
    # If month was December add a year
    year <- year + as.integer(mon == 13)
    mon[mon == 13] <- 1
    iso <- ISOdate(1900 + year, mon, 1, hour = 0, tz = tz(date))
    result <- as.POSIXct(iso) - 86400 # subtract one day
    result <- result + (as.POSIXlt(iso)$isdst - as.POSIXlt(result)$isdst)*3600
    result <- as.Date(result)
    return(result)
}



#' End of Quarter Date
#' 
#' Returns the last day of the quarter from a Date.
#' 
#' Given a date (such as 2019-01-02), will return the last date within that
#' year quarter.
#' 
#' @param date a Date vector
#' @return Date object
#' @export
#' @rdname calendar_eoq
#' @author Alexios Galanos
calendar_eoq <- function(date) {
    if (!is(date, "Date")) date <- as.Date(date)
    date.lt <- as.POSIXlt(date)
    # Date character string containing POSIXct date
    date.lt$mon <- (date.lt$mon %/% 3 + 1) * 3 %% 12
    date.lt$mday <- 1
    date.lt$year <- date.lt$year + as.integer(date.lt$mon == 0)
    result <- as.Date(date.lt) - 1
    return(result)
}



#' End of Year Date
#' 
#' Returns the last day of the year from a Date.
#' 
#' Given a date (such as 2019-01-02), will return the last date within that
#' year.
#' 
#' @param date a Date vector
#' @return Date object
#' @export
#' @rdname calendar_eoy
#' @author Alexios Galanos
calendar_eoy <- function(date) {
    if (!is(date, "Date")) date <- as.Date(date)
    result <- calendar_eom(as.Date(paste0(year(date),"-12-01")))
    return(result)
}
