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

calendar_eoy <- function(date) {
    if (!is(date, "Date")) date <- as.Date(date)
    result <- calendar_eom(as.Date(paste0(year(date),"-12-01")))
    return(result)
}
