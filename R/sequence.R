#' Generate Regular Interval Future Dates
#' 
#' Generates regular interval future dates for use in forecast routine.
#' 
#' 
#' @param start a Date string for the start date.
#' @param frequency frequency of the interval (daily, weekly, monthly or
#' yearly).
#' @param n number of future periods to generate dates for.
#' @return A Date vector
#' @export
#' @rdname future_dates
#' @author Alexios Galanos
future_dates <- function(start, frequency, n = 1)
{
    if (frequency %in% c("days", "weeks", "months","years")) {
        switch(frequency,
               "days"   = as.Date(start) %m+% days(1:n),
               "weeks"  = as.Date(start) %m+% weeks(1:n),
               "months" = calendar_eom(as.Date(start) %m+% months(1:n)),
               "years"  = as.Date(start) %m+% years(1:n))
    } else if(grepl("secs|mins|hours|",frequency)) {
        # Add one extra point and eliminate first one
        seq(as.POSIXct(start), length.out = n + 1, by = frequency)[-1]
    } else{
        as.Date(start) + (1:n)
    }
}
