#' Generate Regular Interval Future Dates
#'
#' Generates regular interval future dates for use in forecast routine.
#'
#'
#' @param start a Date string for the start date.
#' @param frequency frequency of the interval (daily, weekly, monthly or
#' yearly).
#' @param n number of future periods to generate dates for.
#' @returns A Date vector
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
    } else if (grepl("secs|mins|hours|",frequency)) {
        # Add one extra point and eliminate first one
        seq(as.POSIXct(start), length.out = n + 1, by = frequency)[-1]
    } else{
        as.Date(start) + (1:n)
    }
}

#' Generate Train/Test Splits
#'
#' Generates train/test splits given a vector of dates and other options
#'
#' @param x a vector of timestamps (POSIXct) or dates (Date) in the dataset
#' @param start starting date (first estimation/train date)
#' @param test_length type of calendar period to split on
#' @param by every how many periods to split on
#' @param window_size the size of the training set (for moving window). If NULL
#' will use an expanding window.
#' @param calendar_end an optional function to use for the period ending split,
#' such as \code{\link{calendar_eow}}, applied to x. This should be greater in
#' frequency than the underlying frequency of x (i.e. do not use calendar_eow on
#' monthly indices). This overwrites the use of window_size.
#' @param complete_index whether to return the full indices for train and test else
#' just the start and end indices.
#' @param \dots any additional parameters passed to the calendar_end function. For
#' example, the \dQuote{day} argument when using the calendar_eow function.
#' @returns A list with each slot having the training dates and test dates
#' @note For months, quarters and years this will split into the end date of
#' these. For splitting into mins or hours, x must also have this resolution else
#' will throw an error. Additionally, the strict requirement of regularly spaced
#' time is required (no gaps).
#' @rdname time_splits
#' @author Alexios Galanos
#' @export
#'
time_splits <- function(x, start = x[1], test_length = 1, by = test_length,
                        window_size = NULL, calendar_end = NULL,
                        complete_index = TRUE, ...)
{
    if (missing(x)) stop("\nx cannot be missing")
    period_end <- NULL
    if (!is.null(calendar_end)) {
        d <- data.table(index = x)
        d[,period_end := calendar_end(index, ...)]
        d[, list(end = index[which.max(index)]), by = "period_end"]
    } else {
        n <- length(x)
        start_n <- which(x == start)
        end_n <- n - test_length + 1
        estimation_idx <- seq.int(start_n, end_n, by = by)
        if (any(!(estimation_idx %in% (start_n:end_n)))) {
            exc <- which(!(estimation_idx %in% (start_n:end_n)))
            estimation_idx <- estimation_idx[-exc]
        }
        out <- lapply(1:length(estimation_idx), function(i){
            train <- x[1:estimation_idx[i]]
            if (!is.null(window_size)) {
                train <- tail(train, pmin(window_size, length(train)))
            }
            test <- x[(estimation_idx[i] + 1):pmin(n, (estimation_idx[i] + test_length))]
            if (complete_index) {
                list(train = train, test = test)
            } else {
                list(train = c(train[1], tail(train,1)), test = c(test[1], tail(test, 1)))
            }
        })
    }
    return(out)
}
