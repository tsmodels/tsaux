\name{future_dates}
\alias{future_dates}
\title{
Generate Regular Interval Future Dates
}
\description{
Generates regular intervaled future dates for use in forecast routine.
}
\usage{
future_dates(start, frequency, n = 1)
}
\arguments{
\item{start}{
A Date string for the start date.
}
\item{frequency}{
The frequency of the interval (daily, weekly, monthly or yearly).
}
\item{n}{
The number of future periods to generate dates for.
}
}
\value{
A Date vector
}
\author{
Alexios Galanos
}
