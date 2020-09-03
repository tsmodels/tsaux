\name{calendar_eow}
\alias{calendar_eow}
\title{
End of Week Date
}
\description{
Returns the last day of the week from a Date given a choice of week days.
}
\usage{
calendar_eow(date, day = 7)
}
\arguments{
\item{date}{
A Date vector
}
\item{day}{
A value between 1 (Monday) and 7 (Sunday).
}
}
\details{
Given a Date (such as 2019-01-02) and a day of 7, will return the Date for the Sunday
at or immediately after that. The week starting day is Monday (1). A simple use case is
when one wants to aggregate daily data to a regular weekly sequence.
}
\value{
Date object
}
\author{
Alexios Galanos
}
