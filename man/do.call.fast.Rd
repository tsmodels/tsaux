\name{do.call.fast}
\alias{do.call.fast}
\title{
Execute a Function Call Faster
}
\description{
do.call.fast is just a faster version of do.call
}
\usage{
do.call.fast(what, args, quote = FALSE, envir = parent.frame())
}
\arguments{
\item{what}{
either a function or a non-empty character string naming the function to be called.
}
\item{args}{
a list of arguments to the function call. The names attribute of args gives the argument names.
}
\item{quote}{
a logical value indicating whether to quote the arguments.
}
\item{envir}{
an environment within which to evaluate the call. This will be most useful
if what is a character string and the arguments are symbols or quoted expressions.
}
}
\details{
See do.call function documentation
}
\value{
The result of the (evaluated) function call.
}
\references{
\url{https://stackoverflow.com/questions/11054208/lapply-and-do-call-running-very-slowly}
}
