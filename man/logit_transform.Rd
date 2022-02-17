\name{logit_transform}
\alias{logit_transform}
\alias{logit}
\alias{logit_inverse}
\title{
The logit transformation
}
\description{
The logit transformation as an alternative to the Box Cox for bounded outcomes.
}
\usage{
logit(xmin = 0, xmax = 1.0)
logit_inverse(x, xmin = 0, xmax = 1.0)
logit_transform(x, xmin = 0, xmax = 1.0)
}
\arguments{
\item{x}{
A numeric vector (or object which can be coerced to such).
}
\item{xmin}{
The lower bound of the variable.
}
\item{xmax}{
The upper bound of the variable.
}
}
\details{
The logit function without a data argument returns a list of the 2 functions
(transform and inverse) with hard coded min and max bounds to be used in other
functions.
}
\value{
The transformed value or list of transformation functions.
}
\author{
Alexios Galanos
}