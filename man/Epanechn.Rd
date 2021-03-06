% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/start2.R
\name{Epanechn}
\alias{Epanechn}
\title{Epanechnikov kernel calculation.}
\usage{
Epanechn(x)
}
\arguments{
\item{x}{a vector of points at which the kernel value is calculated.}
}
\value{
vector of the kernel values at given points.
}
\description{
\code{Epanechn} calculates the value of Epanechnikov kernel at a given point.
}
\examples{
curve(Epanechn,-2,2)
}
