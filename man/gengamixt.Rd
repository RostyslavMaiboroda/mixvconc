% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/start2.R
\name{gengamixt}
\alias{gengamixt}
\title{Sample from
mixture of gamma distributions with varying concentrations.}
\usage{
gengamixt(p, shape, rate)
}
\arguments{
\item{p}{matrix (or data frame) of mixing probabilities.
with rows corresponding to subjects.
and columns coresponding to the mixture components.}

\item{shape}{vector of shape parameters for gamma distributions of
components.}

\item{rate}{vector of rate parameters for gamma distributions of
components.}
}
\value{
Vector  with the sample.
The sample size equals the rows number of p.
}
\description{
\code{gengamixt} generates a sample from the mixture
of gamma distributions
with mixing probabilities of components given by
the matrix p.
}
\examples{
  x <- genormixt(genunifp(10,2),c(2,1),c(0.5,0.5))
  plot(x)
}
