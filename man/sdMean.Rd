% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/start2.R
\name{sdMean}
\alias{sdMean}
\title{Estimates for standard deviations and CI
for weighted means by observations from the mixture.}
\usage{
sdMean(x, p, comp = 1:ncol(p), means = FALSE, CI = FALSE,
  alpha = 0.05)
}
\arguments{
\item{x}{numeric vector with the observed sample or a
\code{wtsamp} object.}

\item{p}{matrix (or data frame) of mixing probabilities
with rows corresponding to subjects
and columns coresponding to the mixture components.}

\item{comp}{a numeric vector with numbers of components
for which the standard deviations are estimated.}

\item{means}{logical, if \code{TRUE} then the estimates for
components' means are included in the function value.}

\item{CI}{logical, if \code{TRUE} then confidence bounds for
components' means are inculded in the function value.}

\item{alpha}{confidense level for the confidence interval.}
}
\value{
if \code{CI & means =FALSE} the function returns a vector
of the estimated standard deviations
with NA for the components which were not estimated.
Else a data frame is returned in which there can be variables:
\code{sd} are standard deviations of estimates;
\code{means} are the estimates of means;
\code{lower} and \code{upper} are lower and upper bounds
of the confidence intervals for means.
}
\description{
\code{sdMean} calculates estimates of standard deviatitions
 and confidence intervals
 for  weighted means with minimax weights by observations
 from the mixture with varying concentrations.
}
\details{
If \code{CI=TRUE} then the function calculates
confidence intervals for the components' means
with covering probability \code{1-alpha}.

If \code{x} is a vector then the weights for components' means and variances
are calculated as \code{lsweight(p)}. If \code{x} is a \code{wtsamp}
object than its own weights are used.
}
\examples{
set.seed(3)
M<-3 # number of mixture components
p <- genunifp(1000,M) # create mixing probabilities
m<-c(0,1,2)      # true means of components
sd<-c(1,1,0.5)   # true sd of components
x<-genormixt(p,m,sd) # sample generation
# Calculate sd only:
sdMean(x,p)
# the same:
sdMean(wtsamp(x,indiv=lsweight(p)),p)
# Calculate confidence intervals:
sdMean(x,p,means=TRUE,CI=TRUE)
# Plot confidence intervals:
CI<-sdMean(x,p,means=TRUE,CI=TRUE)
library(plotrix)
plotCI(1:M,CI$means,ui=CI$upper,li=CI$lower,
       xlab=" ",ylab="means",xaxt="n")
axis(1,at=1:M,labels=row.names(CI))
}
\seealso{
Maiboroda R.  and  Kubaichuk O.
Asymptotic normality of improved weighted empirical distribution functions.
Theor. Probability and Math. Statist. 69 (2004), 95-102
}
