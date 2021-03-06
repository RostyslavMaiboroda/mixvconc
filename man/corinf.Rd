% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/start2.R
\name{corinf}
\alias{corinf}
\title{Inf-correction of cummulative weights of a sample}
\usage{
corinf(xs)
}
\arguments{
\item{xs}{object from class \code{wtsamp} containing the sample
and weights for components' distributions.}
}
\value{
object from class \code{wtsamp} whose cummulative weights are
inf-corrected and individual weights are set to NULL.
}
\description{
\code{corinf} calculates cummulative weights for
the lower nondecreasing envelope of the
CDFs of mixture components (inf-corrected weights). The weights are
truncated by 0.
}
\details{
If cummulative weights are NULL they will be calculated
from the individual ones by \code{indiv2cum}.
}
\examples{
p <- genunifp(15,2) # create mixing probabilities
a <- lsweight(p) # calculate minimax weights
# create a weighted sample:
xs <- wtsamp(genormixt(p,c(0,1),c(1,1)),indiv=a)
F1 <- edfgen(xs,1) # minimax the estimate for 1-st component
xs_inf <- corinf(xs)
F1_inf <- edfgen(xs_inf,1) # inf-corrected estimate for 1-st component
F1(0)
F1_inf(0)
plot(F1,-3,3)
curve(F1_inf,col="red",,lty="dashed",add=TRUE)
}
\seealso{
Maiboroda R.  and  Kubaichuk O.
Asymptotic normality of improved weighted empirical distribution functions.
Theor. Probability and Math. Statist. 69 (2004), 95-102
}
