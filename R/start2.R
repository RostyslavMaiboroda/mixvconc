#__________________________________________________
# DATA GENERATORS
#__________________________________________________
#' Generates matrix of mixing probabilities.
#'
#' \code{genunifp} returns a matrix whose rows are independent
#' random vectors uniformly distributed over the M-simplex.
#'
#' @param n  number of subjects (rows).
#' @param M  number of mixing components (columns).
#' @return Matrix p of mixing probabilities with n rows and M columns.
#' The names of components (columns) are A,B,C...
#' The names of cases (rows) are 1,2,3,...
#'
#' @examples
#' p <- genunifp(10,2)
#' p
#'
#' p <- genunifp(1000,3)
#' plot(p[,3],p[,1],cex=0.2)
#' @seealso
#' Maiboroda R., Sugakova O. "Statistics of mixtures with varying concentrations 
#' with application to DNA microarray data analysis". 
#' Nonparametric statistics (2012) v.24:1, p. 201 - 215.
#' @export 
#__________________________________________________
genunifp <- function(n, M) {
    p <- matrix(stats::rexp(n * M), nrow = n,dimnames=list(1:n,LETTERS[1:M]))
    p <- t(apply(p, 1, function(x) x/sum(x)))
}
#__________________________________________________
#' Sample form mixture of normals with warying concentrations.
#'
#' \code{genormixt} generates a sample from the mixture
#' of normal distributions
#' with mixing probabilities of components given by
#' the matrix p.
#'
#' @param p  matrix (or data frame) of mixing probabilities
#'  with rows corresponding to subjects
#'  and columns coresponding to the mixture components.
#' @param mean  vector of components' means.
#' @param sd  vector of components' standard deviations.
#' @return   Vector  with the sample.
#' The sample size equals the rows number of p.
#' @examples
#'   x <- genormixt(genunifp(10,2),c(0,1),c(0.5,0.5))
#'   plot(x)
#' @seealso
#' Maiboroda R., Sugakova O. "Statistics of mixtures with varying concentrations 
#' with application to DNA microarray data analysis". 
#' Nonparametric statistics (2012) v.24:1, p. 201 - 215
#' @export 
# __________________________________________________
genormixt <- function(p, mean, sd) {
    nc <- 1:ncol(p)
    obj <- apply(p, 1, function(x) sample(nc, 1, prob = x))
    smp <- stats::rnorm(nrow(p), mean[obj], sd[obj])
}
#_______________________________________________
#' Sample from
#' mixture of gamma distributions with varying concentrations.
#'
#' \code{gengamixt} generates a sample from the mixture
#' of gamma distributions
#' with mixing probabilities of components given by
#' the matrix p.
#'
#' @param p  matrix (or data frame) of mixing probabilities.
#'  with rows corresponding to subjects.
#'  and columns coresponding to the mixture components.
#' @param shape  vector of shape parameters for gamma distributions of
#' components.
#' @param rate   vector of rate parameters for gamma distributions of
#' components.
#' @return   Vector  with the sample.
#' The sample size equals the rows number of p.
#' @examples
#'   x <- genormixt(genunifp(10,2),c(2,1),c(0.5,0.5))
#'   plot(x)
#' @export 
# __________________________________________________
gengamixt <- function(p, shape, rate) {
    nc <- 1:ncol(p)
    obj <- apply(p, 1, function(x) sample(nc, 1, prob = x))
    smp <- stats::rgamma(nrow(p), shape[obj], rate[obj])
}
# __________________________________________________
#' Calculates minimax weights for components' distributions estimation.
#'
#' \code{lsweight} returns a matrix of individual weights
#' which correspond to minimax unbiased estimates for CDFs of
#' mixture components.
#'
#' @param p  matrix (or data frame) of mixing probabilities
#'  with rows corresponding to subjects
#'  and columns coresponding to the mixture components.
#' @return  matrix (or data frame) of minimax weights of the same structure
#' as p
#' @examples
#' set.seed(3)
#' p <- genunifp(10,3)
#' a <- lsweight(p)
#' t(a)%*%p # the result is a unit matrix
#' @seealso
#' Maiboroda R., Sugakova O. "Statistics of mixtures with varying concentrations 
#' with application to DNA microarray data analysis". 
#' Nonparametric statistics (2012) v.24:1, p. 201 - 215.
#' @export 
#____________________________________________________
lsweight <- function(p) {
    p <- as.matrix(p)
    a <- as.data.frame(p %*% MASS::ginv(t(p) %*% p))
    names(a)<-colnames(p)
    a
}
# ______________________________________________
#' Constructor for class \code{wtsamp}
#'
#' \code{wtsamp} returns an object of S3 class \code{wtsamp}
#' containing a sorted sample and a set of individual
#' and/or cummulative weights representing distributions
#' of different components.
#'
#' @param x numeric vector containing the  sample values.
#' @param cumm matrix (or data frame) of cummulative weights
#' of components.
#' @param indiv matrix (or data frame) of individual weights
#' of components.
#' @return object of class \code{wtsamp} which contains the
#' following attributes:
#'
#' \describe{
#' \item{\code{xo}}{vector of sample values sorted
#' in the ascending order with \code{-Inf} as the first element.}
#' \item{\code{cumm}}{matrix of cummulative weigts reordered
#' at the same order as xo with 0 at the first row.}
#' \item{\code{indiv}}{matrix of individual weigts reordered
#' at the same order as xo.}
#' }
#'
#' set.seed(3)
#' p <- genunifp(10,2) # create mixing probabilities
#' a <- lsweight(p) # calculate minimax weights
#' # create a weighted sample:
#' xs <- wtsamp(genormixt(p,c(0,1),c(1,1)),indiv=a)
#' @seealso
#' Maiboroda R., Sugakova O. "Statistics of mixtures with varying concentrations 
#' with application to DNA microarray data analysis". 
#' Nonparametric statistics (2012) v.24:1, p. 201 - 215.
#' @export 
#_______________________________________________
wtsamp <- function(x, cumm = NULL, indiv = NULL) {
    o <- order(x)
    xo <- c(-Inf, x[o])
    if (!is.null(cumm))
        cumm <- rbind(rep(0, ncol(cumm)), cumm[o, ])
    if (!is.null(indiv))
        indiv <- indiv[o, ]
    structure(list(xo = xo, cumm = cumm, indiv = indiv), class = "wtsamp")
}
# ______________________________________________
#' Constructor for class \code{wtcens}
#'
#' \code{wtcens} returns an object of S3 class \code{wtcens}
#' containing the vector of sorted sample values, vector
#' of indicators of non-censoring
#' and a set of individual
#' and/or cummulative weights representing distributions
#' of different components.
#'
#' @param x numeric vector containing the sorted sample values.
#' @param delta logical vector of non-censoring indicators
#' @param cumm matrix (or data frame) of cummulative weights
#' of components.
#' @param indiv matrix (or data frame) of individual weights
#' of components.
#' @return object of class \code{wtcens}which contains the
#' following attributes:
#'
#' \describe{
#' \item{\code{xo}}{vector of sample values sorted
#' in the ascending order with \code{-Inf} as the first element.}
#' \item{\code{deltao}}{vector of non-censoring indicators reordered
#' at the same order as xo.}
#' \item{\code{cumm}}{matrix of cummulative weigts reordered
#' at the same order as xo with 0 at the first row.}
#' \item{\code{indiv}}{matrix of individual weigts reordered
#' at the same order as xo.}
#' }
#' @export 
#_______________________________________________
wtcens <- function(x, delta, cumm = NULL, indiv = NULL) {
    o <- order(x)
    xo <- c(-Inf, x[o])
    deltao <- delta[o]
    if (!is.null(cumm))
        cumm <- rbind(rep(0, ncol(cumm)), cumm[o, ])
    if (!is.null(indiv))
      indiv <- indiv[o, ]
    structure(list(xo = xo, deltao = deltao, cumm = cumm, indiv = indiv), class = "wtcens")
}
# __________________________________________________
#' Calculates cummulative weights.
#'
#' \code{indiv2cumm} calculates cummulative sums of
#' individual weights and put them into cummulative weights
#' attribute.
#'
#' @param xs a \code{wtsamp} or \code{wtcens} object representing
#' a sample with distributions of different components.
#'
#' @return  an object of the same class as \code{xs} with recalculated
#' cummulative weights.
#'
#' @examples
#'set.seed(3)
#' p <- genunifp(10,2) # create mixing probabilities
#' a <- lsweight(p) # calculate minimax weights
#' # create a weighted sample:
#' xs <- wtsamp(genormixt(p,c(0,1),c(1,1)),indiv=a)
#' ys <- indiv2cumm(xs) # create cummulative weights
#' xs1 <- cumm2indiv(ys) #xs1 is the same as xs
#' @export 
#____________________________________________________
indiv2cumm <- function(xs) {
    xs$cumm <- apply(xs$indiv, 2, cumsum)
    xs$cumm <- rbind(rep(0, ncol(xs$cumm)), xs$cumm)
    xs
}
# __________________________________________________
#' Calculates individual weights.
#'
#' \code{cumm2indiv} calculates differences with lag 1 from
#' cummulative weights and put them into individual weights
#' attribute.
#'
#' @param xs a \code{wtsamp} or \code{wtcens} object representing
#' a sample with distributions of different components.
#' @return  an object of the same class as \code{xs} with recalculated
#' individual weights.
#'
#' @examples
#' set.seed(3)
#' p <- genunifp(10,2) # create mixing probabilities
#' a <- lsweight(p) # calculate minimax weights
#' # create a weighted sample:
#' xs <- wtsamp(genormixt(p,c(0,1),c(1,1)),indiv=a)
#' ys <- indiv2cumm(xs) # create cummulative weights
#' xs1 <- cumm2indiv(ys) #xs1 is the same as xs
#' @export 
#____________________________________________________
cumm2indiv <- function(xs) {
    xs$indiv <- apply(xs$cumm, 2, diff)
    xs
}
# __________________________________________________
#' Generator of weighted empirical CDS.
#'
#' @param xs a \code{wtsamp} or \code{wtcens} object representing
#' a sample with distributions of different components.
#'
#' @param m  number of the component whose CDF is estimated.
#'
#' @return  a function with the call \code{f(t)} where \code{t} is
#' the vector of points at which the estimate is calculated.  \code{f(t)}
#' returns the vector of estimates.
#' @examples
#' set.seed(3)
#' p <- genunifp(1000,2) # create mixing probabilities
#' a <- lsweight(p) # calculate minimax weights
#' # create a weighted sample:
#' xs <- wtsamp(genormixt(p,c(0,1),c(1,1)),indiv=a)
#' F1<-edfgen(xs,1) # generate the estimate for 1-st component
#' F1(0)            # 0.5289866 approximately 0.5
#' plot(F1,-3,3) # plot the estimate (approx. standard normal CDF )
#' @export 
#____________________________________________________
edfgen <- function(xs, m) {
    if (is.null(xs$cumm))
        xs <- indiv2cumm(xs)
    a <- xs$cumm[, m]
    x <- xs$xo
    f <- function(t) {
        a[findInterval(t, x)]
    }
}
#' Sup-correction of cummulative weights of a sample
#'
#' \code{corsup} calculates cummulative weights for
#' the upper nondecreasing envelope of the
#' CDFs of mixture components (sup-corrected weights). The weights are
#' truncated by 1.
#'
#' @param xs object from class \code{wtsamp} containing the sample
#' and weights for components' distributions.
#'
#' @return object from class \code{wtsamp} whose cummulative weights are
#' sup-corrected and individual weights are set to NULL.
#'
#' @details If cummulative weights are NULL they will be calculated
#' from the individual ones by \code{indiv2cum}.
#'
#' @examples
#' p <- genunifp(15,2) # create mixing probabilities
#' a <- lsweight(p) # calculate minimax weights
#' # create a weighted sample:
#' xs <- wtsamp(genormixt(p,c(0,1),c(1,1)),indiv=a)
#' F1 <- edfgen(xs,1) # minimax the estimate for 1-st component
#' xs_sup <- corsup(xs)
#' F1_sup <- edfgen(xs_sup,1) # sup-corrected estimate for 1-st component
#' F1(0)
#' F1_sup(0)
#' plot(F1,-3,3)
#' curve(F1_sup,col="red",,lty="dashed",add=TRUE)
#' @seealso
#' Maiboroda R., Sugakova O. "Statistics of mixtures with varying concentrations 
#' with application to DNA microarray data analysis". 
#' Nonparametric statistics (2012) v.24:1, p. 201 - 215.
#' @export 
# ________________________________________________________
corsup <- function(xs) {
    if (is.null(xs$cumm))
        xs <- indiv2cumm(xs)
    xs$cumm <- pmin(apply(xs$cumm, 2, cummax), 1)
    xs$indiv <- NULL
    xs
}
#' Inf-correction of cummulative weights of a sample
#'
#' \code{corinf} calculates cummulative weights for
#' the lower nondecreasing envelope of the
#' CDFs of mixture components (inf-corrected weights). The weights are
#' truncated by 0.
#'
#' @param xs object from class \code{wtsamp} containing the sample
#' and weights for components' distributions.
#'
#' @return object from class \code{wtsamp} whose cummulative weights are
#' inf-corrected and individual weights are set to NULL.
#'
#' @details If cummulative weights are NULL they will be calculated
#' from the individual ones by \code{indiv2cum}.
#'
#' @examples
#' p <- genunifp(15,2) # create mixing probabilities
#' a <- lsweight(p) # calculate minimax weights
#' # create a weighted sample:
#' xs <- wtsamp(genormixt(p,c(0,1),c(1,1)),indiv=a)
#' F1 <- edfgen(xs,1) # minimax the estimate for 1-st component
#' xs_inf <- corinf(xs)
#' F1_inf <- edfgen(xs_inf,1) # inf-corrected estimate for 1-st component
#' F1(0)
#' F1_inf(0)
#' plot(F1,-3,3)
#' curve(F1_inf,col="red",,lty="dashed",add=TRUE)
#' @seealso
#' Maiboroda R.  and  Kubaichuk O.
#' Asymptotic normality of improved weighted empirical distribution functions.
#' Theor. Probability and Math. Statist. 69 (2004), 95-102 
#' @export 
# ________________________________________________________
corinf <- function(xs) {
    if (is.null(xs$cumm))
        xs <- indiv2cumm(xs)
    n <- length(xs$x) - 1
    xs$cumm <- xs$cumm[-1, ]
    xs$cumm <- rbind(rep(0, ncol(xs$cumm)), pmax(apply(xs$cumm[n:1, ], 2, cummin), 0)[n:1,
        ])
    xs$indiv <- NULL
    xs
}
#' Mid-correction of cummulative weights of a sample
#'
#' \code{cormid} calculates cummulative weights as the mean
#' of sup- and inf-corrected weights
#'   (mid-corrected weights).
#'
#' @param xs object from class \code{wtsamp} containing the sample
#' and weights for components' distributions.
#'
#' @return object from class \code{wtsamp} whose cummulative weights are
#' mid-corrected and individual weights are set to NULL.
#'
#' @details If cummulative weights are NULL they will be calculated
#' from the individual ones by \code{indiv2cum}.
#'
#' @examples
#' p <- genunifp(15,2) # create mixing probabilities
#' a <- lsweight(p) # calculate minimax weights
#' # create a weighted sample:
#' xs <- wtsamp(genormixt(p,c(0,1),c(1,1)),indiv=a)
#' F1 <- edfgen(xs,1) # minimax the estimate for 1-st component
#' xs_mid <- cormid(xs)
#' F1_mid <- edfgen(xs_mid,1) # mid-corrected estimate for 1-st component
#' F1(0)
#' F1_mid(0)
#' plot(F1,-3,3)
#' curve(F1_mid,col="red",,lty="dashed",add=TRUE)
#' @seealso
#' Maiboroda R.  and  Kubaichuk O.
#' Asymptotic normality of improved weighted empirical distribution functions.
#' Theor. Probability and Math. Statist. 69 (2004), 95-102 
#' @export 
# ________________________________________________________
cormid <- function(xs) {
    if (is.null(xs$cumm))
        xs <- indiv2cumm(xs)
    cummin <- pmin(apply(xs$cumm, 2, cummax), 1)
    n <- length(xs$x) - 1
    xs$cumm <- xs$cumm[-1, ]
    xs$cumm <- (rbind(rep(0, ncol(xs$cumm)), pmax(apply(xs$cumm[n:1, ], 2, cummin), 0)[n:1,
        ]) + cummin)/2
    xs$indiv <- NULL
    xs
}
#_______________________________________
# RESAMPLING FUNCTIONS
#_______________________________________
#__________________________________________________
#' Resample form mixture  with warying concentrations.
#'
#' \code{gensampmixt} generates a sample from the mixture
#'  with mixing probabilities \code{p} and distributions of components given by
#' the weighted sample \code{xs}.
#'
#' @param p  matrix (or data frame) of mixing probabilities
#'  with rows corresponding to subjects
#'  and columns coresponding to the mixture components.
#' @param xs  object of class \code{wtsamp} with the weighted
#' sample. Defines the mixture components' distribution.
#'
#' @return a vector with the resampled sample.
#'
#' @details \code{gensampmixt} uses \code{randwtgen} for sampling
#' from components with the default value of the option
#' \code{delta0}.
#' @examples
#' set.seed(3)
#' p <- genunifp(10,2) # create mixing probabilities
#' a <- lsweight(p) # calculate minimax weights
#' # create a weighted sample:
#' xs <- wtsamp(genormixt(p,c(0,1),c(1,1)),indiv=a)
#' xs<-cormid(xs) # correct the weights
#' x_resampled<-gensampmixt(p,xs) # resampled data
#' @seealso
#' Maiboroda R., Sugakova O. "Statistics of mixtures with varying concentrations 
#' with application to DNA microarray data analysis". 
#' Nonparametric statistics (2012) v.24:1, p. 201 - 215.
#' @seealso
#' Maiboroda R.  and  Kubaichuk O.
#' Asymptotic normality of improved weighted empirical distribution functions.
#' Theor. Probability and Math. Statist. 69 (2004), 95-102 
#' @export 
# __________________________________________________
gensampmixt <- function(p, xs) {
  randwt <- randwtgen(xs)
  crange <- 1:ncol(p)
  components <- apply(p, 1, function(pr) sample(crange, 1, replace = TRUE, prob = pr))
  sapply(components, randwt, n = 1)
}
#' Creates a random number generator according to a weighted sample distribution
#'
#' \code{randwtgen} returns a function wchich produces random samples with the
#' distribution of a prescribed mixture component. The distribution is estimated
#' by a weighted sample.
#'
#' @param xs object from class \code{wtsamp} containing the sample
#' and weights for components' distributions.
#'
#' @return a function \code{f(m,n,delta)}.
#' The function generates a sample of size \code{n} with the distribution
#' of the \code{m}-th component of the mixture. \code{delta} is the blurring
#' parameter.
#'
#' If \code{delta=0} the generated sample contains values from
#' \code{xs$x0} sampled with probabilities given by the \code{m}-th
#' column of \code{xs$indiv}.
#' If  \code{delta>0} a random variable uniform on [-\code{delta},\code{delta}]
#' is added to each sampled value.
#'
#'  The default value of \code{delta} is a half of the minimal distance
#'  between \code{xs$x0} points.
#' @seealso
#' Maiboroda R., Sugakova O. "Statistics of mixtures with varying concentrations 
#' with application to DNA microarray data analysis". 
#' Nonparametric statistics (2012) v.24:1, p. 201 - 215.
#' @export 
#__________________________________________________________________
randwtgen <- function(xs) {
    x <- xs$x[-1]
    if (is.null(xs$indiv))
        xs <- cumm2indiv(xs)
    prob <- xs$indiv
    delta0 = min(diff(x))/2
    randwt <- function(m, n, delta = delta0) {
        r <- sample(x, n, prob = prob[, m], replace = TRUE)
        if (delta > 0)
            r <- r + stats::runif(n, -delta, delta)
        r
    }
}
#________________________________________
#' Calculates the means of all components for a weighted sample.
#'
#'  \code{meanw} calculates weighted means of a sample
#'  using the individual weights of all components.
#'
#' @param xs object from class \code{wtsamp} containing the sample
#' and weights for components' distributions.
#'
#' @return a vector of components' means
#'
#' @examples
#' set.seed(3)
#' p <- genunifp(1000,2) # create mixing probabilities
#' a <- lsweight(p) # calculate minimax weights
#' # create a weighted sample:
#' xs <- wtsamp(genormixt(p,c(0,1),c(1,1)),indiv=a)
#' meanw(xs)
#' @export 
# _______________________________________
meanw <- function(xs) {
    if (is.null(xs$indiv))
        xs <- cumm2indiv(xs)
    mx<-as.vector(xs$xo[-1] %*% as.matrix(xs$indiv))
    names(mx)<-colnames(xs$indiv)
    return(mx)
}
#________________________________________
#' Calculates the variances of all components for a weighted sample.
#'
#'  \code{varw} calculates weighted variances of a sample
#'  using the individual weights of all components.
#'
#' @param xs object from class \code{wtsamp} containing the sample
#' and weights for components' distributions.
#'
#' @return a vector of components' variances
#'
#' @examples
#' set.seed(3)
#' p <- genunifp(1000,2) # create mixing probabilities
#' a <- lsweight(p) # calculate minimax weights
#' # create a weighted sample:
#' xs <- wtsamp(genormixt(p,c(0,1),c(1,1)),indiv=a)
#' varw(xs)
#' @seealso
#' Maiboroda R., Sugakova O. "Statistics of mixtures with varying concentrations 
#' with application to DNA microarray data analysis". 
#' Nonparametric statistics (2012) v.24:1, p. 201 - 215.
#' @export 
# _______________________________________
varw <- function(xs){
  if (is.null(xs$indiv))
    xs <- cumm2indiv(xs)
  sx<-as.vector((xs$xo[-1])^2 %*% as.matrix(xs$indiv)) -
    (meanw(xs))^2
  names(sx)<-colnames(xs$indiv)
  return(sx)
}
#________________________________________
#' Calculates the standard deviations of all components for a weighted sample.
#'
#'  \code{sdw} calculates weighted standard deviations of a sample
#'  using the individual weights of all components.
#'
#' @param xs object from class \code{wtsamp} containing the sample
#' and weights for components' distributions.
#'
#' @param corr function used for correction of weights before sd caclucation.
#' (no correction if code{corr=NULL}).
#'
#' @return a vector of components' standard deviations
#'
#' @examples
#' set.seed(3)
#' p <- genunifp(1000,2) # create mixing probabilities
#' a <- lsweight(p) # calculate minimax weights
#' # create a weighted sample:
#' xs <- wtsamp(genormixt(p,c(0,1),c(1,1)),indiv=a)
#' sdw(xs)
#' @seealso
#' Maiboroda R., Sugakova O. "Statistics of mixtures with varying concentrations 
#' with application to DNA microarray data analysis". 
#' Nonparametric statistics (2012) v.24:1, p. 201 - 215.
#' @export 
# _______________________________________
sdw <- function(xs,corr=cormid){
  if(is.null(corr))
    {if (is.null(xs$indiv)) xs <- cumm2indiv(xs)}
  else 
   {xs <- cumm2indiv(corr(xs))}
  sx<-sqrt(as.vector((xs$xo[-1])^2 %*% as.matrix(xs$indiv)) -
    (meanw(xs))^2)
  names(sx)<-colnames(xs$indiv)
  return(sx)
}
#________________________________________
#' Calculates the quantiles of all components for a weighted sample.
#'
#'  \code{quantilew} calculates weighted quantiles of a sample
#'  using the individual weights of all components.
#'
#' @param xs object from class \code{wtsamp} containing the sample
#' and weights for components' distributions.
#'
#' @param prob the level of the quantiles (one number).
#'
#' @return a vector of components' quantiles.
#'
#' @examples
#' set.seed(3)
#' p <- genunifp(1000,2) # create mixing probabilities
#' a <- lsweight(p) # calculate minimax weights
#' # create a weighted sample:
#' xs <- wtsamp(genormixt(p,c(0,1),c(1,1)),indiv=a)
#' quantilew(xs,1/3)
#' @seealso
#' Maiboroda R., Sugakova O. "Statistics of mixtures with varying concentrations 
#' with application to DNA microarray data analysis". 
#' Nonparametric statistics (2012) v.24:1, p. 201 - 215.
#' @export 
# _______________________________________
quantilew <- function(xs,prob){
  if (is.null(xs$cumm))
    xs <- indiv2cumm(xs)
  n <- nrow(xs$cumm)
  M <- ncol(xs$cumm)
  q <- numeric(M)
  for( m in 1:M ){
    j <- 1
    while( xs$cumm[j,m]<prob ) j <- j+1
    q_left <- xs$xo[j]
    j <- n
    while( xs$cumm[j,m]>prob ) j <- j-1
    q_right <- xs$xo[j]
    q[m] <- ( q_left + q_right )/2
  }
  names(q)<-colnames(xs$cumm)
  return(q)
}
#________________________________________
#' Calculates the medians of all components for a weighted sample.
#'
#'  \code{medianw} calculates weighted medians of a sample
#'  using the individual weights of all components.
#'
#' @param xs object from class \code{wtsamp} containing the sample
#' and weights for components' distributions.
#'
#' @return a vector of components' medians
#'
#' @examples
#' set.seed(3)
#' p <- genunifp(1000,2) # create mixing probabilities
#' a <- lsweight(p) # calculate minimax weights
#' # create a weighted sample:
#' xs <- wtsamp(genormixt(p,c(0,1),c(1,1)),indiv=a)
#' medianw(xs)
#' @export 
# _______________________________________
medianw <- function(xs)quantilew(xs,0.5)
#________________________________________
#' Calculates the interquartile ranges of all components for a weighted sample.
#'
#'  \code{IQRw} calculates weighted interquartile ranges of a sample
#'  using the individual weights of all components.
#'
#' @param xs object from class \code{wtsamp} containing the sample
#' and weights for components' distributions.
#'
#' @return a vector of components' interquartile ranges
#'
#' @examples
#' set.seed(3)
#' p <- genunifp(1000,2) # create mixing probabilities
#' a <- lsweight(p) # calculate minimax weights
#' # create a weighted sample:
#' xs <- wtsamp(genormixt(p,c(0,1),c(1,1)),indiv=a)
#' IQRw(xs)
#' @export 
# _______________________________________
IQRw <- function(xs){
  quantilew(xs,0.75) - quantilew(xs,0.25)
}
#_________________________________________
#' Epanechnikov kernel calculation.
#'
#' \code{Epanechn} calculates the value of Epanechnikov kernel at a given point.
#'
#' @param x a vector of points at which the kernel value is calculated.
#'
#' @return vector of the kernel values at given points.
#'
#' @examples
#' curve(Epanechn,-2,2)
#' @export
#_________________________________________
Epanechn<-function(x)ifelse(abs(x)<1,0.75*(1-x^2),0)
#_________________________________________
#' Generator of kernel density estimator for mixture components.
#'
#' \code{densgen} generates a function which calculates
#' a kernel density estimate for a prescribed mixture component
#' at a given set of points.
#'
#' @param xs a \code{wtsamp}  object representing
#' a sample with distributions of different components.
#'
#' @param m a number of component whose density is estimated.
#'
#' @param Kern a function which calculates the kernel values
#' at given points (must be vectorized).
#'
#' @return a function
#' \code{f(x,h)} which calculates the estimate at points given in the
#' vector \code{x} with the bandwidth  \code{h}.
#'
#' @examples
#' set.seed(3)
#' p <- genunifp(1000,2) # create mixing probabilities
#' a <- lsweight(p) # calculate minimax weights
#' # create a weighted sample:
#' xs <- wtsamp(genormixt(p,c(0,1),c(1,1)),indiv=a)
#' f<-densgen(xs,1) # create the estimator
#' f(c(0,1),1)      # calculate the estimate at two points
#' curve(f(x,1),-3,3) # plot the graph (estimates N(0,1) density)
#'
#' @seealso
#' Maiboroda R., Sugakova O. "Statistics of mixtures with varying concentrations 
#' with application to DNA microarray data analysis". 
#' Nonparametric statistics (2012) v.24:1, p. 201 - 215.
#' @export 
#_________________________________________
densgen<-function(xs,m,Kern=Epanechn){
  if (is.null(xs$indiv))
    xs <- cumm2indiv(xs)
  a <- xs$indiv[, m]
  x <- xs$xo[-1]
  f <- Vectorize(
    function(t,h) {
      sum(a*Kern((t-x)/h))/h
    },
    vectorize.args ="t"
  )
}
#_________________________________________
#' Silverman's rule of thumb for bandwidth selection
#'
#' \code{silvbw} calculates quasi-optimal bandwidth
#' for kernel density estimate by weighted sample
#' based on the mixture with varying concentrations approach
#'
#' @param xs a \code{wtsamp}  object representing
#' a sample with distributions of different components.
#'
#' @param m a number of component for which the density is
#' estimated.
#'
#' @param delta the canonical bandwidth of the kernel
#' used in the estimate.
#'@details The  default value of \code{delta}=1.7188 corresponds to
#' the Epanechnikov kernel.
#'
#' @return a numeric value of the Silverman's optimal bandwidth.
#'
#' @examples
#' #' @examples
#' set.seed(3)
#' p <- genunifp(1000,2) # create mixing probabilities
#' a <- lsweight(p) # calculate minimax weights
#' # create a weighted sample:
#' xs <- wtsamp(genormixt(p,c(0,1),c(1,1)),indiv=a)
#' f<-densgen(xs,1) # create the estimator
#' h<-silvbw(xs,1) # calculates the bandwidth by the Silverman's rule
#' curve(f(x,h),-3,3) # plot the graph (estimates N(0,1) density)
#'
#' @seealso
#' Maiboroda R., Sugakova O. "Statistics of mixtures with varying concentrations 
#' with application to DNA microarray data analysis". 
#' Nonparametric statistics (2012) v.24:1, p. 201 - 215.
#' @export 
#_________________________________________
silvbw<-function(xs,m,delta=1.7188){
  if (is.null(xs$indiv))
    xs <- cumm2indiv(xs)
  const<-1.011354 # (8*sqrt(pi)/3)^(1/5)/(stats::qnorm(0.75)-stats::qnorm(0.25))
  const*delta*(sum((xs$indiv[,m])^2))^(1/5)*min(IQRw(xs)[m],sdw(xs)[m])
}

#________________________________________
# INFERENCE
#________________________________________
#
#' Estimates for standard deviations and CI
#' for weighted means by observations from the mixture.
#'
#'  \code{sdMean} calculates estimates of standard deviatitions
#'  and confidence intervals
#'  for  weighted means with minimax weights by observations
#'  from the mixture with varying concentrations.
#'
#' @param x numeric vector with the observed sample or a
#' \code{wtsamp} object.
#'
#' @param p matrix (or data frame) of mixing probabilities
#'  with rows corresponding to subjects
#'  and columns coresponding to the mixture components.
#' @param comp a numeric vector with numbers of components
#' for which the standard deviations are estimated.
#'
#' @param means logical, if \code{TRUE} then the estimates for
#' components' means are included in the function value.
#'
#' @param CI logical, if \code{TRUE} then confidence bounds for
#' components' means are inculded in the function value.
#'
#' @param alpha confidense level for the confidence interval.
#'
#' @details If \code{CI=TRUE} then the function calculates
#' confidence intervals for the components' means
#' with covering probability \code{1-alpha}.
#'
#' If \code{x} is a vector then the weights for components' means and variances
#' are calculated as \code{lsweight(p)}. If \code{x} is a \code{wtsamp}
#' object than its own weights are used.
#'
#' @return if \code{CI & means =FALSE} the function returns a vector
#' of the estimated standard deviations
#' with NA for the components which were not estimated.
#' Else a data frame is returned in which there can be variables:
#' \code{sd} are standard deviations of estimates;
#' \code{means} are the estimates of means;
#' \code{lower} and \code{upper} are lower and upper bounds
#' of the confidence intervals for means.
#'
#'
#' @examples
#' set.seed(3)
#' M<-3 # number of mixture components
#' p <- genunifp(1000,M) # create mixing probabilities
#' m<-c(0,1,2)      # true means of components
#' sd<-c(1,1,0.5)   # true sd of components
#' x<-genormixt(p,m,sd) # sample generation
#' # Calculate sd only:
#' sdMean(x,p)
#' # the same:
#' sdMean(wtsamp(x,indiv=lsweight(p)),p)
#' # Calculate confidence intervals:
#' sdMean(x,p,means=TRUE,CI=TRUE)
#' # Plot confidence intervals:
#' CI<-sdMean(x,p,means=TRUE,CI=TRUE)
#' library(plotrix)
#' plotCI(1:M,CI$means,ui=CI$upper,li=CI$lower,
#'        xlab=" ",ylab="means",xaxt="n")
#' axis(1,at=1:M,labels=row.names(CI))
#' @seealso
#' Maiboroda R.  and  Kubaichuk O.
#' Asymptotic normality of improved weighted empirical distribution functions.
#' Theor. Probability and Math. Statist. 69 (2004), 95-102 
#' @export 
# _______________________________________
sdMean<-function(x,p,comp=1:ncol(p),
                 means=FALSE,CI=FALSE,alpha=0.05){
  M<-ncol(p)
  D<-rep(NA,M)
  if(is.vector(x)&is.numeric(x))
  sx<-wtsamp(x,indiv=lsweight(p))
  else{
    if(class(x)=="wtsamp")
      sx<-x
    else{
      warning("x must be vector or wtsamp")
    return(NA)
    }
  }
  m<-meanw(sx)
  m2<-varw(sx)+m^2
  for(k in comp){
    app<-matrix(ncol=M,nrow=M)
    for(i in 1:M){
      for(l in 1:i){
        app[i,l]<-sum(sx$indiv[,k]^2*p[,i]*p[,l])
        app[l,i]<-app[i,l]
      }
    }
    ap<-apply(app,1,sum)
    D[k]=sum(ap*m2)-m%*%app%*%m
    if(D[k]<0){
      warning("Negative estimate of variance is obtained",D[k])
      D[k]=NA
    }
  }
  D=sqrt(D)
  if(!(means|CI)){
    names(D)<-colnames(p)
    return(D)
  }else{
    R<-data.frame(sd=D)
    if(means)R$means<-m
    if(CI){
      lambda=stats::qnorm(1-alpha/2)
      R$lower<-m-lambda*D
      R$upper<-m+lambda*D
    }
    row.names(R)<-colnames(p)
    return(R)
  }
}
#______________________________________________________
#' Estimates for standard deviations and CI
#' for weighted medians by observations from the mixture.
#'
#'  \code{sdMedian} calculates estimates of standard deviatitions
#'  and confidence intervals
#'  for  weighted medians with minimax weights by observations
#'  from the mixture with varying concentrations.
#'
#' @param x numeric vector with the observed sample or a
#' \code{wtsamp} object.
#'
#' @param p matrix (or data frame) of mixing probabilities
#'  with rows corresponding to subjects
#'  and columns coresponding to the mixture components.
#' @param comp a numeric vector with numbers of components
#' for which the standard deviations are estimated.
#'
#' @param medians logical, if \code{TRUE} then the estimates for
#' components' medians are included in the function value.
#'
#' @param CI logical, if \code{TRUE} then confidence bounds for
#' components' means are inculded in the function value.
#'
#' If \code{x} is a vector then the weights for components' medians and variances
#' are calculated as \code{lsweight(p)}. If \code{x} is a \code{wtsamp}
#' object than its own weights are used.
#'
#' @param alpha confidense level for the confidence interval.
#'
#' @details if \code{CI=TRUE} then the function calculates
#' confidence intervals for the components' medians
#' with covering probability \code{1-alpha}.
#'
#' @return if \code{CI & medians =FALSE} the function returns a vector
#' of the estimated standard deviations
#' with NA for the components which were not estimated.
#' Else a data frame is returned in which there can be variables:
#' \code{sd} are standard deviations of estimates;
#' \code{medians} are the estimates of medians;
#' \code{lower} and  b    n\code{upper} are lower and upper bounds
#' of the confidence intervals for medians.
#'c
#'
#' @examples
#' set.seed(3)
#' M<-3 # number of mixture components
#' p <- genunifp(1000,M) # create mixing probabilities
#' m<-c(0,1,2)      # true means of components
#' sd<-c(1,1,0.5)   # true sd of components
#' x<-genormixt(p,m,sd) # sample generation
#' # Calculate sd only:
#' sdMedian(x,p)
#' # the same result:
#' sdMedian(wtsamp(x,indiv=lsweight(p)),p)
#' # Calculate confidence intervals:
#' sdMedian(x,p,medians=TRUE,CI=TRUE)
#' # Plot confidence intervals:
#' CI<-sdMedian(x,p,medians=TRUE,CI=TRUE)
#' library(plotrix)
#' plotCI(1:M,CI$medians,ui=CI$upper,li=CI$lower,
#'        xlab=" ",ylab="medians",xaxt="n")
#' axis(1,at=1:M,labels=row.names(CI))
#' @export 
#________________________________________
sdMedian<-function(x,p,comp=1:ncol(p),
                   medians=FALSE,CI=FALSE,alpha=0.05){
  M<-ncol(p)
  D<-rep(NA,M)
   if(is.vector(x)&is.numeric(x)){
  sx<-wtsamp(x,indiv=lsweight(p))
  sx <- indiv2cumm(sx)
   }
  else{
    if(class(x)=="wtsamp")
      sx<-x
	  if(is.null(sx$cumm)) sx <- indiv2cumm(sx)
    else{
      warning("x must be vector or wtsamp")
    return(NA)
    }
  }
  med<-medianw(sx)
  for(k in comp){
    Fm<-sx$cumm[findInterval(med[k], sx$xo),]
    app<-matrix(ncol=M,nrow=M)
    for(i in 1:M){
      for(l in 1:i){
        app[i,l]<-sum(sx$indiv[,k]^2*p[,i]*p[,l])
        app[l,i]<-app[i,l]
      }
    }
    ap<-apply(app,1,sum)
    f<-densgen(sx,k)
    zz<-sum(ap*Fm)-Fm%*%app%*%Fm
    if(zz<0){
      warning("Negative estimate of variance is obtained",zz)
      zz=NA
      }
    else
      D[k]<-sqrt(zz)/f(med[k],silvbw(sx,k))
  }
  if(!(medians|CI)){
    names(D)<-colnames(p)
    return(D)
  }else{
    R<-data.frame(sd=D)
    if(medians)R$medians<-med
    if(CI){
      lambda=stats::qnorm(1-alpha/2)
      R$lower<-med-lambda*D
      R$upper<-med+lambda*D
    }
    row.names(R)<-colnames(p)
    return(R)
  }
}
# _______________________________________
# CENSORED DATA
# _______________________________________
#' Censored sample from
#' mixture of gamma distributions with varying concentrations.
#'
#' \code{gencensg} generates a censored sample from the mixture
#' of gamma distributions
#' with mixing probabilities of components given by
#' the matrix p.
#'
#' @param p  matrix (or data frame) of mixing probabilities.
#'  with rows corresponding to subjects.
#'  and columns coresponding to the mixture components.
#' @param shape  vector of shape parameters for gamma distributions of
#' non-censored components.
#' @param rate   vector of rate parameters for gamma distributions of
#' non-censored components.
#' @param shapec  vector of shape parameters for gamma distributions of
#' components' censors.
#' @param ratec   vector of rate parameters for gamma distributions of
#' components' censors.
#' @return   a data frame  with the censored sample. The attribute
#' \code{$x} contains the sample values, the attribute \code{$delta}
#' contains the indicators of non-censoring.
#'
#' The sample size equals the rows number of p.
#' @examples
#' set.seed(3)
#' cs <- gencensg(genunifp(50,2),shape=c(2,1),rate=c(0.5,0.5),
#' shapec=c(1,1),ratec=c(0.3,0.3))
#' plot(cs$x,col=cs$delta+1)
#' @seealso
#' Maiboroda R. Khizanov V.
#' A modified Kaplan–Meier estimator for a model of mixtures with varying concentrations
#' Theor. Probability and Math. Statist. 92 (2016), 109-116 
#' @export 
# __________________________________________________
gencensg <- function(p, shape, rate, shapec, ratec) {
    nc <- 1:ncol(p)
    obj <- apply(p, 1, function(x) sample(nc, 1, prob = x))
    x <- stats::rgamma(nrow(p), shape = shape[obj], rate = rate[obj])
    c <- stats::rgamma(nrow(p), shape = shapec[obj], rate = ratec[obj])
    delta <- x < c
    x <- pmin(x, c)
    cs <- data.frame(x, delta)
}
# ______________________________________________________
#' Kaplan-Meyer's estimator for CDFs of mixture components
#'
#' \code{KMcdf} calculates cummulative weigts for the Kaplan-Meyer
#' estimator of CDFs of mixture components by a censored sample
#' of class \code{wtcens} and returns an object of class
#' \code{wtsamp} with the sample values and the calculated weights.
#'
#' @param cs object of class \code{wtcens} censored sample with
#' weights representing the distributions
#' of the censored components.
#'
#' @return object of class \code{wtsamp} containing the sample values
#' and the cummulative weights for the distributions of the components.
#'
#' @examples
#' set.seed(3)
#' n<-5000
#' p<-genunifp(n,2)
#' xs<-gencensg(p,shape=c(1.5,4),rate=c(1,1),shapec=c(1,1),ratec=c(0.3,0.3))
#' cs<-wtcens(x=xs$x,delta=xs$delta,indiv=lsweight(p))
#' KMest<-KMcdf(cs)
#' FKM<-edfgen(KMest,2)
#' curve(FKM,0,8)
#' curve(pgamma(x,shape=4,rate=1),col="red",add=TRUE)
#' @seealso
#' Maiboroda R. Khizanov V.
#' A modified Kaplan–Meier estimator for a model of mixtures with varying concentrations
#' Theor. Probability and Math. Statist. 92 (2016), 109-116 
#' @export 
# ______________________________________________________
KMcdf <- function(cs) {
    if (is.null(cs$cumm))
        cs <- indiv2cumm(cs)
    if (is.null(cs$indiv))
        cs <- cumm2indiv(cs)
    wt <- 1 - apply(1 - cs$indiv * cs$deltao/(1 - utils::head(cs$cumm, n = -1)), 2, cumprod)
    wt[is.nan(wt)] <- 1
    wtsamp(cs$xo[-1], cumm = wt, indiv = NULL)
}
# ___________ The end.
