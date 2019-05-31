library(mixvconc)
context("Descriptives")
tol=0.001 # numerical tolerance
set.seed(3)
pp <- genunifp(10,2) # create mixing probabilities
aa <- lsweight(pp) # calculate minimax weights
# create a weighted sample:
xxs <- wtsamp(genormixt(pp,c(0,1),c(1,1)),indiv=aa)
x_mean<-meanw(xxs)
t_mean<-setNames(c( 0.3994509, 0.4432229),c("A","B"))
x_q<-quantilew(xxs,1/3)
t_q<-setNames(c( -0.9409584,  0.2232919),c("A","B"))
x_med<-medianw(xxs)
t_med<-setNames(c( -0.1984272,  0.2232919),c("A","B"))
x_var<-varw(xxs)
t_var<-setNames(c(2.468498, 0.807258),c("A","B"))
x_sd<-sdw(xxs)
t_sd<-setNames(c(1.5516570, 0.8305492),c("A","B"))
x_IQR<-IQRw(xxs)
t_IQR<-setNames(c(Inf, 0.1383226),c("A","B"))
test_that("Normal mixture",{
  expect_equal(x_mean,t_mean,tolerance = tol,info="means")
  expect_equal(x_q,t_q,tolerance = tol,info="quantiles")
  expect_equal(x_med,t_med,tolerance = tol,info="median")
  expect_equal(x_var,t_var,tolerance = tol,info="var")
  expect_equal(x_sd,t_sd,tolerance = tol,info="sd")
  expect_equal(x_IQR,t_IQR,tolerance = tol,info="IQR")
})

