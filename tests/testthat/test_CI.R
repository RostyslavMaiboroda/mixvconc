library(mixvconc)
context("Confidence intervals")
tol=0.001 # numerical tolerance
set.seed(3)
n<-1000
pp <- genunifp(n,2) # create mixing probabilities
aa <- lsweight(pp) # calculate minimax weights
# create a weighted sample:
xx <- genormixt(pp,c(0,1),c(1,1))
CImean<-sdMean(xx,pp,means=TRUE,CI=TRUE,alpha=0.01)
CImedian<-sdMedian(xx,pp,medians=TRUE,CI=TRUE,alpha=0.01)
test_that("Normal mixture",{
  expect_equal(CImean["A","means"],-0.0008397801,tolerance = tol,info="mean_mean")
  expect_equal(CImean["B","sd"],0.06890969,tolerance = tol,info="mean_sd")
  expect_equal(CImean["B","lower"],0.8039885,tolerance = tol,info="mean_lower")
  expect_equal(CImean["A","upper"],0.1757424,tolerance = tol,info="mean_upper")
  expect_equal(CImedian["A","medians"],-0.07682753,tolerance = tol,info="mean_mean")
  expect_equal(CImedian["B","sd"],0.07789651,tolerance = tol,info="mean_sd")
  expect_equal(CImedian["B","lower"],0.7622856,tolerance = tol,info="mean_lower")
  expect_equal(CImedian["A","upper"],0.1271079,tolerance = tol,info="mean_upper")
})
