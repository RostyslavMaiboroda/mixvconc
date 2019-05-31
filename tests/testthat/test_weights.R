library(mixvconc)
context("Weights calculation and improvement")
tol=0.001 # numerical tolerance
set.seed(3)
M<-3
pp <- genunifp(10,M) # create mixing probabilities
aa <- lsweight(pp) # calculate minimax weights
# test the lsweight
unit_matrix<-t(pp)%*%as.matrix(aa) # must be a unit matrix
#________________________________
# create a weighted sample:
xxs <- wtsamp(gengamixt(pp,c(3,2,1),c(1,2,3)),indiv=aa)
xxsc<-indiv2cumm(xxs)
xxsi<-cumm2indiv(xxsc)
indivdiff<-xxs$indiv-xxsi$indiv # must be zero
#________________________________
# testing corsup
xxsup<-corsup(xxs)
corsup01<-abs(max(xxsup$cumm)-1)+abs(min(xxsup$cumm)) #> 0
corsupC<-min(diff(xxsup$cumm[,2]))>=0 #> TRUE
#________________________________
# testing corinf
xxinf<-corinf(xxs)
corinf01<-abs(max(xxinf$cumm)-1)+abs(min(xxinf$cumm)) #> 0
corinfC<-min(diff(xxinf$cumm[,2]))>=0 #> TRUE
#________________________________
# sup correction is greater then the inf one:
supinfC<-min(xxsup$cumm-xxinf$cumm)>=0
#________________________________
# testing cormid:
xxmid<-cormid(xxs)
cormidC<-(min(xxsup$cumm-xxmid$cumm)>=0)&
        (min(xxmid$cumm-xxinf$cumm)>=0)
#_________________________________
test_that("Gamma mixture",{
  expect_equal(sum(abs(unit_matrix-diag(rep(1,M)))),
               0,tolerance = tol,info="lsweight")
  expect_equal(sum(abs(indivdiff)),0,tolerance = tol,info="indiv|cumm")
  expect_equal(corsup01,0,tolerance = tol,info="corsup between 0 and 1")
  expect_true(corsupC,info="corsup increasing")
  expect_equal(corinf01,0,tolerance = tol,info="corinf between 0 and 1")
  expect_true(corinfC,info="corinf increasing")
  expect_true(supinfC,info="sup > inf")
  expect_true(cormidC,info="sup > mid > inf")
})
