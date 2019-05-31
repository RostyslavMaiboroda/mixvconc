library(mixvconc)
context("Kaplan - Meyer estimate")
tol=0.001 # numerical tolerance
set.seed(3)
n<-1000
pp <- genunifp(n,2) # create mixing probabilities
aa <- lsweight(pp) # calculate minimax weights
xc<-gencensg(pp,c(1,2),c(1,0.5),c(1,1),c(1,0.3))
wxc<-wtcens(xc$x,xc$delta,indiv=aa)
kms<-KMcdf(wxc)
F1<-edfgen(kms,1)
F2<-edfgen(kms,2)
# for plot:
# curve(F1(x),0,4,col="red")
# curve(pgamma(x,shape=1,rate=1),add=T)
# curve(F2(x),0,10,col="red")
# curve(pgamma(x,shape=2,rate=0.5),add=T)
FKM1<-F1(2)
names(FKM1)<-NULL
FKM2<-F2(2)
names(FKM2)<-NULL
test_that("Normal mixture",{
  expect_equal(FKM1,0.8616484,tolerance = tol,info="means")
  expect_equal(FKM2,0.270229,tolerance = tol,info="means")
})
