library(mixvconc)
context("CDF and density estimates")
tol=0.001 # numerical tolerance
set.seed(3)
n<-10 # for plot: 1000
pp <- genunifp(n,2) # create mixing probabilities
aa <- lsweight(pp) # calculate minimax weights
# create a weighted sample:
xxs <- wtsamp(genormixt(pp,c(0,1),c(1,1)),indiv=aa)
# testing empirical CDF
F1<-edfgen(xxs,1)
F2<-edfgen(xxs,2)
# plot with n<-1000 at the begining of the test:
# curve(F1(x),-3,3,col="red")
# curve(pnorm(x),add=T)
# curve(F2(x),-2,4,col="red")
# curve(pnorm(x,mean=1),add=T)
#__________________________________
Femp1<-F1(0)
names(Femp1)<-NULL
t_Femp1<-0.4053883
Femp2<-F2(1)
names(Femp2)<-NULL
t_Femp2<-0.9725569
#_________________________________
# testing kernel density estimator
Epan<-Epanechn(0)
f1<-densgen(xxs,1)
f2<-densgen(xxs,2)
h1<-silvbw(xxs,1)
t_h1<-2.190939
h2<-silvbw(xxs,2)
t_h2<-0.1829349
fkern1<-f1(0,h1)
t_fk1<-0.2234312
fkern2<-f2(1,h2)
t_fk2<-0.2086134
# plot with n<-1000 at the begining of the test:
# curve(f1(x,h1),-3,3,col="red")
# curve(dnorm(x),add=T)
# curve(f2(x,h2),-2,4,col="red")
# curve(dnorm(x,mean=1),add=T)
#__________________________________
test_that("Normal mixture",{
  expect_equal(c(Femp1,Femp2),c(t_Femp1,t_Femp2),tolerance = tol,info="empCDF")
  expect_equal(Epan,0.75,tolerance=tol,info="Epanechnikov")
  expect_equal(c(h1,h2),c(t_h1,t_h2),
               tolerance=tol,info="Silverman")
  expect_equal(c(fkern1,fkern2),c(t_fk1,t_fk2),
               tolerance=tol,info="kernelPDF")
})
