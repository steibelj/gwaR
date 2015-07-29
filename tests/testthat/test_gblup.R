

rm(list=ls())
library(gwaR)
library(regress)
library(testthat)
data("MSUPRP_sample")
data("Z_G")
data("test")

design_G<-~sex+car_wt


test_that("summary gblup",{
  expect_that(gblup(rsp="temp_24h",data=MSUPRP_sample,design=design_G,A,pos=c(T,T)),equals(tst))
})


test_that("testing gwa ",{
  expect_that(GWAS(x=t(Z_thin),gblup=tst),equals(gwtst))
})

test_that("testing gwa summary",{
  expect_that(summary(gwtst,correction = "none"),equals(smr))
})

test_that("testing plots",{
  expect_that(plot(gwas=gwtst,correction="none",gpdata=MSUPRP_sample,
                   plotlog10=F,pvalue=0.05,q.qplot=F),equals(pvals))
})

sm<-summary(gblup(rsp="temp_24h",data=MSUPRP_sample,design=design_G,A,pos=c(T,T)),fe=T,sigma=T,ehat=T)
anova(gblup(rsp="temp_24h",data=MSUPRP_sample,design=design_G,A,pos=c(T,T)))
plot(GWAS(x=t(Z_thin),gblup=tst),pvalue=0.0001)
GWAS

