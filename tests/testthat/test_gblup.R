

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
  expect_that(gwas(x=t(Z_thin),gblup=tst),equals(gwtst))
})

test_that("testing gwa summary",{
  expect_that(summary(gwtst,correction = "none"),equals(smr))
})

test_that("testing plots",{
  expect_that(plot(gwas=gwtst,correction="none",gpdata=MSUPRP_sample,
                   plotlog10=F,pvalue=0.05,q.qplot=F),equals(pvals))
})

design_G<-c(~sex+car_wt,~slgdt_cd)
design_G<-c(~sex+car_wt+slgdt_cd)

system.time({
gb<-gblup(rsp="ph_24h",data=MSUPRP_sample,design=design_G,G_autosome,pos=c(T,T))
#sm<-summary(gb)
#anova(gblup(rsp="temp_24h",data=MSUPRP_sample,design=design_G,A,pos=c(T,T)))

gw<-gwas(x=t(Z_thin),gblup=gb)
})
gw

system.time({gw1<-sapply(colnames(MSUPRP_sample$pheno),run.gwa, data=MSUPRP_sample, design=design_G, G=G_autosome, x=t(Z_thin), 
             LRT = F, threshold = 0.01, 
             returnz = T, 
             saveblup = F, basename = "", pos=c(T,T))})

pv<-getpvalue(gwas = gw)
snpeak<-which.max(pv)
Zpeak<-Z_thin[,(snpeak-20):(snpeak+20)]
Gpeak<-Zpeak%*%t(Zpeak)
Gpeak<-Gpeak/mean(diag(Gpeak))*mean(diag(G_autosome))
designp<-c(~sex+car_wt,~G_autosome)
gp<-gblup(rsp="temp_24h",data=MSUPRP_sample,design=designp,vdata=list(G_autosome=G_autosome),G=Gpeak,pos=c(T,T,T))


test.peak(gb,t(Z_thin),pks)

plot(gwas(x=t(Z_thin),gblup=tst),pvalue=0.0001)
tst


