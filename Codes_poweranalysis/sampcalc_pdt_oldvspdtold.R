rm(list=ls())

##functions
library(lme4)
library(ggplot2); library(ggpubr)
source('SWSamp_binomial.R')

##calculate OR based on p0 and p1
or.calc <- function(p0, p1) {
  or <- p1/(1-p1)/(p0/(1-p0))
  or
}
##calculate p1 based on p0 and OR
p1.calc <- function(p0, OR) {
  odds1 <- OR*p0/(1-p0)
  p1 <- odds1/(1+odds1)
  p1
}
##we are comparing (CS+venous) vs PDT+(CS+venous) 
##for CS+venous, observed % of under-threshold donations are 5.4% (men) and 7.9% (women)
##for PDT+(CS+venous), simulated % of under-threshold donations are 4.9% (men) and 7.8% (women)
or.calc(p0=0.054, p1=0.049)##OR=0.90
or.calc(p0=0.079, p1=0.078)##OR=0.99
p1.calc(p0=0.054, OR=0.9)
p1.calc(p0=0.079, OR=0.99)

##final sample size calculation, vary the ss candidate
##I: num of clusters; J: number of phases (excluding baseline)
##rho=0.5 used in the original calculation
ss_calc <- function(ss_candidate, p0, OR, n.sims, ni.margin, I, J, rho=NULL, sigma.a=NULL, formula=NULL) {
  set.seed(123)  
  start <- Sys.time()
  res <- sim.power.ni.simplified(I=I, J=J, K=ss_candidate, design='cross-sec', mu=p0, b.trt=OR, 
                                 b.time=0,##timeeff=0
                                 rho=rho, rho.ind=NULL, family='binomial', natural.scale = T,
                                 sigma.a=sigma.a,
                                 ni.margin=ni.margin,
                                 formula=formula,
                                 n.sims = n.sims)
  runtime <- Sys.time()-start
  print(paste0('Simulation finished for sample size of ', ss_candidate))
  ##return all the results
  result<-data.frame(ni.margin=ni.margin, ss=ss_candidate, power=res, n.sims=n.sims, 
                     runtime=round(as.numeric(runtime, units='mins'), 2))
  result
}

##men: 0.054 vs 0.049, or=0.90
##3 centres, 4 phases
##n=(175+203+190)/3*(30/7)*4=~3250, number of donations per center per 4 months for men
##ni margin: 0.01; power: 0.904; time: 28min, 3 centres, 4 months per phase
lapply(3250, ss_calc, p0=0.054, OR=0.9, n.sims=500, ni.margin=0.01, I=3, J=3, rho=0.5)

##women: 0.079 vs 0.078, or=0.99
##3 centres, 4 phases
##n=(140+154+185)/3*(30/7)*4=~2740, number of donations per center per 4 months for women
##ni margin: 0.01; power: 0.42; time: 21min, 3 centres, 4 months per phase
lapply(2740, ss_calc, p0=0.079, OR=0.99, n.sims=500, ni.margin=0.01, I=3, J=3, rho=0.5)
##ni margin: 0.015; power: 0.732; time: 25min, 3 centres, 4 months per phase
lapply(2740, ss_calc, p0=0.079, OR=0.99, n.sims=500, ni.margin=0.015, I=3, J=3, rho=0.5)
##ni margin: 0.02; power: 0.922; time: 21min, 3 centres, 4 months per phase
lapply(2740, ss_calc, p0=0.079, OR=0.99, n.sims=500, ni.margin=0.02, I=3, J=3, rho=0.5)

##formal calculation
power_curve <- function(sex, ss_candidate, p0, OR, n.sims, ni.margin, I, J, rho=0.5) {
  ss <- lapply(ss_candidate, ss_calc, p0=p0, OR=OR, n.sims=n.sims, ni.margin=ni.margin, 
               I=I, J=J, rho=rho)
  ss1 <- do.call(rbind, ss)
  ss1$sex <- sex
  write.csv(ss1, paste0('ss_', sex, '_', ni.margin, 
                        '_', min(ss_candidate), '-', max(ss_candidate), '_', I,'centres.csv'), row.names = F)
}
##men
power_curve(sex='m', ss_candidate = seq(1200, 3200, 400), p0=0.054, OR=0.9, n.sims = 1000, ni.margin=0.01, 
            I=3, J=3)
##women
power_curve(sex='f', ss_candidate = seq(800, 2800, 400), p0=0.079, OR=0.99, n.sims = 1000, ni.margin=0.02, 
            I=3, J=3)
