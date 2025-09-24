rm(list=ls())

##Packages
library(MASS); library(flexsurv); library(simsurv); library(lme4); library(splines); library(ggplot2)

setwd('C:/Research/Cambridge/Post-donation-Hb')
source("./Analysis/Simulation/onsite hemocue simulation/strides_donor_des_current_hemocue.R")

#############
#### MEN ####
#############
df_ana_lmm_male <- df_ana_lmm[df_ana_lmm$sex==1&df_ana_lmm$time<=52*1.5,]
quantiles <- c(0.05, 0.1, 0.15, 0.25, 0.5, 0.75, 0.85, 0.9, 0.95)
quantile(df_ana_lmm_male$hb0[!duplicated(df_ana_lmm_male$identifier)], quantiles)

##Spaghetti plot
id_lowhb <- df_ana_lmm_male$identifier[df_ana_lmm_male$outcome==2]
id_lowhb <- id_lowhb[!duplicated(id_lowhb)]
ggplot(data = df_ana_lmm_male, 
       aes(x = time, y = hb, group = identifier)) + geom_line()
ggplot(data = df_ana_lmm_male[df_ana_lmm_male$identifier%in%id_lowhb,], 
       aes(x = time, y = hb, group = identifier)) + geom_line()
ggplot(data = df_ana_lmm_male[!(df_ana_lmm_male$identifier%in%id_lowhb),], 
       aes(x = time, y = hb, group = identifier)) + geom_line()

##Predictors
plot(df_ana_lmm_male$exage, df_ana_lmm_male$hb)##linear
##time is almost linear
plot(df_ana_lmm_male$hb0, df_ana_lmm_male$hb)
plot(df_ana_lmm_male$eth, df_ana_lmm_male$hb)
##Model fitting: alternative covariates: visit, categorical age, nonlinear time effect
lmm_male <- lmer(hb ~ 
                   exage+
                   blood1+blood2+blood3+blood4+blood5+blood6+blood7+
                   eth + time +
                   # ndon1 + ndon2 + nlow1 + nlow2 +
                   hb0+
                   (1|identifier), 
                 data = df_ana_lmm_male[df_ana_lmm_male$visit>0,])
lmm_male
pvo <- data.frame(predicted=predict(lmm_male), 
                  observed=df_ana_lmm_male$hb[!is.na(df_ana_lmm_male$hb)&df_ana_lmm_male$visit>0])
threshold <- 13.5
pvo$agreement[pvo$predicted>=threshold&pvo$observed>=threshold] <- 'agree'
pvo$agreement[pvo$predicted<threshold&pvo$observed<threshold] <- 'agree'
pvo$agreement[pvo$predicted<threshold&pvo$observed>=threshold] <- 'p under, o over'
pvo$agreement[pvo$predicted>=threshold&pvo$observed<threshold] <- 'p over, o under'
table(pvo$agreement)

summary(lm(observed ~ predicted, data=pvo)) ##reasonable calibration

ggplot(pvo, aes(x=predicted, y=observed, col=agreement))+geom_point()+lims(x=c(10, 20),y=c(10,20))+
  geom_abline(intercept=0, slope=1)

##time to return model
nknots <- 4
fpmod_male <- flexsurv::flexsurvspline(Surv(timediff, return) ~ eth+exage+
                                    blood1+blood2+blood3+blood4+blood5+blood6+blood7+
                                    nlow1+nlow2+ndon1+ndon2+hb0
                                    , 
                                  data = df_ana_surv[df_ana_surv$sex==1&
                                                       df_ana_surv$timediff<24*7*(78-12),], k = nknots)


## MODEL

###########################################################################################################

## STRATEGY PARAMS:
## current: indicator for implementing current strategy (1=yes, 0=no)
## optimal: indicator for implementing alternative strategy
## atrisk: indicator for strategy with on-site testing for those with medium certainty of being over the threshold to donate
## indemand: indicator for strategy preferentially recalling donors with in-demand blood groups
## atrisk.thresh: threshold for those with on-site testing (only used if atrisk=1)
## indemand.thresh: threshold for recalling donors with in-demand blood groups (only used if indemand=1)
## lowhbdef.recall: wait time after low hb deferral
## vlowhbdef.recall: wait time after very low hb deferral
## othdef.recall: wait time after other deferral
## prob.thresh: threshold of probability over threshold to donate used for recalling 
## min.recall: earliest permitted recall time
## curr.recall: current earliest permitted recall time (only used if current=1)
## thresh: hb threshold for donation
## hb.min: min hb in observed data (modelled data rounded to this value if lower, for purpose of drawing BLUPs)
## hb.max: max hb in observed data (modelled data rounded to this value if higher, for purpose of drawing BLUPs)

## CURRENT STRATEGY 
## NOTE: with 12 week minimum recall

set.seed(1234)
params_male <- list("dropout"=0.073,##no need to compare number of dropouts
                    "lowhbdefer"=0.273,
                    'highhbdefer'=0,
                    "othdefer"=0.1396,##no problem
                    "strategy"=c(current=1,
                         lowhbdef.recall=12, vlowhbdef.recall=52,
                         othdef.recall=4, 
                         prob.thresh=0.9, min.recall=12, curr.recall=12,
                         length=52*1.5,
                         ##modify the threshold
                         thresh=13.5))
##Run the simulation
start <- Sys.time()
nsim <- 100
res <- data.frame(donation=NA, inapp=NA, lowhb=NA,
                  highhb=NA, lowhb_new=NA,
                  p_don=NA, p_inapp=NA, p_lowhb=NA,
                  p_highhb=NA, p_lowhb_new=NA)
for (i in 1:nsim) {
  res_current_male<-des_current(N=N, sex = 1, params=params_male, lmm=lmm_male, df_baseline = df_baseline,
                                recall=recall.curr, attend=attend.dist, fpmod=fpmod_male, logcumhaz=logcumhaz)
  #Events during the follow-up
  res[i,] <- events_desc(res_current_male$eventhistory, N=N)
  if (i%%20==0) {print(paste('nsim:', i))}
}
Sys.time()-start
##save results
write.csv(res, 'C:/Research/Cambridge/PDT pilot study/New simulation_CSvenous/res_male_current.csv', 
          row.names = F)
res <- read.csv('C:/Research/Cambridge/PDT pilot study/New simulation_CSvenous/res_male_current.csv')
##summary of the results
res_summary <- data.frame(matrix(NA, ncol=ncol(res)))
names(res_summary) <- c('don', 'inapp', 'lowhb', 
                        'highhb', 'lowhb_new',
                        'p_don', 'p_inapp', 'p_lowhb', 
                        'p_highhb', 'p_lowhb_new')
for (i in 1:ncol(res)) {
  res_summary[,i] <- paste0(round(mean(res[,i]), 1),
                            ' (', round(quantile(res[,i], 0.025), 1),
                            ', ', 
                            round(quantile(res[,i], 0.975), 1), ')')
}
res_summary
write.csv(res_summary, 
          'C:/Research/Cambridge/PDT pilot study/New simulation_CSvenous/summary_male_current.csv', 
          row.names = F)
