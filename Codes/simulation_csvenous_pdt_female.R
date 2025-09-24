rm(list=ls())

##Packages
library(MASS); library(flexsurv); library(simsurv); library(lme4); library(splines); library(ggplot2)
library(tidyr); library(dplyr)

setwd('C:/Research/Cambridge/Post-donation-Hb')
source("./Analysis/Simulation/onsite hemocue simulation/strides_donor_des_onsite_hemocue.R")
des.donor_onsite
des_onsite

#############
### WOMEN ###
#############
df_ana_lmm_female <- df_ana_lmm[df_ana_lmm$sex==2&df_ana_lmm$time<=52*1.5,]

##Spaghetti plot
ggplot(data = df_ana_lmm_female[df_ana_lmm_female$time>0,], 
       aes(x = time, y = hb, group = identifier)) + geom_line()
id_lowhb <- df_ana_lmm_female$identifier[df_ana_lmm_female$outcome==2]
id_lowhb <- id_lowhb[!duplicated(id_lowhb)]
ggplot(data = df_ana_lmm_female[df_ana_lmm_female$identifier%in%id_lowhb,], 
       aes(x = time, y = hb, group = identifier)) + geom_line()
ggplot(data = df_ana_lmm_female[!(df_ana_lmm_female$identifier%in%id_lowhb),], 
       aes(x = time, y = hb, group = identifier)) + geom_line()

##Predictors
plot(df_ana_lmm_female$exage, df_ana_lmm_female$hb)##linear
##time is almost linear
plot(df_ana_lmm_female$hb0, df_ana_lmm_female$hb)
plot(df_ana_lmm_female$eth, df_ana_lmm_female$hb)
##Model fitting
lmm_female <- lmer(hb ~ 
                     exage+
                     blood1+blood2+blood3+blood4+blood5+blood6+blood7+
                     eth+ time + 
                     # ndon1 + ndon2 + nlow1 + nlow2 +
                     hb0+
                     (1|identifier), 
                   data = df_ana_lmm_female[df_ana_lmm_female$visit>0,])
lmm_female
pvo <- data.frame(predicted=predict(lmm_female), 
                  observed=df_ana_lmm_female$hb[!is.na(df_ana_lmm_female$hb)&df_ana_lmm_female$visit>0])
threshold <- 12.5
pvo$agreement[pvo$predicted>=threshold&pvo$observed>=threshold] <- 'agree'
pvo$agreement[pvo$predicted<threshold&pvo$observed<threshold] <- 'agree'
pvo$agreement[pvo$predicted<threshold&pvo$observed>=threshold] <- 'p under, o over'
pvo$agreement[pvo$predicted>=threshold&pvo$observed<threshold] <- 'p over, o under'
table(pvo$agreement)

lm(observed ~ predicted, data=pvo)

ggplot(pvo, aes(x=observed, y=predicted, col=agreement))+geom_point()+lims(x=c(10, 20),y=c(10,20))+
  geom_abline(intercept=0, slope=1)

##time to return model
nknots <- 4
fpmod_female <- flexsurv::flexsurvspline(Surv(timediff, return) ~ eth+exage+
                                           blood1+blood2+blood3+blood4+blood5+blood6+blood7+
                                           nlow1+nlow2+ndon1+ndon2+hb0
                                         , 
                                         data = df_ana_surv[df_ana_surv$sex==2&
                                                              df_ana_surv$timediff<24*7*(78-16),], k = nknots)


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
## NOTE: with 16 week minimum recall
params_female <- list("dropout"=0.16,##no need to compare number of dropouts
                      "lowhbdefer"=0.375,
                      'highhbdefer'=0,
                      "othdefer"=0.1544,##no problem
                      "strategy"=c(lowhbdef.recall=12, vlowhbdef.recall=52,
                                   othdef.recall=4, 
                                   prob.thresh=0.9, min.recall=16, curr.recall=16,
                                   length=52*1.5, 
                                   ##modify the threshold
                                   thresh=12.5))

recall_low <- 26
# if (recall_low==52) {recall_interval <- seq(16, 52, 4)}
# if (recall_low==26) {recall_interval <- c(16, 20, 24)}
# if (recall_low==16) {recall_interval <- 16}
##recall interval for medium group
recall_interval <- 16
onsite <- 'medium'
res_summary0 <- list()
res_onsite_female <- list()
N <- 1000
ul_mediumhb <- 13.5 ##candidate values: 13.5/14
nsim <- 100

start <- Sys.time()
{
  for (i in 1:length(recall_interval)) {
    res_summary0_interval <- data.frame(recall_interval=rep(recall_interval[i], nsim),
                                        don=NA, inapp=NA, lowhb=NA, highhbdef=NA, otherdef=NA,
                                        p_don=NA, p_inapp=NA, p_lowhb=NA, p_highhbdef=NA, p_otherdef=NA)
    for (j in 1:nsim) {
      set.seed(seeds[j])
      res_onsite_female<-des_onsite(N=N, sex = 2, params=params_female, lmm=lmm_female,df_baseline = df_baseline,
                                    attend=attend.dist,fpmod=fpmod_female,logcumhaz=logcumhaz,
                                    recall_interval = c(recall_low, recall_interval[i], 16),
                                    recall_thres = c(12.5, ul_mediumhb),
                                    onsite=onsite)
      res_summary0_interval[j, -1] <- events_desc(res_onsite_female$eventhistory, N=N) 
      ##simulation progress
      if (j%%20==0) {print(paste0('Female, recall low: ', recall_low, 
                                  ', medium Hb recall interval: ', recall_interval[i], 
                                  ', onsite: ', onsite, ', nsim: ', j))}
    }
    res_summary0[[i]] <- res_summary0_interval
  }
}
Sys.time()-start
# res_summary0
write.csv(do.call(rbind, res_summary0), paste0('C:/Research/Cambridge/PDT pilot study/New simulation_CSvenous/res_female_onsite_', onsite, '_', 
                                               recall_low, 'w.csv'), row.names = F)
