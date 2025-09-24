rm(list=ls())

##working dir
setwd('C:/Research/Cambridge/PDT pilot study/New simulation_CSvenous')

##Summary tables
summary_table <- function(sex, onsite, recall_low) {
  # if (sex=='male'&recall_low==12) {recall_interval <- 12}
  # if (sex=='male'&recall_low==26) {recall_interval <- c(12, 16, 20, 24)}
  # if (sex=='male'&recall_low==52) {recall_interval <- c(12, 32, 52)}
  # if (sex=='female'&recall_low==16) {recall_interval <- 16}
  # if (sex=='female'&recall_low==26) {recall_interval <- c(16, 20, 24)}
  # if (sex=='female'&recall_low==52) {recall_interval <- c(16, 32, 52)}
  if (sex=='male') {recall_interval <- 12}
  if (sex=='female') {recall_interval <- 16}
  
  ##Read data
  res_summary_all <- read.csv(paste0('res_', sex, '_onsite_hemocue_',onsite,'_', 
                                     recall_low,'w.csv'))
  nrow(res_summary_all)
  
  ##new percentages
  res_summary_all <- res_summary_all %>% mutate(total=don+inapp+lowhb+highhbdef, 
                             lowhb_new=lowhb+highhbdef,
                             p_don=don/total, p_inapp=inapp/total, p_lowhb_new=lowhb_new/total)
  
  ##Summary table
  summary_table <- data.frame(recall_interval=recall_interval)
  
  ##number of don
  summary_table$don <- paste0(
    round(as.vector(tapply(res_summary_all$don, res_summary_all$recall_interval, mean))),
    ' (', round(as.vector(tapply(res_summary_all$don, res_summary_all$recall_interval, 
                                 function(x)quantile(x, 0.025)))),
    ', ', round(as.vector(tapply(res_summary_all$don, res_summary_all$recall_interval, 
                                 function(x)quantile(x, 0.975)))), ')')
  ##number of inapp
  summary_table$inapp <- paste0(
    round(as.vector(tapply(res_summary_all$inapp, res_summary_all$recall_interval, mean))),
    ' (', round(as.vector(tapply(res_summary_all$inapp, res_summary_all$recall_interval, 
                                 function(x)quantile(x, 0.025)))),
    ', ', round(as.vector(tapply(res_summary_all$inapp, res_summary_all$recall_interval, 
                                 function(x)quantile(x, 0.975)))), ')')
  ##number of lowhb
  summary_table$lowhb_new <- paste0(
    round(as.vector(tapply(res_summary_all$lowhb_new, res_summary_all$recall_interval, mean))),
    ' (', round(as.vector(tapply(res_summary_all$lowhb_new, res_summary_all$recall_interval, 
                                 function(x)quantile(x, 0.025)))),
    ', ', round(as.vector(tapply(res_summary_all$lowhb_new, res_summary_all$recall_interval, 
                                 function(x)quantile(x, 0.975)))), ')')
  ##perc don
  summary_table$p_don <- paste0(
    round(as.vector(tapply(res_summary_all$p_don, res_summary_all$recall_interval, mean))*100, 1),
    ' (', round(as.vector(tapply(res_summary_all$p_don, res_summary_all$recall_interval, 
                                 function(x)quantile(x, 0.025)))*100, 1),
    ', ', round(as.vector(tapply(res_summary_all$p_don, res_summary_all$recall_interval, 
                                 function(x)quantile(x, 0.975)))*100, 1), ')')
  ##perc inapp
  summary_table$p_inapp <- paste0(
    round(as.vector(tapply(res_summary_all$p_inapp, res_summary_all$recall_interval, mean))*100, 1),
    ' (', round(as.vector(tapply(res_summary_all$p_inapp, res_summary_all$recall_interval, 
                                 function(x)quantile(x, 0.025)))*100, 1),
    ', ', round(as.vector(tapply(res_summary_all$p_inapp, res_summary_all$recall_interval, 
                                 function(x)quantile(x, 0.975)))*100, 1), ')')
  ##perc lowhb
  summary_table$p_lowhb_new <- paste0(
    round(as.vector(tapply(res_summary_all$p_lowhb_new, res_summary_all$recall_interval, mean))*100, 1),
    ' (', round(as.vector(tapply(res_summary_all$p_lowhb_new, res_summary_all$recall_interval, 
                                 function(x)quantile(x, 0.025)))*100, 1),
    ', ', round(as.vector(tapply(res_summary_all$p_lowhb_new, res_summary_all$recall_interval, 
                                 function(x)quantile(x, 0.975)))*100, 1), ')')
  
  write.csv(summary_table, paste0('summary_table_', sex, '_onsite_', onsite,'_',
                                  recall_low,'w.csv'), 
            row.names = F)
  
  summary_table
}
##onsite for medium group
summary_table('male', 'medium', 26)
summary_table('female', 'medium', 26)
