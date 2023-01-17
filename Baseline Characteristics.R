#### -- Load libraries -- ####

library(plyr)
library(dplyr)
library(EnvStats)
library(ggplot2)
library(triangle)
library(binom)
library(parallel)
library(tidyr)

sourcelocation <- "~/Documents/JHU/IPT Model/"


#load parameter list, set cohort size, set number of simulations 

parameter_list = read.csv(paste0(sourcelocation, "parameter_list.csv"), header = FALSE, fill = TRUE, row.names = NULL)
ncohort = 1000
nsims = 1000

#create triangle distribution for each parameter value 
parameter_dist = sapply(1:nrow(parameter_list), function(i){
  rtriangle(nsims,a = parameter_list[i,3],b = parameter_list[i,4],c = parameter_list[i,2])
})


parameter_dist = as.data.frame(parameter_dist)
colnames(parameter_dist) = parameter_list[,1]

cohort_sims = vector("list", length = nsims) 

for(j in 1:nsims){
  
  #initialize cohort 
  cohort = matrix(nrow = ncohort, ncol = 5)
  cohort = as.data.frame(cohort)
  
  #baseline charachteristics
  colnames(cohort) = c("CD4","TB","Xpert","Sx","CRP")
  
  #CD4 Status 
  
  cohort[,"CD4"] = sample(c('100','100_200','200_350','350_500'), size = ncohort,
                          c(((parameter_dist$CD4_l100[j])*parameter_dist$CD4_l200[j]),
                            ((1-parameter_dist$CD4_l100[j])*parameter_dist$CD4_l200[j]),
                            ((1-parameter_dist$CD4_l200[j])*(1-parameter_dist$CD4_g350[j])),
                            ((1-parameter_dist$CD4_l200[j])*(parameter_dist$CD4_g350[j]))), replace = TRUE) 
  
  #TB Status 
  
  cohort[which(cohort[,"CD4"] == "100"),"TB"] = sample(c('TB','No_TB'), size = length(which(cohort[,"CD4"] == "100")),
                                                       c(parameter_dist$TB_CD41[j]*parameter_dist$TB_uncertainty[j],
                                                         1-parameter_dist$TB_CD41[j]*parameter_dist$TB_uncertainty[j]), replace = TRUE)
  
  cohort[which(cohort[,"CD4"] == "100_200"),"TB"] = sample(c('TB','No_TB'), size = length(which(cohort[,"CD4"] == "100_200")),
                                                           c(parameter_dist$TB_CD42[j]*parameter_dist$TB_uncertainty[j],
                                                             1-parameter_dist$TB_CD42[j]*parameter_dist$TB_uncertainty[j]), replace = TRUE)
  
  cohort[which(cohort[,"CD4"] == "200_350"),"TB"] = sample(c('TB','No_TB'), size = length(which(cohort[,"CD4"] == "200_350")),
                                                           c(parameter_dist$TB_CD43[j]*parameter_dist$TB_uncertainty[j],
                                                             1-parameter_dist$TB_CD43[j]*parameter_dist$TB_uncertainty[j]), replace = TRUE)
  
  cohort[which(cohort[,"CD4"] == "350_500"),"TB"] = sample(c('TB','No_TB'), size = length(which(cohort[,"CD4"] == "350_500")),
                                                           c(parameter_dist$TB_CD44[j]*parameter_dist$TB_uncertainty[j],
                                                             1-parameter_dist$TB_CD44[j]*parameter_dist$TB_uncertainty[j]), replace = TRUE)
  
  
  #Xpert status 
  
  cohort_TB = cohort[which(cohort[,"TB"] == "TB"),]
  cohort_TB[,"Xpert"] = sample(c("Xpert+","Xpert-","no_sputum"), size = nrow(cohort_TB), 
                               c((1-parameter_dist$No_Sputum_TB[j])*parameter_dist$Xpert_Sens[j],
                                 (1-parameter_dist$No_Sputum_TB[j])*(1-parameter_dist$Xpert_Sens[j]), 
                                 parameter_dist$No_Sputum_TB[j]), replace = TRUE)
  
  cohort_No_TB = cohort[which(cohort[,"TB"] == "No_TB"),]      
  cohort_No_TB[,"Xpert"] = sample(c("Xpert+","Xpert-","no_sputum"), size = nrow(cohort_No_TB), 
                                  c((1-parameter_dist$No_Sputum_No_TB[j])*(1-parameter_dist$Xpert_Spec[j]),
                                    (1-parameter_dist$No_Sputum_No_TB[j])*(parameter_dist$Xpert_Spec[j]), 
                                    parameter_dist$No_Sputum_No_TB[j]), replace = TRUE)
  
  #Symptom status, TB
  
  cohort_xp_pos = cohort_TB[which(cohort_TB[,"Xpert"] == "Xpert+"),]
  cohort_xp_neg = cohort_TB[which(cohort_TB[,"Xpert"] == "Xpert-"),]
  cohort_no_spu = cohort_TB[which(cohort_TB[,"Xpert"] == "no_sputum"),]
  
  cohort_xp_pos[,"Sx"] = sample(c("Sx+","Sx-"), size = nrow(cohort_xp_pos), 
                                c(parameter_dist$Symp_TB_xp_pos[j],
                                  1-parameter_dist$Symp_TB_xp_pos[j]), replace = TRUE)
  
  cohort_xp_neg[,"Sx"] = sample(c("Sx+","Sx-"), size = nrow(cohort_xp_neg), 
                                c(parameter_dist$Symp_TB_xp_neg[j],
                                  1-parameter_dist$Symp_TB_xp_neg[j]), replace = TRUE)
  
  cohort_no_spu[,"Sx"] = sample(c("Sx+","Sx-"), size = nrow(cohort_no_spu), 
                                c(parameter_dist$Symp_TB_spu[j],
                                  1-parameter_dist$Symp_TB_spu[j]), replace = TRUE)
  
  #Symptom status, no TB 
  
  cohort_No_TB[,"Sx"] = sample(c("Sx+","Sx-"), size = nrow(cohort_No_TB), 
                               c(parameter_dist$Symp_No_TB[j],
                                 1-parameter_dist$Symp_No_TB[j]), replace = TRUE)
  
  
  #CRP Status, by Xpert status TB
  
  cohort_xp_pos[,"CRP"] = sample(c("CRP+","CRP-"), size = nrow(cohort_xp_pos), 
                                 c(parameter_dist$CRP_TB_pos[j],
                                   1-parameter_dist$CRP_TB_pos[j]), replace = TRUE)
  
  cohort_xp_neg[,"CRP"] = sample(c("CRP+","CRP-"), size = nrow(cohort_xp_neg),
                                 c(parameter_dist$CRP_TB_neg[j],
                                   1-parameter_dist$CRP_TB_neg[j]), replace = TRUE)
  
  cohort_no_spu[,"CRP"] = sample(c("CRP+","CRP-"), size = nrow(cohort_no_spu),
                                 c(parameter_dist$CRP_TB_No_Spu[j],
                                   1-parameter_dist$CRP_TB_No_Spu[j]), replace = TRUE)
  
  cohort_TB = rbind(cohort_xp_pos,cohort_xp_neg,cohort_no_spu)
  
  #CRP Status, no TB
  
  cohort_No_TB[,"CRP"] = sample(c("CRP+","CRP-"), size = nrow(cohort_No_TB), 
                                c(1-parameter_dist$CRP_Spec[j],
                                  parameter_dist$CRP_Spec[j]), replace = TRUE)
  
  
  #CXR Status 
  
  #consolidate split dfs
  cohort = rbind(cohort_TB, cohort_No_TB)     
  
  cohort_sims[[j]] = cohort
}


return(cohort_sims)


