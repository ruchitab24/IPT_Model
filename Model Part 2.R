#-- Load libraries -- 

library(plyr)
library(dplyr)
library(EnvStats)
library(ggplot2)

#load parameter list, set cohort size, set number of simulations 

parameter_list = read.table("~/Documents/JHU/IPT Model/parameter_list.txt", header = TRUE, fill = TRUE, row.names = NULL)cohort = 1000
ncohort = 1000
nsims = 1000

#create triangle distribution for each parameter value 
parameter_dist = sapply(1:nrow(parameter_list), function(i){
  rtri(nsims,min = parameter_list[i,3],max = parameter_list[i,4],mode = parameter_list[i,2])
})

#sample each parameter from triangle distribution for 1000 simulations 
colnames(parameter_dist) = parameter_list$Parameter

#initialize cohort 
cohort = matrix(nrow = ncohort, ncol = 6)

#baseline charachteristics
colnames(cohort) = c("CD4","TB","Xpert","Sx","CRP")


#Create list of variable names from this for setup 
#Create cohort for each simulation 

cohort_sims = sapply(1:nsims, function(j){
 
   #CD4 Status 
  
  cohort[,1] = sample(c('50','50_100','100_200','200_350','350'),size = ncohort,c(parameter_dist$CD41[j],parameter_dist$CD42[j],parameter_dist$CD43[j],parameter_dist$CD44[j],parameter_dist$CD45[j]),replace = TRUE)
  
  #TB Status (independent of CD4 count?)
  
  cohort[,2] = sample(c('TB','No_TB'),size = ncohort,c(parameter_dist$TB[j],1-parameter_dist$TB),replace = TRUE)
  
  #Xpert status 
  
  cohort_TB = cohort[which(cohort[,2] == "TB"),]
  cohort_TB[,3] = sample(c("Xpert+","Xpert-","no_sputum"), size = nrow(cohort_TB), c(parameter_dist$Xpert_Sens[j],1-parameter_dist$Xpert_Sens[j], parameter_dist$No_Sputum_TB[j]), replace = TRUE)
  
  
  cohort_No_TB = cohort[which(cohort[,2] == "No_TB"),]      
  cohort_No_TB[,3] = sample(c("Xpert+","Xpert-","no_sputum"), size = nrow(cohort_No_TB), c(parameter_dist$Xpert_Spec[j],1-parameter_dist$Xpert_Spec[j], parameter_dist$No_Sputum_No_TB[j]), replace = TRUE)
  
  #Symptom status, TB
  
  cohort_xp_pos = cohort_TB[which(cohort_TB[,3] == "Xpert+"),]
  cohort_xp_neg = cohort_TB[which(cohort_TB[,3] == "Xpert-"),]
  cohort_no_spu = cohort_TB[which(cohort_TB[,3] == "no_sputum"),]
  
  cohort_xp_pos[,4] = sample(c("Sx+","Sx-"), size = nrow(cohort_xp_pos), c(parameter_dist$Symp_TB_xp_pos[j],1-parameter_dist$Symp_TB_xp_pos[j]), replace = TRUE)
  cohort_xp_neg[,4] = sample(c("Sx+","Sx-"), size = nrow(cohort_xp_neg), c(parameter_dist$Symp_TB_xp_neg[j],1-parameter_dist$Symp_TB_xp_neg[j]), replace = TRUE)
  cohort_no_spu[,4] = sample(c("Sx+","Sx-"), size = nrow(cohort_no_spu), c(parameter_dist$Symp_TB_spu[j],1-parameter_dist$Symp_TB_spu[j]), replace = TRUE)
  
  #Symptom status, no TB 
  
  cohort_No_TB[,4] = sample(c("Sx+","Sx-"), size = nrow(cohort_No_TB), c(parameter_dist$Symp_No_TB[j],1-parameter_dist$Symp_No_TB[j]), replace = TRUE)
  
  
  #CRP Status, by Xpert status TB
  
  cohort_xp_pos[,5] = sample(c("CRP+","CRP-"), size = nrow(cohort_xp_pos), c(parameter_dist$CRP_TB_pos[j],1-parameter_dist$Symp_CRP_TB_pos[j]), replace = TRUE)
  cohort_xp_neg[,5] = sample(c("CRP+","CRP-"), size = nrow(cohort_xp_neg),c(parameter_dist$CRP_TB_neg[j],1-parameter_dist$Symp_CRP_TB_neg[j]), replace = TRUE)
  
  cohort_TB = rbind(cohort_xp_pos,cohort_xp_neg,cohort_no_spu)
  
  #CRP Status, no TB
  
  cohort_No_TB[,5] = sample(c("CRP+","CRP-"), size = nrow(cohort_No_TB), c(parameter_dist$CRP_Spec[j],1-parameter_dist$CRP_Spec[j]), replace = TRUE)
  
  #consolidate split dfs
  cohort = rbind(cohort_TB, cohort_No_TB)     
  
  return(cohort)
})













