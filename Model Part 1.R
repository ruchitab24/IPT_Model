
#### Load libraries ####

library(tidyr)
library(ggplot2)
library(triangle)
library(ggpubr)

# make the source location modifiable, for running on different machines
if(missing(sourcelocation)) sourcelocation <- "~/Documents/JHU/IPT Model/"
# "~/Library/CloudStorage/OneDrive-JohnsHopkins/Research/TPT and inital Xpert in PWH model/Model/"

#source the baseline_cohort; 
#  creates cohort_sims (a list of nsims dataframes, each with ncohort observations of 10 patient characteristic variables)
source(paste0(sourcelocation,"Baseline Characteristics.R"))


#load parameter list for part 1 of model (all values must be as percentages)
#for primary analysis, use parameter_list_part1_main.csv, where ltf is 0. For secondary analysis, use original csv parameter_list_part1.csv where ltf is not 0. For tertiary analysis, assume ltf is 0 and there is no TPT discontinuation. 
parameter_list_part1 = read.csv(paste0(sourcelocation, "parameter_list_part1_main.csv"), header = FALSE, fill = TRUE, row.names = NULL)
#create triangle distribution for each parameter

parameter_dist_1 = sapply(1:nrow(parameter_list_part1), function(i){
  
  rtriangle(nsims,a = parameter_list_part1[i,3],b = parameter_list_part1[i,4],c = parameter_list_part1[i,2])
  
})



parameter_dist_1 = as.data.frame(parameter_dist_1)
colnames(parameter_dist_1) = parameter_list_part1[,1]

#library(rriskDistributions)
#pars <- get.gamma.par(p=c(0.28, 0.56), q=c(28/30, 90/30), plot = T, tol = 0.001)

visit2 = read.csv(paste0(sourcelocation, "visit2.csv"), header = FALSE, fill = TRUE, row.names = NULL)

#Assign visit 2 day to all cohorts, repeat for 1000 simulations 
t_value = vector("list", length = nsims)
for (a in 1:nsims){
  ltf = parameter_dist_1$ltf[j]/100 #sample the loss to follow up rate 
  t_prob = append(visit2[,2]*(1-ltf),ltf)
  t = append(visit2[,1],NA)
  
  t_value[[a]] = sample(t, size = ncohort,
             t_prob, replace = TRUE)
}


#Assign ART and TPT stop day to all cohorts, this will be modified as appropriate in later parts of the model  

ART_Rate = rtriangle(nsims,a = .01,b =.05,c = .035)

TPT_Rate = rtriangle(nsims,a = .04,b =.10,c = .07)

Stop_Days = vector("list", length = nsims)
for(a in 1:nsims){
  
  #modify such that everyone completes TPT 
  TPT_Stop_Day = round(rexp(ncohort, rate = (TPT_Rate[a]/30)),0)
  TPT_Stop_Day[which(TPT_Stop_Day>730)]<-730
  #TPT_Stop_Day = rep(730,ncohort)
  ART_Stop_Day = round(rexp(ncohort, rate = (ART_Rate[a]/30)),0)
  ART_Stop_Day[which(ART_Stop_Day>730)]<-730
  Stop_Days[[a]] = cbind(TPT_Stop_Day,ART_Stop_Day)
  colnames(Stop_Days[[a]]) = c("TPT_Stop_Day","ART_Stop_Day")
  
}



####Algorithm A####

outcome_cohort_A = vector("list", length = nsims)

for (j in 1:nsims){
  
  #select cohort for simulation 
  cohort = cohort_sims[[j]] 
  
  #initialize outcomes of cohort matrix
  outcome_cohort = matrix(nrow = ncohort, ncol = 10)
  
  #outcomes to track 
  colnames(outcome_cohort) = c("Sx_Screen","CRP_Screen","Xpert_Screen", # True if these tests are performed at visit 1
                               "TPT","Day_TPT", # TPT, treatment, and ART are true if these are ever started; Day_X is day of initiation
                               "Treatment","Day_Treatment",
                               "ART","Day_ART",
                               "t") # t is days to visit 2; will be NA if visit 2 does not occur 
  
  #Sample time here 
  
  #pars not working
  #returntimes <- rgamma(n = ncohort, shape = pars[1], rate =pars[2])*30
  #returntimes_observed <- returntimes
  #returntimes_observed[returntimes>90]<-NA
  #returntimes_observed[returntimes<7]<-7
  #summary(returntimes_observed)
  #hist(returntimes_observed, breaks = 1000)
  # sd(returntimes_observed, na.rm=T)
  
  


  
  outcome_cohort[,"t"] = t_value[[j]]
  
  #All recieve Symptom Screen, Sx- get ART at first visit
  outcome_cohort[which(cohort$Sx == "Sx-"),"ART"] = TRUE
  outcome_cohort[which(cohort$Sx == "Sx-"),"Day_ART"] = 0
  
  
  
  
  for(i in 1:ncohort){
    
    
    t = outcome_cohort[i,"t"] 

    #50% Sx+ get ART at first visit 
    if(cohort[i,"Sx"] == "Sx+"){
      
      if(runif(1)<.5){
        outcome_cohort[i,"ART"] = 
          outcome_cohort[i,"Day_ART"] = 0
      }
      
    }
    
    #Probability of symptom screen 
    if(runif(1)*100<parameter_dist_1$Pss[j]){
      outcome_cohort[i,"Sx_Screen"] = TRUE
      
      if(cohort[i,"Sx"] == "Sx+"){
        
        #Probability of confirmatory test
        if(runif(1)*100<parameter_dist_1$Pxc[j]){
          outcome_cohort[i,"Xpert_Screen"] = TRUE
          
          if(!is.na(t)){
            
            if(cohort[i,"Xpert"] == "Xpert+"){
              #Probability of starting treatment at visit 2
              if(runif(1)*100<parameter_dist_1$Pt2p[j]){
                outcome_cohort[i,"Treatment"] = TRUE
                outcome_cohort[i,"Day_Treatment"] = t
              }
              
            }
            else if(cohort[i,"Xpert"] == "Xpert-"){
              #Probability of considering TPT at visit 2
              if(runif(1)*100<((parameter_dist_1$Pp1n[j]/100)*(parameter_dist_1$RRp2n[j]/100)*100)){
                outcome_cohort[i,"TPT"] = TRUE
                outcome_cohort[i,"Day_TPT"] = t
              
              }
              
            }
            
          }
          
          
        }
        else{
       
          
          if(!is.na(t)){
            #Probability of considering TPT at visit 2 
            if(runif(1)*100<((parameter_dist_1$Pp1n[j]/100)*(parameter_dist_1$RRp2n[j]/100)*100)){
              outcome_cohort[i,"TPT"] = TRUE
              outcome_cohort[i,"Day_TPT"] = t 
              
            }
            
          }
        }
        
      }
      
      
      
      else if (cohort[i,"Sx"] == "Sx-"){
     
        
        #Probability of starting TPT
        if(runif(1)*100<parameter_dist_1$Pp1n[j]){
          outcome_cohort[i,"TPT"] = TRUE
          outcome_cohort[i,"Day_TPT"] = 0
          
          #Probability of continuing TPT at Visit 2 
           if(!is.na(t)){
             
             if(t>31){
               if(runif(1)*100>parameter_dist_1$Ptcp[j]){
                 outcome_cohort[i,"Day_TPT"] = "discontinued" 
               }
               
             }
             else{
               if(runif(1)*100<(parameter_dist_1$dt[j])){
                 outcome_cohort[i,"Day_TPT"] = "discontinued" 
               }
             }
           }
          
        }
        else{
          
          if(!is.na(t)){
            if(runif(1)*100<((parameter_dist_1$Pp1n[j]/100)*(parameter_dist_1$RRp2n[j]/100)*100)){
              outcome_cohort[i,"TPT"] = TRUE
              outcome_cohort[i,"Day_TPT"] = t 
              
            }
          }
        }
      }
      
      
    }
    
    else{

     
      
      if(!is.na(t)){
        
        #Probability of considering TPT at visit 2
        if(runif(1)*100<((parameter_dist_1$Pp1n[j]/100)*(parameter_dist_1$RRp2n[j]/100)*100)){
          outcome_cohort[i,"TPT"] = TRUE
          outcome_cohort[i,"Day_TPT"] = t
        }
      }
    }
    
    
    
  }
  
  outcome_cohort[is.na(outcome_cohort[,"TPT"]),"TPT"] <- 0
  outcome_cohort[is.na(outcome_cohort[,"Treatment"]),"Treatment"] <- 0
  outcome_cohort[is.na(outcome_cohort[,"ART"]),"ART"] <- 0
  outcome_cohort[is.na(outcome_cohort[,"Sx_Screen"]),"Sx_Screen"] <- 0
  outcome_cohort[is.na(outcome_cohort[,"Xpert_Screen"]),"Xpert_Screen"] <- 0
  outcome_cohort[is.na(outcome_cohort[,"CRP_Screen"]),"CRP_Screen"] <- 0
  
  outcome_cohort[which(outcome_cohort[,"TPT"] == TRUE),"TPT"] <- 1
  outcome_cohort[which(outcome_cohort[,"Treatment"] == TRUE),"Treatment"] <- 1
  outcome_cohort[which(outcome_cohort[,"ART"] == TRUE),"ART"] <- 1
  outcome_cohort[which(outcome_cohort[,"Sx_Screen"] == TRUE),"Sx_Screen"] <- 1
  outcome_cohort[which(outcome_cohort[,"Xpert_Screen"] == TRUE),"Xpert_Screen"] <- 1
  outcome_cohort[which(outcome_cohort[,"CRP_Screen"] == TRUE),"CRP_Screen"] <- 1
  
  #return outcomes for cohort j 
  outcome_cohort_A[[j]] = outcome_cohort
  outcome_cohort_A[[j]] = cbind(cohort_sims[[j]],outcome_cohort)
  outcome_cohort_A[[j]] = cbind(outcome_cohort_A[[j]],Stop_Days[[j]])
  
  
}


#General outcomes
outcome_discon_a = vector("list",length = nsims)
for(a in 1:nsims){
  x = outcome_cohort_A[[a]]
  x$Day_TPT[which(x$Day_TPT == "discontinued")] <- 0
  x$Day_TPT = as.numeric(as.character(x$Day_TPT))
  outcome_discon_a[[a]] = x
  
}

Treatment = do.call(rbind,lapply(outcome_cohort_A, function(x) sum(x[,"Treatment"]==1 & x[,"TB"]=="TB")/sum(x[,"TB"]=="TB"))) 
Treatment_60 = do.call(rbind,lapply(outcome_cohort_A, function(x) sum(x[,"Treatment"]==1 & as.numeric(as.character(x[,"Day_Treatment"]))<=60 & x[,"TB"]=="TB")/sum(x[,"TB"]=="TB")))
ART = do.call(rbind,lapply(outcome_cohort_A, function(x) sum(x[,"ART"]==1 & x[,"TB"]=="No_TB")/sum(x[,"TB"]=="No_TB")))
ART_Visit1 = do.call(rbind,lapply(outcome_cohort_A, function(x) sum(x[,"Day_ART"]==0 & x[,"TB"]=="No_TB")/sum(x[,"TB"]=="No_TB")))
ART_60 = do.call(rbind,lapply(outcome_cohort_A, function(x) sum(x[,"ART"]==1 & x[,"Day_ART"]<=60 & x[,"TB"]=="No_TB")/sum(x[,"TB"]=="No_TB")))
TPT = do.call(rbind,lapply(outcome_cohort_A, function(x) sum(x[,"TPT"]==1 & x[,"TB"]=="No_TB")/sum(x[,"TB"]=="No_TB")))
TPT_60 = do.call(rbind,lapply(outcome_discon_a, function(x) sum(x[,"TPT"]==1 & x[,"Day_TPT"]<=60 & x[,"TB"]=="No_TB")/sum(x[,"TB"]=="No_TB")))
Diagnosed = do.call(rbind,lapply(outcome_cohort_A, function(x) length(which(x[,"Xpert_Screen"]==1 & x[,"TB"]=="TB" & x[,"Xpert"] == "Xpert+"))/sum(x[,"TB"]=="TB")))
TB_ART_30 =  0#do.call(rbind,lapply(outcome_cohort_A, function(x) sum(x[,"TPT"]==1 & x[,"ART"]==1 & abs(x[,"Day_TPT"]- x[,"Day_ART"])<=30 & x[,"TB"]=="No_TB")/sum(x[,"TB"]=="No_TB")))
Treatment_Total_60 = do.call(rbind,lapply(outcome_cohort_A, function(x) sum(x[,"Treatment"]==1 & as.numeric(as.character(x[,"Day_Treatment"]))<=60 & x[,"TB"]=="TB")))
TPT_Total_60 = do.call(rbind,lapply(outcome_discon_a, function(x) sum(x[,"TPT"]==1 & x[,"Day_TPT"]<=60 & x[,"TB"]=="No_TB")))
Diagnosed_Total = do.call(rbind,lapply(outcome_cohort_A, function(x) length(which(x[,"Xpert_Screen"]==1 & x[,"TB"]=="TB" & x[,"Xpert"] == "Xpert+"))))
eval_outcomes_A = as.data.frame(cbind(Treatment,Treatment_60,ART,ART_60,TPT,TPT_60,Diagnosed,Treatment_Total_60 ,TPT_Total_60,Diagnosed_Total))
colnames(eval_outcomes_A) = c("Treatment for TB","Treatment <60d for TB","ART for No TB","ART <60d for No TB","TPT for No TB","TPT <60d for No TB","Diagnosed TB","Total Treatment","Total TPT","Total Diagnosed")

summary(do.call(rbind,lapply(outcome_cohort_D, function(x) sum(x[,"Treatment"]==1 & as.numeric(as.character(x[,"Day_Treatment"]))<=60 & x[,"TB"]=="TB"))) - (do.call(rbind,lapply(outcome_cohort_A, function(x) sum(x[,"Treatment"]==1 & as.numeric(as.character(x[,"Day_Treatment"]))<=60 & x[,"TB"]=="TB")))))
summary(do.call(rbind,lapply(outcome_cohort_D, function(x) sum(x[,"Treatment"]==1 & as.numeric(as.character(x[,"Day_Treatment"]))<=60 & x[,"TB"]=="TB"))) - (do.call(rbind,lapply(outcome_cohort_B, function(x) sum(x[,"Treatment"]==1 & as.numeric(as.character(x[,"Day_Treatment"]))<=60 & x[,"TB"]=="TB")))))
summary(do.call(rbind,lapply(outcome_cohort_D, function(x) sum(x[,"Treatment"]==1 & as.numeric(as.character(x[,"Day_Treatment"]))<=60 & x[,"TB"]=="TB"))) - (do.call(rbind,lapply(outcome_cohort_C, function(x) sum(x[,"Treatment"]==1 & as.numeric(as.character(x[,"Day_Treatment"]))<=60 & x[,"TB"]=="TB")))))
summary(do.call(rbind,lapply(outcome_cohort_C, function(x) sum(x[,"Treatment"]==1 & as.numeric(as.character(x[,"Day_Treatment"]))<=60 & x[,"TB"]=="TB"))) - (do.call(rbind,lapply(outcome_cohort_A, function(x) sum(x[,"Treatment"]==1 & as.numeric(as.character(x[,"Day_Treatment"]))<=60 & x[,"TB"]=="TB")))))


summary(do.call(rbind,lapply(outcome_discon_d, function(x) sum(x[,"TPT"]==1 & x[,"Day_TPT"]<=60 & x[,"TB"]=="No_TB"))) - do.call(rbind,lapply(outcome_discon_a, function(x) sum(x[,"TPT"]==1 & x[,"Day_TPT"]<=60 & x[,"TB"]=="No_TB"))))
summary(do.call(rbind,lapply(outcome_discon_d, function(x) sum(x[,"TPT"]==1 & x[,"Day_TPT"]<=60 & x[,"TB"]=="No_TB"))) - do.call(rbind,lapply(outcome_discon_b, function(x) sum(x[,"TPT"]==1 & x[,"Day_TPT"]<=60 & x[,"TB"]=="No_TB"))))
summary(do.call(rbind,lapply(outcome_discon_d, function(x) sum(x[,"TPT"]==1 & x[,"Day_TPT"]<=60 & x[,"TB"]=="No_TB"))) - do.call(rbind,lapply(outcome_discon_c, function(x) sum(x[,"TPT"]==1 & x[,"Day_TPT"]<=60 & x[,"TB"]=="No_TB"))))
summary(do.call(rbind,lapply(outcome_discon_c, function(x) sum(x[,"TPT"]==1 & x[,"Day_TPT"]<=60 & x[,"TB"]=="No_TB"))) - do.call(rbind,lapply(outcome_discon_b, function(x) sum(x[,"TPT"]==1 & x[,"Day_TPT"]<=60 & x[,"TB"]=="No_TB"))))
#Diagnostic evaluations 

symp = do.call(rbind, lapply(outcome_cohort_A, function(x) length(which(x[,"Sx_Screen"] == 1))))
CRP = do.call(rbind, lapply(outcome_cohort_A, function(x) length(which(x[,"CRP_Screen"] == 1))))
Xperts = do.call(rbind, lapply(outcome_cohort_A, function(x) length(which(x[,"Xpert_Screen"] == 1))))

Diagnostic_eval_A = as.data.frame(cbind(symp,CRP, Xperts))
colnames(Diagnostic_eval_A) =  c("Symptom Screen","CRP Screen","Xpert Screen")


#Unecessary or inappropriate therapy
Wrong_Treatment = do.call(rbind,lapply(outcome_cohort_A, function(x) sum(x[,"Treatment"]==1 & x[,"TB"]=="No_TB")/sum(x[,"TB"]=="No_TB")) ) 
Wrong_TPT = do.call(rbind,lapply(outcome_cohort_A, function(x) sum(x[,"TPT"]==1 & x[,"TB"]=="TB")/sum(x[,"TB"]=="TB")) ) 
Wrong_Treatment_Total =do.call(rbind,lapply(outcome_cohort_A, function(x) sum(x[,"Treatment"]==1 & x[,"TB"]=="No_TB")) ) 
Wrong_TPT_Total = do.call(rbind,lapply(outcome_cohort_A, function(x) sum(x[,"TPT"]==1 & x[,"TB"]=="TB")) ) 
Inappropriate_Therapy_A = as.data.frame(cbind(Wrong_Treatment,Wrong_TPT,Wrong_Treatment_Total,Wrong_TPT_Total))
colnames(Inappropriate_Therapy_A) = c("Treatment Incorrectly Given Without TB","TPT Incorrectly Given With TB","Total Incorrect Treatment","Total Incorrect TPT")

TPT_30_A = do.call(rbind, lapply(outcome_cohort_A, function(x) sum(x[,"TPT"]==1 & x[,"TB"]=="TB" & x[,"Day_TPT"] != "discontinued" & as.numeric(as.character(x[,"TPT_Stop_Day"])) > 30)/sum(x[,"TB"]=="TB")))

#Plot values

#Evaluation outcomes
Eval_A = ggplot(gather(eval_outcomes_A, cols, value), aes(x = value)) + 
  geom_histogram() + xlim(0,1.1) + ylim(0,1000) + facet_grid(.~cols)+                                                                
  theme(strip.text.x = element_text(size = 12))+
  theme(axis.text.x = element_text(size = 10))

#Diagnostic evaluations 
Diag_A = ggplot(gather(Diagnostic_eval_A, cols, value), aes(x = value)) + 
  geom_histogram() + ylim(0,1000) + facet_grid(.~cols)+                                                              
  theme(strip.text.x = element_text(size = 12))+
  theme(axis.text.x = element_text(size = 10))

#Unecessary or inappropriate therapy
Unecessary_A = ggplot(gather(Inappropriate_Therapy_A, cols, value), aes(x = value)) + 
  geom_histogram() +xlim(0,1) + ylim(0,1000) + facet_grid(.~cols)+                                                                
  theme(strip.text.x = element_text(size = 12))+
  theme(axis.text.x = element_text(size = 10))


####Algorithm A2####

# outcome_cohort_A2 = vector("list", length = nsims)
# 
# for (j in 1:nsims){
#   
#   #select cohort for simulation 
#   cohort = cohort_sims[[j]] 
#   
#   #initialize outcomes of cohort matrix
#   outcome_cohort = matrix(nrow = ncohort, ncol = 11)
#   
#   #outcomes to track 
#   colnames(outcome_cohort) = c("Sx_Screen","CRP_Screen","Xpert_Screen","CXR_Screen", # True if these tests are performed at visit 1
#                                "TPT","Day_TPT", # TPT, treatment, and ART are true if these are ever started; Day_X is day of initiation
#                                "Treatment","Day_Treatment",
#                                "ART","Day_ART",
#                                "t") # t is days to visit 2; will be NA if visit 2 does not occur 
#   
#   
#   
#    outcome_cohort[,"t"] = t_value[[j]]

          #All recieve Symptom Screen, Sx- get ART at first visit
          #outcome_cohort[which(cohort$Sx == "Sx-"),"ART"] = TRUE
          #outcome_cohort[which(cohort$Sx == "Sx-"),"Day_ART"] = 0

#   
#   for(i in 1:ncohort){
#     
#     
#     t = outcome_cohort[i,"t"] 
#     

#50% Sx+ get ART at first visit 
#if(cohort[i,"Sx"] == "Sx+"){
  
#  if(runif(1)<.5){
#    outcome_cohort[i,"ART"] = 
#      outcome_cohort[i,"Day_ART"] = 0
#  }
  
#}
#     
#     #Probability of symptom screen 
#     if(runif(1)*100<parameter_dist_1$Pscxr[j]){
#       outcome_cohort[i,"Sx_Screen"] = TRUE
#       outcome_cohort[i,"CXR_Screen"] = TRUE
#       
#       if(cohort[i,"Sx"] == "Sx+" || cohort[i,"CXR+"] == "CXR+"){
#         
#         #Probability of confirmatory test
#         if(runif(1)*100<parameter_dist_1$Pxc[j]){
#           outcome_cohort[i,"Xpert_Screen"] = TRUE
#           
#           if(!is.na(t)){
#             
#             if(cohort[i,"Xpert"] == "Xpert+"){
#               #Probability of starting treatment at visit 2
#               if(runif(1)*100<parameter_dist_1$Pt2p[j]){
#                 outcome_cohort[i,"Treatment"] = TRUE
#                 outcome_cohort[i,"Day_Treatment"] = t
#               }
#               
#             }
#             else if(cohort[i,"Xpert"] == "Xpert-"){
#               #Probability of considering TPT at visit 2
#               if(runif(1)*100<((parameter_dist_1$Pp1n[j]/100)*(parameter_dist_1$RRp2n[j]/100)*100)){
#                 outcome_cohort[i,"TPT"] = TRUE
#                 outcome_cohort[i,"Day_TPT"] = t
#                 outcome_cohort[i,"ART"] = TRUE
#                 outcome_cohort[i,"Day_ART"] = t
#               }
#               
#             }
#             
#           }
#           
#           
#         }
#         else{
#           outcome_cohort[i,"ART"] = TRUE 
#           outcome_cohort[i,"Day_ART"] = 0
#           
#           if(!is.na(t)){
#             #Probability of considering TPT at visit 2 
#             if(runif(1)*100<((parameter_dist_1$Pp1n[j]/100)*(parameter_dist_1$RRp2n[j]/100)*100)){
#               outcome_cohort[i,"TPT"] = TRUE
#               outcome_cohort[i,"Day_TPT"] = t 
#               
#             }
#             
#           }
#         }
#         
#       }
#       
#       
#       
#       else{
#         
#         outcome_cohort[i,"ART"] = TRUE
#         outcome_cohort[i,"Day_ART"] = 0
#         
#         #Probability of starting TPT
#         if(runif(1)*100<parameter_dist_1$Pp1n[j]){
#           outcome_cohort[i,"TPT"] = TRUE
#           outcome_cohort[i,"Day_TPT"] = 0
#         }
#         else{
#           
#           if(!is.na(t)){
#             if(runif(1)*100<((parameter_dist_1$Pp1n[j]/100)*(parameter_dist_1$RRp2n[j]/100)*100)){
#               outcome_cohort[i,"TPT"] = TRUE
#               outcome_cohort[i,"Day_TPT"] = t 
#               
#             }
#           }
#         }
#       }
#       
#       
#     }
#     
#     else{
#       outcome_cohort[i,"Sx_Screen"] = FALSE
#       outcome_cohort[i,"CXR_Screen"] = FALSE
#       outcome_cohort[i,"ART"] = TRUE
#       outcome_cohort[i,"Day_ART"] = 0
#       
#       if(!is.na(t)){
#         
#         #Probability of considering TPT at visit 2
#         if(runif(1)*100<((parameter_dist_1$Pp1n[j]/100)*(parameter_dist_1$RRp2n[j]/100)*100)){
#           outcome_cohort[i,"TPT"] = TRUE
#           outcome_cohort[i,"Day_TPT"] = t
#         }
#       }
#     }
#     
#     
#     
#   }
#   
#   outcome_cohort[is.na(outcome_cohort[,"TPT"]),"TPT"] <- 0
#   outcome_cohort[is.na(outcome_cohort[,"Treatment"]),"Treatment"] <- 0
#   outcome_cohort[is.na(outcome_cohort[,"ART"]),"ART"] <- 0
#   outcome_cohort[is.na(outcome_cohort[,"Sx_Screen"]),"Sx_Screen"] <- 0
#   outcome_cohort[is.na(outcome_cohort[,"Xpert_Screen"]),"Xpert_Screen"] <- 0
#   outcome_cohort[is.na(outcome_cohort[,"CRP_Screen"]),"CRP_Screen"] <- 0
#   
#   outcome_cohort[which(outcome_cohort[,"TPT"] == TRUE),"TPT"] <- 1
#   outcome_cohort[which(outcome_cohort[,"Treatment"] == TRUE),"Treatment"] <- 1
#   outcome_cohort[which(outcome_cohort[,"ART"] == TRUE),"ART"] <- 1
#   outcome_cohort[which(outcome_cohort[,"Sx_Screen"] == TRUE),"Sx_Screen"] <- 1
#   outcome_cohort[which(outcome_cohort[,"Xpert_Screen"] == TRUE),"Xpert_Screen"] <- 1
#   outcome_cohort[which(outcome_cohort[,"CRP_Screen"] == TRUE),"CRP_Screen"] <- 1
#   
#   #return outcomes for cohort j 
#   outcome_cohort_A2[[j]] = outcome_cohort
#   outcome_cohort_A2[[j]] = cbind(cohort_sims[[j]],outcome_cohort)
#   outcome_cohort_A2[[j]] = cbind(outcome_cohort_A2[[j]],Stop_Days[[j]])
#   
#   
#   
# }
# 
# 
# #General outcomes
# 
# Treatment = do.call(rbind,lapply(outcome_cohort_A2, function(x) sum(x[,"Treatment"]==1 & x[,"TB"]=="TB")/sum(x[,"TB"]=="TB")) ) 
# Treatment_60 = do.call(rbind,lapply(outcome_cohort_A2, function(x) sum(x[,"Treatment"]==1 & x[,"Day_Treatment"]<=60 & x[,"TB"]=="TB")/sum(x[,"TB"]=="TB")))
# ART = do.call(rbind,lapply(outcome_cohort_A2, function(x) sum(x[,"ART"]==1 & x[,"TB"]=="No_TB")/sum(x[,"TB"]=="No_TB")))
# ART_Visit1 = do.call(rbind,lapply(outcome_cohort_A2, function(x) sum(x[,"Day_ART"]==0 & x[,"TB"]=="No_TB")/sum(x[,"TB"]=="No_TB")))
# ART_60 = do.call(rbind,lapply(outcome_cohort_A2, function(x) sum(x[,"ART"]==1 & x[,"Day_ART"]<=60 & x[,"TB"]=="No_TB")/sum(x[,"TB"]=="No_TB")))
# TPT = do.call(rbind,lapply(outcome_cohort_A2, function(x) sum(x[,"TPT"]==1 & x[,"TB"]=="No_TB")/sum(x[,"TB"]=="No_TB")))
# TPT_60 = do.call(rbind,lapply(outcome_cohort_A2, function(x) sum(x[,"TPT"]==1 & x[,"Day_TPT"]<=60 & x[,"TB"]=="No_TB")/sum(x[,"TB"]=="No_TB")))
# Diagnosed = do.call(rbind,lapply(outcome_cohort_A2, function(x) length(which(x[,"Xpert_Screen"]==1 & x[,"TB"]=="TB" & x[,"Xpert"] == "Xpert+"))/sum(x[,"TB"]=="TB")))
# TB_ART_30 =  0#do.call(rbind,lapply(outcome_cohort_A, function(x) sum(x[,"TPT"]==1 & x[,"ART"]==1 & abs(x[,"Day_TPT"]- x[,"Day_ART"])<=30 & x[,"TB"]=="No_TB")/sum(x[,"TB"]=="No_TB")))
# 
# eval_outcomes_A2 = as.data.frame(cbind(Treatment,Treatment_60,ART,ART_60,TPT,TPT_60,Diagnosed))
# colnames(eval_outcomes_A2) = c("Treatment for TB","Treatment <60d for TB","ART for No TB","ART <60d for No TB","TPT for No TB","TPT <60d for No TB","Diagnosed TB")
# 
# 
# #Diagnostic evaluations 
# 
# symp = do.call(rbind, lapply(outcome_cohort_A2, function(x) length(which(x[,"Sx_Screen"] == 1))))
# CRP = do.call(rbind, lapply(outcome_cohort_A2, function(x) length(which(x[,"CRP_Screen"] == 1))))
# Xperts = do.call(rbind, lapply(outcome_cohort_A2, function(x) length(which(x[,"Xpert_Screen"] == 1))))
# CXR = do.call(rbind, lapply(outcome_cohort_A2, function(x) length(which(x[,"CXR_Screen"] == 1))))
# 
# Diagnostic_eval_A2 = as.data.frame(cbind(symp,CRP, Xperts,CXR))
# colnames(Diagnostic_eval_A2) =  c("Symptom Screen","CRP Screen","Xpert Screen","CXR Screen")
# 
# 
# #Unecessary or inappropriate therapy
# Wrong_Treatment = do.call(rbind,lapply(outcome_cohort_A2, function(x) sum(x[,"Treatment"]==1 & x[,"TB"]=="No_TB")/sum(x[,"TB"]=="No_TB")) ) 
# Wrong_TPT = do.call(rbind,lapply(outcome_cohort_A2, function(x) sum(x[,"TPT"]==1 & x[,"TB"]=="TB")/sum(x[,"TB"]=="TB")) ) 
# 
# Inappropriate_Therapy_A2 = as.data.frame(cbind(Wrong_Treatment,Wrong_TPT))
# colnames(Inappropriate_Therapy_A2) = c("Treatment Incorrectly Given Without TB","TPT Incorrectly Given With TB")
# 

#Plot values

# #Evaluation outcomes
# Eval_A2 = ggplot(gather(eval_outcomes_A2, cols, value), aes(x = value)) + 
#   geom_histogram() + xlim(0,1.1) + ylim(0,1000) + facet_grid(.~cols)+                                                               
#   theme(strip.text.x = element_text(size = 12))+
#   theme(axis.text.x = element_text(size = 10))
# 
# #Diagnostic evaluations 
# Diag_A2 = ggplot(gather(Diagnostic_eval_A2, cols, value), aes(x = value)) + 
#   geom_histogram() + ylim(0,1000) + facet_grid(.~cols)+                                                                
#   theme(strip.text.x = element_text(size = 12))+
#   theme(axis.text.x = element_text(size = 10))
# 
# #Unecessary or inappropriate therapy
# Unecessary_A2 = ggplot(gather(Inappropriate_Therapy_A2, cols, value), aes(x = value)) + 
#   geom_histogram() +xlim(0,1) + ylim(0,1000) + facet_grid(.~cols)+                                                                
#   theme(strip.text.x = element_text(size = 12))+
#   theme(axis.text.x = element_text(size = 10))



####Algorithm B####

outcome_cohort_B = vector("list", length = nsims)

for(j in 1:nsims){
  #select cohort for simulation 
  cohort = cohort_sims[[j]] 
  
  #initialize outcomes of cohort matrix
  outcome_cohort = matrix(nrow = ncohort, ncol = 10)
  
  #outcomes to track 
  colnames(outcome_cohort) = c("Sx_Screen","CRP_Screen","Xpert_Screen","TPT","Day_TPT","Treatment","Day_Treatment","ART","Day_ART","t")

  
  outcome_cohort[,"t"] = t_value[[j]]
  
  for(i in 1:ncohort){
    #sample visit2
    t = outcome_cohort[i,"t"]
    

    #Probability of CRP screen 
    if(runif(1)*100<(parameter_dist_1$Pss[j]/100)*(parameter_dist_1$RRcs[j]/100)*100){
      outcome_cohort[i,"CRP_Screen"] = TRUE
      
      if(cohort[i,"CRP"] == "CRP+"){
        
        
        
        #Probability of confirmatory test
        if(runif(1)*100<parameter_dist_1$Pxc[j]){
          outcome_cohort[i,"Xpert_Screen"] = TRUE
          
          if(!is.na(t)){
            
            if(cohort[i,"Xpert"] == "Xpert+"){
              #Probability of starting treatment at visit 2
              if(runif(1)*100<parameter_dist_1$Pt2p[j]){
                outcome_cohort[i,"Treatment"] = TRUE
                outcome_cohort[i,"Day_Treatment"] = t
              }
            }
            else if(cohort[i,"Xpert"] == "Xpert-"){
              #Probability of considering TPT at visit 2
              if(runif(1)*100<((parameter_dist_1$Pp1n[j]/100)*(parameter_dist_1$RRp2n[j]/100)*100)){
                outcome_cohort[i,"TPT"] = TRUE
                outcome_cohort[i,"Day_TPT"] = t
                
              }
            }
          }
          
        }
        else{
         
          
          if(!is.na(t)){
            #Probability of considering TPT at visit 2 
            if(runif(1)*100<((parameter_dist_1$Pp1n[j]/100)*(parameter_dist_1$RRp2n[j]/100)*100)){
              outcome_cohort[i,"TPT"] = TRUE
              outcome_cohort[i,"Day_TPT"] = t 
            }
          }
          
          
        }
      }
      else if (cohort[i,"CRP"] == "CRP-"){
        
        outcome_cohort[i,"ART"] = TRUE
        outcome_cohort[i,"Day_ART"] = 0
        
        #Probability of starting TPT
        if(runif(1)*100<parameter_dist_1$Pp1n[j]){
          outcome_cohort[i,"TPT"] = TRUE
          outcome_cohort[i,"Day_TPT"] = 0
          
          #Probability of continuing TPT at Visit 2 
          if(!is.na(t)){
            
            if(t>31){
              if(runif(1)*100>parameter_dist_1$Ptcp[j]){
                outcome_cohort[i,"Day_TPT"] = "discontinued" 
              }
              
            }
            else{
              if(runif(1)*100<(parameter_dist_1$dt[j])){
                outcome_cohort[i,"Day_TPT"] = "discontinued" 
              }
            }
          }
        }
        else{
          
          if(!is.na(t)){
            #Probability of considering TPT at visit 2 
            if(runif(1)*100<((parameter_dist_1$Pp1n[j]/100)*(parameter_dist_1$RRp2n[j]/100)*100)){
              outcome_cohort[i,"TPT"] = TRUE
              outcome_cohort[i,"Day_TPT"] = t 
            }
          }
        }
        
      }
      
      
      else{
        outcome_cohort[i,"CRP_Screen"] = FALSE
       
        if(!is.na(t)){
          
          #Probability of considering TPT at visit 2
          if(runif(1)*100<((parameter_dist_1$Pp1n[j]/100)*(parameter_dist_1$RRp2n[j]/100)*100)){
            outcome_cohort[i,"TPT"] = TRUE
            outcome_cohort[i,"Day_TPT"] = t
          }
        }
      }
      
    }
  }
  
  outcome_cohort[is.na(outcome_cohort[,"TPT"]),"TPT"] <- 0
  outcome_cohort[is.na(outcome_cohort[,"Treatment"]),"Treatment"] <- 0
  outcome_cohort[is.na(outcome_cohort[,"ART"]),"ART"] <- 0
  outcome_cohort[is.na(outcome_cohort[,"Sx_Screen"]),"Sx_Screen"] <- 0
  outcome_cohort[is.na(outcome_cohort[,"Xpert_Screen"]),"Xpert_Screen"] <- 0
  outcome_cohort[is.na(outcome_cohort[,"CRP_Screen"]),"CRP_Screen"] <- 0
  
  outcome_cohort[which(outcome_cohort[,"TPT"] == TRUE),"TPT"] <- 1
  outcome_cohort[which(outcome_cohort[,"Treatment"] == TRUE),"Treatment"] <- 1
  outcome_cohort[which(outcome_cohort[,"ART"] == TRUE),"ART"] <- 1
  outcome_cohort[which(outcome_cohort[,"Sx_Screen"] == TRUE),"Sx_Screen"] <- 1
  outcome_cohort[which(outcome_cohort[,"Xpert_Screen"] == TRUE),"Xpert_Screen"] <- 1
  outcome_cohort[which(outcome_cohort[,"CRP_Screen"] == TRUE),"CRP_Screen"] <- 1
  
  #return outcomes for cohort j 
  outcome_cohort_B[[j]] = outcome_cohort
  outcome_cohort_B[[j]] = cbind(cohort_sims[[j]],outcome_cohort)
  outcome_cohort_B[[j]] = cbind(outcome_cohort_B[[j]],Stop_Days[[j]])
  
}


#General outcomes
outcome_discon_b = vector("list",length = nsims)
for(a in 1:nsims){
  x = outcome_cohort_B[[a]]
  x$Day_TPT[which(x$Day_TPT == "discontinued")] <- 0
  x$Day_TPT = as.numeric(as.character(x$Day_TPT))
  outcome_discon_b[[a]] = x
  
}

Treatment = do.call(rbind,lapply(outcome_cohort_B, function(x) sum(x[,"Treatment"]==1 & x[,"TB"]=="TB")/sum(x[,"TB"]=="TB")) ) 
Treatment_60 = do.call(rbind,lapply(outcome_cohort_B, function(x) sum(x[,"Treatment"]==1 & as.numeric(as.character(x[,"Day_Treatment"]))<=60 & x[,"TB"]=="TB")/sum(x[,"TB"]=="TB")))
ART = do.call(rbind,lapply(outcome_cohort_B, function(x) sum(x[,"ART"]==1 & x[,"TB"]=="No_TB")/sum(x[,"TB"]=="No_TB")))
ART_Visit1 = do.call(rbind,lapply(outcome_cohort_B, function(x) sum(x[,"Day_ART"]==0 & x[,"TB"]=="No_TB")/sum(x[,"TB"]=="No_TB")))
ART_60 = do.call(rbind,lapply(outcome_cohort_B, function(x) sum(x[,"ART"]==1 & x[,"Day_ART"]<=60 & x[,"TB"]=="No_TB")/sum(x[,"TB"]=="No_TB")))
TPT = do.call(rbind,lapply(outcome_cohort_B, function(x) sum(x[,"TPT"]==1 & x[,"TB"]=="No_TB")/sum(x[,"TB"]=="No_TB")))
TPT_60 = do.call(rbind,lapply(outcome_discon_b, function(x) sum(x[,"TPT"]==1 & x[,"Day_TPT"]<=60 & x[,"TB"]=="No_TB")/sum(x[,"TB"]=="No_TB")))
Diagnosed = do.call(rbind,lapply(outcome_cohort_B, function(x) length(which(x[,"Xpert_Screen"]==1 & x[,"TB"]=="TB" & x[,"Xpert"] == "Xpert+"))/sum(x[,"TB"]=="TB")))
TB_ART_30 =  0#do.call(rbind,lapply(outcome_cohort_B, function(x) sum(x[,"TPT"]==1 & x[,"ART"]==1 & abs(x[,"Day_TPT"]- x[,"Day_ART"])<=30 & x[,"TB"]=="No_TB")/sum(x[,"TB"]=="No_TB")))
Treatment_Total_60 = do.call(rbind,lapply(outcome_cohort_B, function(x) sum(x[,"Treatment"]==1 & as.numeric(as.character(x[,"Day_Treatment"]))<=60 & x[,"TB"]=="TB")))
TPT_Total_60 = do.call(rbind,lapply(outcome_discon_b, function(x) sum(x[,"TPT"]==1 & x[,"Day_TPT"]<=60 & x[,"TB"]=="No_TB")))
Diagnosed_Total = do.call(rbind,lapply(outcome_cohort_B, function(x) length(which(x[,"Xpert_Screen"]==1 & x[,"TB"]=="TB" & x[,"Xpert"] == "Xpert+"))))
eval_outcomes_B = as.data.frame(cbind(Treatment,Treatment_60,ART,ART_60,TPT,TPT_60,Diagnosed,Treatment_Total_60 ,TPT_Total_60,Diagnosed_Total))
colnames(eval_outcomes_B) = c("Treatment for TB","Treatment <60d for TB","ART for No TB","ART <60d for No TB","TPT for No TB","TPT <60d for No TB","Diagnosed TB","Total Treatment","Total TPT","Total Diagnosed")


#Diagnostic evaluations 

symp = do.call(rbind, lapply(outcome_cohort_B, function(x) length(which(x[,"Sx_Screen"] == 1))))
CRP = do.call(rbind, lapply(outcome_cohort_B, function(x) length(which(x[,"CRP_Screen"] == 1))))
Xperts = do.call(rbind, lapply(outcome_cohort_B, function(x) length(which(x[,"Xpert_Screen"] == 1))))

Diagnostic_eval_B = as.data.frame(cbind(symp,CRP, Xperts))
colnames(Diagnostic_eval_B) =  c("Symptom Screen","CRP Screen","Xpert Screen")


#Unecessary or inappropriate therapy
Wrong_Treatment = do.call(rbind,lapply(outcome_cohort_B, function(x) sum(x[,"Treatment"]==1 & x[,"TB"]=="No_TB")/sum(x[,"TB"]=="No_TB")) ) 
Wrong_TPT = do.call(rbind,lapply(outcome_cohort_B, function(x) sum(x[,"TPT"]==1 & x[,"TB"]=="TB")/sum(x[,"TB"]=="TB")) ) 
Wrong_Treatment_Total =do.call(rbind,lapply(outcome_cohort_B, function(x) sum(x[,"Treatment"]==1 & x[,"TB"]=="No_TB")) ) 
Wrong_TPT_Total = do.call(rbind,lapply(outcome_cohort_B, function(x) sum(x[,"TPT"]==1 & x[,"TB"]=="TB")) ) 
Inappropriate_Therapy_B = as.data.frame(cbind(Wrong_Treatment,Wrong_TPT,Wrong_Treatment_Total,Wrong_TPT_Total))
colnames(Inappropriate_Therapy_B) = c("Treatment Incorrectly Given Without TB","TPT Incorrectly Given With TB","Total Incorrect Treatment","Total Incorrect TPT")

TPT_30_B = do.call(rbind, lapply(outcome_cohort_B, function(x) sum(x[,"TPT"]==1 & x[,"TB"]=="TB" & x[,"Day_TPT"] != "discontinued" & as.numeric(as.character(x[,"TPT_Stop_Day"])) > 30)/sum(x[,"TB"]=="TB")))
#Plot values

#Evaluation outcomes
Eval_B = ggplot(gather(eval_outcomes_B, cols, value), aes(x = value)) + 
  geom_histogram() + xlim(0,1.1) + ylim(0,1000) + facet_grid(.~cols)+                                                                
  theme(strip.text.x = element_text(size = 12))+
  theme(axis.text.x = element_text(size = 10))

#Diagnostic evaluations 
Diag_B = ggplot(gather(Diagnostic_eval_B, cols, value), aes(x = value)) + 
  geom_histogram() + ylim(0,1000) + facet_grid(.~cols)+                                                                
  theme(strip.text.x = element_text(size = 12))+
  theme(axis.text.x = element_text(size = 10))

#Unecessary or inappropriate therapy
Unecessary_B = ggplot(gather(Inappropriate_Therapy_B, cols, value), aes(x = value)) + 
  geom_histogram() +xlim(0,1) + ylim(0,1000) + facet_grid(.~cols)+                                                                
  theme(strip.text.x = element_text(size = 12))+
  theme(axis.text.x = element_text(size = 10))

####Algorithm C (delayed TPT)#### 

outcome_cohort_C = vector("list", length = nsims)
for(j in 1:nsims){
  #select cohort for simulation 
  cohort = cohort_sims[[j]] 
  
  #initialize outcomes of cohort matrix
  outcome_cohort = matrix(nrow = ncohort, ncol = 10)
  
  #outcomes to track 
  colnames(outcome_cohort) = c("Sx_Screen","CRP_Screen","Xpert_Screen","TPT","Day_TPT","Treatment","Day_Treatment","ART","Day_ART","t")

  
  
  outcome_cohort[,"t"] = t_value[[j]]
  
  #All recieve Symptom Screen, Sx- get ART at first visit, 50% Sx+ get ART 
  outcome_cohort[which(cohort$Sx == "Sx-"),"ART"] = TRUE
  outcome_cohort[which(cohort$Sx == "Sx-"),"Day_ART"] = 0


  
  for(i in 1:ncohort){
    
    #sample visit2
    t = outcome_cohort[i,"t"]
    
    #50% Sx+ get ART at first visit 
    if(cohort[i,"Sx"] == "Sx+"){
      
      if(runif(1)<.5){
        outcome_cohort[i,"ART"] = 
          outcome_cohort[i,"Day_ART"] = 0
      }
      
    }

    #Probability of universal xpert screen 
    
    if(runif(1)*100<((parameter_dist_1$Pss[j]/100)*(parameter_dist_1$RRxs[j]/100)*100)){
      
      outcome_cohort[i,"Xpert_Screen"] = TRUE
  
      
      #Probability of following Universal Algorithm 
      if(runif(1)*100<parameter_dist_1$Puc[j]){
        
        #All delay TPT 
        
        if(!is.na(t)){
          if(cohort[i,"Xpert"] == "Xpert+"){
            if(runif(1)*100<parameter_dist_1$Pt2p[j]){
              outcome_cohort[i,"Treatment"] = TRUE
              outcome_cohort[i,"Day_Treatment"] = t
            }
          }
          else if(cohort[i,"Xpert"] == "Xpert-"){
            
            if(runif(1)*100<((parameter_dist_1$Pp1n[j]/100)*(parameter_dist_1$RRp2n[j]/100)*100)){
              outcome_cohort[i,"TPT"] = TRUE
              outcome_cohort[i,"Day_TPT"] = t
             
            }
          }
        }
        
      }
      else{
        outcome_cohort[i,"Sx_Screen"] = TRUE
        
        if(cohort[i,"Sx"] == "Sx+"){
          
          if(!is.na(t)){
            if(cohort[i,"Xpert"] == "Xpert+"){
              if(runif(1)*100<parameter_dist_1$Pt2p[j]){
                outcome_cohort[i,"Treatment"] = TRUE
                outcome_cohort[i,"Day_Treatment"] = t
              }
            }
            else if(cohort[i,"Xpert"] == "Xpert-"){
              
              if(runif(1)*100<((parameter_dist_1$Pp1n[j]/100)*(parameter_dist_1$RRp2n[j]/100)*100)){
                outcome_cohort[i,"TPT"] = TRUE
                outcome_cohort[i,"Day_TPT"] = t
              }
            }
          }
          
        } else if(cohort[i,"Sx"] == "Sx-"){
          
          #Probability of starting TPT simultaneously
          if(runif(1)*100<parameter_dist_1$Pp1n[j]){
            outcome_cohort[i,"TPT"] = TRUE
            outcome_cohort[i,"Day_TPT"] = 0
            
            if(!is.na(t)){
              if(cohort[i,"Xpert"] == "Xpert+"){
                if(runif(1)*100<parameter_dist_1$Pt2p[j]){
                  outcome_cohort[i,"Treatment"] = TRUE
                  outcome_cohort[i,"Day_Treatment"] = t
                  outcome_cohort[i,"Day_TPT"] = "discontinued"
                }
              }
              else if(cohort[i,"Xpert"] == "Xpert-"){
                #Probability of continuing TPT?
                if(t>31){
                  if(runif(1)*100>parameter_dist_1$Ptcp[j]){
                    outcome_cohort[i,"Day_TPT"] = "discontinued" 
                  }
                  
                }
                else{
                  if(runif(1)*100<(parameter_dist_1$dt[j])){
                    outcome_cohort[i,"Day_TPT"] = "discontinued" 
                  }
                }
                
              }
            }
            
          }
          else{
            if(!is.na(t)){
              if(cohort[i,"Xpert"] == "Xpert+"){
                if(runif(1)*100<parameter_dist_1$Pt2p[j]){
                  outcome_cohort[i,"Treatment"] = TRUE
                  outcome_cohort[i,"Day_Treatment"] = t
                }
              }
              else if(cohort[i,"Xpert"] == "Xpert-"){
                
                if(runif(1)*100<((parameter_dist_1$Pp1n[j]/100)*(parameter_dist_1$RRp2n[j]/100)*100)){
                  outcome_cohort[i,"TPT"] = TRUE
                  outcome_cohort[i,"Day_TPT"] = t
                }
                
              }
            }
          }
          
        }
        
        
      }
      
    }
    
    else{
      
      
      if(!is.na(t)){
        
        #Probability of considering TPT at visit 2
        if(runif(1)*100<((parameter_dist_1$Pp1n[j]/100)*(parameter_dist_1$RRp2n[j]/100)*100)){
          outcome_cohort[i,"TPT"] = TRUE
          outcome_cohort[i,"Day_TPT"] = t
        }
      }
      
      
    }
    
    
    
    
  }
  
  outcome_cohort[is.na(outcome_cohort[,"TPT"]),"TPT"] <- 0
  outcome_cohort[is.na(outcome_cohort[,"Treatment"]),"Treatment"] <- 0
  outcome_cohort[is.na(outcome_cohort[,"ART"]),"ART"] <- 0
  outcome_cohort[is.na(outcome_cohort[,"Sx_Screen"]),"Sx_Screen"] <- 0
  outcome_cohort[is.na(outcome_cohort[,"Xpert_Screen"]),"Xpert_Screen"] <- 0
  outcome_cohort[is.na(outcome_cohort[,"CRP_Screen"]),"CRP_Screen"] <- 0
   
  outcome_cohort[which(outcome_cohort[,"TPT"] == TRUE),"TPT"] <- 1
  outcome_cohort[which(outcome_cohort[,"Treatment"] == TRUE),"Treatment"] <- 1
  outcome_cohort[which(outcome_cohort[,"ART"] == TRUE),"ART"] <- 1
  outcome_cohort[which(outcome_cohort[,"Sx_Screen"] == TRUE),"Sx_Screen"] <- 1
  outcome_cohort[which(outcome_cohort[,"Xpert_Screen"] == TRUE),"Xpert_Screen"] <- 1
  outcome_cohort[which(outcome_cohort[,"CRP_Screen"] == TRUE),"CRP_Screen"] <- 1
  
  
  #return outcomes for cohort j 
  outcome_cohort_C[[j]] = outcome_cohort
  outcome_cohort_C[[j]] = cbind(cohort_sims[[j]],outcome_cohort)
  outcome_cohort_C[[j]] = cbind(outcome_cohort_C[[j]],Stop_Days[[j]])

  
  
}

#General outcomes

outcome_discon_c = vector("list",length = nsims)
for(a in 1:nsims){
  x = outcome_cohort_C[[a]]
  x$Day_TPT[which(x$Day_TPT == "discontinued")] <- 0
  x$Day_TPT = as.numeric(as.character(x$Day_TPT))
  outcome_discon_c[[a]] = x
  
}

Treatment = do.call(rbind,lapply(outcome_cohort_C, function(x) sum(x[,"Treatment"]==1 & x[,"TB"]=="TB")/sum(x[,"TB"]=="TB")) )
Treatment_60 = do.call(rbind,lapply(outcome_cohort_C, function(x) sum(x[,"Treatment"]==1 & as.numeric(as.character(x[,"Day_Treatment"]))<=60 & x[,"TB"]=="TB")/sum(x[,"TB"]=="TB")))
ART = do.call(rbind,lapply(outcome_cohort_C, function(x) sum(x[,"ART"]==1 & x[,"TB"]=="No_TB")/sum(x[,"TB"]=="No_TB")))
ART_Visit1 = do.call(rbind,lapply(outcome_cohort_C, function(x) sum(x[,"Day_ART"]==0 & x[,"TB"]=="No_TB")/sum(x[,"TB"]=="No_TB")))
ART_60 = do.call(rbind,lapply(outcome_cohort_C, function(x) sum(x[,"ART"]==1 & as.numeric(as.character(x[,"Day_ART"]))<=60 & x[,"TB"]=="No_TB")/sum(x[,"TB"]=="No_TB")))
TPT = do.call(rbind,lapply(outcome_cohort_C, function(x) sum(x[,"TPT"]==1 & x[,"TB"]=="No_TB")/sum(x[,"TB"]=="No_TB")))
TPT_60 = do.call(rbind,lapply(outcome_discon_c, function(x) sum(x[,"TPT"]==1 & as.numeric(as.character(x[,"Day_TPT"]))<=60 & x[,"TB"]=="No_TB")/sum(x[,"TB"]=="No_TB")))
Diagnosed = do.call(rbind,lapply(outcome_cohort_C, function(x) length(which(x[,"Xpert_Screen"]==1 & x[,"TB"]=="TB" & x[,"Xpert"] == "Xpert+"))/sum(x[,"TB"]=="TB")))
TB_ART_30 =  0#do.call(rbind,lapply(outcome_cohort_C, function(x) sum(x[,"TPT"]==1 & x[,"ART"]==1 & abs(x[,"Day_TPT"]- x[,"Day_ART"])<=30 & x[,"TB"]=="No_TB")/sum(x[,"TB"]=="No_TB")))
Treatment_Total_60 = do.call(rbind,lapply(outcome_cohort_C, function(x) sum(x[,"Treatment"]==1 & as.numeric(as.character(x[,"Day_Treatment"]))<=60 & x[,"TB"]=="TB")))
TPT_Total_60 = do.call(rbind,lapply(outcome_discon_c, function(x) sum(x[,"TPT"]==1 & as.numeric(as.character(x[,"Day_TPT"]))<=60 & x[,"TB"]=="No_TB")))
Diagnosed_Total = do.call(rbind,lapply(outcome_cohort_C, function(x) length(which(x[,"Xpert_Screen"]==1 & x[,"TB"]=="TB" & x[,"Xpert"] == "Xpert+"))))
eval_outcomes_C = as.data.frame(cbind(Treatment,Treatment_60,ART,ART_60,TPT,TPT_60,Diagnosed,Treatment_Total_60 ,TPT_Total_60,Diagnosed_Total))
colnames(eval_outcomes_C) = c("Treatment for TB","Treatment <60d for TB","ART for No TB","ART <60d for No TB","TPT for No TB","TPT <60d for No TB","Diagnosed TB","Total Treatment","Total TPT","Total Diagnosed")


#Diagnostic evaluations 

symp = do.call(rbind, lapply(outcome_cohort_C, function(x) length(which(x[,"Sx_Screen"] == 1))))
CRP = do.call(rbind, lapply(outcome_cohort_C, function(x) length(which(x[,"CRP_Screen"] == 1))))
Xperts = do.call(rbind, lapply(outcome_cohort_C, function(x) length(which(x[,"Xpert_Screen"] == 1))))

Diagnostic_eval_C = as.data.frame(cbind(symp,CRP, Xperts))
colnames(Diagnostic_eval_C) = c("Symptom Screen","CRP Screen","Xpert Screen")


#Unecessary or inappropriate therapy
Wrong_Treatment = do.call(rbind,lapply(outcome_cohort_C, function(x) sum(x[,"Treatment"]==1 & x[,"TB"]=="No_TB")/sum(x[,"TB"]=="No_TB")) ) 
Wrong_TPT = do.call(rbind,lapply(outcome_cohort_C, function(x) sum(x[,"TPT"]==1 & x[,"TB"]=="TB")/sum(x[,"TB"]=="TB")) ) 
Wrong_Treatment_Total =do.call(rbind,lapply(outcome_cohort_C, function(x) sum(x[,"Treatment"]==1 & x[,"TB"]=="No_TB")) ) 
Wrong_TPT_Total = do.call(rbind,lapply(outcome_cohort_C, function(x) sum(x[,"TPT"]==1 & x[,"TB"]=="TB")) ) 
Inappropriate_Therapy_C = as.data.frame(cbind(Wrong_Treatment,Wrong_TPT,Wrong_Treatment_Total,Wrong_TPT_Total))
colnames(Inappropriate_Therapy_C) = c("Treatment Incorrectly Given Without TB","TPT Incorrectly Given With TB","Total Incorrect Treatment","Total Incorrect TPT")

TPT_30_C = do.call(rbind, lapply(outcome_cohort_C, function(x) sum(x[,"TPT"]==1 & x[,"TB"]=="TB" & x[,"Day_TPT"] != "discontinued" & as.numeric(as.character(x[,"TPT_Stop_Day"])) > 30)/sum(x[,"TB"]=="TB")))

#Plot values

#Evaluation outcomes
Eval_C = ggplot(gather(eval_outcomes_C, cols, value), aes(x = value)) + 
  geom_histogram() + xlim(0,1.1) + ylim(0,1000) + facet_grid(.~cols)+                                                                # Change font size
  theme(strip.text.x = element_text(size = 12))+
  theme(axis.text.x = element_text(size = 10))

#Diagnostic evaluations 
Diag_C = ggplot(gather(Diagnostic_eval_C, cols, value), aes(x = value)) + 
  geom_histogram() + ylim(0,1000) + facet_grid(.~cols)+                                                                # Change font size
  theme(strip.text.x = element_text(size = 12))+
  theme(axis.text.x = element_text(size = 10))

#Unecessary or inappropriate therapy
Unecessary_C = ggplot(gather(Inappropriate_Therapy_C, cols, value), aes(x = value)) + 
  geom_histogram() +xlim(0,1) + ylim(0,1000) + facet_grid(.~cols)+                                                                # Change font size
  theme(strip.text.x = element_text(size = 12))+
  theme(axis.text.x = element_text(size = 10))


####Algorithm D (Universal TPT)#####

outcome_cohort_D = vector("list", length = nsims)

for (j in 1:nsims){
  #select cohort for simulation 
  cohort = cohort_sims[[j]] 
  
  #initialize outcomes of cohort matrix
  outcome_cohort = matrix(nrow = ncohort, ncol = 10)
  
  #outcomes to track 
  colnames(outcome_cohort) = c("Sx_Screen","CRP_Screen","Xpert_Screen","TPT","Day_TPT","Treatment","Day_Treatment","ART","Day_ART","t")
  
  
  
  outcome_cohort[,"t"] = t_value[[j]]
  
  #All recieve Symptom Screen, Sx- get ART at first visit, 50% Sx+ get ART at first visit 
  outcome_cohort[which(cohort$Sx == "Sx-"),"ART"] = TRUE
  outcome_cohort[which(cohort$Sx == "Sx-"),"Day_ART"] = 0
  
  for(i in 1:ncohort){
  
    #sample visit2
    t = outcome_cohort[i,"t"]
    
    #50% Sx+ get ART at first visit 
    if(cohort[i,"Sx"] == "Sx+"){
      
      if(runif(1)<.5){
        outcome_cohort[i,"ART"] = 
          outcome_cohort[i,"Day_ART"] = 0
      }
      
    }
    
    #Probability of universal xpert screen 
    
    if(runif(1)*100<(parameter_dist_1$Pss[j]/100)*(parameter_dist_1$RRxs[j]/100)*100){
      outcome_cohort[i,"Xpert_Screen"] = TRUE
      
      #Probability of following universal algorithm 
      
      if(runif(1)*100<parameter_dist_1$Pud[j]){
        outcome_cohort[i,"TPT"] = TRUE
        outcome_cohort[i,"Day_TPT"] = 0
        
        if(!is.na(t)){
          if(cohort[i,"Xpert"] == "Xpert+"){
            if(runif(1)*100<parameter_dist_1$Pt2p[j]){
              outcome_cohort[i,"Treatment"] = TRUE
              outcome_cohort[i,"Day_Treatment"] = t
              outcome_cohort[i,"Day_TPT"] = "discontinued"
            }
          }
          else if(cohort[i,"Xpert"] == "Xpert-"){
            #Probability of continuing TPT? 
            
            if(t>31){
              if(runif(1)*100>parameter_dist_1$Ptcp[j]){
                outcome_cohort[i,"Day_TPT"] = "discontinued" 
              }
              
            }
            else{
              if(runif(1)*100<(parameter_dist_1$dt[j])){
                outcome_cohort[i,"Day_TPT"] = "discontinued" 
              }
            }
          }
        }
        
        
      }
      else{
        outcome_cohort[i,"Sx_Screen"] = TRUE
        if(cohort[i,"Sx"] == "Sx+"){
          
          if(!is.na(t)){
            if(cohort[i,"Xpert"] == "Xpert+"){
              if(runif(1)*100<parameter_dist_1$Pt2p[j]){
                outcome_cohort[i,"Treatment"] = TRUE
                outcome_cohort[i,"Day_Treatment"] = t
              }
            }
            else if(cohort[i,"Xpert"] == "Xpert-"){
              
              if(runif(1)*100<((parameter_dist_1$Pp1n[j]/100)*(parameter_dist_1$RRp2n[j]/100)*100)){
                outcome_cohort[i,"TPT"] = TRUE
                outcome_cohort[i,"Day_TPT"] = t
            
              }
            }
          }
          
        }
        else if(cohort[i,"Sx"] == "Sx-"){
          
          #Probability of starting TPT simultaneously
          if(runif(1)*100<parameter_dist_1$Pp1n[j]){
            outcome_cohort[i,"TPT"] = TRUE
            outcome_cohort[i,"Day_TPT"] = 0
            
            if(!is.na(t)){
              if(cohort[i,"Xpert"] == "Xpert+"){
                if(runif(1)*100<parameter_dist_1$Pt2p[j]){
                  outcome_cohort[i,"Treatment"] = TRUE
                  outcome_cohort[i,"Day_Treatment"] = t
                  outcome_cohort[i,"Day_TPT"] = "discontinued"
                }
              }
              else if(cohort[i,"Xpert"] == "Xpert-"){
                #Probability of continuing TPT? remember the <31 >31 issue with PTCP and 1-dt? 
                
               if(t>31){
                  if(runif(1)*100>parameter_dist_1$Ptcp[j]){
                    outcome_cohort[i,"Day_TPT"] = "discontinued" 
              }
                  
               }
               else{
                if(runif(1)*100<(parameter_dist_1$dt[j])){
                   outcome_cohort[i,"Day_TPT"] = "discontinued" 
                 }
               }
                
              }
            }
            
          }
          else{
            if(!is.na(t)){
              if(cohort[i,"Xpert"] == "Xpert+"){
                if(runif(1)*100<parameter_dist_1$Pt2p[j]){
                  outcome_cohort[i,"Treatment"] = TRUE
                  outcome_cohort[i,"Day_Treatment"] = t
                }
              }
              else if(cohort[i,"Xpert"] == "Xpert-"){
                
                if(runif(1)*100<((parameter_dist_1$Pp1n[j]/100)*(parameter_dist_1$RRp2n[j]/100)*100)){
                  outcome_cohort[i,"TPT"] = TRUE
                  outcome_cohort[i,"Day_TPT"] = t
                }
                
              }
            }
          }
          
        }
        
      }
      
    }
    else{
      
      if(!is.na(t)){
        
        #Probability of considering TPT at visit 2
        if(runif(1)*100<((parameter_dist_1$Pp1n[j]/100)*(parameter_dist_1$RRp2n[j]/100)*100)){
          outcome_cohort[i,"TPT"] = TRUE
          outcome_cohort[i,"Day_TPT"] = t
        }
      }
      
    }
    
    
  }
  outcome_cohort[is.na(outcome_cohort[,"TPT"]),"TPT"] <- 0
  outcome_cohort[is.na(outcome_cohort[,"Treatment"]),"Treatment"] <- 0
  outcome_cohort[is.na(outcome_cohort[,"ART"]),"ART"] <- 0
  outcome_cohort[is.na(outcome_cohort[,"Sx_Screen"]),"Sx_Screen"] <- 0
  outcome_cohort[is.na(outcome_cohort[,"Xpert_Screen"]),"Xpert_Screen"] <- 0
  outcome_cohort[is.na(outcome_cohort[,"CRP_Screen"]),"CRP_Screen"] <- 0
  
  outcome_cohort[which(outcome_cohort[,"TPT"] == TRUE),"TPT"] <- 1
  outcome_cohort[which(outcome_cohort[,"Treatment"] == TRUE),"Treatment"] <- 1
  outcome_cohort[which(outcome_cohort[,"ART"] == TRUE),"ART"] <- 1
  outcome_cohort[which(outcome_cohort[,"Sx_Screen"] == TRUE),"Sx_Screen"] <- 1
  outcome_cohort[which(outcome_cohort[,"Xpert_Screen"] == TRUE),"Xpert_Screen"] <- 1
  outcome_cohort[which(outcome_cohort[,"CRP_Screen"] == TRUE),"CRP_Screen"] <- 1
  
  #return outcomes for cohort j 
  outcome_cohort_D[[j]] = outcome_cohort
  outcome_cohort_D[[j]] = cbind(cohort_sims[[j]],outcome_cohort)
  outcome_cohort_D[[j]] = cbind(outcome_cohort_D[[j]],Stop_Days[[j]])
}

#General outcomes
outcome_discon_d = vector("list",length = nsims)
for(a in 1:nsims){
  x = outcome_cohort_D[[a]]
  x$Day_TPT[which(x$Day_TPT == "discontinued")] <- 0
  x$Day_TPT = as.numeric(as.character(x$Day_TPT))
  outcome_discon_d[[a]] = x
  
}

Treatment = do.call(rbind,lapply(outcome_cohort_D, function(x) sum(x[,"Treatment"]==1 & x[,"TB"]=="TB")/sum(x[,"TB"]=="TB")) ) 
Treatment_60 = do.call(rbind,lapply(outcome_cohort_D, function(x) sum(x[,"Treatment"]==1 & as.numeric(as.character(x[,"Day_Treatment"]))<=60 & x[,"TB"]=="TB")/sum(x[,"TB"]=="TB")))
ART = do.call(rbind,lapply(outcome_cohort_D, function(x) sum(x[,"ART"]==1 & x[,"TB"]=="No_TB")/sum(x[,"TB"]=="No_TB")))
ART_Visit1 = do.call(rbind,lapply(outcome_cohort_D, function(x) sum(x[,"Day_ART"]==0 & x[,"TB"]=="No_TB")/sum(x[,"TB"]=="No_TB")))
ART_60 = do.call(rbind,lapply(outcome_cohort_D, function(x) sum(x[,"ART"]==1 & as.numeric(as.character(x[,"Day_ART"]))<=60 & x[,"TB"]=="No_TB")/sum(x[,"TB"]=="No_TB")))
TPT = do.call(rbind,lapply(outcome_cohort_D, function(x) sum(x[,"TPT"]==1 & x[,"TB"]=="No_TB")/sum(x[,"TB"]=="No_TB")))
TPT_60 = do.call(rbind,lapply(outcome_discon_d, function(x) sum(x[,"TPT"]==1 & as.numeric(as.character(x[,"Day_TPT"]))<=60 & x[,"TB"]=="No_TB")/sum(x[,"TB"]=="No_TB")))
Diagnosed = do.call(rbind,lapply(outcome_cohort_D, function(x) length(which(x[,"Xpert_Screen"]==1 & x[,"TB"]=="TB" & x[,"Xpert"] == "Xpert+"))/sum(x[,"TB"]=="TB")))
TB_ART_30 =  0#do.call(rbind,lapply(outcome_cohort_D, function(x) sum(x[,"TPT"]==1 & x[,"ART"]==1 & abs(x[,"Day_TPT"]- x[,"Day_ART"])<=30 & x[,"TB"]=="No_TB")/sum(x[,"TB"]=="No_TB")))
Treatment_Total_60 = do.call(rbind,lapply(outcome_cohort_D, function(x) sum(x[,"Treatment"]==1 & as.numeric(as.character(x[,"Day_Treatment"]))<=60 & x[,"TB"]=="TB")))
TPT_Total_60 = do.call(rbind,lapply(outcome_discon_d, function(x) sum(x[,"TPT"]==1 & as.numeric(as.character(x[,"Day_TPT"]))<=60 & x[,"TB"]=="No_TB")))
Diagnosed_Total = do.call(rbind,lapply(outcome_cohort_D, function(x) length(which(x[,"Xpert_Screen"]==1 & x[,"TB"]=="TB" & x[,"Xpert"] == "Xpert+"))))
eval_outcomes_D = as.data.frame(cbind(Treatment,Treatment_60,ART,ART_60,TPT,TPT_60,Diagnosed,Treatment_Total_60 ,TPT_Total_60,Diagnosed_Total))
colnames(eval_outcomes_D) = c("Treatment for TB","Treatment <60d for TB","ART for No TB","ART <60d for No TB","TPT for No TB","TPT <60d for No TB","Diagnosed TB","Total Treatment","Total TPT","Total Diagnosed")


#Diagnostic evaluations 

symp = do.call(rbind, lapply(outcome_cohort_D, function(x) length(which(x[,"Sx_Screen"] == 1))))
CRP = do.call(rbind, lapply(outcome_cohort_D, function(x) length(which(x[,"CRP_Screen"] == 1))))
Xperts = do.call(rbind, lapply(outcome_cohort_D, function(x) length(which(x[,"Xpert_Screen"] == 1))))

Diagnostic_eval_D = as.data.frame(cbind(symp,CRP, Xperts))
colnames(Diagnostic_eval_D) =  c("Symptom Screen","CRP Screen","Xpert Screen")


#Unecessary or inappropriate therapy
Wrong_Treatment = do.call(rbind,lapply(outcome_cohort_D, function(x) sum(x[,"Treatment"]==1 & x[,"TB"]=="No_TB")/sum(x[,"TB"]=="No_TB")) ) 
Wrong_TPT = do.call(rbind,lapply(outcome_cohort_D, function(x) sum(x[,"TPT"]==1 & x[,"TB"]=="TB")/sum(x[,"TB"]=="TB")) ) 
Wrong_Treatment_Total =do.call(rbind,lapply(outcome_cohort_D, function(x) sum(x[,"Treatment"]==1 & x[,"TB"]=="No_TB")) ) 
Wrong_TPT_Total = do.call(rbind,lapply(outcome_cohort_D, function(x) sum(x[,"TPT"]==1 & x[,"TB"]=="TB")) ) 
Inappropriate_Therapy_D = as.data.frame(cbind(Wrong_Treatment,Wrong_TPT,Wrong_Treatment_Total,Wrong_TPT_Total))
colnames(Inappropriate_Therapy_D) = c("Treatment Incorrectly Given Without TB","TPT Incorrectly Given With TB","Total Incorrect Treatment","Total Incorrect TPT")

TPT_30_D = do.call(rbind, lapply(outcome_cohort_D, function(x) sum(x[,"TPT"]==1 & x[,"TB"]=="TB" & x[,"Day_TPT"] != "discontinued" & as.numeric(as.character(x[,"TPT_Stop_Day"])) > 30)/sum(x[,"TB"]=="TB")))

#Plot values

#Evaluation outcomes
Eval_D = ggplot(gather(eval_outcomes_D, cols, value), aes(x = value)) + 
  geom_histogram() + xlim(0,1.1) + ylim(0,1000) + facet_grid(.~cols)+                                                               
  theme(strip.text.x = element_text(size = 12))+
  theme(axis.text.x = element_text(size = 10))

#Diagnostic evaluations 
Diag_D = ggplot(gather(Diagnostic_eval_D, cols, value), aes(x = value)) + 
  geom_histogram() + ylim(0,1000) + facet_grid(.~cols)+                                                                
  theme(strip.text.x = element_text(size = 12))+
  theme(axis.text.x = element_text(size = 10))

#Unecessary or inappropriate therapy
Unecessary_D = ggplot(gather(Inappropriate_Therapy_D, cols, value), aes(x = value)) + 
  geom_histogram() +xlim(0,1) + ylim(0,1000) + facet_grid(.~cols)+                                                                
  theme(strip.text.x = element_text(size = 12))+
  theme(axis.text.x = element_text(size = 10))



#Combine Plots


figure_eval <- ggarrange(Eval_A, Eval_B, Eval_C,Eval_D,
                    ncol = 1, nrow = 4)

figure_diag <- ggarrange(Diag_A, Diag_B, Diag_C,Diag_D,
                         ncol = 1, nrow = 4)

figure_unecessary <-ggarrange(Unecessary_A, Unecessary_B, Unecessary_C,Unecessary_D,
                              ncol = 1, nrow = 4)

#Number treatment treated/TPT within 60 days
summary(do.call(rbind,lapply(outcome_cohort_D, function(x) sum(x[,"Treatment"]==1 & as.numeric(as.character(x[,"Day_Treatment"]))<=60 & x[,"TB"]=="TB")))-do.call(rbind,lapply(outcome_cohort_A, function(x) sum(x[,"Treatment"]==1 & as.numeric(as.character(x[,"Day_Treatment"]))<=60 & x[,"TB"]=="TB"))))
summary(do.call(rbind,lapply(outcome_cohort_D, function(x) sum(x[,"Treatment"]==1 & as.numeric(as.character(x[,"Day_Treatment"]))<=60 & x[,"TB"]=="TB")))-do.call(rbind,lapply(outcome_cohort_B, function(x) sum(x[,"Treatment"]==1 & as.numeric(as.character(x[,"Day_Treatment"]))<=60 & x[,"TB"]=="TB"))))
summary(do.call(rbind,lapply(outcome_cohort_D, function(x) sum(x[,"Treatment"]==1 & as.numeric(as.character(x[,"Day_Treatment"]))<=60 & x[,"TB"]=="TB")))-do.call(rbind,lapply(outcome_cohort_C, function(x) sum(x[,"Treatment"]==1 & as.numeric(as.character(x[,"Day_Treatment"]))<=60 & x[,"TB"]=="TB"))))
summary(do.call(rbind,lapply(outcome_cohort_C, function(x) sum(x[,"Treatment"]==1 & as.numeric(as.character(x[,"Day_Treatment"]))<=60 & x[,"TB"]=="TB")))-do.call(rbind,lapply(outcome_cohort_A, function(x) sum(x[,"Treatment"]==1 & as.numeric(as.character(x[,"Day_Treatment"]))<=60 & x[,"TB"]=="TB"))))


summary(do.call(rbind,lapply(outcome_discon_d, function(x) sum(x[,"TPT"]==1 & as.numeric(as.character(x[,"Day_TPT"]))<=60 & x[,"TB"]=="No_TB")))-do.call(rbind,lapply(outcome_discon_a, function(x) sum(x[,"TPT"]==1 & as.numeric(as.character(x[,"Day_TPT"]))<=60 & x[,"TB"]=="No_TB"))))
summary(do.call(rbind,lapply(outcome_discon_d, function(x) sum(x[,"TPT"]==1 & as.numeric(as.character(x[,"Day_TPT"]))<=60 & x[,"TB"]=="No_TB")))-do.call(rbind,lapply(outcome_discon_b, function(x) sum(x[,"TPT"]==1 & as.numeric(as.character(x[,"Day_TPT"]))<=60 & x[,"TB"]=="No_TB"))))
summary(do.call(rbind,lapply(outcome_discon_d, function(x) sum(x[,"TPT"]==1 & as.numeric(as.character(x[,"Day_TPT"]))<=60 & x[,"TB"]=="No_TB")))-do.call(rbind,lapply(outcome_discon_c, function(x) sum(x[,"TPT"]==1 & as.numeric(as.character(x[,"Day_TPT"]))<=60 & x[,"TB"]=="No_TB"))))
summary(do.call(rbind,lapply(outcome_discon_c, function(x) sum(x[,"TPT"]==1 & as.numeric(as.character(x[,"Day_TPT"]))<=60 & x[,"TB"]=="No_TB")))-do.call(rbind,lapply(outcome_discon_a, function(x) sum(x[,"TPT"]==1 & as.numeric(as.character(x[,"Day_TPT"]))<=60 & x[,"TB"]=="No_TB"))))

summary(do.call(rbind, lapply(outcome_cohort_A, function(x) sum(x[,"TPT"]==1 & x[,"TB"]=="No_TB" & x[,"Day_TPT"] != "discontinued" & as.numeric(as.character(x[,"TPT_Stop_Day"])) > 30)/sum(x[,"TB"]=="No_TB"))))
summary(do.call(rbind, lapply(outcome_cohort_B, function(x) sum(x[,"TPT"]==1 & x[,"TB"]=="No_TB" & x[,"Day_TPT"] != "discontinued" & as.numeric(as.character(x[,"TPT_Stop_Day"])) > 30)/sum(x[,"TB"]=="No_TB"))))
summary(do.call(rbind, lapply(outcome_cohort_C, function(x) sum(x[,"TPT"]==1 & x[,"TB"]=="No_TB" & x[,"Day_TPT"] != "discontinued" & as.numeric(as.character(x[,"TPT_Stop_Day"])) > 30)/sum(x[,"TB"]=="No_TB"))))
summary(do.call(rbind, lapply(outcome_cohort_D, function(x) sum(x[,"TPT"]==1 & x[,"TB"]=="No_TB" & x[,"Day_TPT"] != "discontinued" & as.numeric(as.character(x[,"TPT_Stop_Day"])) > 30)/sum(x[,"TB"]=="No_TB"))))

        