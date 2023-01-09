
sourcelocation <- "~/Documents/JHU/IPT Model/"
library(reshape2)
#Costing outline sketch 
source(paste0(sourcelocation,"Model Part 3.R"))

#CNR = .61 (.43-.95) 
#CNR = rtriangle(nsims,a = .61,b =.43,c = .95)
#Proportion that gets treatment 
CNR = .61
#TB susceptible treatment cost per month
TB_S_cost = 130/12 #Double check
#Xpert test cost
Xpert_cost = 33.5
#TB MDR treatment cost total
TB_MDR_cost = 9700.57 #Still need to update
#CRP cost per person
CRP_cost = 2.5
#Monthly cost of TPT
TPT_Cost = 34.11/12
#Monthly cost of ART
ART_Cost = 249.15
#Proprotion MDR 
MDR_p = .021

#Cost of incident TB 


Cost_Incident_TB <-function(long_outcome){
  
  incident_tb = sum(long_outcome$expected_incidence_total)
  treated = CNR*incident_tb #deterministic fraction that gets TB 
  total_cost = treated*(Xpert_cost+6*TB_S_cost) #Assume that all incident TB get full 6 month course

  #return total cost
  return(total_cost)

 }

#Cost of baseline TB 

Cost_Baseline_TB <-function(outcome){
  
  #randomly assign some baseline prevalent TB as MDR 
  outcome$MDR = NA
  outcome_TB = outcome[which(outcome$TB == "TB"),]
  outcome_No_TB = outcome[which(outcome$TB == "No_TB"),]
  MDR_index= sample(1:length(which(outcome$TB == "TB")),round(MDR_p*1000,0))
  outcome_No_TB$MDR = NA
  outcome_TB$MDR[MDR_index] <- 1
  outcome = rbind(outcome_TB,outcome_No_TB)
  
  cost = matrix(nrow = nsims, ncol = 1)
  
  for(a in 1:nsims){
    if(outcome$Treatment[a] == 1){
      if(!is.na(outcome$MDR[a])){
        cost[a,] = TB_MDR_cost = 9700.57
      } else{
        cost[a,] = TB_S_cost*6
      }
    } else{
      cost[a,] =  0
    }
  
  }

  total_cost = sum(cost)
  
  #return total cost
  return(total_cost)
  
}

#Cost Basline

Cost_Baseline_A = unlist(lapply(outcome_cohort_A,Cost_Baseline_TB))
Cost_Baseline_B = unlist(lapply(outcome_cohort_B,Cost_Baseline_TB))
Cost_Baseline_C = unlist(lapply(outcome_cohort_C,Cost_Baseline_TB))
Cost_Baseline_D = unlist(lapply(outcome_cohort_D,Cost_Baseline_TB))

Baseline_Cost = cbind(Cost_Baseline_A,Cost_Baseline_B,Cost_Baseline_C,Cost_Baseline_D)
colnames(Baseline_Cost) = c("Baseline Cost A","Baseline Cost B","Baseline Cost C","Baseline Cost D")
Baseline_Cost = melt(Baseline_Cost)
Baseline_Cost = Baseline_Cost[,-1]
colnames(Baseline_Cost) = c("Algorithm","Cost")

#Cost Incident

Cost_Incident_A = unlist(lapply(long_outcome_A_TB,Cost_Incident_TB))
Cost_Incident_B = unlist(lapply(long_outcome_B_TB,Cost_Incident_TB))
Cost_Incident_C = unlist(lapply(long_outcome_C_TB,Cost_Incident_TB))
Cost_Incident_D = unlist(lapply(long_outcome_D_TB,Cost_Incident_TB))

Incident_Cost = cbind(Cost_Incident_A,Cost_Incident_B,Cost_Incident_C,Cost_Incident_D)
colnames(Incident_Cost) = c("Incident Cost A","Incident Cost B","Incident Cost C","Incident Cost D")
Incident_Cost = melt(Incident_Cost)
Incident_Cost = Incident_Cost[,-1]
colnames(Incident_Cost) = c("Algorithm","Cost")

#TPT 
#Cannot make assumptions about whether TPT stops after incident TB because we don't know the day it occurs

Cost_TPT <-function(outcome,long_outcome){
  
  outcome_TB = outcome[which(outcome$TB == "TB"),]
 
  #Inappropriate therapy, up to 6 months if continued
  TB_Cost = matrix(0,nrow = dim(outcome_TB)[1],ncol = 1)
  for(a in 1:dim(outcome_TB)[1]){
    if(outcome_TB$TPT[a] == 1){
      if(outcome_TB$Day_TPT[a] == "discontinued"){
        TB_Cost[a,1] = TPT_Cost
      } else {
        if(outcome_TB$TPT_Stop_Day[a]>360){
          TB_Cost[a,1] = TPT_Cost*12
        } else{
          TB_Cost[a,1] = TPT_Cost*(outcome_TB$TPT_Stop_Day[a])/30
        }
      }
    } else{
      TB_Cost[a,1] = 0
    }
  }
  
  TB_TPT_Cost = sum(TB_Cost)
  
  No_TB_TPT_Cost = matrix(0,nrow = dim(long_outcome)[1],ncol = 1)
  for(a in 1:dim(long_outcome)[1]){
  
    if(long_outcome$TPT_Stop_Day[a]>360){
      No_TB_TPT_Cost[a,1]  = TPT_Cost*12
    } else{
      No_TB_TPT_Cost[a,1]  = TPT_Cost*(long_outcome$TPT_Stop_Day[a])/30
    }
    
  }
  
  total_cost = sum(No_TB_TPT_Cost) + TB_TPT_Cost
  
  return(total_cost)
  
  
}

TPT_Cost_A = matrix(nrow = 1000,ncol= 1)
TPT_Cost_B = matrix(nrow = 1000,ncol= 1)
TPT_Cost_C = matrix(nrow = 1000,ncol= 1)
TPT_Cost_D = matrix(nrow = 1000,ncol= 1)
for(a in 1:1000){
  TPT_Cost_A[a,] = Cost_TPT(outcome_cohort_A[[a]],long_outcome_A_TB[[a]])
  TPT_Cost_B[a,] = Cost_TPT(outcome_cohort_B[[a]],long_outcome_B_TB[[a]])
  TPT_Cost_C[a,] = Cost_TPT(outcome_cohort_C[[a]],long_outcome_C_TB[[a]])
  TPT_Cost_D[a,] = Cost_TPT(outcome_cohort_D[[a]],long_outcome_D_TB[[a]])
}

TPT_cost = cbind(TPT_Cost_A,TPT_Cost_B,TPT_Cost_C,TPT_Cost_D)
colnames(TPT_cost) = c("TPT Cost A","TPT Cost B","TPT Cost C","TPT Cost D")
TPT_cost = melt(TPT_cost)
TPT_cost = TPT_cost[,-1]
colnames(TPT_cost) = c("Algorithm","Cost")

#Xpert 

Xpert_Cost_A = matrix(nrow = 1000,ncol= 1)
Xpert_Cost_B = matrix(nrow = 1000,ncol= 1)
Xpert_Cost_C = matrix(nrow = 1000,ncol= 1)
Xpert_Cost_D = matrix(nrow = 1000,ncol= 1)

Xpert_Cost_A[,1] = Diagnostic_eval_A$`Xpert Screen`*Xpert_cost
Xpert_Cost_B[,1] = Diagnostic_eval_B$`Xpert Screen`*Xpert_cost
Xpert_Cost_C[,1] = Diagnostic_eval_C$`Xpert Screen`*Xpert_cost
Xpert_Cost_D[,1] = Diagnostic_eval_D$`Xpert Screen`*Xpert_cost

Xpert_Cost = cbind(Xpert_Cost_A,Xpert_Cost_B, Xpert_Cost_C, Xpert_Cost_D)
colnames(Xpert_Cost) = c("Xpert Cost A","Xpert Cost B","Xpert Cost C","Xpert Cost D")
Xpert_Cost = melt(Xpert_Cost)
Xpert_Cost = Xpert_Cost[,-1]
colnames(Xpert_Cost) = c("Algorithm","Cost")

#CRP 

Cost_CRP <-function(outcome_cohort){
  
  CRP_cost = sum(as.numeric(as.character(outcome_cohort$CRP_Screen)))*CRP_cost
  return(CRP_cost)
}

CRP_Cost_A =do.call("rbind",lapply(outcome_cohort_A,Cost_CRP))
CRP_Cost_B =do.call("rbind",lapply(outcome_cohort_B,Cost_CRP))
CRP_Cost_C =do.call("rbind",lapply(outcome_cohort_C,Cost_CRP))
CRP_Cost_D =do.call("rbind",lapply(outcome_cohort_D,Cost_CRP))



#Visualize outputs

#Total cost of Xpert across algorithms 
ggplot(Xpert_Cost, aes(x = Cost)) + 
  geom_histogram() + facet_wrap(~Algorithm)                                                       
  theme(strip.text.x = element_text(size = 12))+
  theme(axis.text.x = element_text(size = 10))


#Total cost of TPT across algorithms 
ggplot(TPT_cost, aes(x = Cost)) + 
    geom_histogram() + facet_wrap(~Algorithm)                                                       
  theme(strip.text.x = element_text(size = 12))+
    theme(axis.text.x = element_text(size = 10))
  

#Total cost of Treatment Baseline TB
  ggplot(Baseline_Cost, aes(x = Cost)) + 
    geom_histogram() + facet_wrap(~Algorithm)                                                       
  theme(strip.text.x = element_text(size = 12))+
    theme(axis.text.x = element_text(size = 10))
  

#Total cost of Treatment Incident TB

  ggplot(Incident_Cost, aes(x = Cost)) + 
    geom_histogram() + facet_wrap(~Algorithm)                                                       
  theme(strip.text.x = element_text(size = 12))+
    theme(axis.text.x = element_text(size = 10))
  
  
#Total cost - excluding cost of Baseline TB that was diagnosed and treated post visit 2 
  
  Total_Cost_A = cbind(TPT_Cost_A,Xpert_Cost_A,Cost_Incident_A,Cost_Baseline_A,CRP_Cost_A)
  Total_Cost_B = cbind(TPT_Cost_B,Xpert_Cost_B,Cost_Incident_B,Cost_Baseline_B,CRP_Cost_B)
  Total_Cost_C = cbind(TPT_Cost_C,Xpert_Cost_C,Cost_Incident_C,Cost_Baseline_C,CRP_Cost_C)
  Total_Cost_D = cbind(TPT_Cost_D,Xpert_Cost_D,Cost_Incident_D,Cost_Baseline_D,CRP_Cost_D)
  
  t_cost_A = rowSums(Total_Cost_A)
  t_cost_B = rowSums(Total_Cost_B)
  t_cost_C = rowSums(Total_Cost_C)
  t_cost_D = rowSums(Total_Cost_D)

##### All Outcomes #####
  
### Costing outcomes ###
  
#Outcome 1 
#proportion recieveing treatment less than 60 days/total cost
treated_per_cost_A = t_cost_A/eval_outcomes_A$`Total Treatment`
treated_per_cost_B = t_cost_B/eval_outcomes_B$`Total Treatment`
treated_per_cost_C = t_cost_C/eval_outcomes_C$`Total Treatment`
treated_per_cost_D = t_cost_D/eval_outcomes_D$`Total Treatment`
  
par(mfrow = c(1,4))
hist(treated_per_cost_A, xlim = c(0,0.00005))
hist(treated_per_cost_B, xlim = c(0,0.00005))
hist(treated_per_cost_C, xlim = c(0,0.00005))
hist(treated_per_cost_D, xlim = c(0,0.00005))

#Outcome 2
#number of incident TB cases/total cost 
incident_TB_A = do.call(rbind,lapply(long_outcome_A_TB, function(x) length(which(x$incident_tb_day<=730))))
incident_TB_B = do.call(rbind,lapply(long_outcome_B_TB, function(x) length(which(x$incident_tb_day<=730))))
incident_TB_C = do.call(rbind,lapply(long_outcome_C_TB, function(x) length(which(x$incident_tb_day<=730))))
incident_TB_D = do.call(rbind,lapply(long_outcome_D_TB, function(x) length(which(x$incident_tb_day<=730))))

incident_per_cost_A = t_cost_A/incident_TB_A
incident_per_cost_B = t_cost_B/incident_TB_B
incident_per_cost_C = t_cost_C/incident_TB_C
incident_per_cost_D = t_cost_D/incident_TB_D


par(mfrow = c(1,4))
hist(incident_per_cost_A,xlim = c(0,0.0015))
hist(incident_per_cost_B,xlim = c(0,0.0015))
hist(incident_per_cost_C,xlim = c(0,0.0015))
hist(incident_per_cost_D,xlim = c(0,0.0015))

#Incident TB cases 
par(mfrow = c(1,4))
hist(incident_TB_A)
hist(incident_TB_B)
hist(incident_TB_C)
hist(incident_TB_D)

#Outcome 3
#proportion diagnosed TB/total cost
diagnosed_per_cost_A = t_cost_A/eval_outcomes_A$`Total Diagnosed`
diagnosed_per_cost_B = t_cost_B/eval_outcomes_B$`Total Diagnosed`
diagnosed_per_cost_C = t_cost_C/eval_outcomes_C$`Total Diagnosed`
diagnosed_per_cost_D = t_cost_D/eval_outcomes_D$`Total Diagnosed`

par(mfrow = c(1,4))
hist(diagnosed_per_cost_A,xlim = c(0,0.00005))
hist(diagnosed_per_cost_B,xlim = c(0,0.00005))
hist(diagnosed_per_cost_C,xlim = c(0,0.00005))
hist(diagnosed_per_cost_D,xlim = c(0,0.00005))

par(mfrow = c(1,4))
hist(eval_outcomes_A$`Total Diagnosed`)
hist(eval_outcomes_B$`Total Diagnosed`)
hist(eval_outcomes_C$`Total Diagnosed`)
hist(eval_outcomes_D$`Total Diagnosed`)


#Outcome 4
#proportion TPT <60 days/total cost
TPT_per_cost_A = t_cost_A/eval_outcomes_A$`Total TPT`
TPT_per_cost_B = t_cost_B/eval_outcomes_B$`Total TPT`
TPT_per_cost_C = t_cost_C/ eval_outcomes_C$`Total TPT`
TPT_per_cost_D = t_cost_D/eval_outcomes_D$`Total TPT`


par(mfrow = c(1,4))
hist(eval_outcomes_A$`Total TPT`)
hist(eval_outcomes_B$`Total TPT`)
hist(eval_outcomes_C$`Total TPT`)
hist(eval_outcomes_D$`Total TPT`)

par(mfrow = c(1,4))
hist(TPT_per_cost_A,xlim = c(0,0.00005))
hist(TPT_per_cost_B,xlim = c(0,0.00005))
hist(TPT_per_cost_C,xlim = c(0,0.00005))
hist(TPT_per_cost_D,xlim = c(0,0.00005))

#Outcome comparisons cost
#Alg A-D
par(mfrow = c(1,4))
hist((eval_outcomes_D$`Total Treatment`-eval_outcomes_A$`Total Treatment`)/(t_cost_D-t_cost_A))
hist((incident_TB_D-incident_TB_A)/(t_cost_D-t_cost_A))
hist((eval_outcomes_D$`Total TPT`-eval_outcomes_A$`Total TPT`)/(t_cost_D-t_cost_A))
hist((eval_outcomes_D$`Total Diagnosed`-eval_outcomes_A$`Total Diagnosed`)/(t_cost_D-t_cost_A))

#Alg B-D
hist((eval_outcomes_D$`Total Treatment`-eval_outcomes_B$`Total Treatment`)/(t_cost_D-t_cost_B))
hist((incident_TB_D-incident_TB_B)/(t_cost_D-t_cost_B))
hist((eval_outcomes_D$`Total TPT`-eval_outcomes_B$`Total TPT`)/(t_cost_D-t_cost_B))
hist((eval_outcomes_D$`Total Diagnosed`-eval_outcomes_B$`Total Diagnosed`)/(t_cost_D-t_cost_B))

#Alg C-D
par(mfrow = c(1,4))
hist((eval_outcomes_D$`Total Treatment`-eval_outcomes_C$`Total Treatment`)/(t_cost_D-t_cost_C))
hist((incident_TB_D-incident_TB_C)/(t_cost_D-t_cost_C))
hist((eval_outcomes_D$`Total TPT`-eval_outcomes_C$`Total TPT`)/(t_cost_D-t_cost_C))
hist((eval_outcomes_D$`Total Diagnosed`-eval_outcomes_C$`Total Diagnosed`)/(t_cost_D-t_cost_C))

#Alg A-C 
par(mfrow = c(1,4))
hist((eval_outcomes_C$`Total Treatment`-eval_outcomes_A$`Total Treatment`)/(t_cost_C-t_cost_A))
hist((incident_TB_C-incident_TB_A)/(t_cost_C-t_cost_A))
hist((eval_outcomes_C$`Total TPT`-eval_outcomes_A$`Total TPT`)/(t_cost_C-t_cost_A))
hist((eval_outcomes_C$`Total Diagnosed`-eval_outcomes_A$`Total Diagnosed`)/(t_cost_C-t_cost_A))

### Clinical outcomes ###

#Outcome 1 - Treatment compared
par(mfrow = c(1,4))
hist(eval_outcomes_D$`Total Treatment`-eval_outcomes_A$`Total Treatment`,main = "Algorithm D-A")
hist(eval_outcomes_D$`Total Treatment`-eval_outcomes_B$`Total Treatment`,main = "Algorithm D-B")
hist(eval_outcomes_D$`Total Treatment`-eval_outcomes_C$`Total Treatment`,main = "Algorithm D-C")
hist(eval_outcomes_C$`Total Treatment`-eval_outcomes_A$`Total Treatment`,main = "Algorithm C-A")

summary(eval_outcomes_D$`Total Treatment`-eval_outcomes_A$`Total Treatment`)
summary(eval_outcomes_D$`Total Treatment`-eval_outcomes_B$`Total Treatment`)
summary(eval_outcomes_D$`Total Treatment`-eval_outcomes_C$`Total Treatment`)
summary(eval_outcomes_C$`Total Treatment`-eval_outcomes_A$`Total Treatment`)
                    
#Outcome 2 - Incident TB compared
par(mfrow = c(1,4))
hist(incident_TB_D-incident_TB_A,main = "Algorithm D-A")
hist(incident_TB_D-incident_TB_B,main = "Algorithm D-B")
hist(incident_TB_D-incident_TB_C,main = "Algorithm D-C")
hist(incident_TB_C-incident_TB_A,main = "Algorithm C-A")

summary(unlist(mapply(FUN = function(a,b) sum(b$expected_incidence_total) - sum(a$expected_incidence_total), a=long_outcome_A_TB, b=long_outcome_D_TB)))
summary(unlist(mapply(FUN = function(a,b) sum(b$expected_incidence_total) - sum(a$expected_incidence_total), a=long_outcome_B_TB, b=long_outcome_D_TB)))
summary(unlist(mapply(FUN = function(a,b) sum(b$expected_incidence_total) - sum(a$expected_incidence_total), a=long_outcome_C_TB, b=long_outcome_D_TB)))
summary(unlist(mapply(FUN = function(a,b) sum(b$expected_incidence_total) - sum(a$expected_incidence_total), a=long_outcome_A_TB, b=long_outcome_C_TB)))

#Outcome 3 - Total TPT compared

par(mfrow = c(1,4))
hist(eval_outcomes_D$`Total TPT`-eval_outcomes_A$`Total TPT`,main = "Algorithm D-A")
hist(eval_outcomes_D$`Total TPT`-eval_outcomes_B$`Total TPT`,main = "Algorithm D-B")
hist(eval_outcomes_D$`Total TPT`-eval_outcomes_C$`Total TPT`,main = "Algorithm D-C")
hist(eval_outcomes_C$`Total TPT`-eval_outcomes_A$`Total TPT`,main = "Algorithm C-A")

summary(eval_outcomes_D$`Total TPT`-eval_outcomes_A$`Total TPT`)
summary(eval_outcomes_D$`Total TPT`-eval_outcomes_B$`Total TPT`)
summary(eval_outcomes_D$`Total TPT`-eval_outcomes_C$`Total TPT`)
summary(eval_outcomes_C$`Total TPT`-eval_outcomes_A$`Total TPT`)

#Outcome 4 - Total Diagnosed TB compared

par(mfrow = c(1,4))
hist(eval_outcomes_D$`Total Diagnosed`-eval_outcomes_A$`Total Diagnosed`,main = "Algorithm D-A")
hist(eval_outcomes_D$`Total Diagnosed`-eval_outcomes_B$`Total Diagnosed`,main = "Algorithm D-B")
hist(eval_outcomes_D$`Total Diagnosed`-eval_outcomes_C$`Total Diagnosed`,main = "Algorithm D-C")
hist(eval_outcomes_C$`Total Diagnosed`-eval_outcomes_A$`Total Diagnosed`,main = "Algorithm C-A")

#Outcome 5 - Unecessary Treatment Compared 

par(mfrow = c(1,4))
hist(Inappropriate_Therapy_D$`Total Incorrect Treatment`-Inappropriate_Therapy_A$`Total Incorrect Treatment`,main = "Algorithm D-A")
hist(Inappropriate_Therapy_D$`Total Incorrect Treatment`-Inappropriate_Therapy_B$`Total Incorrect Treatment`, main = "Algorithm D-B")
hist(Inappropriate_Therapy_D$`Total Incorrect Treatment`-Inappropriate_Therapy_C$`Total Incorrect Treatment`,main = "Algorithm D-C")
hist(Inappropriate_Therapy_C$`Total Incorrect Treatment`-Inappropriate_Therapy_A$`Total Incorrect Treatment`,main = "Algorithm C-A")

summary(Inappropriate_Therapy_D$`Total Incorrect Treatment`-Inappropriate_Therapy_A$`Total Incorrect Treatment`)
summary(Inappropriate_Therapy_D$`Total Incorrect Treatment`-Inappropriate_Therapy_B$`Total Incorrect Treatment`)
summary(Inappropriate_Therapy_D$`Total Incorrect Treatment`-Inappropriate_Therapy_C$`Total Incorrect Treatment`)
summary(Inappropriate_Therapy_C$`Total Incorrect Treatment`-Inappropriate_Therapy_A$`Total Incorrect Treatment`)

#Outcome 6 - Unecessary TPT Compared 
par(mfrow = c(1,4))
hist(Inappropriate_Therapy_D$`Total Incorrect TPT`-Inappropriate_Therapy_A$`Total Incorrect TPT`,main = "Algorithm D-A")
hist(Inappropriate_Therapy_D$`Total Incorrect TPT`-Inappropriate_Therapy_B$`Total Incorrect TPT`,main = "Algorithm D-B")
hist(Inappropriate_Therapy_D$`Total Incorrect TPT`-Inappropriate_Therapy_C$`Total Incorrect TPT`,main = "Algorithm D-C")
hist(Inappropriate_Therapy_C$`Total Incorrect TPT`-Inappropriate_Therapy_A$`Total Incorrect TPT`,main = "Algorithm C-A")


summary(Inappropriate_Therapy_D$`Total Incorrect TPT`-Inappropriate_Therapy_A$`Total Incorrect TPT`)
summary(Inappropriate_Therapy_D$`Total Incorrect TPT`-Inappropriate_Therapy_B$`Total Incorrect TPT`)
summary(Inappropriate_Therapy_D$`Total Incorrect TPT`-Inappropriate_Therapy_C$`Total Incorrect TPT`)
summary(Inappropriate_Therapy_C$`Total Incorrect TPT`-Inappropriate_Therapy_A$`Total Incorrect TPT`)


#Outcome 7- Averted TB Mortality 
summary(unlist(mapply(FUN = function(a,b) (sum(b$expected_incidence_total) - sum(a$expected_incidence_total))*(1/1.5), a=long_outcome_A_TB, b=long_outcome_D_TB)))
summary(unlist(mapply(FUN = function(a,b) (sum(b$expected_incidence_total) - sum(a$expected_incidence_total))*(1/1.5), a=long_outcome_B_TB, b=long_outcome_D_TB)))
summary(unlist(mapply(FUN = function(a,b) (sum(b$expected_incidence_total) - sum(a$expected_incidence_total))*(1/1.5), a=long_outcome_C_TB, b=long_outcome_D_TB)))
summary(unlist(mapply(FUN = function(a,b) (sum(b$expected_incidence_total) - sum(a$expected_incidence_total))*(1/1.5), a=long_outcome_A_TB, b=long_outcome_C_TB)))

hist(unlist(mapply(FUN = function(a,b) (sum(b$expected_incidence_total) - sum(a$expected_incidence_total))*(1/1.5), a=long_outcome_A_TB, b=long_outcome_D_TB)))
hist(unlist(mapply(FUN = function(a,b) (sum(b$expected_incidence_total) - sum(a$expected_incidence_total))*(1/1.5), a=long_outcome_B_TB, b=long_outcome_D_TB)))
hist(unlist(mapply(FUN = function(a,b) (sum(b$expected_incidence_total) - sum(a$expected_incidence_total))*(1/1.5), a=long_outcome_C_TB, b=long_outcome_D_TB)))
hist(unlist(mapply(FUN = function(a,b) (sum(b$expected_incidence_total) - sum(a$expected_incidence_total))*(1/1.5), a=long_outcome_A_TB, b=long_outcome_C_TB)))


#Individual clinical outcomes


#summary TB 
summary(do.call("rbind",lapply(outcome_cohort_A, function(x) length(which(x$TB == "TB")))))

#summary No TB 
summary(do.call("rbind",lapply(outcome_cohort_A, function(x) length(which(x$TB == "No_TB")))))

#Outcome 1 - Treatment compared
par(mfrow = c(1,4))
hist(eval_outcomes_A$`Total Treatment`,main = "Algorithm A")
hist(eval_outcomes_B$`Total Treatment`,main = "Algorithm B")
hist(eval_outcomes_C$`Total Treatment`,main = "Algorithm C")
hist(eval_outcomes_D$`Total Treatment`,main = "Algorithm D")

summary(eval_outcomes_A$`Total Treatment`/do.call("rbind",lapply(outcome_cohort_A, function(x) length(which(x$TB == "TB")))))
summary(eval_outcomes_B$`Total Treatment`/do.call("rbind",lapply(outcome_cohort_B, function(x) length(which(x$TB == "TB")))))
summary(eval_outcomes_C$`Total Treatment`/do.call("rbind",lapply(outcome_cohort_C, function(x) length(which(x$TB == "TB")))))
summary(eval_outcomes_D$`Total Treatment`/do.call("rbind",lapply(outcome_cohort_D, function(x) length(which(x$TB == "TB")))))


#Outcome 2 - Incident TB compared
#total absolute incidence with individual algorithms
summary(unlist(lapply(long_outcome_A_TB, function(x) sum(x$expected_incidence_total)))/do.call("rbind",lapply(outcome_cohort_A, function(x) length(which(x$TB == "No_TB")))))
summary(unlist(lapply(long_outcome_B_TB, function(x) sum(x$expected_incidence_total)))/do.call("rbind",lapply(outcome_cohort_B, function(x) length(which(x$TB == "No_TB")))))
summary(unlist(lapply(long_outcome_C_TB, function(x) sum(x$expected_incidence_total)))/do.call("rbind",lapply(outcome_cohort_C, function(x) length(which(x$TB == "No_TB")))))
summary(unlist(lapply(long_outcome_D_TB, function(x) sum(x$expected_incidence_total)))/do.call("rbind",lapply(outcome_cohort_D, function(x) length(which(x$TB == "No_TB")))))

par(mfrow=c(1,4))
hist(unlist(lapply(long_outcome_A_TB, function(x) sum(x$expected_incidence_total))), main = "Algorithm A")
hist(unlist(lapply(long_outcome_B_TB, function(x) sum(x$expected_incidence_total))), main = "Algorithm B")
hist(unlist(lapply(long_outcome_C_TB, function(x) sum(x$expected_incidence_total))), main = "Algorithm C")
hist(unlist(lapply(long_outcome_D_TB, function(x) sum(x$expected_incidence_total))), main = "Algorithm D")


#Outcome 3 - Total TPT compared

par(mfrow = c(1,4))
hist(eval_outcomes_A$`Total TPT`,main = "Algorithm A")
hist(eval_outcomes_B$`Total TPT`,main = "Algorithm B")
hist(eval_outcomes_C$`Total TPT`,main = "Algorithm C")
hist(eval_outcomes_D$`Total TPT`,main = "Algorithm D")

summary(eval_outcomes_A$`Total TPT`/do.call("rbind",lapply(outcome_cohort_A, function(x) length(which(x$TB == "No_TB")))))
summary(eval_outcomes_B$`Total TPT`/do.call("rbind",lapply(outcome_cohort_B, function(x) length(which(x$TB == "No_TB")))))
summary(eval_outcomes_C$`Total TPT`/do.call("rbind",lapply(outcome_cohort_C, function(x) length(which(x$TB == "No_TB")))))
summary(eval_outcomes_D$`Total TPT`/do.call("rbind",lapply(outcome_cohort_D, function(x) length(which(x$TB == "No_TB")))))


#Outcome 4 - Total Diagnosed TB compared

par(mfrow = c(1,4))
hist(eval_outcomes_A$`Total Diagnosed`,main = "Algorithm A")
hist(eval_outcomes_B$`Total Diagnosed`,main = "Algorithm B")
hist(eval_outcomes_C$`Total Diagnosed`,main = "Algorithm C")
hist(eval_outcomes_D$`Total Diagnosed`,main = "Algorithm D")

#Outcome 5 - Unecessary Treatment Compared 

par(mfrow = c(1,4))
hist(Inappropriate_Therapy_A$`Total Incorrect Treatment`,main = "Algorithm A"/do.call("rbind",lapply(outcome_cohort_A, function(x) length(which(x$TB == "No_TB")))))
hist(Inappropriate_Therapy_B$`Total Incorrect Treatment`,main = "Algorithm B"/do.call("rbind",lapply(outcome_cohort_B, function(x) length(which(x$TB == "No_TB")))))
hist(Inappropriate_Therapy_C$`Total Incorrect Treatment`,main = "Algorithm C"/do.call("rbind",lapply(outcome_cohort_C, function(x) length(which(x$TB == "No_TB")))))
hist(Inappropriate_Therapy_D$`Total Incorrect Treatment`,main = "Algorithm D"/do.call("rbind",lapply(outcome_cohort_D, function(x) length(which(x$TB == "No_TB")))))

summary(Inappropriate_Therapy_A$`Total Incorrect Treatment`/do.call("rbind",lapply(outcome_cohort_A, function(x) length(which(x$TB == "No_TB")))))
summary(Inappropriate_Therapy_B$`Total Incorrect Treatment`/do.call("rbind",lapply(outcome_cohort_B, function(x) length(which(x$TB == "No_TB")))))
summary(Inappropriate_Therapy_C$`Total Incorrect Treatment`/do.call("rbind",lapply(outcome_cohort_C, function(x) length(which(x$TB == "No_TB")))))
summary(Inappropriate_Therapy_D$`Total Incorrect Treatment`/do.call("rbind",lapply(outcome_cohort_D, function(x) length(which(x$TB == "No_TB")))))

#Outcome 6 - Unecessary TPT Compared 
par(mfrow = c(1,4))
hist(Inappropriate_Therapy_A$`Total Incorrect TPT`,main = "Algorithm A")
hist(Inappropriate_Therapy_B$`Total Incorrect TPT`,main = "Algorithm B")
hist(Inappropriate_Therapy_C$`Total Incorrect TPT`,main = "Algorithm C")
hist(Inappropriate_Therapy_D$`Total Incorrect TPT`,main = "Algorithm D")

summary(Inappropriate_Therapy_A$`Total Incorrect TPT`/do.call("rbind",lapply(outcome_cohort_A, function(x) length(which(x$TB == "TB")))))
summary(Inappropriate_Therapy_B$`Total Incorrect TPT`/do.call("rbind",lapply(outcome_cohort_A, function(x) length(which(x$TB == "TB")))))
summary(Inappropriate_Therapy_C$`Total Incorrect TPT`/do.call("rbind",lapply(outcome_cohort_A, function(x) length(which(x$TB == "TB")))))
summary(Inappropriate_Therapy_D$`Total Incorrect TPT`/do.call("rbind",lapply(outcome_cohort_A, function(x) length(which(x$TB == "TB")))))

#Cost effectiveness fixed 

#Arrange treated per cost
par(mfrow = c(1,4))

D_A = (t_cost_D-t_cost_A)/(do.call(rbind,lapply(outcome_cohort_D, function(x) sum(x[,"Treatment"]==1 & as.numeric(as.character(x[,"Day_Treatment"]))<=60 & x[,"TB"]=="TB")))-do.call(rbind,lapply(outcome_cohort_A, function(x) sum(x[,"Treatment"]==1 & as.numeric(as.character(x[,"Day_Treatment"]))<=60 & x[,"TB"]=="TB"))))
D_B = (t_cost_D-t_cost_B)/(do.call(rbind,lapply(outcome_cohort_D, function(x) sum(x[,"Treatment"]==1 & as.numeric(as.character(x[,"Day_Treatment"]))<=60 & x[,"TB"]=="TB")))-do.call(rbind,lapply(outcome_cohort_B, function(x) sum(x[,"Treatment"]==1 & as.numeric(as.character(x[,"Day_Treatment"]))<=60 & x[,"TB"]=="TB"))))
D_C = (t_cost_D-t_cost_C)/(do.call(rbind,lapply(outcome_cohort_D, function(x) sum(x[,"Treatment"]==1 & as.numeric(as.character(x[,"Day_Treatment"]))<=60 & x[,"TB"]=="TB")))-do.call(rbind,lapply(outcome_cohort_C, function(x) sum(x[,"Treatment"]==1 & as.numeric(as.character(x[,"Day_Treatment"]))<=60 & x[,"TB"]=="TB"))))
C_A = (t_cost_C-t_cost_A)/(do.call(rbind,lapply(outcome_cohort_C, function(x) sum(x[,"Treatment"]==1 & as.numeric(as.character(x[,"Day_Treatment"]))<=60 & x[,"TB"]=="TB")))-do.call(rbind,lapply(outcome_cohort_A, function(x) sum(x[,"Treatment"]==1 & as.numeric(as.character(x[,"Day_Treatment"]))<=60 & x[,"TB"]=="TB"))))

D_A[which(D_A<=0)]<-1000000
D_B[which(D_B<=0)]<-1000000
D_C[which(D_C<=0)]<-1000000
C_A[which(C_A<=0)]<-1000000

summary(D_A)
summary(D_B)
summary(D_C)
summary(C_A)

#Arrange Incident TB per cost

#comparing  algorithms 
D_A =  (t_cost_D-t_cost_A)/unlist(mapply(FUN = function(a,b) sum(b$expected_incidence_total) - sum(a$expected_incidence_total), a=long_outcome_D_TB, b=long_outcome_A_TB))
D_B =  (t_cost_D-t_cost_B)/unlist(mapply(FUN = function(a,b) sum(b$expected_incidence_total) - sum(a$expected_incidence_total), a=long_outcome_D_TB, b=long_outcome_B_TB))
D_C =  (t_cost_D-t_cost_C)/unlist(mapply(FUN = function(a,b) sum(b$expected_incidence_total) - sum(a$expected_incidence_total), a=long_outcome_D_TB, b=long_outcome_C_TB))
C_A =  (t_cost_C-t_cost_A)/unlist(mapply(FUN = function(a,b) sum(b$expected_incidence_total) - sum(a$expected_incidence_total), a=long_outcome_C_TB, b=long_outcome_A_TB))

D_A[which(D_A<=0)]<-1000000
D_B[which(D_B<=0)]<-1000000
D_C[which(D_C<=0)]<-1000000
C_A[which(C_A<=0)]<-1000000

summary(D_A)
summary(D_B)
summary(D_C)
summary(C_A)

#Arrange TPT per cost

TPT_A_6_months = do.call("rbind",lapply(long_outcome_A_TB,function(x) length(which(x$TPT_Stop_Day>180))))
TPT_B_6_months = do.call("rbind",lapply(long_outcome_B_TB,function(x) length(which(x$TPT_Stop_Day>180))))
TPT_C_6_months = do.call("rbind",lapply(long_outcome_C_TB,function(x) length(which(x$TPT_Stop_Day>180))))
TPT_D_6_months = do.call("rbind",lapply(long_outcome_D_TB,function(x) length(which(x$TPT_Stop_Day>180))))

D_A = (t_cost_D-t_cost_A)/(TPT_D_6_months-TPT_A_6_months)
D_B = (t_cost_D-t_cost_B)/(TPT_D_6_months-TPT_B_6_months)
D_C = (t_cost_D-t_cost_C)/(TPT_D_6_months-TPT_C_6_months)
C_A = (t_cost_C-t_cost_A)/(TPT_C_6_months-TPT_A_6_months)

D_A[which(D_A<=0)]<-10000000
D_B[which(D_B<=0)]<-10000000
D_C[which(D_C<=0)]<-10000000
C_A[which(C_A<=0)]<-10000000

summary(D_A)
summary(D_B)
summary(D_C)
summary(C_A)


#Averted TB mortality 

D_A = (t_cost_D-t_cost_A)/(unlist(mapply(FUN = function(a,b) (sum(b$expected_incidence_total) - sum(a$expected_incidence_total))*(1/1.5), a=long_outcome_D_TB, b=long_outcome_A_TB)))
D_B = (t_cost_D-t_cost_B)/(unlist(mapply(FUN = function(a,b) (sum(b$expected_incidence_total) - sum(a$expected_incidence_total))*(1/1.5), a=long_outcome_D_TB, b=long_outcome_B_TB)))
D_C = (t_cost_D-t_cost_C)/(unlist(mapply(FUN = function(a,b) (sum(b$expected_incidence_total) - sum(a$expected_incidence_total))*(1/1.5), a=long_outcome_D_TB, b=long_outcome_C_TB)))
C_A = (t_cost_C-t_cost_A)/(unlist(mapply(FUN = function(a,b) (sum(b$expected_incidence_total) - sum(a$expected_incidence_total))*(1/1.5), a=long_outcome_C_TB, b=long_outcome_A_TB)))

D_A[which(D_A<=0)]<-1000000 
D_B[which(D_B<=0)]<-1000000 
D_C[which(D_C<=0)]<-1000000 
C_A[which(C_A<=0)]<-1000000 

summary(D_A)
summary(D_B)
summary(D_C)
summary(C_A)

#Cascade of Care TPT

#Proportion Starting TPT at visit 1, of no Baseline TB 
TPT_A_Visit1 = do.call("rbind",lapply(long_outcome_A_TB,function(x) sum(x$TPT1)/dim(x)[1]))
TPT_B_Visit1 = do.call("rbind",lapply(long_outcome_B_TB,function(x) sum(x$TPT1)/dim(x)[1]))
TPT_C_Visit1 = do.call("rbind",lapply(long_outcome_C_TB,function(x) sum(x$TPT1)/dim(x)[1]))
TPT_D_Visit1 = do.call("rbind",lapply(long_outcome_D_TB,function(x) sum(x$TPT1)/dim(x)[1]))

#Proportion starting TPT at visit 1 or 2, of no Baseline TB 

TPT_A_Visit1_2 = do.call("rbind", lapply(long_outcome_A_TB,function(x) length(which(x$TPT1 == 1 | x$TPT2 == 1))/dim(x)[1]))
TPT_B_Visit1_2 = do.call("rbind", lapply(long_outcome_B_TB,function(x) length(which(x$TPT1 == 1 | x$TPT2 == 1))/dim(x)[1]))
TPT_C_Visit1_2 = do.call("rbind", lapply(long_outcome_C_TB,function(x) length(which(x$TPT1 == 1 | x$TPT2 == 1))/dim(x)[1]))
TPT_D_Visit1_2 = do.call("rbind", lapply(long_outcome_D_TB,function(x) length(which(x$TPT1 == 1 | x$TPT2 == 1))/dim(x)[1]))

#Proportion starting TPT at visit 1 or 2, and goes to follow up visit, of no Baseline TB 

TPT_A_Visit1_2_Follow = do.call("rbind", lapply(long_outcome_A_TB, function(x) length(which((x$TPT1 == 1 | x$TPT2 == 1) & !is.na(x$Visit2)))/dim(x)[1]))
TPT_B_Visit1_2_Follow = do.call("rbind", lapply(long_outcome_B_TB, function(x) length(which((x$TPT1 == 1 | x$TPT2 == 1) & !is.na(x$Visit2)))/dim(x)[1]))
TPT_C_Visit1_2_Follow = do.call("rbind", lapply(long_outcome_C_TB, function(x) length(which((x$TPT1 == 1 | x$TPT2 == 1) & !is.na(x$Visit2)))/dim(x)[1]))
TPT_D_Visit1_2_Follow = do.call("rbind", lapply(long_outcome_D_TB, function(x) length(which((x$TPT1 == 1 | x$TPT2 == 1) & !is.na(x$Visit2)))/dim(x)[1]))

#Proportion continuing on TPT beyond visit 2
TPT_A_Visit2 = do.call("rbind",lapply(long_outcome_A_TB, function(x) sum(x$TPT2)/dim(x)[1]))
TPT_B_Visit2 = do.call("rbind",lapply(long_outcome_B_TB, function(x) sum(x$TPT2)/dim(x)[1]))
TPT_C_Visit2 = do.call("rbind",lapply(long_outcome_C_TB, function(x) sum(x$TPT2)/dim(x)[1]))
TPT_D_Visit2 = do.call("rbind",lapply(long_outcome_D_TB, function(x) sum(x$TPT2)/dim(x)[1]))

#Proportion completing 6 months of TPT

TPT_A_6_months = do.call("rbind",lapply(long_outcome_A_TB,function(x) length(which(x$TPT_Stop_Day>180))/dim(x)[1]))
TPT_B_6_months = do.call("rbind",lapply(long_outcome_B_TB,function(x) length(which(x$TPT_Stop_Day>180))/dim(x)[1]))
TPT_C_6_months = do.call("rbind",lapply(long_outcome_C_TB,function(x) length(which(x$TPT_Stop_Day>180))/dim(x)[1]))
TPT_D_6_months = do.call("rbind",lapply(long_outcome_D_TB,function(x) length(which(x$TPT_Stop_Day>180))/dim(x)[1]))

summary(do.call("rbind",lapply(long_outcome_D_TB,function(x) length(which(x$TPT_Stop_Day>180)))) - do.call("rbind",lapply(long_outcome_A_TB,function(x) length(which(x$TPT_Stop_Day>180)))))
summary(do.call("rbind",lapply(long_outcome_D_TB,function(x) length(which(x$TPT_Stop_Day>180)))) - do.call("rbind",lapply(long_outcome_B_TB,function(x) length(which(x$TPT_Stop_Day>180)))))
summary(do.call("rbind",lapply(long_outcome_D_TB,function(x) length(which(x$TPT_Stop_Day>180)))) - do.call("rbind",lapply(long_outcome_C_TB,function(x) length(which(x$TPT_Stop_Day>180)))))
summary(do.call("rbind",lapply(long_outcome_C_TB,function(x) length(which(x$TPT_Stop_Day>180)))) - do.call("rbind",lapply(long_outcome_A_TB,function(x) length(which(x$TPT_Stop_Day>180)))))

mean(1-(do.call("rbind",lapply(long_outcome_A_TB,function(x) length(which(x$TPT_Stop_Day>180))))/do.call("rbind", lapply(long_outcome_A_TB,function(x) length(which(x$TPT1 == 1 | x$TPT2 == 1))))))
mean(1-(do.call("rbind",lapply(long_outcome_B_TB,function(x) length(which(x$TPT_Stop_Day>180))))/do.call("rbind", lapply(long_outcome_B_TB,function(x) length(which(x$TPT1 == 1 | x$TPT2 == 1))))))
mean(1-(do.call("rbind",lapply(long_outcome_C_TB,function(x) length(which(x$TPT_Stop_Day>180))))/do.call("rbind", lapply(long_outcome_C_TB,function(x) length(which(x$TPT1 == 1 | x$TPT2 == 1))))))
mean(1-(do.call("rbind",lapply(long_outcome_D_TB,function(x) length(which(x$TPT_Stop_Day>180))))/do.call("rbind", lapply(long_outcome_D_TB,function(x) length(which(x$TPT1 == 1 | x$TPT2 == 1))))))


#Proportion completing 12 months of TPT 

TPT_A_12_months = do.call("rbind",lapply(long_outcome_A_TB,function(x) length(which(x$TPT_Stop_Day>360))/dim(x)[1]))
TPT_B_12_months = do.call("rbind",lapply(long_outcome_B_TB,function(x) length(which(x$TPT_Stop_Day>360))/dim(x)[1]))
TPT_C_12_months = do.call("rbind",lapply(long_outcome_C_TB,function(x) length(which(x$TPT_Stop_Day>360))/dim(x)[1]))
TPT_D_12_months = do.call("rbind",lapply(long_outcome_D_TB,function(x) length(which(x$TPT_Stop_Day>360))/dim(x)[1]))

mean(1- (do.call("rbind",lapply(long_outcome_A_TB,function(x) length(which(x$TPT_Stop_Day>180))))/do.call("rbind", lapply(long_outcome_A_TB, function(x) length(which((x$TPT1 == 1 | x$TPT2 == 1) & !is.na(x$Visit2)))))))
mean(1- (do.call("rbind",lapply(long_outcome_B_TB,function(x) length(which(x$TPT_Stop_Day>180))))/do.call("rbind", lapply(long_outcome_B_TB, function(x) length(which((x$TPT1 == 1 | x$TPT2 == 1) & !is.na(x$Visit2)))))))
mean(1- (do.call("rbind",lapply(long_outcome_C_TB,function(x) length(which(x$TPT_Stop_Day>180))))/do.call("rbind", lapply(long_outcome_C_TB, function(x) length(which((x$TPT1 == 1 | x$TPT2 == 1) & !is.na(x$Visit2)))))))
mean(1- (do.call("rbind",lapply(long_outcome_D_TB,function(x) length(which(x$TPT_Stop_Day>180))))/do.call("rbind", lapply(long_outcome_D_TB, function(x) length(which((x$TPT1 == 1 | x$TPT2 == 1) & !is.na(x$Visit2)))))))
#Cascades of Care Treatment 

#Proportion tested for TB at visit 1, of those with baseline TB

TB_Test_Visit1_A = do.call("rbind", lapply(outcome_cohort_A, function(x) length(which(x$TB == "TB" & x$`Xpert_Screen` == 1))/length(which(x$TB == "TB"))))
TB_Test_Visit1_B = do.call("rbind", lapply(outcome_cohort_B, function(x) length(which(x$TB == "TB" & x$`Xpert_Screen` == 1))/length(which(x$TB == "TB"))))
TB_Test_Visit1_C = do.call("rbind", lapply(outcome_cohort_C, function(x) length(which(x$TB == "TB" & x$`Xpert_Screen` == 1))/length(which(x$TB == "TB"))))
TB_Test_Visit1_D = do.call("rbind", lapply(outcome_cohort_D, function(x) length(which(x$TB == "TB" & x$`Xpert_Screen` == 1))/length(which(x$TB == "TB"))))

#Proportion with positive Xpert test, of those with baseline TB

TB_Xp_Pos_Visit1_A = do.call("rbind", lapply(outcome_cohort_A, function(x) length(which(x$TB == "TB" & x$`Xpert_Screen` == 1 & x$Xpert == "Xpert+"))/length(which(x$TB == "TB"))))
TB_Xp_Pos_Visit1_B = do.call("rbind", lapply(outcome_cohort_B, function(x) length(which(x$TB == "TB" & x$`Xpert_Screen` == 1 & x$Xpert == "Xpert+"))/length(which(x$TB == "TB"))))
TB_Xp_Pos_Visit1_C = do.call("rbind", lapply(outcome_cohort_C, function(x) length(which(x$TB == "TB" & x$`Xpert_Screen` == 1 & x$Xpert == "Xpert+"))/length(which(x$TB == "TB"))))
TB_Xp_Pos_Visit1_D = do.call("rbind", lapply(outcome_cohort_D, function(x) length(which(x$TB == "TB" & x$`Xpert_Screen` == 1 & x$Xpert == "Xpert+"))/length(which(x$TB == "TB"))))

#Proportion returning to Visit 2, of those with baseline TB

TB_Visit_2_A = do.call("rbind",lapply(outcome_cohort_A, function(x) length(which(x$TB == "TB" & x$`Xpert_Screen` == 1 & x$Xpert == "Xpert+" & !is.na(x$t)))/length(which(x$TB == "TB"))))
TB_Visit_2_B = do.call("rbind",lapply(outcome_cohort_B, function(x) length(which(x$TB == "TB" & x$`Xpert_Screen` == 1 & x$Xpert == "Xpert+" & !is.na(x$t)))/length(which(x$TB == "TB"))))
TB_Visit_2_C = do.call("rbind",lapply(outcome_cohort_C, function(x) length(which(x$TB == "TB" & x$`Xpert_Screen` == 1 & x$Xpert == "Xpert+" & !is.na(x$t)))/length(which(x$TB == "TB"))))
TB_Visit_2_D = do.call("rbind",lapply(outcome_cohort_D, function(x) length(which(x$TB == "TB" & x$`Xpert_Screen` == 1 & x$Xpert == "Xpert+" & !is.na(x$t)))/length(which(x$TB == "TB"))))


#Proportion starting treatment at Visit 2 

TB_Treatment_Visit_2_A = do.call("rbind",lapply(outcome_cohort_A, function(x) length(which(x$TB == "TB" & x$`Xpert_Screen` == 1 & x$Xpert == "Xpert+" & !is.na(x$t) & !is.na(x$Day_Treatment)))/length(which(x$TB == "TB"))))
TB_Treatment_Visit_2_B = do.call("rbind",lapply(outcome_cohort_B, function(x) length(which(x$TB == "TB" & x$`Xpert_Screen` == 1 & x$Xpert == "Xpert+" & !is.na(x$t) & !is.na(x$Day_Treatment)))/length(which(x$TB == "TB"))))
TB_Treatment_Visit_2_C = do.call("rbind",lapply(outcome_cohort_C, function(x) length(which(x$TB == "TB" & x$`Xpert_Screen` == 1 & x$Xpert == "Xpert+" & !is.na(x$t) & !is.na(x$Day_Treatment)))/length(which(x$TB == "TB"))))
TB_Treatment_Visit_2_D = do.call("rbind",lapply(outcome_cohort_D, function(x) length(which(x$TB == "TB" & x$`Xpert_Screen` == 1 & x$Xpert == "Xpert+" & !is.na(x$t) & !is.na(x$Day_Treatment)))/length(which(x$TB == "TB"))))

#TPT Cascade:
  # % TPT Visit 1 of those without baseline TB
  # % completing visit 2
  # % on TPT after visit 2
  # % completing at least 6 months of TPT
  # % completing 12 months of TPT

TPT_Cascade_V1 = rbind(cbind(mean(TPT_A_Visit1),"Symptom Screen","TPT Started at Visit 1"), cbind(mean(TPT_B_Visit1),"CRP Screen","TPT Started at Visit 1"), cbind(mean(TPT_C_Visit1),"Universal Xpert + Delayed TPT","TPT Started at Visit 1"), cbind(mean(TPT_D_Visit1),"Universal Xpert + TPT", "TPT Started at Visit 1"))
TPT_Cascade_V12 = rbind(cbind(mean(TPT_A_Visit1_2),"Symptom Screen","TPT Started at Anytime (Visit 1 or 2)"),cbind(mean(TPT_B_Visit1_2),"CRP Screen","TPT Started at Anytime (Visit 1 or 2)"),cbind(mean(TPT_C_Visit1_2),"Universal Xpert + Delayed TPT","TPT Started at Anytime (Visit 1 or 2)"),cbind(mean(TPT_D_Visit1_2),"Universal Xpert + TPT","TPT Started at Anytime (Visit 1 or 2)"))
TPT_Cascade_Follow = rbind(cbind(mean(TPT_A_Visit1_2_Follow),"Symptom Screen","TPT Started at Anytime, and Attends Visit 2"),cbind(mean(TPT_B_Visit1_2_Follow),"CRP Screen","TPT Started at Anytime, and Attends Visit 2"),cbind(mean(TPT_C_Visit1_2_Follow),"Universal Xpert + Delayed TPT","TPT Started at Anytime, and Attends Visit 2"),cbind(mean(TPT_D_Visit1_2_Follow),"Universal Xpert + TPT","TPT Started at Anytime, and Attends Visit 2"))
TPT_Cascade_Cont = rbind(cbind(mean(TPT_A_Visit2),"Symptom Screen","TPT Started or Continued at Visit 2"),cbind(mean(TPT_B_Visit2),"CRP Screen","TPT Started or Continued at Visit 2"),cbind(mean(TPT_C_Visit2),"Universal Xpert + Delayed TPT","TPT Started or Continued at Visit 2"),cbind(mean(TPT_D_Visit2),"Universal Xpert + TPT","TPT Started or Continued at Visit 2"))
TPT_Cascade_6 = rbind(cbind(mean(TPT_A_6_months),"Symptom Screen", "Completes at least 6 Months TPT"),cbind(mean(TPT_B_6_months),"CRP Screen", "Completes at least 6 Months TPT"),cbind(mean(TPT_C_6_months),"Universal Xpert + Delayed TPT", "Completes at least 6 Months TPT"),cbind(mean(TPT_D_6_months),"Universal Xpert + TPT", "Completes at least 6 Months TPT"))
TPT_Cascade_12 = rbind(cbind(mean(TPT_A_12_months), "Symptom Screen", "Completes at least 12 Months TPT"),cbind(mean(TPT_B_12_months), "CRP Screen", "Completes at least 12 Months TPT"),cbind(mean(TPT_C_12_months), "Universal Xpert + Delayed TPT", "Completes at least 12 Months TPT"),cbind(mean(TPT_D_12_months), "Universal Xpert + TPT", "Completes at least 12 Months TPT"))

#Remove TPT_Cascade_Follow
TPT_Cascade = rbind(TPT_Cascade_V1,TPT_Cascade_V12,TPT_Cascade_Cont,TPT_Cascade_6,TPT_Cascade_12)
colnames(TPT_Cascade) = c("Value","Algorithm","Care Cascade")
TPT_Cascade = as.data.frame(TPT_Cascade)
TPT_Cascade$Value = as.numeric(as.character(TPT_Cascade$Value))
TPT_Cascade$`Care Cascade` = as.character(TPT_Cascade$`Care Cascade`)
TPT_Cascade$`Care Cascade` = factor(TPT_Cascade$`Care Cascade`, levels = unique(TPT_Cascade$`Care Cascade`))
TPT_Cascade$`Algorithm` = factor(TPT_Cascade$`Algorithm`, levels = c("Symptom Screen", "CRP Screen", "Universal Xpert + Delayed TPT", "Universal Xpert + TPT"))
ggplot(TPT_Cascade, aes(fill=`Algorithm`, y=`Value`, x = `Care Cascade`))  +
  geom_bar(position = "dodge",stat="identity") + ggtitle("TPT Care Cascade") + xlab("Care Cascade") + ylab("Average Proportion")+ theme(axis.text.x = element_text(angle = 90),axis.text = element_text(size = 10))


#Treatent Cascade:
  # % tested for TB at visit 1
  # % with positive Xpert result
  # % returning to visit 2
  # % starting treatment at visit 2

Treatment_Cascade_V1 = rbind(cbind(mean(TB_Test_Visit1_A),"Symptom Screen","Receives TB Test at Visit 1"),cbind(mean(TB_Test_Visit1_B),"CRP Screen","Receives TB Test at Visit 1"),cbind(mean(TB_Test_Visit1_C),"Universal Xpert + Delayed TPT","Receives TB Test at Visit 1"),cbind(mean(TB_Test_Visit1_D),"Universal Xpert + TPT","Receives TB Test at Visit 1") )
Treatment_Cascade_Xp_Pos = rbind(cbind(mean(TB_Xp_Pos_Visit1_A),"Symptom Screen","Has Positive Xpert Test"),cbind(mean(TB_Xp_Pos_Visit1_B),"CRP Screen","Has Positive Xpert Test"),cbind(mean(TB_Xp_Pos_Visit1_C),"Universal Xpert + Delayed TPT","Has Positive Xpert Test"),cbind(mean(TB_Xp_Pos_Visit1_D),"Universal Xpert + TPT","Has Positive Xpert Test"))
Treatment_Cascade_V2 = rbind(cbind(mean(TB_Visit_2_A),"Symptom Screen","Has Positive Xpert and Returns to Visit 2"),cbind(mean(TB_Visit_2_B),"CRP Screen","Has Positive Xpert and Returns to Visit 2"),cbind(mean(TB_Visit_2_C),"Universal Xpert + Delayed TPT","Has Positive Xpert and Returns to Visit 2"),cbind(mean(TB_Visit_2_D),"Universal Xpert + TPT","Has Positive Xpert and Returns to Visit 2"))
Treatment_Cascade_V2_Treatment = rbind(cbind(mean(TB_Treatment_Visit_2_A),"Symptom Screen", "Starts Treatment Visit 2"),cbind(mean(TB_Treatment_Visit_2_B),"CRP Screen", "Starts Treatment Visit 2" ),cbind(mean(TB_Treatment_Visit_2_C),"Universal Xpert + Delayed TPT", "Starts Treatment Visit 2"),cbind(mean(TB_Treatment_Visit_2_D),"Universal Xpert + TPT", "Starts Treatment Visit 2"))

#Remove Treatment_Cascade_V2
Treatment_Cascade = rbind(Treatment_Cascade_V1,Treatment_Cascade_Xp_Pos,Treatment_Cascade_V2_Treatment)
colnames(Treatment_Cascade) = c("Value","Algorithm","Care Cascade")
Treatment_Cascade = as.data.frame(Treatment_Cascade)
Treatment_Cascade$Value = as.numeric(as.character(Treatment_Cascade$Value))
Treatment_Cascade$`Care Cascade` = as.character(Treatment_Cascade$`Care Cascade`)
Treatment_Cascade$`Care Cascade` = factor(Treatment_Cascade$`Care Cascade`, levels = unique(Treatment_Cascade$`Care Cascade`))
Treatment_Cascade$`Algorithm` <- factor(Treatment_Cascade$`Algorithm`, levels = c("Symptom Screen", "CRP Screen", "Universal Xpert + Delayed TPT", "Universal Xpert + TPT"))


ggplot(Treatment_Cascade, aes(fill=`Algorithm`, y=`Value`, x=`Care Cascade`)) + 
  geom_bar(position = "dodge",stat="identity") + ggtitle("Treatment Care Cascade") + xlab("Care Cascade") + ylab("Average Proportion") + theme(axis.text.x = element_text(angle = 90), axis.text = element_text(size = 10))

