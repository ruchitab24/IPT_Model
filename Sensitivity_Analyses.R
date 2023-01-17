sourcelocation <- "~/Documents/JHU/IPT Model/"
source(paste0(sourcelocation,"Model Part 4.R"))

library(gridExtra)
library(grid)
library(tornado)
library(ggplot2)
library(sensitivity)
library(ppcor)


 outcomes_incidence = unlist(mapply(FUN = function(a,b) sum(b$expected_incidence_total) - sum(a$expected_incidence_total), a=long_outcome_C_TB, b=long_outcome_D_TB))
 outcomes_treatment = eval_outcomes_D$`Total Treatment`-eval_outcomes_A$`Total Treatment`
 outcomes = cbind(outcomes_incidence, outcomes_treatment)
 
 parameters = cbind(parameter_dist_1[,-c(1,2,7,8,13)],parameter_dist[,-c(1,4,5,6,7,8)],Incidence_Multiplier,TPT_Multiplier) #remove ltf, add incidence uncertainity, tpt uncertainity, baseline charachteristics 

 #Incidence Difference
 
 # PRCC - Outcomes Incidence - None are greater than .2 so include all, rank in terms of PRCC 
 
 prcc <- pcor(cbind(parameters,outcomes_incidence),method = "spearman")$estimate[ncol(parameters)+1,-(ncol(parameters)+1)]

 #Reorder parameters
 parameters_incidence = matrix(nrow = dim(parameters)[1], ncol = dim(parameters)[2])
 for(a in 1:length(parameters)){
   parameters_incidence[,a] = parameters[,order(abs(prcc))[a]]
 }
 colnames(parameters_incidence) = colnames(parameters)[order(abs(prcc))]
 parameters_incidence = as.data.frame(parameters_incidence)

 b <- boxplot(outcomes_incidence)
 par(mfrow=c(1,1))
 a <- 55
 b$stats <- array(quantile(outcomes_incidence, c(0.025,0.25,0.5,0.75,0.975)), dim=c(5,1))
 bxp(b, outline = F, horizontal = T, xlim=c(0,a+2),  at=a, ylim = c(-10,0), varwidth=T,  boxfill="orange")
 abline(v=0, col="gray"); abline(v=median(outcomes_incidence), col="gray")
 bxp(b, outline = F, horizontal = T, xlim=c(0,a+2), at=a, ylim = c(-10,0), varwidth=T, add=T,
      xaxt='n', boxfill="orange")
 mtext("Sensitivity of Difference in Incident TB Cases between the Universal Xpert and TPT Algorithm and the Universal Xpert and Delayed TPT Algorithm",side=1,line=4); 
 a <- a-2
 
 for (param in 1:length(parameters_incidence))
 {
   b$stats <- array(quantile(subset(outcomes_incidence, parameters_incidence[,param] > quantile(parameters_incidence[,param],0.8)), c(0.025,0.25,0.5,0.75,0.975)), dim=c(5,1))
   bxp(b, horizontal=T, xaxt='n', add=T, at=a-0.25, outline=F, boxfill="darkorange3", varwidth=T)
   b$stats <- array(quantile(subset(outcomes_incidence, parameters_incidence[,param] < quantile(parameters_incidence[,param],0.2)), c(0.025,0.25,0.5,0.75,0.975)), dim=c(5,1))
   bxp(b, horizontal=T, xaxt='n', add=T, at=a-0.75, outline=F, boxfill="yellow", varwidth=T)
   a <- a-2
 }
 
 axis(side=2, at=c(1 + 1 + 2*length(parameters_incidence), 2*(length(parameters_incidence):1)-0.5), 
      labels=c("All Simulations", colnames(parameters_incidence)),
      las=1, cex.axis=0.5)
 
 dev.off()

 #Treatment Differnece
 
 # PRCC - Outcomes Treatment - None are greater than .2 so include all, rank in terms of PRCC 
 
 prcc <- pcor(cbind(parameters,outcomes_treatment),method = "spearman")$estimate[ncol(parameters)+1,-(ncol(parameters)+1)]
 
 #Reorder parameters
 parameters_treatment = matrix(nrow = dim(parameters)[1], ncol = dim(parameters)[2])
 for(a in 1:length(parameters)){
   parameters_treatment[,a] = parameters[,order(abs(prcc))[a]]
 }
 colnames(parameters_treatment) = colnames(parameters)[order(abs(prcc))]
 parameters_treatment = as.data.frame(parameters_treatment)
 
 b <- boxplot(outcomes_treatment)
 par(mfrow=c(1,1))
 a <- 55
 b$stats <- array(quantile(outcomes_treatment, c(0.025,0.25,0.5,0.75,0.975)), dim=c(5,1))
 bxp(b, outline = F, horizontal = T, xlim=c(0,a+2),  at=a, ylim = c(0,50), varwidth=T,  boxfill="orange")
 abline(v=0, col="gray"); abline(v=median(outcomes_treatment), col="gray")
 bxp(b, outline = F, horizontal = T, xlim=c(0,a+2), at=a, ylim = c(0,50), varwidth=T, add=T,
     xaxt='n', boxfill="orange")
 mtext("Sensitivity of Difference in Treatment Initiation between the Universal Xpert and TPT Algorithm and the Symptom Screen Algorithm",side=1,line=4); 
 a <- a-2
 
 for (param in 1:length(parameters_treatment))
 {
   b$stats <- array(quantile(subset(outcomes_treatment, parameters_treatment[,param] > quantile(parameters_treatment[,param],0.8)), c(0.025,0.25,0.5,0.75,0.975)), dim=c(5,1))
   bxp(b, horizontal=T, xaxt='n', add=T, at=a-0.25, outline=F, boxfill="darkorange3", varwidth=T)
   b$stats <- array(quantile(subset(outcomes_treatment, parameters_treatment[,param] < quantile(parameters_treatment[,param],0.2)), c(0.025,0.25,0.5,0.75,0.975)), dim=c(5,1))
   bxp(b, horizontal=T, xaxt='n', add=T, at=a-0.75, outline=F, boxfill="yellow", varwidth=T)
   a <- a-2
 }
 
 axis(side=2, at=c(1 + 1 + 2*length(parameters_treatment), 2*(length(parameters_treatment):1)-0.5), 
      labels=c("All Simulations", colnames(parameters_treatment)),
      las=1, cex.axis=0.5)
 
 
 
 dev.off()
 
 
  
