# make the source location modifiable, for running on different machines

sourcelocation <- "~/Documents/JHU/IPT Model/"
#"~/OneDrive - Johns Hopkins/Research/TPT and inital Xpert in PWH model/Model/" #Emily's location

#source the baseline_cohort; 
#  creates cohort_sims (a list of nsims dataframes, each with ncohort observations of 5 patient characteristic variables)
source(paste0(sourcelocation,"Model Part 1.R"))

ncohort = 5000
nsims = 1000

ncores=7

#Subset out no TB individuals for each outcome cohort Algs A,B,C,D

#Round all values 

outcome_A = lapply(outcome_cohort_A, function(x) x[which(x$TB == "No_TB"),])
outcome_B = lapply(outcome_cohort_B, function(x) x[which(x$TB == "No_TB"),])
outcome_C = lapply(outcome_cohort_C, function(x) x[which(x$TB == "No_TB"),])
outcome_D = lapply(outcome_cohort_D, function(x) x[which(x$TB == "No_TB"),])


#Create a subsetted dataframe with longitudonal outcomes we need to track 

create_long_outcome <-function(outcome_cohort){
  
  long_outcome = matrix(nrow = dim(outcome_cohort)[1], ncol = 8)
  colnames(long_outcome) = c("CD4","TPT1","ART1","TPT2","ART2","Visit2","ART_Stop_Day","TPT_Stop_Day")
  long_outcome = as.data.frame(long_outcome)
  
  #Add discontinued column to outcome cohort
  outcome_cohort[,"Discon"] <-NA
  
  #CD4?
  long_outcome$CD4 = outcome_cohort$CD4
  
  #Time of visit 2 or NA? 
  long_outcome$Visit2 = round(as.numeric(as.character(outcome_cohort$t)),0)
  
  
  #TPT at visit 1 or 2 or never?
  
  
  outcome_cohort$Discon[which(outcome_cohort$Day_TPT=="discontinued")] = "discontinued"
  outcome_cohort$Day_TPT[which(outcome_cohort$Day_TPT=="discontinued")] <- 0
  outcome_cohort$Day_TPT = round(as.numeric(as.character(outcome_cohort$Day_TPT)),0)
  
  long_outcome$TPT1[which(outcome_cohort$Day_TPT==0)] <- 1 
  long_outcome$TPT1[which(outcome_cohort$Day_TPT>0)] <- 0
  long_outcome$TPT1[is.na(outcome_cohort$Day_TPT)] <- 0
  
  
  long_outcome$TPT2[which(outcome_cohort$Day_TPT >= 0 )] <- 1 
  long_outcome$TPT2[which(outcome_cohort$Discon == "discontinued")] <-0
  long_outcome$TPT2[is.na(outcome_cohort$Day_TPT)] <- 0
  long_outcome$TPT2[is.na(outcome_cohort$t)] <- 0
  
  #Overwrite with Treatment, assume that treatment is the same as TPT, Treatment can only be administered post visit 1
  
  outcome_cohort$Day_Treatment = as.numeric(as.character(outcome_cohort$Day_Treatment))
  
  long_outcome$TPT2[which(outcome_cohort$Day_Treatment > 0 )] <- 1
  
  
  #ART at visit 1 or 2 or never?
  
  outcome_cohort$Day_ART = round(as.numeric(as.character(outcome_cohort$Day_ART)),0)
  
  long_outcome$ART1[which(outcome_cohort$Day_ART ==0)] <- 1
  long_outcome$ART1[which(outcome_cohort$Day_ART !=0)] <- 0
  long_outcome$ART1[is.na(outcome_cohort$Day_ART)] <- 0
  
  
  long_outcome$ART2[which(outcome_cohort$Day_ART >= 0)] <- 1
  long_outcome$ART2[is.na(outcome_cohort$Day_ART)] <- 0
  long_outcome$ART2[is.na(outcome_cohort$t)] <- 0
  
  #Time of TPT stop?
  
  long_outcome$TPT_Stop_Day = outcome_cohort$TPT_Stop_Day 
  long_outcome$TPT_Stop_Day = as.numeric(as.character(long_outcome$TPT_Stop_Day))
  
  #Double check this? potentially remove
  long_outcome$TPT_Stop_Day[which(long_outcome$TPT2 == 1 & long_outcome$TPT1 == 0)] <- long_outcome$TPT_Stop_Day[which(long_outcome$TPT2 == 1 & long_outcome$TPT1 == 0)] - long_outcome$Visit2[which(long_outcome$TPT2 == 1 & long_outcome$TPT1 == 0)] #add to visit 2 date 
  long_outcome$TPT_Stop_Day[which(long_outcome$TPT2 == 1 & long_outcome$TPT1 == 1)] <- long_outcome$TPT_Stop_Day[which(long_outcome$TPT2 == 1 & long_outcome$TPT1 == 1)] - long_outcome$Visit2[which(long_outcome$TPT2 == 1 & long_outcome$TPT1 == 0)] + min(30, long_outcome$Visit2[which(long_outcome$TPT2 == 1 & long_outcome$TPT1 == 0)]) #add to visit 2 date 
  long_outcome$TPT_Stop_Day[which(long_outcome$TPT2 == 0 & long_outcome$TPT1 == 1)] <-30 #Made this 30 days 
  long_outcome$TPT_Stop_Day[which(long_outcome$TPT2 == 0 & long_outcome$TPT1 == 0)] <-0 
  long_outcome$TPT_Stop_Day[which(long_outcome$TPT_Stop_Day > 730)] <- 730
  
  #Time of ART stop?
  
  long_outcome$ART_Stop_Day = outcome_cohort$ART_Stop_Day
  long_outcome$ART_Stop_Day = as.numeric(as.character(long_outcome$ART_Stop_Day))
  
  long_outcome$ART_Stop_Day[which(long_outcome$ART2 == 1)] <- long_outcome$ART_Stop_Day[which(long_outcome$ART2 == 1)]+long_outcome$Visit2[which(long_outcome$ART2 == 1)] 
  long_outcome$ART_Stop_Day[which(long_outcome$ART2 == 0 & long_outcome$ART1 == 1)] <-30 #Made 30 days of ART
  long_outcome$ART_Stop_Day[which(long_outcome$ART2 == 0  & long_outcome$ART1 == 0)] <-0
  long_outcome$ART_Stop_Day[which(long_outcome$ART_Stop_Day > 730)] <- 730
  
  #If ART discontinues before TPT, make TPT discontinuation day = ART discontinuation day
  
  long_outcome$TPT_Stop_Day[which(long_outcome$TPT_Stop_Day > long_outcome$ART_Stop_Day)] <- long_outcome$ART_Stop_Day[which(long_outcome$TPT_Stop_Day > long_outcome$ART_Stop_Day)]
  
  return(long_outcome)
  
  
}

long_outcome_A = mclapply(outcome_A,create_long_outcome, mc.cores=ncores)
long_outcome_B = mclapply(outcome_B,create_long_outcome, mc.cores=ncores)
long_outcome_C = mclapply(outcome_C,create_long_outcome, mc.cores=ncores)
long_outcome_D = mclapply(outcome_D,create_long_outcome, mc.cores=ncores)




#Visit 1 to 2 

#Create incidence matrix between visit 1 and 2
#Read in incidence rates
#Sensitvity Analysis, change Daily_Incidence_Rates csv to "Incidence_Rates_labeled_sensitivity.csv" and long_outcome_cohort$AL = 60  + long_outcome_cohort$Visit2 - long_outcome_cohort$ART1*pmin(long_outcome_cohort$Visit2,30) 
Daily_Incidence_Rates = read.csv(paste0(sourcelocation, "Incidence_Rates_labeled_sensitivity.csv"), header = TRUE, fill = TRUE, row.names = 1)
colnames(Daily_Incidence_Rates) = c("100","100_200","200_350","350_500")
Incidence_Multiplier = rtriangle(nsims,a = .8,b =1.2,c = 1)
TPT_Multiplier = rtriangle(nsims, a = .8, b =1.2,c=1)
Daily_Incidence_Rates_Unc = vector("list", length = nsims)

for(a in 1:nsims){
  Rates = Daily_Incidence_Rates
  Rates[c("AI","AC","LI","LC","I","C"),]= Rates[c("AI","AC","LI","LC","I","C"),]*TPT_Multiplier[a]
  Daily_Incidence_Rates_Unc[[a]] = Rates*Incidence_Multiplier[a] 
  
}

#Create end dates between visit 1 and 2 


visit12_end_date <- function(long_outcome_cohort){
  visit_12_end_date = matrix(0, nrow = dim(long_outcome_cohort)[1], ncol = 4)
  colnames(visit_12_end_date) =c("A1I1","A1","I1","null1")
  visit_12_end_date = as.data.frame(visit_12_end_date)
  
  for(a in 1:dim(long_outcome_cohort)[1]){
    
    if(!is.na(long_outcome_cohort$Visit2[a])){
      if(long_outcome_cohort$ART1[a] == 1 && long_outcome_cohort$TPT1[a] == 1){
        visit_12_end_date$A1I1[a] = min(30, long_outcome_cohort$TPT_Stop_Day[a],long_outcome_cohort$Visit2[a])
      }
      else if(long_outcome_cohort$ART1[a] == 1 && long_outcome_cohort$TPT1[a] == 0){
        visit_12_end_date$A1[a] = min(30, long_outcome_cohort$ART_Stop_Day[a],long_outcome_cohort$Visit2[a]) 
      }
      else if(long_outcome_cohort$ART1[a] == 0 && long_outcome_cohort$TPT1[a] == 1){
        visit_12_end_date$I1[a] = min(30, long_outcome_cohort$TPT_Stop_Day[a],long_outcome_cohort$Visit2[a]) 
      }
      visit_12_end_date$A1[a] = visit_12_end_date$A1I1[a] + visit_12_end_date$A1[a]
      visit_12_end_date$I1[a] = visit_12_end_date$I1[a] + visit_12_end_date$A1[a] 
      visit_12_end_date$null1[a] = long_outcome_cohort$Visit2[a]
    } else{
      if(long_outcome_cohort$ART1[a] == 1 && long_outcome_cohort$TPT1[a] == 1){
        visit_12_end_date$A1I1[a] = min(30, long_outcome_cohort$TPT_Stop_Day[a])
      }
      else if(long_outcome_cohort$ART1[a] == 1 && long_outcome_cohort$TPT1[a] == 0){
        visit_12_end_date$A1[a] = min(30, long_outcome_cohort$ART_Stop_Day[a]) 
      }
      else if(long_outcome_cohort$ART1[a] == 0 && long_outcome_cohort$TPT1[a] == 1){
        visit_12_end_date$I1[a] = min(30, long_outcome_cohort$TPT_Stop_Day[a]) 
      }
      visit_12_end_date$A1[a] = visit_12_end_date$A1I1[a] + visit_12_end_date$A1[a]
      visit_12_end_date$I1[a] = visit_12_end_date$I1[a] + visit_12_end_date$A1[a] 
      visit_12_end_date$null1[a] = 730
      
    }
    
  }
  
  
  return(visit_12_end_date)
  
  
}


V12_end_dates_A = mclapply(long_outcome_A, function(x) visit12_end_date(x), mc.cores=ncores)
V12_end_dates_B = mclapply(long_outcome_B, function(x) visit12_end_date(x), mc.cores=ncores)
V12_end_dates_C = mclapply(long_outcome_C, function(x) visit12_end_date(x), mc.cores=ncores)
V12_end_dates_D = mclapply(long_outcome_D, function(x) visit12_end_date(x), mc.cores=ncores)

# # for checking: 
# cbind(long_outcome_A[[1]], V12_end_dates_A[[1]])


#** This version calculates each patient's risk of incident TB in the first year, and total over two years.
#** We no longer determine the day of incident TB, and we'll have to make some assumptions to assign treatment costs that span the two years.
#*
#old 12 incident_tb_deterministic function
#visit12_incident_tb_deterministic <- function(outcome_cohort,visit_end_date, Daily_Incidence_Rates) {
  
  outcome_cohort$expected_incidence_y1_12 <- 0
  
  CD4columns <- match(outcome_cohort$CD4,colnames(Daily_Incidence_Rates))
  rate_a1i1 = unlist(Daily_Incidence_Rates["AI", CD4columns])
  rate_a1 = unlist(Daily_Incidence_Rates["A", CD4columns])
  rate_i1 = unlist(Daily_Incidence_Rates["I", CD4columns])
  rate_01 = unlist(Daily_Incidence_Rates["null", CD4columns])
  
  outcome_cohort$expected_incidence_y1_12 = outcome_cohort$expected_incidence_y1_12 + 
                (rate_a1i1)*pmin(365, visit_end_date$A1I1) + 
                (rate_a1)*(pmin(365, visit_end_date$A1) - pmin(365, visit_end_date$A1I1)) +
                (rate_i1)*(pmin(365, visit_end_date$I1) - pmin(365, visit_end_date$A1)) +
                (rate_01)*(pmin(365, visit_end_date$null1) - pmin(365, visit_end_date$I1))
  
  outcome_cohort$expected_incidence_total_12 = 
    (rate_a1i1)*(visit_end_date$A1I1) + 
    (rate_a1)*(visit_end_date$A1 - visit_end_date$A1I1) +
    (rate_i1)*(visit_end_date$I1 - visit_end_date$A1) +
    (rate_01)*(visit_end_date$null1 - visit_end_date$I1)
  
  return(outcome_cohort)
}

visit12_incident_tb_deterministic <- function(outcome_cohort,visit_end_date, Daily_Incidence_Rates) {
  
  outcome_cohort$expected_incidence_y1_12_art <- 0
  outcome_cohort$expected_incidence_y1_12_noart <- 0
  
  CD4columns <- match(outcome_cohort$CD4,colnames(Daily_Incidence_Rates))
  rate_a1i1 = unlist(Daily_Incidence_Rates["AI", CD4columns])
  rate_a1 = unlist(Daily_Incidence_Rates["A", CD4columns])
  rate_i1 = unlist(Daily_Incidence_Rates["I", CD4columns])
  rate_01 = unlist(Daily_Incidence_Rates["null", CD4columns])
  
  outcome_cohort$expected_incidence_y1_12_art = outcome_cohort$expected_incidence_y1_12_art + 
    (rate_a1i1)*pmin(365, visit_end_date$A1I1) + 
    (rate_a1)*(pmin(365, visit_end_date$A1) - pmin(365, visit_end_date$A1I1))
  outcome_cohort$expected_incidence_y1_12_noart = outcome_cohort$expected_incidence_y1_12_noart + 
    (rate_i1)*(pmin(365, visit_end_date$I1) - pmin(365, visit_end_date$A1)) +
    (rate_01)*(pmin(365, visit_end_date$null1) - pmin(365, visit_end_date$I1))
  
  outcome_cohort$expected_incidence_total_12_art = 
    (rate_a1i1)*(visit_end_date$A1I1) + 
    (rate_a1)*(visit_end_date$A1 - visit_end_date$A1I1)
  outcome_cohort$expected_incidence_total_12_noart = 
    (rate_i1)*(visit_end_date$I1 - visit_end_date$A1) +
    (rate_01)*(visit_end_date$null1 - visit_end_date$I1)
  
  outcome_cohort$expected_incidence_y1_12 <- 
    outcome_cohort$expected_incidence_y1_12_art + 
    outcome_cohort$expected_incidence_y1_12_noart
  outcome_cohort$expected_incidence_total_12 <- 
    outcome_cohort$expected_incidence_total_12_art + 
    outcome_cohort$expected_incidence_total_12_noart
  
  return(outcome_cohort)
}

#Update incidence
for(a in 1:nsims){
  long_outcome_A[[a]] = visit12_incident_tb_deterministic(long_outcome_A[[a]],V12_end_dates_A[[a]],Daily_Incidence_Rates_Unc[[a]])
  long_outcome_B[[a]] = visit12_incident_tb_deterministic(long_outcome_B[[a]],V12_end_dates_B[[a]],Daily_Incidence_Rates_Unc[[a]])
  long_outcome_C[[a]] = visit12_incident_tb_deterministic(long_outcome_C[[a]],V12_end_dates_C[[a]],Daily_Incidence_Rates_Unc[[a]])
  long_outcome_D[[a]] = visit12_incident_tb_deterministic(long_outcome_D[[a]],V12_end_dates_D[[a]],Daily_Incidence_Rates_Unc[[a]])
}


# check that cumulative incidence is reasonable (this is only between visits 1 and 2)
#long_outcome_A[[1]] %>% group_by(CD4) %>% summarise(median(expected_incidence_total_12))
#long_outcome_B[[1]] %>% group_by(CD4) %>% summarise(median(expected_incidence_total_12))
#long_outcome_C[[1]] %>% group_by(CD4) %>% summarise(median(expected_incidence_total_12))
#long_outcome_D[[1]] %>% group_by(CD4) %>% summarise(median(expected_incidence_total_12))
#long_outcome_A[[1]] %>% group_by(CD4) %>% summarise(mean(expected_incidence_total_12>0.05), mean(expected_incidence_total_12>0.05 & is.na(Visit2))) # almost of the high incidence is among those with no visit 2 before day 730, as it shouldbe

#** I fixed visit2 -> Visit2 and min -> pmin (since you're working with vectors) in this initial section
#** And then I vectorized the for loops over 'a', to try to speed this up and make it cleaner for debugging.

#sensitivity analysis: 
#long_outcome_cohort$AL = 60  + long_outcome_cohort$Visit2 - long_outcome_cohort$ART1*pmin(long_outcome_cohort$Visit2,30) 
#Original analysis: 
#180  + long_outcome_cohort$Visit2 - long_outcome_cohort$ART1*pmin(long_outcome_cohort$Visit2,30)

Vlong_end_date <-function(long_outcome_cohort){
 long_outcome_cohort$AL = 60  + long_outcome_cohort$Visit2 - long_outcome_cohort$ART1*pmin(long_outcome_cohort$Visit2,30)   #long term ART
 long_outcome_cohort$IE = 360 + long_outcome_cohort$Visit2 - long_outcome_cohort$TPT1*pmin(long_outcome_cohort$Visit2,30) #threshold for completed IPT
 long_outcome_cohort$IP = 180  + long_outcome_cohort$Visit2 - long_outcome_cohort$TPT1*pmin(long_outcome_cohort$Visit2,30)#point at which you are still completed 
 long_outcome_cohort$Tmax = 730 # total period of observation? 
 long_outcome_cohort$AS =  long_outcome_cohort$ART_Stop_Day #Stop day ART? -should be counted from day 0
 long_outcome_cohort$IS =  long_outcome_cohort$TPT_Stop_Day #Stop day TPT? -should be counted from day 0
 
 Vlong_end_date <- long_outcome_cohort %>% mutate(AI = case_when(is.na(Visit2) ~ 730,
                                                                 ART2==1 & TPT2==1 ~ pmin(AL, AS, IE, IS) + Visit2,
                                                                 TRUE ~ Visit2),
                                                  AC = case_when(is.na(Visit2) ~ 730,
                                                                 ART2==1 & TPT2==1 & IP < pmin(IS,IE) & min(IS,IE)<AL~ pmin(IS, IE, AL) + AI,
                                                                 TRUE ~ AI),
                                                  LI = case_when(is.na(Visit2) ~ 730,
                                                                 ART2==1 & TPT2==1 & AL < pmin(IS,IE) & AL<AS ~ pmin(IS, IE, AS) + AC,
                                                                 TRUE ~ AC),
                                                  LC = case_when(is.na(Visit2) ~ 730,
                                                                 ART2==1 & TPT2==1 & AL<AS & pmin(IS,IE)>AL & IP < pmin(IS,IE) ~ AS + LI,
                                                                 TRUE ~ LI),
                                                  A = case_when(is.na(Visit2) ~ 730,
                                                                (ART2==1 & TPT2==0) | (ART2==1 & TPT2==1 & IS<IP) ~ pmin(AS,AL) + LC,
                                                                TRUE ~ LC),
                                                  L = case_when(is.na(Visit2) ~ 730,
                                                                (ART2==1 & TPT2==0 & AL<AS) | 
                                                                 (ART2==1 & TPT2==1 & IS<IP & AL<AS) ~ AS + A,
                                                                TRUE ~ A),
                                                  I = case_when(is.na(Visit2) ~ 730,
                                                                (ART2==0 & TPT2==1) | 
                                                                 (ART2==1 & TPT2==1 & AS<IS) ~ IS + L,
                                                                TRUE ~ L),
                                                  C = case_when(is.na(Visit2) ~ 730,
                                                                (TPT2==1 & IP<pmin(IS,IE)) ~ Tmax,
                                                                TRUE ~ I),
                                                  nullL = case_when(is.na(Visit2) ~ 730,
                                                                    TRUE ~ Tmax)) %>%
  dplyr::select(AI, AC, LI, LC, A, L, I, C, nullL)
 
 Vlong_end_date[Vlong_end_date >730]<- 730
 
 
 return(Vlong_end_date)
 
}

Vlong_end_dates_A = mclapply(long_outcome_A, function(x) Vlong_end_date(x), mc.cores=ncores)
Vlong_end_dates_B = mclapply(long_outcome_B, function(x) Vlong_end_date(x), mc.cores=ncores)
Vlong_end_dates_C = mclapply(long_outcome_C, function(x) Vlong_end_date(x), mc.cores=ncores)
Vlong_end_dates_D = mclapply(long_outcome_D, function(x) Vlong_end_date(x), mc.cores=ncores)


#old Vlong incident tb deterministic function
#Vlong_incident_tb_deterministic <- function(outcome_cohort,Vlong_end_date,Daily_Incidence_Rates) {
  
  CD4columns <- match(outcome_cohort$CD4,colnames(Daily_Incidence_Rates))
  
  rate_ai = unlist(Daily_Incidence_Rates["AI",CD4columns])
  rate_ac = unlist(Daily_Incidence_Rates["AC",CD4columns])
  rate_li = unlist(Daily_Incidence_Rates["LI",CD4columns])
  rate_lc = unlist(Daily_Incidence_Rates["LC",CD4columns])
  rate_a = unlist(Daily_Incidence_Rates["A",CD4columns])
  rate_l = unlist(Daily_Incidence_Rates["L",CD4columns])
  rate_i = unlist(Daily_Incidence_Rates["I",CD4columns])
  rate_c = unlist(Daily_Incidence_Rates["C",CD4columns])
  rate_nullL = unlist(Daily_Incidence_Rates["null",CD4columns])
  
  outcome_cohort$expected_incidence_y1 <- outcome_cohort$expected_incidence_y1_12
  outcome_cohort$expected_incidence_total <- outcome_cohort$expected_incidence_total_12
  
    outcome_cohort$expected_incidence_y1[outcome_cohort$Visit2 < 365 & !is.na(outcome_cohort$Visit2)] = 
      outcome_cohort$expected_incidence_y1[outcome_cohort$Visit2 < 365 & !is.na(outcome_cohort$Visit2)] + 
      ((rate_ai)*(pmin(365,Vlong_end_date$AI) - outcome_cohort$Visit2) + 
      (rate_ac)*(pmin(365, Vlong_end_date$AC) - pmin(365, Vlong_end_date$AI)) +
      (rate_li)*(pmin(365, Vlong_end_date$LI) - pmin(365, Vlong_end_date$AC)) +
      (rate_lc)*(pmin(365, Vlong_end_date$LC) - pmin(365, Vlong_end_date$LI)) +
      (rate_a)*(pmin(365, Vlong_end_date$A) - pmin(365, Vlong_end_date$LC)) +
      (rate_l)*(pmin(365, Vlong_end_date$L) - pmin(365, Vlong_end_date$A)) +
      (rate_i)*(pmin(365, Vlong_end_date$I) - pmin(365, Vlong_end_date$L)) +
      (rate_c)*(pmin(365, Vlong_end_date$C) - pmin(365, Vlong_end_date$I)) +
      (rate_nullL)*(pmin(365, Vlong_end_date$nullL) - pmin(365, Vlong_end_date$C)))[outcome_cohort$Visit2 < 365 & !is.na(outcome_cohort$Visit2)]
  
  outcome_cohort$expected_incidence_total[!is.na(outcome_cohort$Visit2)] = 
    outcome_cohort$expected_incidence_total[!is.na(outcome_cohort$Visit2)] + 
    (
      (rate_ai)*(Vlong_end_date$AI - outcome_cohort$Visit2) + 
      (rate_ac)*(Vlong_end_date$AC - Vlong_end_date$AI) +
      (rate_li)*(Vlong_end_date$LI - Vlong_end_date$AC) +
      (rate_lc)*(Vlong_end_date$LC - Vlong_end_date$LI) +
      (rate_a)*(Vlong_end_date$A - Vlong_end_date$LC) +
      (rate_l)*(Vlong_end_date$L - Vlong_end_date$A) +
      (rate_i)*(Vlong_end_date$I - Vlong_end_date$L) +
      (rate_c)*(Vlong_end_date$C - Vlong_end_date$I) +
      (rate_nullL)*(Vlong_end_date$nullL - Vlong_end_date$C))[!is.na(outcome_cohort$Visit2)]
    
  
  return(outcome_cohort)
  
}

Vlong_incident_tb_deterministic <- function(outcome_cohort,Vlong_end_date,Daily_Incidence_Rates) {
  
  CD4columns <- match(outcome_cohort$CD4,colnames(Daily_Incidence_Rates))
  
  rate_ai = unlist(Daily_Incidence_Rates["AI",CD4columns])
  rate_ac = unlist(Daily_Incidence_Rates["AC",CD4columns])
  rate_li = unlist(Daily_Incidence_Rates["LI",CD4columns])
  rate_lc = unlist(Daily_Incidence_Rates["LC",CD4columns])
  rate_a = unlist(Daily_Incidence_Rates["A",CD4columns])
  rate_l = unlist(Daily_Incidence_Rates["L",CD4columns])
  rate_i = unlist(Daily_Incidence_Rates["I",CD4columns])
  rate_c = unlist(Daily_Incidence_Rates["C",CD4columns])
  rate_nullL = unlist(Daily_Incidence_Rates["null",CD4columns])
  
  outcome_cohort$expected_incidence_y1_art <- outcome_cohort$expected_incidence_y1_12_art
  outcome_cohort$expected_incidence_y1_noart <- outcome_cohort$expected_incidence_y1_12_noart
  outcome_cohort$expected_incidence_total_art <- outcome_cohort$expected_incidence_total_12_art
  outcome_cohort$expected_incidence_total_noart <- outcome_cohort$expected_incidence_total_12_noart
  
  outcome_cohort$expected_incidence_y1_art[outcome_cohort$Visit2 < 365 & !is.na(outcome_cohort$Visit2)] = 
    outcome_cohort$expected_incidence_y1_art[outcome_cohort$Visit2 < 365 & !is.na(outcome_cohort$Visit2)] + 
    ((rate_ai)*(pmin(365,Vlong_end_date$AI) - outcome_cohort$Visit2) + 
       (rate_ac)*(pmin(365, Vlong_end_date$AC) - pmin(365, Vlong_end_date$AI)) +
       (rate_li)*(pmin(365, Vlong_end_date$LI) - pmin(365, Vlong_end_date$AC)) +
       (rate_lc)*(pmin(365, Vlong_end_date$LC) - pmin(365, Vlong_end_date$LI)) +
       (rate_a)*(pmin(365, Vlong_end_date$A) - pmin(365, Vlong_end_date$LC)) +
       (rate_l)*(pmin(365, Vlong_end_date$L) - pmin(365, Vlong_end_date$A)))[outcome_cohort$Visit2 < 365 & !is.na(outcome_cohort$Visit2)]
  outcome_cohort$expected_incidence_y1_noart[outcome_cohort$Visit2 < 365 & !is.na(outcome_cohort$Visit2)] = 
    outcome_cohort$expected_incidence_y1_noart[outcome_cohort$Visit2 < 365 & !is.na(outcome_cohort$Visit2)] + 
    ((rate_i)*(pmin(365, Vlong_end_date$I) - pmin(365, Vlong_end_date$L)) +
       (rate_c)*(pmin(365, Vlong_end_date$C) - pmin(365, Vlong_end_date$I)) +
       (rate_nullL)*(pmin(365, Vlong_end_date$nullL) - pmin(365, Vlong_end_date$C)))[outcome_cohort$Visit2 < 365 & !is.na(outcome_cohort$Visit2)]
  
  outcome_cohort$expected_incidence_total_art[!is.na(outcome_cohort$Visit2)] = 
    outcome_cohort$expected_incidence_total_art[!is.na(outcome_cohort$Visit2)] + 
    ( (rate_ai)*(Vlong_end_date$AI - outcome_cohort$Visit2) + 
        (rate_ac)*(Vlong_end_date$AC - Vlong_end_date$AI) +
        (rate_li)*(Vlong_end_date$LI - Vlong_end_date$AC) +
        (rate_lc)*(Vlong_end_date$LC - Vlong_end_date$LI) +
        (rate_a)*(Vlong_end_date$A - Vlong_end_date$LC) +
        (rate_l)*(Vlong_end_date$L - Vlong_end_date$A))[!is.na(outcome_cohort$Visit2)]
  outcome_cohort$expected_incidence_total_noart[!is.na(outcome_cohort$Visit2)] = 
    outcome_cohort$expected_incidence_total_noart[!is.na(outcome_cohort$Visit2)] + 
    ( (rate_i)*(Vlong_end_date$I - Vlong_end_date$L) +
        (rate_c)*(Vlong_end_date$C - Vlong_end_date$I) +
        (rate_nullL)*(Vlong_end_date$nullL - Vlong_end_date$C))[!is.na(outcome_cohort$Visit2)]
  
  outcome_cohort$expected_incidence_y1 <- outcome_cohort$expected_incidence_y1_art + outcome_cohort$expected_incidence_y1_noart
  outcome_cohort$expected_incidence_total <- outcome_cohort$expected_incidence_total_art + outcome_cohort$expected_incidence_total_noart
  
  
  return(outcome_cohort)
  
}

#Update incidence 
long_outcome_A_TB = vector("list",length = nsims)
long_outcome_B_TB = vector("list",length = nsims)
long_outcome_C_TB = vector("list",length = nsims)
long_outcome_D_TB = vector("list",length = nsims)

for(a in 1:nsims){
  long_outcome_A_TB[[a]] = Vlong_incident_tb_deterministic(long_outcome_A[[a]],Vlong_end_dates_A[[a]],Daily_Incidence_Rates_Unc[[a]])
  long_outcome_B_TB[[a]] = Vlong_incident_tb_deterministic(long_outcome_B[[a]],Vlong_end_dates_B[[a]],Daily_Incidence_Rates_Unc[[a]])
  long_outcome_C_TB[[a]] = Vlong_incident_tb_deterministic(long_outcome_C[[a]],Vlong_end_dates_C[[a]],Daily_Incidence_Rates_Unc[[a]])
  long_outcome_D_TB[[a]] = Vlong_incident_tb_deterministic(long_outcome_D[[a]],Vlong_end_dates_D[[a]],Daily_Incidence_Rates_Unc[[a]])
}



# As above check that incidence rates are reasonable
long_outcome_A_TB[[1]] %>% group_by(CD4) %>% summarise(mean(expected_incidence_total), mean(expected_incidence_y1))
long_outcome_B_TB[[2]] %>% group_by(CD4) %>% summarise(mean(expected_incidence_total), mean(expected_incidence_y1))
long_outcome_C_TB[[3]] %>% group_by(CD4) %>% summarise(mean(expected_incidence_total), mean(expected_incidence_y1))
long_outcome_D_TB[[4]] %>% group_by(CD4) %>% summarise(mean(expected_incidence_total), mean(expected_incidence_y1))

# calculate total incidence in each year:
#y1, absolute incidence with individual algorithms
summary(unlist(lapply(long_outcome_A_TB, function(x) sum(x$expected_incidence_y1)/dim(x)[1])))
summary(unlist(lapply(long_outcome_B_TB, function(x) sum(x$expected_incidence_y1)/dim(x)[1])))
summary(unlist(lapply(long_outcome_C_TB, function(x) sum(x$expected_incidence_y1)/dim(x)[1])))
summary(unlist(lapply(long_outcome_D_TB, function(x) sum(x$expected_incidence_y1)/dim(x)[1])))

par(mfrow=c(1,4))
hist(unlist(lapply(long_outcome_A_TB, function(x) sum(x$expected_incidence_y1))))
hist(unlist(lapply(long_outcome_B_TB, function(x) sum(x$expected_incidence_y1))))
hist(unlist(lapply(long_outcome_C_TB, function(x) sum(x$expected_incidence_y1))))
hist(unlist(lapply(long_outcome_D_TB, function(x) sum(x$expected_incidence_y1))))


#total absolute incidence with individual algorithms
summary(unlist(lapply(long_outcome_A_TB, function(x) sum(x$expected_incidence_total)))/unlist(lapply(long_outcome_A_TB, function(X) dim(x)[1])))
summary(unlist(lapply(long_outcome_B_TB, function(x) sum(x$expected_incidence_total)))/unlist(lapply(long_outcome_B_TB, function(X) dim(x)[1])))
summary(unlist(lapply(long_outcome_C_TB, function(x) sum(x$expected_incidence_total)))/unlist(lapply(long_outcome_C_TB, function(X) dim(x)[1])))
summary(unlist(lapply(long_outcome_D_TB, function(x) sum(x$expected_incidence_total)))/unlist(lapply(long_outcome_D_TB, function(X) dim(x)[1])))

par(mfrow=c(1,4))
hist(unlist(lapply(long_outcome_A_TB, function(x) sum(x$expected_incidence_total))), main = "Algorithm A")
hist(unlist(lapply(long_outcome_B_TB, function(x) sum(x$expected_incidence_total))), main = "Algorithm B")
hist(unlist(lapply(long_outcome_C_TB, function(x) sum(x$expected_incidence_total))), main = "Algorithm C")
hist(unlist(lapply(long_outcome_D_TB, function(x) sum(x$expected_incidence_total))), main = "Algorithm D")

#comparing  algorithms 
summary(unlist(mapply(FUN = function(a,b) sum(b$expected_incidence_total) - sum(a$expected_incidence_total), a=long_outcome_A_TB, b=long_outcome_D_TB)))
summary(unlist(mapply(FUN = function(a,b) sum(b$expected_incidence_total) - sum(a$expected_incidence_total), a=long_outcome_B_TB, b=long_outcome_D_TB)))
summary(unlist(mapply(FUN = function(a,b) sum(b$expected_incidence_total) - sum(a$expected_incidence_total), a=long_outcome_C_TB, b=long_outcome_D_TB)))
summary(unlist(mapply(FUN = function(a,b) sum(b$expected_incidence_total) - sum(a$expected_incidence_total), a=long_outcome_A_TB, b=long_outcome_C_TB)))

par(mfrow=c(1,4))
hist(unlist(mapply(FUN = function(a,b) sum(b$expected_incidence_total) - sum(a$expected_incidence_total), a=long_outcome_A_TB, b=long_outcome_D_TB)), main = "Algorithm D-A")
hist(unlist(mapply(FUN = function(a,b) sum(b$expected_incidence_total) - sum(a$expected_incidence_total), a=long_outcome_B_TB, b=long_outcome_D_TB)), main = "Algorithm D-B")
hist(unlist(mapply(FUN = function(a,b) sum(b$expected_incidence_total) - sum(a$expected_incidence_total), a=long_outcome_C_TB, b=long_outcome_D_TB)), main = "Algorithm D-C")
hist(unlist(mapply(FUN = function(a,b) sum(b$expected_incidence_total) - sum(a$expected_incidence_total), a=long_outcome_A_TB, b=long_outcome_C_TB)), main = "Algorithm C-A")



# y1, comparing D vs A:
summary(unlist(mapply(FUN = function(a,b) sum(b$expected_incidence_y1) - sum(a$expected_incidence_y1), a=long_outcome_A_TB, b=long_outcome_D_TB)))

#y2
summary(unlist(lapply(long_outcome_A_TB, function(x) sum(x$expected_incidence_total - x$expected_incidence_y1))))

# Difference in TPT beyond visit 2:
summary(unlist(lapply(long_outcome_A_TB, function(x) sum(x$TPT2))))
summary(unlist(lapply(long_outcome_B_TB, function(x) sum(x$TPT2))))
summary(unlist(lapply(long_outcome_C_TB, function(x) sum(x$TPT2))))
summary(unlist(lapply(long_outcome_D_TB, function(x) sum(x$TPT2))))

# Advantages of B and D over A and C may have more to do with ART coverage than with TPT:
summary(unlist(lapply(long_outcome_A_TB, function(x) sum(x$ART2))))
summary(unlist(lapply(long_outcome_B_TB, function(x) sum(x$ART2)))) 
summary(unlist(lapply(long_outcome_C_TB, function(x) sum(x$ART2)))) 
summary(unlist(lapply(long_outcome_D_TB, function(x) sum(x$ART2))))

# So ~150 more people getting TPT (1 in 6 of those eligible), and their incidence is reduced by about 1/2, means 

#Incidence while on/off ART 
summary(unlist(lapply(long_outcome_A_TB, function(x) sum(x$expected_incidence_total_noart)))/unlist(lapply(long_outcome_A_TB, function(x) sum(x$expected_incidence_total))))
summary(unlist(lapply(long_outcome_B_TB, function(x) sum(x$expected_incidence_total_noart)))/unlist(lapply(long_outcome_B_TB, function(x) sum(x$expected_incidence_total))))
summary(unlist(lapply(long_outcome_C_TB, function(x) sum(x$expected_incidence_total_noart)))/unlist(lapply(long_outcome_C_TB, function(x) sum(x$expected_incidence_total))))
summary(unlist(lapply(long_outcome_D_TB, function(x) sum(x$expected_incidence_total_noart)))/unlist(lapply(long_outcome_D_TB, function(x) sum(x$expected_incidence_total))))



















