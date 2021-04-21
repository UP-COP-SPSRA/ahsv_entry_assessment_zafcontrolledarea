# Contents ####
# A: Setup
  # Libraries and functions
  # Setup output directory
  # Model parameters
  # Retrieval of data used in analysis
  # Dataset processing
      # 1. Complete AHS cases dataset
      # 2. Create series of sub-clinical rates
      # 3. Establish underlying p_inf
      # 4. Movement dataset
          # 4.1 Summary stats for Table 1

# B: Analysis code
  # 5. probability of entry
    # 5.1 Standard movements
    # 5.2 Standard SOQ movements
    # 5.3 VPQO movements
    # 5.4 VPSOQ movements
    # 5.5 Zebra movements
    # 5.6 What-if analysis - p_ent with no control
  # 6. Sensitivity analysis
    #6.1 Standard movements
    #6.2 Standard SOQ movements

# C: Model output code 
  # 7. Movement map per municipality
  # 8. P_inf with map output
  # 9. P_entry with table and graphics outputs
  # 10. Sensitivity analysis outputs - tornado plots
  # 11. What-if scenario with tablular outputs

# A: Setup ####
rm(list = ls())

#Libraries and functions####
library(mc2d)
library(readr)
library (plyr)
library(dplyr)
library(tidyr)
library(fitdistrplus)
library(httr)
require(ggplot2)
require(sf)

`%notin%` = function(x,y) !(x %in% y)

round_df <- function(x, digits) {
  numeric_columns <- sapply(x, mode) == 'numeric'
  x[numeric_columns] <-  round(x[numeric_columns], digits)
  x
}

# Setup output directory ####
subDir <- "modeloutputs"
subDirTables<-"tables"
subDirGIS<-"gis"
subDirfigures<-"figures"

output_dir <- file.path(getwd(), subDir)
output_dirtables <- file.path(getwd(), subDir, subDirTables)
output_dirgis <- file.path(getwd(), subDir, subDirGIS)
output_dirfigures <- file.path(getwd(), subDir, subDirfigures)

if (!dir.exists(output_dir)){
  dir.create(output_dir)
} else {
  print("Dir already exists!")
}

if (!dir.exists(output_dirtables)){
  dir.create(output_dirtables)
} else {
  print("Tables Dir already exists!")
}

if (!dir.exists(output_dirgis)){
  dir.create(output_dirgis)
} else {
  print("GIS Dir already exists!")
}


if (!dir.exists(output_dirfigures)){
  dir.create(output_dirfigures)
} else {
  print("Figures Dir already exists!")
}

# Model parameters ####
iterations<-10000
evaluationperiods<-1:12 # number of months from 1 Jan 2019 to 31 Dec 2019

# Infectious period (IP) distributions ####
# Horses (H)
# based on published parameter estimate: Backer et al. https://doi.org/10.1371/journal.pone.0023066
hip.alpha<-29.75
hip.beta<-4.95

# Zebra (Z)
# Based on experimental data: Barnard BJ et al. PMID: 7501371
dpi<-c(20,13,17,24,24,24,40) #days post infection where virus isolation still occurred
fit=fitdist(dpi, "gamma", method="mle")# Computes MLE parameters for gamma distribution using dpi data
zip.alpha<-as.numeric(fit$estimate[1])
zip.beta<-as.numeric(fit$estimate[2])

#Retrieval of data used in analysis####
#Online from GitHUB
# Import and unzip repository into working directory
url <- "https://github.com/UP-COP-SPSRA/ahsv_entry_assessment_zafcontrolledarea/archive/refs/heads/main.zip"
GET(url, write_disk("manuscriptdata.zip", overwrite = TRUE))
#unzip file directly into working directory - ensure that the outcome should be a folder with associated data called 
# ahsv_entry_assessment_zafcontrolledarea-main in the working directory

dfoutbreakdata <-read.csv("./ahsv_entry_assessment_zafcontrolledarea-main/datafiles/data_outbreaks_ca.csv")
dfmovement <-read.csv("./ahsv_entry_assessment_zafcontrolledarea-main/datafiles/data_2019movements.csv")
dfcases <-read.csv("./ahsv_entry_assessment_zafcontrolledarea-main/datafiles/data_2019cases_ecod.csv")
dfcensus <-read.csv("./ahsv_entry_assessment_zafcontrolledarea-main/datafiles/data_census_rsa.csv")

# Data set processing ####
# 1: Complete AHS cases dataset ####
dfcases.analyse <- dfcases %>% dplyr::select(c(lmgid, casemonth, ahs_cases))

# Ensure that all movement areas of origin are included in the areas that cases occurred in with replacement of zero when not
# Initial allocation of zero to arbitrary month 1 when the area does not contain a movement
# This will be extrapolated to fill the dataset for the remaining months in future step

dfcases.analyse<-rbind(dfcases.analyse,
               tibble(
                      lmgid = unique(dfmovement$lmgid)[which(unique(dfmovement$lmgid) %in% unique(dfcases.analyse$lmgid)==FALSE)],
                      casemonth = 1,
                      ahs_cases= 0))

# Ensure case months that do not occur within the cases dataset are allocated with zero cases
# allocate to an arbitrary minimum area (local municipality) GID that is in the cases dataset
# This will be extrapolated to fill the dataset for the remaining areas in future step

dfcases.analyse<-rbind(dfcases.analyse,
                       tibble(lmgid = min(unique(dfcases.analyse$lmgid)),
                              casemonth = (evaluationperiods[which(evaluationperiods %in% unique(dfcases.analyse$casemonth)==FALSE)]),
                              ahs_cases= 0))

# Fill the data set for any case:area permutations that do not exist
dfcases.analyse<- dfcases.analyse %>% tidyr::complete(casemonth, 
                              nesting(lmgid), 
                              fill = list(ahs_cases = 0))

#error checking 
nrow(dfcases.analyse) == length(unique(dfcases.analyse$lmgid))*length(evaluationperiods) # Must be TRUE

#create vector of local municipalities where cases occurred from to make est_pinf more realistic
# i.e. horses moving through SOQ and VPQO likely to come from areas where cases were reported
dfcases.lmorigins<-unique(dfcases$lmgid)

#2 - create a series of sub-clinical rates ####
#Base the sub-clinical rate on the outbreaks where this has been evident (and measured) in the controlled area i.e. in 2011, 2014 (x2) and 2016
#See Grewar et al. 2020 DOI: 10.1111/tbed.13566 where a similar approach was taken

dfoutbreakdata.subclin<-dfoutbreakdata %>% dplyr::filter(year %in% c("2011","2014","2016"))

o.2011<-rbeta(iterations,
              dfoutbreakdata.subclin[1,]$cases_subclinical+1,
              dfoutbreakdata.subclin[1,]$cases_all - dfoutbreakdata.subclin[1,]$cases_subclinical+1)
o.2014p<-rbeta(iterations,
              dfoutbreakdata.subclin[2,]$cases_subclinical+1,
              dfoutbreakdata.subclin[2,]$cases_all - dfoutbreakdata.subclin[2,]$cases_subclinical+1)
o.2014r<-rbeta(iterations,
              dfoutbreakdata.subclin[3,]$cases_subclinical+1,
              dfoutbreakdata.subclin[3,]$cases_all - dfoutbreakdata.subclin[3,]$cases_subclinical+1)
o.2016<-rbeta(iterations,
              dfoutbreakdata.subclin[4,]$cases_subclinical+1,
              dfoutbreakdata.subclin[4,]$cases_all - dfoutbreakdata.subclin[4,]$cases_subclinical+1)

o.dat<-cbind(o.2011,o.2014p,o.2014r,o.2016)

sample.subclin.rate=rep(NA, iterations)
n.outbreaks=ncol(o.dat)
for(i in 1:iterations){
  which.outbreak=sample(1:n.outbreaks,1)
  sample.subclin.rate[i] = o.dat[i,which.outbreak]
}

# 3 - establish underlying probability of infection (Pinf)
# including that based on cumulative incidence data ####
# 3.1 Merge case and census data ####
dfprev<-merge(dfcases.analyse, dfcensus) %>% arrange(lmgid, casemonth) #merge will automatically be done on lmgid

# 3.2 Estimate infected zone Pinf for use where case totals are zero and an informed prior is required (to prevent overestimation of risk)
# In this case the approach is similar to de Vos et al (https://doi.org/10.1016/j.prevetmed.2012.01.019) where the maximum 
# cumulative incidence is used as PInf but estimated across the South African population.
# In our case we exclude the AHS controlled area to get realistic infected zone case Pinf using a more realistic population at risk

lm.ahsca<-c(210,212,213,214,215,216,217,218,219,220,221,222,223,224) # Local municipalities overlapping the AHS controlled area

PAR.RSA<-dfcensus %>%  filter(lmgid %notin% lm.ahsca) %>% summarise(sum(totalhorses))
# we know that no cases of AHS were found in the controlled zone in 2019 so no need to subset the cases dataset
cases.RSA<-dfcases.analyse %>% group_by(casemonth) %>% dplyr::summarise(totalcases = sum(ahs_cases))
data.RSA<-as.data.frame(cases.RSA)
data.RSA$par<-PAR.RSA$`sum(totalhorses)`

#note the two months where no cases occurred - in that case we account using an uninformed prior of Beta(1,1) when 
#estimating Pinf for those months acroos the whole PAR (i.e. assume 1 case at least)

p.inf<-list() # list to hold the PInf values for each area for each month
parameter.subclin.pinf<-list() # list to hold the subclin rates associated with Pinf for sensitivity analysis

#fill p.inf with NA
for (i in unique(dfprev$lmgid)) {
  for (j in unique(dfprev$casemonth)) {
    p.inf[[paste(i)]][[paste(j)]]<-rep(NA,iterations)
  }
}

# there are occurrences where Pinf is not possible to estimate when iteration permutations estimate more cases
# of AHS than there are horses in the area. A check for whether movements occurred from these areas is performed,
# in the interim an allocation of 0.5 is allocated as Pinf to these iterations and a note is logged for when it occurs

count.pinfNA<-list() # to keep track of where a Pinf of 0.5 was used

for (i in unique(dfprev$lmgid)) { #unique area
  tempdata.lm <- dfprev[which(dfprev$lmgid == i), ]
  for (j in unique(tempdata.lm$casemonth)) { #unique month - will be 1:12
    tempdata.rsa.month <- data.RSA[which(data.RSA$casemonth == j), ] # for use when case data is zero to establish informed prior
    tempdata.lm.month <- tempdata.lm[which(tempdata.lm$casemonth == j), ] # for use where cases data is > 0
    tempvector <- vector() # temporary vector holding pinfo values
    tempvector.parameter <- vector() #temporary vector for holding subclin parameter values for sensitivity analysis
    if (tempdata.lm.month$ahs_cases > 0) { # cases reported in month
      # cases exist
      for (k in 1:iterations) {
        temp.parameter.subclin.pinf <-
          sample(sample.subclin.rate,
                 size = 1,
                 replace = TRUE)
        tempvector.parameter <-
          c(tempvector.parameter, temp.parameter.subclin.pinf)
        temp <- rbeta(
          1,
          tempdata.lm.month$ahs_cases / (1 - temp.parameter.subclin.pinf) + 1,
          tempdata.lm.month$totalhorses - tempdata.lm.month$ahs_cases/(1 - temp.parameter.subclin.pinf) + 1
        )
        
        if (is.na(temp)) {
          # check for cases>PAR
          tempvector <- c(tempvector, 0.5) # concatenate 0.5 = Pinf
          count.pinfNA[[paste(i)]][[paste(j)]] <-
            1 #maintain log of this occurrence
        } else {
          tempvector <- c(tempvector, temp) # concatenate Pinf onto temporary vector
        }
      }# cases reported in month END
      
    } else { # - cases are not reported
      for (k in 1:iterations) {
        temp.parameter.subclin.pinf <-
          sample(sample.subclin.rate,
                 size = 1,
                 replace = TRUE)
        tempvector.parameter <- c(tempvector.parameter, temp.parameter.subclin.pinf)
        temp <- rbeta(
          1,
          data.RSA$totalcases / (1 - temp.parameter.subclin.pinf) + 1,
          data.RSA$par - data.RSA$totalcases / (1 - temp.parameter.subclin.pinf) + 1
        )
        tempvector <- c(tempvector, temp)
      }
    } # END cases not reported
    p.inf[[paste(i)]][[paste(j)]] <- tempvector #concatenate iterations results for area i and month j into Pinf list
    parameter.subclin.pinf[[paste(i)]][[paste(j)]] <- tempvector.parameter
  } #END unique month 
} #END unique area

# 3.2 Output of P inf ####
p.inf.median<-list() 

for (i in names(p.inf)) { # every lm id
  for (j in names(p.inf[[i]])) { # every month of year
    p.inf.median[[paste0(i,'_',j)]]$median.p.inf<-median(p.inf[[i]][[j]])
    p.inf.median[[paste0(i,'_',j)]]$lowerPI.p.inf<-quantile(p.inf[[i]][[j]], prob=0.025)
    p.inf.median[[paste0(i,'_',j)]]$upperPI.p.inf<-quantile(p.inf[[i]][[j]], prob=0.975)
    p.inf.median[[paste0(i,'_',j)]]$lm<-paste(i)
    p.inf.median[[paste0(i,'_',j)]]$casemonth<-as.numeric(paste(j))
  }
}

p.inf.summarytable <- ldply(p.inf.median, data.frame,.id = NULL) #convert list to table

# 4. Movement dataset ####
dfmovement.analyse<-dfmovement
dfmovement.analyse$movementdate<-as.Date(dfmovement.analyse$movementdate) # ensure class of date is such
dfmovement.analyse$movementmonth<-as.numeric(format(as.Date(dfmovement.analyse$movementdate), "%m")) #create month of movement

#4.1 summary statistics for Table 1 ####
table1summary<-dfmovement.analyse %>% 
  dplyr::group_by(movementtype) %>% 
  dplyr::summarise(moved = sum(totalmoved)) %>% 
  ungroup() %>% 
  dplyr::mutate (totalmoved = sum(moved))  %>%
  dplyr::mutate(percent = moved / sum(moved)) %>% as.data.frame()

write.csv(table1summary, "./modeloutputs/tables/table1_movementsummary.csv")

#4.2 aggregate movement information ####
dfmovement.analyse<-dfmovement.analyse %>% 
  dplyr::group_by(lmgid, movementmonth, movementtype) %>% 
  dplyr::summarise(total_moved_agg = sum(totalmoved)) %>% 
  ungroup() %>% arrange(movementtype, lmgid, movementmonth) %>% 
  as.data.frame()

#create zero rows for zero movement data by movement type and area 
dfmovement.analyse<-dfmovement.analyse %>% tidyr::complete(movementmonth, 
                     nesting(lmgid, movementtype), 
                     fill = list(total_moved_agg = 0)) %>% arrange(movementtype, lmgid, movementmonth)

#create vector of local municipalities where movements occurred from to make pinf_est_origin more realistic
# i.e. equids moving through SOQ, VPQO and ZebraQ likely to come from areas where equids move from in general
dfmovement.lmorigins<-unique(dfmovement.analyse$lmgid)

# B: analysis code ####
# 5. Entry risk - main analysis####
#5.1 Standard movements ##########
dfmovement.analyse.standard<-dfmovement.analyse %>% filter(movementtype == 'standard') %>% dplyr::select(-movementtype) # subset data and extract columns needed

# 5.1.1 Movements per area per month per iteration ####
movementtotals.list.standard<-list() # list to hold the movement totals for each area for each month

for (i in unique(dfmovement.analyse.standard$lmgid)) { # unique area
  movementtotals.list.standard[[paste(i)]] <- list()
  tempdata <- dfmovement.analyse.standard[which(dfmovement.analyse.standard$lmgid == i),]
  for (j in unique(tempdata$movementmonth)) { #unique month
    temp<-tempdata[which(tempdata$movementmonth == j),]
    movementtotals.list.standard[[paste(i)]][[paste(j)]]<- replicate(iterations, temp$total_moved_agg)
  } #END unique month
} #END unique area

# 5.1.2 Subclinical rates to be used in for analysis ####
# limited to areas where movement actually occurred from
psubclin.standard<-list()
for (i in unique(dfmovement.analyse.standard$lmgid)) { #unique area
  for (j in evaluationperiods) { #unique month
    psubclin.standard[[paste(i)]][[paste(j)]]<- sample(sample.subclin.rate, size=iterations, replace =T)
  } # END unique month
}# END unique area

# 5.1.3 P_detection to be used in for analysis ####
pdetect.standard<-list()
for (i in unique(dfmovement.analyse.standard$lmgid)) {# unique area
  for (j in evaluationperiods) { #unique month
    pdetect.standard[[paste(i)]][[paste(j)]]<- rpert(iterations, 0.7,0.9,1)
  } #END unique month
} #END unique area

# 5.1.4 Subset P_inf to limit to areas where movement occurred from ####
# Note that some months will still include 0 movements
p.inf.standard<-p.inf[names(movementtotals.list.standard)]

# 5.1.5 Standard movement outcome ####
p.ent.standard<-list()

for (i in names(movementtotals.list.standard)) { #unique areas
  for (j in  names(movementtotals.list.standard[[1]])) { #unique months
    if (movementtotals.list.standard[[i]][[j]][1] != 0){ #if there are at least some movements in that month
        p.ent.standard[[i]][[j]]<-1-(1-((p.inf.standard[[i]][[j]]*psubclin.standard[[i]][[j]]) + (p.inf.standard[[i]][[j]]*(1-psubclin.standard[[i]][[j]])*(1-pdetect.standard[[i]][[j]]))))^movementtotals.list.standard[[i]][[j]]
    } else {
      p.ent.standard[[i]][[j]]<-replicate(iterations, 0)
    }
  } #END unique months
} #END unique areas

#5.2 Standard Stop over movements ##########
dfmovement.analyse.soq<-dfmovement.analyse %>% filter(movementtype == 'soq') %>% dplyr::select(-movementtype)

# 5.2.1 Movements per area per month per iteration ####
movementtotals.list.soq<-list() # list to hold the movement totals for each area for each month

for (i in unique(dfmovement.analyse.soq$lmgid)) {
  movementtotals.list.soq[[paste(i)]] <- list()
  tempdata <- dfmovement.analyse.soq[which(dfmovement.analyse.soq$lmgid == i),]
  for (j in unique(tempdata$movementmonth)) {
    temp<-tempdata[which(tempdata$movementmonth == j),]
    movementtotals.list.soq[[paste(i)]][[paste(j)]]<- replicate(iterations, temp$total_moved_agg)
  }
}

# 5.2.2 Subclinical rates to be used in for analysis ####
psubclin.soq<-list()
for (i in unique(dfmovement.analyse.soq$lmgid)) {
  for (j in evaluationperiods) {
    psubclin.soq[[paste(i)]][[paste(j)]]<- sample(sample.subclin.rate, size=iterations, replace =T)
  }
}

# 5.2.3 P_detection to be used in for analysis ####
pdetect.soq<-list()
for (i in unique(dfmovement.analyse.soq$lmgid)) {
  for (j in evaluationperiods) {
    pdetect.soq[[paste(i)]][[paste(j)]]<- rpert(iterations, 0.7,0.9,1)
  }
}

# 5.2.4 PCR Sensitivity to be used in for analysis ####
sensPCR.soq<-list()
for (i in unique(dfmovement.analyse.soq$lmgid)) {
  for (j in evaluationperiods) {
    sensPCR.soq[[paste(i)]][[paste(j)]]<- rbeta(iterations, 9.65,1.19)
  }
}

# 5.2.5 Subset P_inf to limit to areas where movement occurred from ####
p.inf.soq<-p.inf[names(movementtotals.list.soq)]

# 5.2.6 Obtain median median Pinf for use as Pinf of origin linked to month ####
p.inf.soq.origin<-p.inf.summarytable %>%  
  filter(lm %in% dfmovement.lmorigins) %>% #realistically horses will originate where other movements have at least taken place from
  filter(lm %in% dfcases.lmorigins) %>% #realistically horses will originate where cases were reported at least some time in the year
  group_by(casemonth) %>%
  dplyr::arrange(desc(median.p.inf), .by_group = TRUE) %>% # order by descending to get maximum P inf where row numbers are even
  slice(ceiling(n()*0.5)) %>% #retrieve row associated with the 50th percentile for median p inf
  dplyr::select(casemonth, lm)

p.inf.est.soq<-list()
for(i in p.inf.soq.origin$casemonth){
  p.inf.est.soq[[paste(i)]]<-p.inf[[paste(p.inf.soq.origin[which(p.inf.soq.origin$casemonth==i),]$lm)]][[paste(i)]]
}

# 5.2.7 Establish clearance periods for quarantine ####
# 5.2.7.1 Incubation period ####
inc.soq<-list()
for (i in unique(dfmovement.analyse.soq$lmgid)) {
  for (j in evaluationperiods) {
    inc.soq[[paste(i)]][[paste(j)]]<- rpert(iterations, 2,6,10)
  }
}

# 5.2.7.2 Infectious (viraemic) period in horses####
inf.soq<-list()
for (i in unique(dfmovement.analyse.soq$lmgid)) {
  for (j in evaluationperiods) {
    inf.soq[[paste(i)]][[paste(j)]]<- rgamma(iterations, hip.alpha, hip.beta)
  }
}

# 5.2.7.3 Aggregation of incubation and infectious period to establish riskperiod ####
riskperiod.soq<-list()
for (i in unique(dfmovement.analyse.soq$lmgid)) {
  for (j in evaluationperiods) {
    riskperiod.soq[[paste(i)]][[paste(j)]]<-inc.soq[[paste(i)]][[paste(j)]]+inf.soq[[paste(i)]][[paste(j)]]
  }
}

# 5.2.7.4 Period to clear - 14 days quarantine + random day of infection prior to quarantine ####
periodtoclear.soq<-list()
for (i in unique(dfmovement.analyse.soq$lmgid)) {
  for (j in evaluationperiods) {
    periodtoclear.soq[[paste(i)]][[paste(j)]]<- 14+sample(iterations,x = 0:30,replace = TRUE)
  }
}

# 5.2.7.5 - Final clearance parameter ####
probclear.soq<-list() #0 = clearance occurred; 1 = clearance did not occur
for (i in unique(dfmovement.analyse.soq$lmgid)) {
  for (j in evaluationperiods) {
    probclear.soq[[paste(i)]][[paste(j)]]<-sapply(periodtoclear.soq[[paste(i)]][[paste(j)]]-riskperiod.soq[[paste(i)]][[paste(j)]], function(x) ifelse(x>=0, 0, 1))
  }
}

# 5.2.8 SOQ movement output ####
p.ent.soq.pathway1<-list()
p.ent.soq.pathway2<-list()
p.ent.soq<-list()

for (i in names(movementtotals.list.soq)) { # for every area where soq movements came from
  for (j in  names(movementtotals.list.soq[[1]])) { # for every month of the year
    if (movementtotals.list.soq[[i]][[j]][1] != 0){ # if there are at least some movements
      p.ent.soq.pathway1[[i]][[j]]<-probclear.soq[[i]][[j]] * ((p.inf.est.soq[[j]] * (1-sensPCR.soq[[i]][[j]]) * psubclin.soq[[i]][[j]]) +  
        (p.inf.est.soq[[j]] * (1-sensPCR.soq[[i]][[j]]) * (1-psubclin.soq[[i]][[j]]) * (1-pdetect.soq[[i]][[j]])))
      p.ent.soq.pathway2[[i]][[j]]<-(p.inf.soq[[i]][[j]] * (1-(sensPCR.soq[[i]][[j]]/2)) * psubclin.soq[[i]][[j]]) +  
        (p.inf.soq[[i]][[j]] * (1-(sensPCR.soq[[i]][[j]]/2)) * (1-psubclin.soq[[i]][[j]]) * (1-pdetect.soq[[i]][[j]]))
    } else { #no movements in that month
      p.ent.soq.pathway1[[i]][[j]]<-replicate(iterations, 0)
      p.ent.soq.pathway2[[i]][[j]]<-replicate(iterations, 0)
    }
  }
}
#summing two non-mutually exclusive probabilities - i.e. horses can be infected in Pathway 1 AND/OR Pathway 2

for (i in names(movementtotals.list.soq)) { # for every lm where soq movements came from
  for (j in  names(movementtotals.list.soq[[1]])) { # for every month of the year
      p.ent.soq[[i]][[j]]<-1-(1-(p.ent.soq.pathway1[[i]][[j]]+
                                     p.ent.soq.pathway2[[i]][[j]]-
                                     (p.ent.soq.pathway1[[i]][[j]]*p.ent.soq.pathway2[[i]][[j]])))^movementtotals.list.soq[[i]][[j]]
    
  }
}

#5.3 Vector Protected Quarantine at Origin - VPQO##########
dfmovement.analyse.vpqo<-dfmovement.analyse %>% filter(movementtype == 'vpqo') %>% dplyr::select(-movementtype)

# 5.3.1 Movements per area per month per iteration ####
movementtotals.list.vpqo<-list()

for (i in unique(dfmovement.analyse.vpqo$lmgid)) { #unique area
  movementtotals.list.vpqo[[paste(i)]] <- list()
  tempdata <- dfmovement.analyse.vpqo[which(dfmovement.analyse.vpqo$lmgid == i),]
  for (j in unique(tempdata$movementmonth)) { #unique month
    temp<-tempdata[which(tempdata$movementmonth == j),]
    movementtotals.list.vpqo[[paste(i)]][[paste(j)]]<- replicate(iterations, temp$total_moved_agg)
  } # END unique area
} # END unique month

# 5.3.2 Subclinical rates to be used in for analysis ####
psubclin.vpqo<-list()
for (i in unique(dfmovement.analyse.vpqo$lmgid)) {
  for (j in evaluationperiods) {
    psubclin.vpqo[[paste(i)]][[paste(j)]]<- sample(sample.subclin.rate, size=iterations, replace =T)
  }
}

# 5.3.3 P_detection to be used in for analysis ####
pdetect.vpqo<-list()
for (i in unique(dfmovement.analyse.vpqo$lmgid)) {
  for (j in evaluationperiods) {
    pdetect.vpqo[[paste(i)]][[paste(j)]]<- rpert(iterations, 0.7,0.9,1)
  }
}

# 5.3.4 PCR Sensitivity to be used in for analysis (entry and exit here) ####
sensPCR.vpqo<-list()
for (i in unique(dfmovement.analyse.vpqo$lmgid)) {
  for (j in evaluationperiods) {
    sensPCR.vpqo[[paste(i)]][[paste(j)]]<- rbeta(iterations, 9.65,1.19)
  }
}

# 5.3.5 Subset P_inf to limit to areas where movement occurred from ####
p.inf.vpqo<-p.inf[names(movementtotals.list.vpqo)]

# 5.3.6 Obtain maximum median Pinf for use as Pinf of origin linked to month ####
p.inf.vpqo.origin<-p.inf.summarytable %>%  
  filter(lm %in% dfmovement.lmorigins) %>% #realistically horses will originate where other movements have at least taken place from
  filter(lm %in% dfcases.lmorigins) %>% #realistically horses will originate where cases were reported
  group_by(casemonth) %>%
  dplyr::arrange(desc(median.p.inf), .by_group = TRUE) %>% # order by descending to get maximum P inf where row numbers are even
  slice(ceiling(n()*0.5)) %>% #retrieve row associated with the 50th percentile for median p inf
  dplyr::select(casemonth, lm)
p.inf.est.vpqo<-list()

for(i in p.inf.vpqo.origin$casemonth){
  p.inf.est.vpqo[[paste(i)]]<-p.inf[[paste(p.inf.vpqo.origin[which(p.inf.vpqo.origin$casemonth==i),]$lm)]][[paste(i)]]
}

# 5.3.7 Establish clearance periods for quarantine ####
# 5.3.7.1 Incubation period ####
inc.vpqo<-list()
for (i in unique(dfmovement.analyse.vpqo$lmgid)) {
  for (j in evaluationperiods) {
    inc.vpqo[[paste(i)]][[paste(j)]]<- rpert(iterations, 2,6,10)
  }
}

# 5.3.7.2 Infectious (viraemic) period in horses####
inf.vpqo<-list()
for (i in unique(dfmovement.analyse.vpqo$lmgid)) {
  for (j in evaluationperiods) {
    inf.vpqo[[paste(i)]][[paste(j)]]<- rgamma(iterations, 29.75,4.9586)
  }
}

# 5.3.7.3 Aggregation of incubation and infectious period to establish riskperiod ####
riskperiod.vpqo<-list()
for (i in unique(dfmovement.analyse.vpqo$lmgid)) {
  for (j in evaluationperiods) {
    riskperiod.vpqo[[paste(i)]][[paste(j)]]<-inc.vpqo[[paste(i)]][[paste(j)]]+inf.vpqo[[paste(i)]][[paste(j)]]
  }
}

# 5.3.7.4 Period to clear - 14 days quarantine + random day of infection prior to quarantine ####
periodtoclear.vpqo<-list()
for (i in unique(dfmovement.analyse.vpqo$lmgid)) {
  for (j in evaluationperiods) {
    periodtoclear.vpqo[[paste(i)]][[paste(j)]]<- 14+sample(iterations,x = 0:30,replace = TRUE)
  }
}

# 5.3.7.5 - Final clearance parameter ####
probclear.vpqo<-list() #0 = clear; 1 = not clear
for (i in unique(dfmovement.analyse.vpqo$lmgid)) {
  for (j in evaluationperiods) {
    probclear.vpqo[[paste(i)]][[paste(j)]]<-sapply(periodtoclear.vpqo[[paste(i)]][[paste(j)]]-riskperiod.vpqo[[paste(i)]][[paste(j)]], function(x) ifelse(x>=0, 0, 1))
  }
}
# 5.3.8 vpqo movement aggregation ####
p.ent.vpqo.pathway1<-list()
p.ent.vpqo.pathway2<-list()
p.ent.vpqo<-list()

for (i in names(movementtotals.list.vpqo)) { # for every lm where vpqo movements came from
  for (j in  names(movementtotals.list.vpqo[[1]])) { # for every month of the year
    if (movementtotals.list.vpqo[[i]][[j]][1] != 0){ #if there are at least some movements
      p.ent.vpqo.pathway1[[i]][[j]]<-probclear.vpqo[[i]][[j]] * ((p.inf.est.vpqo[[j]] * (1-sensPCR.vpqo[[i]][[j]])^2 * psubclin.vpqo[[i]][[j]]) +  
                                                                        (p.inf.est.vpqo[[j]] * (1-sensPCR.vpqo[[i]][[j]])^2 * (1-psubclin.vpqo[[i]][[j]]) * (1-pdetect.vpqo[[i]][[j]])))
      p.ent.vpqo.pathway2[[i]][[j]]<-(p.inf.vpqo[[i]][[j]] * (1-(sensPCR.vpqo[[i]][[j]]/2)) * psubclin.vpqo[[i]][[j]]) +  
        (p.inf.vpqo[[i]][[j]] * (1-(sensPCR.vpqo[[i]][[j]]/2)) * (1-psubclin.vpqo[[i]][[j]]) * (1-pdetect.vpqo[[i]][[j]]))
    } else {
      p.ent.vpqo.pathway1[[i]][[j]]<-replicate(iterations, 0)
      p.ent.vpqo.pathway2[[i]][[j]]<-replicate(iterations, 0)
    }
  }
}

#summing two non-mutually exclusive probabilities - i.e. horses can be infected in Pathway 1 AND/OR Pathway 2

for (i in names(movementtotals.list.vpqo)) { # for every lm where vpqo movements came from
  for (j in  names(movementtotals.list.vpqo[[1]])) { # for every month of the year
    p.ent.vpqo[[i]][[j]]<-1-(1-(p.ent.vpqo.pathway1[[i]][[j]]+
                                   p.ent.vpqo.pathway2[[i]][[j]]-
                                   (p.ent.vpqo.pathway1[[i]][[j]]*p.ent.vpqo.pathway2[[i]][[j]])))^movementtotals.list.vpqo[[i]][[j]]
    
  }
}

#5.4 Controlled area vector Protected Stop over movements VPSOQ##########

dfmovement.analyse.vpsoq<-dfmovement.analyse %>% filter(movementtype == 'vpsoq') %>% dplyr::select(-movementtype)
# 5.4.1 Movements per area per month per iteration ####

movementtotals.list.vpsoq<-list() # list to hold the movement totals for each area for each month

for (i in unique(dfmovement.analyse.vpsoq$lmgid)) { #unique area
  movementtotals.list.vpsoq[[paste(i)]] <- list()
  tempdata <- dfmovement.analyse.vpsoq[which(dfmovement.analyse.vpsoq$lmgid == i),]
  for (j in unique(tempdata$movementmonth)) { #unique month
    temp<-tempdata[which(tempdata$movementmonth == j),]
    movementtotals.list.vpsoq[[paste(i)]][[paste(j)]]<- replicate(iterations, temp$total_moved_agg)
  } #END unique month
} #END unique area

# 5.4.2 Subclinical rates to be used in for analysis ####
psubclin.vpsoq<-list()
for (i in unique(dfmovement.analyse.vpsoq$lmgid)) {
  for (j in evaluationperiods) {
    psubclin.vpsoq[[paste(i)]][[paste(j)]]<- sample(sample.subclin.rate, size=iterations, replace =T)
  }
}

# 5.4.3 P_detection to be used in for analysis ####
pdetect.vpsoq<-list()
for (i in unique(dfmovement.analyse.vpsoq$lmgid)) {
  for (j in evaluationperiods) {
    pdetect.vpsoq[[paste(i)]][[paste(j)]]<- rpert(iterations, 0.7,0.9,1)
  }
}

# 5.4.4 PCR Sensitivity to be used in for analysis (entry and exit here) ####
sensPCR.vpsoq<-list()
for (i in unique(dfmovement.analyse.vpsoq$lmgid)) {
  for (j in evaluationperiods) {
    sensPCR.vpsoq[[paste(i)]][[paste(j)]]<- rbeta(iterations, 9.65,1.19)
  }
}

# 5.4.5 Subset P_inf to limit to areas where movement occurred from ####
p.inf.vpsoq<-p.inf[names(movementtotals.list.vpsoq)]

# 5.4.6 Establish clearance periods for quarantine ####
# 5.4.6.1 Incubation period ####
inc.vpsoq<-list()
for (i in unique(dfmovement.analyse.vpsoq$lmgid)) {
  for (j in evaluationperiods) {
    inc.vpsoq[[paste(i)]][[paste(j)]]<- rpert(iterations, 2,6,10)
  }
}

# 5.4.6.2 Infectious (viraemic) period in horses####
inf.vpsoq<-list()
for (i in unique(dfmovement.analyse.vpsoq$lmgid)) {
  for (j in evaluationperiods) {
    inf.vpsoq[[paste(i)]][[paste(j)]]<- rgamma(iterations, 29.75,4.9586)
  }
}

# 5.4.6.3 Aggregation of incubation and infectious period to establish riskperiod ####
riskperiod.vpsoq<-list()
for (i in unique(dfmovement.analyse.vpsoq$lmgid)) {
  for (j in evaluationperiods) {
    riskperiod.vpsoq[[paste(i)]][[paste(j)]]<-inc.vpsoq[[paste(i)]][[paste(j)]]+inf.vpsoq[[paste(i)]][[paste(j)]]
  }
}

# 5.4.6.4 Period to clear - 14 days quarantine + random day of infection prior to quarantine ####
periodtoclear.vpsoq<-list()
for (i in unique(dfmovement.analyse.vpsoq$lmgid)) {
  for (j in evaluationperiods) {
    periodtoclear.vpsoq[[paste(i)]][[paste(j)]]<- 14+sample(iterations,x = 0:30,replace = TRUE)
  }
}

# 5.4.6.5 - Final clearance parameter ####
probclear.vpsoq<-list() 
for (i in unique(dfmovement.analyse.vpsoq$lmgid)) {
  for (j in evaluationperiods) {
    probclear.vpsoq[[paste(i)]][[paste(j)]]<-sapply(periodtoclear.vpsoq[[paste(i)]][[paste(j)]]-riskperiod.vpsoq[[paste(i)]][[paste(j)]], function(x) ifelse(x>=0, 0, 1))
  }
}

# 5.4.7 VPSOQ movement aggregation ####
# Note no pathway 2 here - refer to the manuscript
p.ent.vpsoq<-list()

for (i in names(movementtotals.list.vpsoq)) { # for every lm where vpsoq movements came from
  for (j in  names(movementtotals.list.vpsoq[[1]])) { # for every month of the year
    if (movementtotals.list.vpsoq[[i]][[j]][1] != 0){ #if there are at least some movements (i.e. the first element of the movement total list for that area and month then the p.ent is NOT 0 and not 1 which is what you get by powering to 0 for the non mutual exclusive adding of probability)
      p.ent.vpsoq[[i]][[j]]<-probclear.vpsoq[[i]][[j]] * ((p.inf.vpsoq[[i]][[j]] * (1-sensPCR.vpsoq[[i]][[j]])^2 * psubclin.vpsoq[[i]][[j]]) +  
                                                                        (p.inf.vpsoq[[i]][[j]] * (1-sensPCR.vpsoq[[i]][[j]])^2 * (1-psubclin.vpsoq[[i]][[j]]) * (1-pdetect.vpsoq[[i]][[j]])))
    } else {
      p.ent.vpsoq[[i]][[j]]<-replicate(iterations, 0)
    }
  }
}

for (i in names(movementtotals.list.vpsoq)) { # for every lm where vpsoq movements came from
  for (j in  names(movementtotals.list.vpsoq[[1]])) { # for every month of the year
    p.ent.vpsoq[[i]][[j]]<-1-(1-(p.ent.vpsoq[[i]][[j]]
                                   ))^movementtotals.list.vpsoq[[i]][[j]]
  }
}

#5.5 Zebra movements ##########
dfmovement.analyse.zebq<-dfmovement.analyse %>% filter(movementtype == 'zebq') %>% dplyr::select(-movementtype)

# 5.5.1 Movements per area per month per iteration ####
movementtotals.list.zebq<-list() # list to hold the movement totals for each area for each month

for (i in unique(dfmovement.analyse.zebq$lmgid)) {
  movementtotals.list.zebq[[paste(i)]] <- list()
  tempdata <- dfmovement.analyse.zebq[which(dfmovement.analyse.zebq$lmgid == i),]
  for (j in unique(tempdata$movementmonth)) {
    temp<-tempdata[which(tempdata$movementmonth == j),]
    movementtotals.list.zebq[[paste(i)]][[paste(j)]]<- replicate(iterations, temp$total_moved_agg)
  }
}

# 5.5.2 PCR Sensitivity to be used in for analysis (entry and exit here) ####
sensPCR.zebq<-list()
for (i in unique(dfmovement.analyse.zebq$lmgid)) {
  for (j in evaluationperiods) {
    sensPCR.zebq[[paste(i)]][[paste(j)]]<- rbeta(iterations, 9.65,1.19)
  }
}

# 5.5.3 Subset P_inf to limit to areas where movement occurred from ####
p.inf.zebq<-p.inf[names(movementtotals.list.zebq)]

# 5.5.4 Obtain maximum median Pinf for use as Pinf of origin linked to month ####
p.inf.zebq.origin<-p.inf.summarytable %>%  
  filter(lm %in% dfmovement.lmorigins) %>% #realistically zebra will originate where other movements have at least taken place from
  # filter(lm %in% dfcases.lmorigins) %>% # zebra can't move from high risk areas so this filter removed for this pathway
  group_by(casemonth) %>%
  dplyr::arrange(desc(median.p.inf), .by_group = TRUE) %>% # order by descending to get maximum P inf where row numbers are even
  slice(ceiling(n()*0.5)) %>% #retrieve row associated with the 50th percentile for median p inf
  dplyr::select(casemonth, lm)
p.inf.est.zebq<-list()

for(i in p.inf.zebq.origin$casemonth){
  p.inf.est.zebq[[paste(i)]]<-p.inf[[paste(p.inf.zebq.origin[which(p.inf.zebq.origin$casemonth==i),]$lm)]][[paste(i)]]
}

# 5.5.5 Establish clearance periods for quarantine ####
# 5.5.5.1 Incubation period ####
inc.zebq<-list()
for (i in unique(dfmovement.analyse.zebq$lmgid)) {
  for (j in evaluationperiods) {
    inc.zebq[[paste(i)]][[paste(j)]]<- rpert(iterations, 2,6,10)
  }
}

# 5.5.5.2 Infectious (viraemic) period in zebra####
inf.zebq<-list()
for (i in unique(dfmovement.analyse.zebq$lmgid)) {
  for (j in evaluationperiods) {
    inf.zebq[[paste(i)]][[paste(j)]]<- rgamma(iterations, zip.alpha,zip.beta)
  }
}

# 5.5.5.3 Aggregation of incubation and infectious period to establish riskperiod ####
riskperiod.zebq<-list()
for (i in unique(dfmovement.analyse.zebq$lmgid)) {
  for (j in evaluationperiods) {
    riskperiod.zebq[[paste(i)]][[paste(j)]]<-inc.zebq[[paste(i)]][[paste(j)]]+inf.zebq[[paste(i)]][[paste(j)]]
  }
}

# 5.5.5.4 Period to clear - 21 days quarantine + random day of infection prior to quarantine ####
periodtoclear.zebq<-list()
for (i in unique(dfmovement.analyse.zebq$lmgid)) {
  for (j in evaluationperiods) {
    periodtoclear.zebq[[paste(i)]][[paste(j)]]<- 21+sample(iterations,x = 0:30,replace = TRUE)
  }
}

# 5.5.5.5 - Final clearance parameter ####
probclear.zebq<-list() 
for (i in unique(dfmovement.analyse.zebq$lmgid)) {
  for (j in evaluationperiods) {
    probclear.zebq[[paste(i)]][[paste(j)]]<-sapply(periodtoclear.zebq[[paste(i)]][[paste(j)]]-riskperiod.zebq[[paste(i)]][[paste(j)]], function(x) ifelse(x>=0, 0, 1))
  }
}

# 5.5.6 Zebra Q movement aggregation ####
p.ent.zebq.pathway1<-list()
p.ent.zebq.pathway2<-list()
p.ent.zebq<-list()


for (i in names(movementtotals.list.zebq)) { # for every lm where zebq movements came from
  for (j in  names(movementtotals.list.zebq[[1]])) { # for every month of the year
    if (movementtotals.list.zebq[[i]][[j]][1] != 0){ #if there are at least some movements (i.e. the first element of the movement total list for that area and month then the p.ent is NOT 0 and not 1 which is what you get by powering to 0 for the non mutual exclusive adding of probability)
      p.ent.zebq.pathway1[[i]][[j]]<-probclear.zebq[[i]][[j]] * (p.inf.est.zebq[[j]] * (1-sensPCR.zebq[[i]][[j]])^2)
      p.ent.zebq.pathway2[[i]][[j]]<-p.inf.zebq[[i]][[j]] * (1-(sensPCR.zebq[[i]][[j]]/2))
    } else {
      p.ent.zebq.pathway1[[i]][[j]]<-replicate(iterations, 0)
      p.ent.zebq.pathway2[[i]][[j]]<-replicate(iterations, 0)
    }
  }
}

#summing two non-mutually exclusive probabilities - i.e. horses can be infected in Pathway 1 AND/OR Pathway 2

for (i in names(movementtotals.list.zebq)) { # for every lm where zebq movements came from
  for (j in  names(movementtotals.list.zebq[[1]])) { # for every month of the year
    p.ent.zebq[[i]][[j]]<-1-(1-(p.ent.zebq.pathway1[[i]][[j]]+
                                     p.ent.zebq.pathway2[[i]][[j]]-
                                     (p.ent.zebq.pathway1[[i]][[j]]*p.ent.zebq.pathway2[[i]][[j]])))^movementtotals.list.zebq[[i]][[j]]
    
  }
}

# 5.6 What if analysis - Uncontrolled movement ##########
# Establish the areas where the 50th percentile Pinf was per month and allocate the SOQ and VPQO and Zebra movements to those LM's

p.inf.est.lm<-p.inf.summarytable %>%  
  filter(lm %in% dfmovement.lmorigins) %>%
  filter(lm %in% dfcases.lmorigins) %>%
  group_by(casemonth) %>%
  dplyr::arrange(desc(median.p.inf), .by_group = TRUE) %>% # order by descending to get maximum P inf where row numbers are even
  slice(ceiling(n()*0.5)) %>% #retrieve row associated with the 50th percentile for median p inf
  dplyr::select(casemonth, lm)

# View(p.inf.summarytable %>% filter(casemonth==1))

dfmovement.analyse.estpinf.totalmoved<-dfmovement.analyse %>% 
  filter(movementtype == 'soq'| movementtype == 'vpqo'| movementtype == 'zebq') %>% 
  group_by(movementmonth) %>% dplyr::summarise(total_moved_agg = sum(total_moved_agg))

dfmovement.analyse.estpinf.totalmoved.agg<-inner_join(x = dfmovement.analyse.estpinf.totalmoved, 
                                                      y = p.inf.est.lm, 
                                                      by = c("movementmonth" = "casemonth")) %>% 
  dplyr::rename("lmgid" = "lm")

dfmovement.analyse.estpinf.totalmoved.agg$movementtype = "estpinf"

# ensure fill of zero movements for months where estpinf movements occurred in
dfmovement.analyse.estpinf.totalmoved.agg<-dfmovement.analyse.estpinf.totalmoved.agg %>% tidyr::complete(movementmonth, 
                                                    nesting(lmgid, movementtype), 
                                                    fill = list(total_moved_agg = 0)) %>% arrange(movementtype, lmgid, movementmonth)

dfmovement.analyse.nocontrol<-dfmovement.analyse %>% 
  filter(movementtype != 'soq'& movementtype != 'vpqo'& movementtype != 'zebq')

dfmovement.analyse.nocontrol<-rbind(dfmovement.analyse.nocontrol, dfmovement.analyse.estpinf.totalmoved.agg)
dfmovement.analyse.nocontrol<-dfmovement.analyse.nocontrol %>% group_by(movementmonth, lmgid) %>% dplyr::summarise(total_moved_agg = sum(total_moved_agg))

# 5.6.1 Movements per area per month per iteration ####
movementtotals.list.nocontrol<-list() # list to hold the movement totals for each area for each month

for (i in unique(dfmovement.analyse.nocontrol$lmgid)) { # unique area
  movementtotals.list.nocontrol[[paste(i)]] <- list()
  tempdata <- dfmovement.analyse.nocontrol[which(dfmovement.analyse.nocontrol$lmgid == i),]
  for (j in unique(tempdata$movementmonth)) { #unique month
    temp<-tempdata[which(tempdata$movementmonth == j),]
    movementtotals.list.nocontrol[[paste(i)]][[paste(j)]]<- replicate(iterations, temp$total_moved_agg)
  } #END unique month
} #END unique area


# 5.6.2 Subset P_inf to limit to areas where movement occurred from ####
# Note that some months will still include 0 movements
p.inf.nocontrol<-p.inf[names(movementtotals.list.nocontrol)]

# 5.6.5 nocontrol movement aggregation ####
p.ent.nocontrol<-list()

for (i in names(movementtotals.list.nocontrol)) { #unique areas
  for (j in  names(movementtotals.list.nocontrol[[1]])) { #unique months
    if (movementtotals.list.nocontrol[[i]][[j]][1] != 0){ #if there are at least some movements in that month
      p.ent.nocontrol[[i]][[j]]<-1-(1-(p.inf.nocontrol[[i]][[j]]))^movementtotals.list.nocontrol[[i]][[j]]
    } else {
      p.ent.nocontrol[[i]][[j]]<-replicate(iterations, 0)
    }
  } #END unique months
} #END unique areas
# 

# 6 Sensitivity analysis ####
# Done prior to aggregation to whole country (i.e. at municipality level)
# After considering the insubstantial risk of entry (and lack of movements) through all but the 
# standard and soq movements, the sensitivity analysis was only performed on standard and standard soq movements,
# and only for those areas and months where movement took place (i.e. zero risk of entry was excluded)
# The SOQ movements were added to get some feeling for incubation, infectious period and PCR sensitivity parameters

# 6.1 Standard movements ####

outputlist.sensitivity.standard<-list()
for (i in names(p.ent.standard)) { # every lm
  for (j in names(p.ent.standard[[i]])) { # every month of year
    # check that movements occurred
    if (sum(movementtotals.list.standard[[i]][[j]])>0) {
      outputlist.sensitivity.standard[[paste0(i,'_',j)]]$outputpintro<-p.ent.standard[[i]][[j]] 
      outputlist.sensitivity.standard[[paste0(i,'_',j)]]$subclin.pinf<-parameter.subclin.pinf[[i]][[j]]
      outputlist.sensitivity.standard[[paste0(i,'_',j)]]$pinf<-p.inf[[i]][[j]]
      outputlist.sensitivity.standard[[paste0(i,'_',j)]]$subclin.vhc<-psubclin.standard[[i]][[j]]
      outputlist.sensitivity.standard[[paste0(i,'_',j)]]$pdetect<-pdetect.standard[[i]][[j]]
    }
  }
}

# 6.2 SOQ movements ####
outputlist.sensitivity.soq<-list()
for (i in names(p.ent.soq)) { # every lm
  for (j in names(p.ent.soq[[i]])) { # every month of year
    if (sum(movementtotals.list.soq[[i]][[j]])>0) {    # check that movements occurred
      outputlist.sensitivity.soq[[paste0(i,'_',j)]]$outputpintro<-p.ent.soq[[i]][[j]] 
      outputlist.sensitivity.soq[[paste0(i,'_',j)]]$subclin.pinf<-parameter.subclin.pinf[[i]][[j]]
      outputlist.sensitivity.soq[[paste0(i,'_',j)]]$pinf_o<-p.inf.est.soq[[j]]
      outputlist.sensitivity.soq[[paste0(i,'_',j)]]$pinf_q<-p.inf.soq[[i]][[j]]
      outputlist.sensitivity.soq[[paste0(i,'_',j)]]$inc<-inc.soq[[i]][[j]]
      outputlist.sensitivity.soq[[paste0(i,'_',j)]]$inf<-inf.soq[[i]][[j]]
      outputlist.sensitivity.soq[[paste0(i,'_',j)]]$pcr_se<-sensPCR.soq[[i]][[j]]   
      outputlist.sensitivity.soq[[paste0(i,'_',j)]]$subclin.vhc<-psubclin.soq[[i]][[j]]
      outputlist.sensitivity.soq[[paste0(i,'_',j)]]$pdetect<-pdetect.soq[[i]][[j]]
    } #end if
  } #end j
} #end i

# C: Model output Code ####

lm_sf <- st_read("./ahsv_entry_assessment_zafcontrolledarea-main/gis/lm.shp")

# 7. Figure 1: Map of movements per local municipality ####
movement.map.annualtotal<-dfmovement.analyse %>% group_by(lmgid) %>% dplyr::summarise(total_annual = sum(total_moved_agg))

movement.map.monthly.max<-dfmovement.analyse %>% 
  group_by(lmgid, movementmonth) %>% 
  dplyr::summarise(total_month = sum(total_moved_agg)) %>% 
  slice_max(total_month, n = 1, with_ties = FALSE) %>% 
  as.data.frame() %>% #establish the maximum per lm
  arrange(as.numeric(lmgid))

movement.map<-inner_join(x = movement.map.annualtotal, y = movement.map.monthly.max, by = c("lmgid" = "lmgid")) 
movement.map.spatial <- merge(lm_sf, movement.map, by.x = "gid", by.y = "lmgid")
st_write(movement.map.spatial, paste0("./modeloutputs/gis/Fig1_movement.shp"))

# 8. Probability of infection ####
# Analyisis 8.1: Establish the maximum p_inf per local municipality and in which month that occurred in ####
p.inf.summarytable.max<-p.inf.summarytable %>% group_by(lm) %>% slice_max(median.p.inf, n = 1, with_ties = FALSE) %>% as.data.frame()

# 8.1.1. Figure 2: Histogram of number of LM's and when they had maximum p_inf during the year ####
#Note that LM's are only included where either movements took place from or where cases were reported 

ggplot(data = p.inf.summarytable.max, aes(x = as.factor(casemonth))) + geom_histogram(stat = "count") +
  theme_bw() +
  xlab("Month of year") +
  ylab("Number of areas where annual maximum \nprobability of infection occurred") +
  theme(
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size = rel(0.5)),
    axis.text = element_text(size = rel(0.4)),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank()
  )

ggsave("./modeloutputs/figures/Figure2.tiff",  width = 6.68, height = 6.68, units = "cm",dpi = 300) # size and units as req by PloS

# 8.1.2 Spatial distribution of maximum p_inf ####
p.inf.summarytable.max.spatial <- merge(lm_sf, p.inf.summarytable.max, by.x = "gid", by.y = "lm")
st_write(p.inf.summarytable.max.spatial, paste0("./modeloutputs/gis/Fig3_MaxPinf.shp"))

# 9: Probability of entry ####

# Bind all p.ent dataframes into a single one after converting list to dataframe
p.ent.overall.lm.agg.analyse<-rbind(
  plyr::ldply (p.ent.standard, data.frame) %>% dplyr::group_by(.id) %>% dplyr::mutate(iteration = row_number(), movementtype = "standard"),
  plyr::ldply (p.ent.soq, data.frame) %>% dplyr::group_by(.id) %>% dplyr::mutate(iteration = row_number(), movementtype = "soq"),
  plyr::ldply (p.ent.vpqo, data.frame) %>% dplyr::group_by(.id) %>% dplyr::mutate(iteration = row_number(), movementtype = "vpqo"),
  plyr::ldply (p.ent.vpsoq, data.frame) %>% dplyr::group_by(.id) %>% dplyr::mutate(iteration = row_number(), movementtype = "vpsoq"),
  plyr::ldply (p.ent.zebq, data.frame) %>% dplyr::group_by(.id) %>% dplyr::mutate(iteration = row_number(), movementtype = "zebq")
)

p.ent.overall.lm.agg.analyse$lmgid<-p.ent.overall.lm.agg.analyse$.id #rename local municipality field
p.ent.overall.lm.agg.analyse<-as.data.frame(p.ent.overall.lm.agg.analyse) #ensure output is dataframe
p.ent.overall.lm.agg.analyse<-p.ent.overall.lm.agg.analyse %>% dplyr::select(-.id) #drop .id field
# every month represented by X_monthnumber

#Analysis 9.1: Monthly overall probability of entry from all areas with all movement types ####
p.ent.overall<-p.ent.overall.lm.agg.analyse  %>% 
  mutate_at(vars(contains("X")),function(x) (1-x)) %>% # probability of freedom for every data point
  dplyr::select(-c(movementtype,lmgid)) %>% group_by(iteration) %>% # start aggregating different movement types
  dplyr::summarise_all(prod) %>% # probability all movement types result in freedom (non-entry)  - across municipalities
  mutate_at(vars(contains("X")),function(x) (1-x)) %>% # for every month (column) probability of entry through at least one movement type
  sapply(., quantile, probs=c(0.025, 0.5, 0.975)) %>% as.data.frame() %>% mutate(estimate = c("2.5%","50%","97.5%")) %>% dplyr::select(-iteration) #summarize outcome with 95% CI

# reorganise data to month, CI range
p.ent.overall.agg.analyse.all.output<-tidyr::gather(p.ent.overall, month, pent, -estimate, 5) %>% 
  pivot_wider(names_from = estimate, values_from = pent) %>% as.data.frame()

#rename columns and rows 
p.ent.overall.agg.analyse.all.output<-p.ent.overall.agg.analyse.all.output %>% mutate_all(~gsub("X", "", .)) # drop X from month names
p.ent.overall.agg.analyse.all.output<-sapply(p.ent.overall.agg.analyse.all.output, as.numeric) %>% as.data.frame()
#final output
round_df(p.ent.overall.agg.analyse.all.output, 5)
write.csv(round_df(p.ent.overall.agg.analyse.all.output, 5), "./modeloutputs/tables/table3_allmovements_monthly.csv")
write.csv(round_df(p.ent.overall.agg.analyse.all.output, 5), "./modeloutputs/tables/table4_allmovements_monthly_controlled.csv")

#Graphing for Plot ####
ggplot(data = p.ent.overall.agg.analyse.all.output,
       aes(x = as.factor(month), y = `50%`)) +
  geom_errorbar(
    aes(
      x = as.factor(month),
      ymin = `2.5%`,
      ymax = `97.5%`
    ),
    color = "dark grey",
    width = rel(0.3),
    size = rel(0.3)
  ) +
  geom_point(size = rel(0.5)) +
  theme_bw() +
  xlab("Month of year") +
  ylab("Median probability of entry \n with 95% CI") +
  theme(
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size = rel(0.5)),
    axis.text = element_text(size = rel(0.4)),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank()
  )

ggsave("./modeloutputs/figures/Figure4.tiff",  width = 6.68, height = 6.68, units = "cm",dpi = 300) # size and units as req by PloS

#Analysis 9.2: Overall probability of entry from all areas with all movement types,  across the year ####
p.ent.overall.annual<-p.ent.overall.lm.agg.analyse  %>% 
  mutate_at(vars(contains("X")),function(x) (1-x)) %>% # probability of freedom for every data point
  dplyr::select(-c(movementtype,lmgid)) %>% group_by(iteration) %>% # start aggregating different movement types
  dplyr::summarise_all(prod) %>% ungroup() %>% as.data.frame() %>% # probability all movement types result in freedom (non-entry)  - across municipalities
  rowwise() %>% 
  dplyr::mutate(pent = 1-prod(c_across(starts_with("X")))) %>% dplyr::select(-starts_with("X") )%>% # for every month (column) probability of entry through at least one movement type
  sapply(., quantile, probs=c(0.025, 0.5, 0.975)) %>% as.data.frame() %>% dplyr::select(-iteration)

#final outcome - year around
write.csv(round_df(p.ent.overall.annual,5), "./modeloutputs/tables/table3_allmovements_annual.csv")
write.csv(round_df(p.ent.overall.annual,5), "./modeloutputs/tables/table4_allmovements_annual_controlled.csv")

# Analysis 9.3: P_ent by movement type - For Table 3####
p.ent.movementtype<-p.ent.overall.lm.agg.analyse  %>% 
  mutate_at(vars(contains("X")),function(x) (1-x)) %>% # probability of freedom for every data point
  dplyr::select(-c(lmgid)) %>% group_by(movementtype, iteration) %>% # start aggregating different movement types
  dplyr::summarise_all(prod) %>% # probability all movement types result in freedom (non-entry) - across municipalities
  mutate_at(vars(contains("X")),function(x) (1-x))

p.ent.movementtype.output<-data.frame()

for(i in unique(p.ent.movementtype$movementtype)){
  subset<-p.ent.movementtype %>% filter(movementtype == i)
  subset.output<-sapply(subset[-1], quantile, probs=c(0.025, 0.5, 0.975)) %>%  #-1 to remove movementtype
    as.data.frame() %>% 
    mutate(estimate = c("2.5%","50%","97.5%")) %>% 
    dplyr::select(-iteration) %>% tidyr::gather(., month, pent, -estimate, 5) %>% 
    pivot_wider(names_from = estimate, values_from = pent) %>% as.data.frame() %>% mutate_all(~gsub("X", "", .)) %>% 
    sapply(., as.numeric) %>% as.data.frame() %>% mutate(movementtype = i)
  
  p.ent.movementtype.output<-rbind(round_df(subset.output,5), p.ent.movementtype.output)
}

write.csv(p.ent.movementtype.output, "./modeloutputs/tables/table3_movementspecific_monthly.csv")

# Check which months would have had any movements to differentiate bewteen zero risk and NULL risk for the table
# This NB for some movement types where 95% CI was (0-0)
dfmovement.analyse %>% group_by(movementtype, movementmonth) %>% dplyr::summarise(totalmoved = sum(total_moved_agg)) %>% as.data.frame()

#9.4 - cumulative overall p_ent by movement type across the year
p.ent.overall.annual.movement<-p.ent.overall.lm.agg.analyse  %>% 
  mutate_at(vars(contains("X")),function(x) (1-x)) %>% # probability of freedom for every data point
  dplyr::select(-c(lmgid)) %>% group_by(movementtype, iteration) %>% # start aggregating different movement types
  dplyr::summarise_all(prod) %>% ungroup() %>% as.data.frame() %>% # probability all movement types result in freedom (non-entry)  - across municipalities
  rowwise() %>% 
  dplyr::mutate(pent = 1-prod(c_across(starts_with("X")))) %>% dplyr::select(-starts_with("X")) %>% 
  group_by(movementtype) %>% 
  summarise_all(funs(list(round(quantile(., probs = c(0.025, 0.5, 0.975)),5)))) %>% as.data.frame() %>% dplyr::select(-iteration) %>% 
  dplyr::mutate(pent_text=as.character(pent)) %>% dplyr::select(-pent)

write.csv(p.ent.overall.annual.movement, "./modeloutputs/tables/table3_movementspecific_annual.csv")

# Analysis 9.5: To get max of the median p_ent by local municipality and retrieve month in which that occurs
p.ent.overall.lm<-p.ent.overall.lm.agg.analyse  %>% 
  mutate_at(vars(contains("X")),function(x) (1-x)) %>% # probability of freedom for every data point
  dplyr::select(-c(movementtype)) %>% group_by(lmgid, iteration) %>%
  dplyr::summarise_all(prod) %>% # probability all movement types result in freedom (non-entry) - within municipalities
  mutate_at(vars(contains("X")),function(x) (1-x)) %>% # for every month (column) probability of entry through at least one movement type
  summarise_all(funs(quantile(., probs = c(0.5)))) %>% as.data.frame() %>% dplyr::select(-iteration) #summarize outcome with 95% CI

colnames(p.ent.overall.lm)<-c("lmgid",c(1:12))
p.ent.overall.lm.long <- tidyr::gather(p.ent.overall.lm, month, median, '1':'12', factor_key=TRUE) # make long dataset from wide dataset resulting in p.ent per month per lm
p.ent.overall.lm.long.max<-p.ent.overall.lm.long %>% group_by(lmgid) %>% 
  slice_max(median, n = 1, with_ties = FALSE) %>% 
  as.data.frame() #establish the maximum per lm

p.ent.overall.lm.long.max.spatial <- merge(lm_sf, p.ent.overall.lm.long.max, by.x = "gid", by.y = "lmgid") #merge with spatial dataset for mapping
st_write(p.ent.overall.lm.long.max.spatial, paste0("./modeloutputs/gis/Fig5_maxPent.shp"))

# Analysis 9.6 - to estabish CI for LM's for paper for LM which is the highest individual per month ####
# Max individual P_ent per month per LM
maxlm<-p.ent.overall.lm.long.max %>% slice_max(median, n = 1, with_ties = FALSE) %>% .$lmgid 

p.ent.overall.lm.CI<-p.ent.overall.lm.agg.analyse  %>% 
  mutate_at(vars(contains("X")),function(x) (1-x)) %>% # probability of freedom for every data point
  dplyr::select(-c(movementtype)) %>% group_by(lmgid, iteration) %>%
  dplyr::summarise_all(funs(prod)) %>% # probability all movement types result in freedom (non-entry) - within municipalities
  mutate_at(vars(contains("X")),function(x) (1-x)) %>% # for every month (column) probability of entry through at least one movement type
  summarise_all(funs(list(round(quantile(., probs = c(0.025, 0.5, 0.975)),5)))) %>% as.data.frame() %>% dplyr::select(-iteration)  %>% tidyr::pivot_longer(cols = starts_with('X'))#summarize outcome with 95% CI

View(p.ent.overall.lm.CI[which(p.ent.overall.lm.CI$lmgid==maxlm),])

# 10 Sensitivity analysis ####
#base inputs
percs.eval<-seq(.01, .99, length.out=3)

#10.1 - Standard movements ####
standard.matrix<-as.matrix(plyr::ldply(outputlist.sensitivity.standard, data.frame,.id = NULL))
standard.parameters=standard.matrix[,-(1)]  #Removes output for calculations - i.e. only parameters left
standard.output=standard.matrix[,1]  #Creates array with only output simulation results 

standard.quants<-apply(standard.parameters,2,function(x) quantile(x,percs.eval,na.rm=T))  #Calculate quantiles from outputs
#Conditional means of the output, for each input
standard.cond.means=standard.quants; standard.cond.means[,]=NA  #Empty dataframe
for(i in 1:ncol(standard.quants)){  #loops through each input variable
  input=standard.parameters[,i]
  for(j in 1:nrow(standard.quants)){  #Loops through each quantile
    if(j==1) {
      standard.cond.means[j,i]=mean(standard.output[input<=standard.quants[j,i]])} else {  #calculation for first quantile
        standard.cond.means[j,i]=mean(standard.output[input<=standard.quants[j,i] & input>standard.quants[j-1,i]])
        if(is.na(standard.cond.means[j,i])) {standard.cond.means[j,i]=standard.cond.means[j-1,i]}
      }  #if
  }  #j
}  #i

standard.size=apply(standard.cond.means,2,function(x) diff(range(x, na.rm=T)))  #Calculates individual variable range (impact)
standard.idx=order(abs(standard.size), decreasing=F)  #Index of input variables by magnitude (regardless of sign)
standard.cond.means.sort=standard.cond.means[,standard.idx]  #Sorted dataset by impact (width of range)
standard.cond.means=apply(standard.cond.means.sort,2, range) # get the range of conditional output result
standard.cond.means.df<-as.data.frame(standard.cond.means)
standard.cond.means.df$level<-c("ymin","ymax")
standard.cond.means.df.plot<- tidyr::gather(standard.cond.means.df, var, value, -level) %>% 
  pivot_wider(names_from = level, values_from = value)

width <- 0.5 # for width of bars in tornado plot

standard.cond.means.df.plot$plotorder <- seq.int(nrow(standard.cond.means.df.plot))
standard.cond.means.df.plot<-standard.cond.means.df.plot %>% 
  mutate(xmin=plotorder-width/2,
         xmax=plotorder+width/2)

ggplot() + 
  geom_rect(data = standard.cond.means.df.plot, 
            aes(ymax=ymax, ymin=ymin, xmax=xmax, xmin=xmin)) +
  scale_x_continuous(breaks = c(1:nrow(standard.cond.means.df.plot)), 
                     labels = unique(standard.cond.means.df.plot$var)) +
  xlab("Input parameters") +
  ylab("Parameter conditional mean outcome of probability of entry\n at local municipality level through standard movements") +
  theme_bw() + 
  theme(
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size = rel(0.5)),
    axis.text = element_text(size = rel(0.4)),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank()
  ) + 
  coord_flip()

ggsave("./modeloutputs/figures/Figure6.tiff",  width = 6.68, height = 6.68, units = "cm",dpi = 300) # size and units as req by PloS

#10.2 - Stop-over movements ####
soq.matrix<-as.matrix(plyr::ldply(outputlist.sensitivity.soq, data.frame,.id = NULL))
soq.parameters=soq.matrix[,-(1)]  #Removes output for calculations - i.e. only parameters left
soq.output=soq.matrix[,1]  #Creates array with only output simulation results 

soq.quants<-apply(soq.parameters,2,function(x) quantile(x,percs.eval,na.rm=T))  #Calculate quantiles from outputs
#Conditional means of the output, for each input
soq.cond.means=soq.quants; soq.cond.means[,]=NA  #Empty dataframe
for(i in 1:ncol(soq.quants)){  #loops through each input variable
  input=soq.parameters[,i]
  for(j in 1:nrow(soq.quants)){  #Loops through each quantile
    if(j==1) {
      soq.cond.means[j,i]=mean(soq.output[input<=soq.quants[j,i]])} else {  #calculation for first quantile
        soq.cond.means[j,i]=mean(soq.output[input<=soq.quants[j,i] & input>soq.quants[j-1,i]])
        if(is.na(soq.cond.means[j,i])) {soq.cond.means[j,i]=soq.cond.means[j-1,i]}
      }  #if
  }  #j
}  #i

soq.size=apply(soq.cond.means,2,function(x) diff(range(x, na.rm=T)))  #Calculates individual variable range (impact)
soq.idx=order(abs(soq.size), decreasing=F)  #Index of input variables by magnitude (regardless of sign)
soq.cond.means.sort=soq.cond.means[,soq.idx]  #Sorted dataset by impact (width of range)
soq.cond.means=apply(soq.cond.means.sort,2, range) # get the range of conditional output result
soq.cond.means.df<-as.data.frame(soq.cond.means)
soq.cond.means.df$level<-c("ymin","ymax")
soq.cond.means.df.plot<- tidyr::gather(soq.cond.means.df, var, value, -level) %>% 
  pivot_wider(names_from = level, values_from = value)

width <- 0.5 # for width of bars in tornado plot

soq.cond.means.df.plot$plotorder <- seq.int(nrow(soq.cond.means.df.plot))
soq.cond.means.df.plot<-soq.cond.means.df.plot %>% 
  mutate(xmin=plotorder-width/2,
         xmax=plotorder+width/2)

ggplot() + 
  geom_rect(data = soq.cond.means.df.plot, 
            aes(ymax=ymax, ymin=ymin, xmax=xmax, xmin=xmin)) +
  scale_x_continuous(breaks = c(1:nrow(soq.cond.means.df.plot)), 
                     labels = unique(soq.cond.means.df.plot$var)) +
  xlab("Input parameters") +
  ylab("Parameter conditional mean outcome of probability of entry\n at local municipality level through standard SOQ movements") +
  theme_bw() + 
  theme(
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size = rel(0.5)),
    axis.text = element_text(size = rel(0.4)),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank()
  ) + 
  coord_flip()

ggsave("./modeloutputs/figures/Figure7.tiff",  width = 6.68, height = 6.68, units = "cm",dpi = 300) # size and units as req by PloS


# 11 - What if analysis (no additional control)####
# 11.1 Monthly probability of entry ####
p.ent.overall.lm.agg.nocontrol.analyse<-
  plyr::ldply(p.ent.nocontrol, data.frame) %>% 
  group_by(.id) %>% 
  dplyr::mutate(iteration = row_number(), movementtype = "nocontrol")

p.ent.overall.lm.agg.nocontrol.analyse$lmgid<-p.ent.overall.lm.agg.nocontrol.analyse$.id #rename local municipality field
p.ent.overall.lm.agg.nocontrol.analyse<-as.data.frame(p.ent.overall.lm.agg.nocontrol.analyse) #ensure output is dataframe
p.ent.overall.lm.agg.nocontrol.analyse<-p.ent.overall.lm.agg.nocontrol.analyse %>% dplyr::select(-.id) #drop .id field

p.ent.nocontrol.analyse<-p.ent.overall.lm.agg.nocontrol.analyse  %>% 
  mutate_at(vars(contains("X")),function(x) (1-x)) %>% # probability of freedom for every data point
  dplyr::select(-c(movementtype,lmgid)) %>% group_by(iteration) %>% # start aggregating different movement types
  dplyr::summarise_all(prod) %>% # probability all movement types result in freedom (non-entry)  - across municipalities
  mutate_at(vars(contains("X")),function(x) (1-x)) %>% # for every month (column) probability of entry through at least one movement type
  sapply(., quantile, probs=c(0.025, 0.5, 0.975)) %>% as.data.frame() %>% mutate(estimate = c("2.5%","50%","97.5%")) %>% dplyr::select(-iteration) #summarize outcome with 95% CI

# reorganise data to month, CI range
p.ent.nocontrol.agg.analyse.all.output<-tidyr::gather(p.ent.nocontrol.analyse, month, pent, -estimate, 5) %>% 
  pivot_wider(names_from = estimate, values_from = pent) %>% as.data.frame() %>% 
  mutate_all(~gsub("X", "", .)) %>%  # drop X from month names
  sapply(., as.numeric) %>% as.data.frame() %>% round_df(., 5)

#final output
write.csv(p.ent.nocontrol.agg.analyse.all.output, "./modeloutputs/tables/table4_whatif_monthly_all_uncontrolled.csv")

#11.2 - P.ent annual aggregation ####
p.ent.nocontrol.aggannual.analyse.all.output<-p.ent.overall.lm.agg.nocontrol.analyse  %>% 
  mutate_at(vars(contains("X")),function(x) (1-x)) %>% # probability of freedom for every data point
  dplyr::select(-c(movementtype,lmgid)) %>% group_by(iteration) %>% # start aggregating across municipalities
  dplyr::summarise_all(prod) %>% ungroup() %>% as.data.frame() %>% # probability freedom across municipalities
  rowwise() %>% 
  dplyr::mutate(pent_annual = 1-prod(c_across(starts_with("X")))) %>% dplyr::select(-starts_with("X") ) %>% # for every month (column) probability of entry through at least one munic
  sapply(., quantile, probs=c(0.025, 0.5, 0.975)) %>% as.data.frame() %>% dplyr::select(-iteration) %>% 
  round_df(.,5)

write.csv(p.ent.nocontrol.aggannual.analyse.all.output, "./modeloutputs/tables/table4_whatif_annual_all_uncontrolled.csv")

#11.3 Risk differentials ####
p.ent.overall.nocontrol.annual.nonagg<-p.ent.overall.lm.agg.nocontrol.analyse  %>% 
  mutate_at(vars(contains("X")),function(x) (1-x)) %>% # probability of freedom for every data point
  dplyr::select(-c(movementtype,lmgid)) %>% group_by(iteration) %>% # start aggregating different movement types
  dplyr::summarise_all(prod) %>% ungroup() %>% as.data.frame() %>% # probability all movement types result in freedom (non-entry)  - across municipalities
  rowwise() %>% 
  dplyr::mutate(pent_nocontrol = 1-prod(c_across(starts_with("X")))) %>% dplyr::select(-starts_with("X"))

p.ent.overall.control.annual<-p.ent.overall.lm.agg.analyse  %>% 
  mutate_at(vars(contains("X")),function(x) (1-x)) %>% # probability of freedom for every data point
  dplyr::select(-c(movementtype,lmgid)) %>% group_by(iteration) %>% # start aggregating different movement types
  dplyr::summarise_all(prod) %>% ungroup() %>% as.data.frame() %>% # probability all movement types result in freedom (non-entry)  - across municipalities
  rowwise() %>% 
  dplyr::mutate(pent_control = 1-prod(c_across(starts_with("X")))) %>% dplyr::select(-starts_with("X") )

pent.riskdiffential<-cbind(p.ent.overall.nocontrol.annual.nonagg %>% dplyr::select(-iteration),p.ent.overall.control.annual)

# #overall risk reduction (factor and percentage change)
p.ent.reduction.annual<-
  pent.riskdiffential %>% 
  mutate(pent_red = as.numeric(pent_nocontrol/pent_control)) %>% 
  mutate(pent_red = as.numeric(pent_red)) %>% 
  mutate(pent_diff = as.numeric((pent_control-pent_nocontrol)/pent_nocontrol)) %>% 
  mutate(pent_diff = as.numeric(pent_diff)) %>% 
  dplyr::select(iteration,  pent_red, pent_diff) %>%
  dplyr::summarise(lower.x.factred = quantile(pent_red, probs = 0.025),
                   mean.x.factred = quantile(pent_red, probs = 0.5),
                   upper.x.factred = quantile(pent_red, probs = 0.975),
                   lower.x.diff = quantile(pent_diff, probs = 0.025),
                   mean.x.diff = quantile(pent_diff, probs = 0.5),
                   upper.x.diff = quantile(pent_diff, probs = 0.975)) %>% 
  as.data.frame()

#final output
write.csv(round_df(p.ent.reduction.annual,5), "./modeloutputs/tables/table4_whatif_annual_riskdiff.csv")

# monthly differential in risk

p.ent.monthly.nocontrol.full<-p.ent.overall.lm.agg.nocontrol.analyse  %>% 
  mutate_at(vars(contains("X")),function(x) (1-x)) %>% # probability of freedom for every data point
  dplyr::select(-c(movementtype,lmgid)) %>% group_by(iteration) %>% # start aggregating different movement types
  dplyr::summarise_all(prod) %>% # probability all movement types result in freedom (non-entry)  - across municipalities
  mutate_at(vars(contains("X")),function(x) (1-x)) %>% tidyr::gather(., month, pent_nocontrol, X1:X12, factor_key=TRUE) %>% 
  as.data.frame()

p.ent.monthly.control.full<-p.ent.overall.lm.agg.analyse  %>% 
  mutate_at(vars(contains("X")),function(x) (1-x)) %>% # probability of freedom for every data point
  dplyr::select(-c(movementtype,lmgid)) %>% group_by(iteration) %>% # start aggregating different movement types
  dplyr::summarise_all(prod) %>% # probability all movement types result in freedom (non-entry)  - across municipalities
  mutate_at(vars(contains("X")),function(x) (1-x)) %>% tidyr::gather(., month, pent_control, X1:X12, factor_key=TRUE) %>% 
  as.data.frame() %>% dplyr::select(pent_control)

#risk differential
p.ent.differential.month<-
  cbind(p.ent.monthly.nocontrol.full, p.ent.monthly.control.full) %>% 
  mutate(pent_diff = as.numeric((pent_control-pent_nocontrol)/pent_nocontrol)) %>% 
  mutate(pent_red = as.numeric(pent_nocontrol/pent_control)) %>% 
  mutate_all(~gsub("X", "", .)) %>% 
  mutate(month = as.numeric(month)) %>% 
  mutate(pent_diff = as.numeric(pent_diff)) %>% 
  mutate(pent_red = as.numeric(pent_red)) %>% 
  dplyr::select(month,  pent_diff, pent_red) %>% group_by(month) %>% 
  dplyr::summarise(lower.x.diff = quantile(pent_diff, probs = 0.025),
                   mean.x.diff = quantile(pent_diff, probs = 0.5),
                   upper.x.diff = quantile(pent_diff, probs = 0.975),
                   lower.x.factred = quantile(pent_red, probs = 0.025),
                   mean.x.factred = quantile(pent_red, probs = 0.5),
                   upper.x.factred = quantile(pent_red, probs = 0.975)) %>% 
  as.data.frame()
#final output
write.csv(round_df(p.ent.differential.month,5), "./modeloutputs/tables/table4_whatif_monthly_riskdiff.csv")

#Jan VPSOQ movements
dfmovement.analyse %>% filter(movementmonth == 1 & movementtype == 'vpsoq')
p.inf.summarytable %>% filter(lm == 160 & casemonth == 1)

# Analysis 11.4: NOT IN MAUSCRIPT: To get max of the median p_ent by local municipality and retrieve month in which that occurs ####
# Here because only one movement type (nocontrol) is present there is no need to aggregate across movement types
# Furthermore becasue the output is LM associated there is no need to aggregate across municipality

p.ent.overall.nocontrol.lm<-p.ent.overall.lm.agg.nocontrol.analyse  %>% 
  dplyr::select(-c(movementtype,iteration)) %>% group_by(lmgid) %>%
  summarise_all(funs(quantile(., probs = c(0.5)))) %>% as.data.frame() #summarize outcome with 95% CI

colnames(p.ent.overall.nocontrol.lm)<-c("lmgid",c(1:12))
p.ent.overall.nocontrol.lm.long <- tidyr::gather(p.ent.overall.lm, month, median, '1':'12', factor_key=TRUE) # make long dataset from wide dataset resulting in p.ent per month per lm
p.ent.overall.nocontrol.lm.long.max<-p.ent.overall.nocontrol.lm.long %>% 
  group_by(lmgid) %>% 
  slice_max(median, n = 1, with_ties = FALSE) %>% 
  as.data.frame() %>% #establish the maximum per lm
  arrange(as.numeric(lmgid))

p.ent.overall.nocontrol.lm.long.max.spatial <- merge(lm_sf, p.ent.overall.nocontrol.lm.long.max, by.x = "gid", by.y = "lmgid") #merge with spatial dataset for mapping
st_write(p.ent.overall.nocontrol.lm.long.max.spatial, paste0("./modeloutputs/gis/NIM_pent_nocontrol.shp"))