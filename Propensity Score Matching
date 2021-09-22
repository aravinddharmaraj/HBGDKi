cat("\f")###to clear console

# LOAD REQUIRED LIBRARIES
library(ggplot2)
library(reshape2)
library(plyr)
library(dplyr)
library(stringr)


rm(list = ls()) # CLEAR ENVIRONMENT

# SET WORK DIRECTORY
setwd("C:\\Aravind work file\\BF_data\\MERGING FILES\\data_clean") 



######## **************** LOAD DATASETS ********************##########
ki_main <- read.csv("main/data_clean/final_data/ki_data_updated_new_5_impute_growth.csv")###for 7 studies


ki <- ki_main
colnames(ki)
table(ki$weand_cat)
ki$weand_cat <- ifelse(ki$age_weand <= 134 & ki$age_weand >= 106, "4 months",
                       ifelse(ki$age_weand <= 194 & ki$age_weand >= 166, "6 months",                                
                              ifelse(ki$age_weand , "others",NA)))
table(ki$weand_cat)


## filter for 6 studies and exclude neovita study for growth variable
table(ki$study_name)
ki<- filter(ki, study_name %in%  c("cmc_bcs", "cmc_irc","cmc_tdc","imnci","vita","zn_sga",
                                   "maled"))


##FILTER 4 AND 6 MONTHS GROUP
ki <- filter(ki, weand_cat %in% c("4 months", "6 months"))
ki1= ki

hist(ki1$age_weand,breaks = 50)
table(ki1$weand_cat)

colnames(ki1)




# SELECT REQUIRED VARIABLES AND EXCLUDE OTHER VARIABLES
ki1<- select(ki, id, study_name,religion, house_type,family_type,sesgrp,cooking_mode,place_cooking,
             no_hh_cat,adult_cat,no_children_cat,no_livechild_cat,hh_occ_cat2,hhedu_cat,mom_age_cat,
             mot_edu_cat,gest_wk_cat,pl_deliv,deliv_mode,abortions,still_birth,sex,
             birth_weight_cat, singleton,weand_cat)
colnames(ki1)
str(ki1)


# CHANGE TO FACTOR 
ki1$weand_cat<- as.factor(ki1$weand_cat)
ki1$religion<- as.factor(ki1$religion)
ki1$house_type<- as.factor(ki1$house_type)
ki1$family_type<- as.factor(ki1$family_type)
ki1$sesgrp<- as.factor(ki1$sesgrp)
ki1$cooking_mode<- as.factor(ki1$cooking_mode)
ki1$place_cooking<- as.factor(ki1$place_cooking)
ki1$no_hh_cat<- as.factor(ki1$no_hh_cat)
ki1$adult_cat<- as.factor(ki1$adult_cat)
ki1$no_children_cat<- as.factor(ki1$no_children_cat)
ki1$no_livechild_cat<- as.factor(ki1$no_livechild_cat)
ki1$hh_occ_cat2<- as.factor(ki1$hh_occ_cat2)
ki1$hhedu_cat<- as.factor(ki1$hhedu_cat)
ki1$mom_age_cat<- as.factor(ki1$mom_age_cat)
ki1$mot_edu_cat<- as.factor(ki1$mot_edu_cat)
ki1$gest_wk_cat<- as.factor(ki1$gest_wk_cat)
ki1$pl_deliv<- as.factor(ki1$pl_deliv)
ki1$deliv_mode<- as.factor(ki1$deliv_mode)
ki1$abortions<- as.factor(ki1$abortions)
ki1$still_birth<- as.factor(ki1$still_birth)
ki1$sex<- as.factor(ki1$sex)
ki1$birth_weight_cat<- as.factor(ki1$birth_weight_cat)
ki1$singleton<-as.factor(ki1$singleton)


# RUN STEPWISE TO SELECT OPTIMAL MODEL WHICH INFLUENCING 4 AND 6 MO. GROUP
####stepwise regression starts
colnames(ki1)
ki2<- select(ki1, religion, house_type,family_type,sesgrp,cooking_mode,place_cooking,
             no_hh_cat,adult_cat,no_children_cat,no_livechild_cat, hh_occ_cat2,hhedu_cat,
             mom_age_cat,mot_edu_cat,gest_wk_cat,pl_deliv,deliv_mode,abortions,
             still_birth,sex,birth_weight_cat,singleton,weand_cat)

# FIT ALL THE VARIABLES
fitall <- glm(weand_cat ~ .,data=ki2, family = binomial(link = logit))
summary(fitall)


step(fitall,direction = "backward")## normal code 

# SAVE THE OPTIMAL MODEL



# RUN PRE-MATCHING ANALYSIS 
colnames(ki1)
str(ki1)
table(ki1$weand_cat)



# CHANGE 4 AND 6 MONTHS TO 0 AND 1
ki1$weand_cat<- as.character(ki1$weand_cat) # FIRST CHANGE TO CHARACTER
ki1$weand_cat[ki1$weand_cat == "4 months"]<- "0"
ki1$weand_cat[ki1$weand_cat == "6 months"]<- "1"
ki1$weand_cat<- as.factor(ki1$weand_cat) # AGAIN CHANGE TO FACTOR


# PUT THE OPTIMAL MODEL VARIABLES FOR PRE MATCHING ANALYSIS AND FOR PSM
ki_predict<- select(ki1,id,study_name,weand_cat,religion , sesgrp , cooking_mode , 
                    place_cooking , no_hh_cat , adult_cat , no_livechild_cat , 
                    hh_occ_cat2 ,hhedu_cat, mom_age_cat , mot_edu_cat , still_birth)



ki2 <- ki_predict

## check a model smd difference without PS MATCHING  
### pre-matching

psmodel<- glm(weand_cat ~ religion + sesgrp + cooking_mode + 
                place_cooking + no_hh_cat + adult_cat + no_livechild_cat + 
                hh_occ_cat2 + hhedu_cat + mom_age_cat + mot_edu_cat + still_birth, 
              family = binomial(link = logit), 
              data = ki_predict)



summary(psmodel)



xvars<- c("religion","sesgrp","cooking_mode",
          "place_cooking","no_hh_cat","adult_cat","no_livechild_cat",
          "hh_occ_cat2","hhedu_cat","mom_age_cat","mot_edu_cat","still_birth")


# CHECK THE FREQUENCY AND SIGNIFICANCE TEST
library(tableone)
matchedtab1<- CreateTableOne(vars = xvars,strata = "weand_cat",
                             data=ki_predict, test = FALSE)
pre<- data.frame(print(matchedtab1, smd = TRUE))
colnames(pre)[1]<- '4 months'
colnames(pre)[2]<-"6 months"
pre # SAVE THIS OUTPUT

# CHECK FOR SMD DIFFERENCE (STANDARD VALUE SHOULD BE LESS THAN 0.1)
addmargins(table(ExtractSmd(matchedtab1) > 0.1))


##########################################################################

############################ pre matching analysis end ####################

##########################################################################






###########################################################################

############################# PROPENSITY SCORE MATCHING ###################

###########################################################################
### REQUIRED LIBRARY FOR PSM
library(MatchIt)
library(optmatch)
library(rcbalance)
library(tableone)


ki_predict1<- ki_predict
colnames(ki_predict1)
#### fit a propensity score matching model based on optimal model (stepwise backward)


psmodel<- glm(weand_cat ~ religion + sesgrp + cooking_mode + 
                place_cooking + no_hh_cat + adult_cat + no_livechild_cat + 
                hh_occ_cat2 + hhedu_cat + mom_age_cat + mot_edu_cat + still_birth, 
              family = binomial(link = logit), 
              data = ki_predict)


## summary
summary(psmodel)

## create ps score and save with data frame
ki_predict$pscore<- psmodel$fitted.values
ki_predict1$pscore<- psmodel$fitted.values





######################################################

# DO SORTING FOR EXPOSURE AND PS SCORE
ki_sort<- ki_predict
ki_sort <- ki_sort[order(ki_sort[,"weand_cat"], ki_sort[,"pscore"]), ] 


# PSM MATCHING THROUGH MATCHIT FUNCTION
m.out <- matchit(weand_cat ~ religion + sesgrp + cooking_mode + 
                   place_cooking + no_hh_cat + adult_cat + no_livechild_cat + 
                   hh_occ_cat2 + hhedu_cat + mom_age_cat + mot_edu_cat + still_birth, 
                 data = ki_sort,distance = "glm",link="probit", caliper = .2,ratio = 1, 
                 method = "nearest", replace = F)


summary(m.out)

# OUTPUT TO CHECK THE MATCHED NUMBERS
dff <- match.data(m.out) 

table(dff$weand_cat)
plot(m.out, type = 'jitter', interactive = FALSE)


## RUN FREQUENCY AND SUMMARY STATISTICS AFTER MATCHING
xvars<- c("religion","sesgrp","cooking_mode",
          "place_cooking","no_hh_cat","adult_cat","no_livechild_cat",
          "hh_occ_cat2","hhedu_cat","mom_age_cat","mot_edu_cat","still_birth")

df <- dff

matchedtab<- CreateTableOne(vars = xvars,strata = "weand_cat",
                            data=df, test = FALSE)
print(matchedtab, smd = TRUE)


# CHECK FOR SMD LESS THAN 0.1 (HERE MAXIMUM VARIABLES ARE LESS THAN 0.1)
addmargins(table(ExtractSmd(matchedtab) > 0.1))
post<- data.frame(print(matchedtab, smd = TRUE))
colnames(post)[1]<- '6 months'
colnames(post)[2]<-"4 months"
post


#############################################################################

############################ PROPENSITY SCORE MATCHING END ##################

#############################################################################

