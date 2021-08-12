cat("\f")###to clear console

# LOAD REQUIRED LIBRARIES
library(ggplot2)
library(reshape2)
library(plyr)
library(dplyr)
library(stringr)
library(caret)
library(grid)


rm(list = ls()) # CLEAR ENVIRONMENT

# SET WORK DIRECTORY
setwd("C:\\Aravind work file\\BF_data\\MERGING FILES\\data_clean")



######## **************** LOAD DATASETS ********************##########
growth_vel <- read.csv("main/data_clean/final_data/growth_outcome_count/post_match_data_height_weight.csv")###for 7 studies

str(growth_vel)

length(unique(growth_vel$idno))


# CHANGE AGE TO MONTHS
growth_vel$age_m <- (growth_vel$age)/30


# VISUALIZE AGE
hist(growth_vel$age)
hist(growth_vel$age_m)

str(growth_vel)

# ROUND AGE WITHOUT INTEGERS
growth_vel$age_m1<- round(growth_vel$age_m,digits=0)

# FILTER AGE 0 TO 24
growth_vel$age_cat <- ifelse(growth_vel$age_m1 <= 24 , "less than 24",
                             ifelse(growth_vel$age_weand >= 25 , "others",NA))

# FILTER
growth_vel <- filter(growth_vel, age_cat %in% c("less than 24"))

# CHANGE TO FACTOR
growth_vel$age_m1<- as.factor(growth_vel$age_m1)
growth_vel$study<- as.factor(growth_vel$study)

table(growth_vel$age_m1)

str(growth_vel)


# FILTER FOR 4 AND 6 MONTHS (EXPOSURE GROUP)
growth_vel$weand_cat <- ifelse(growth_vel$age_weand <= 134 & growth_vel$age_weand >= 106, "4 months",
                               ifelse(growth_vel$age_weand <= 194 & growth_vel$age_weand >= 166, "6 months",                                
                                      ifelse(growth_vel$age_weand , "others",NA)))
table(growth_vel$weand_cat)

table(growth_vel$study)


##filter for multiple categories
growth_vel <- filter(growth_vel, weand_cat %in% c("4 months", "6 months"))

hist(growth_vel$age_weand,breaks = 100)

head(growth_vel)

O<-growth_vel

# keep only first entries based on ID and Age
O <- distinct(O, idno , age_m1 ,.keep_all = TRUE) 


## CHECK OUTLIERS AND REMOVE FOR HEIGHT AND WEIGHT 
boxplot(O$height)

height<- O[-which(O$height %in% boxplot.stats(O$height)$out),]
boxplot(height$height)
summary(height$height)
summary(O$height)

boxplot(O$weight)

weight<- O[-which(O$weight %in% boxplot.stats(O$weight)$out),]
boxplot(weight$weight)
summary(weight$weight)
summary(O$weight)

### removed outliers
## take children with minimum 5 serial measurement
head(height)
d1_sort<- select(height, idno, age_m1, height,weight) # SELECT REQUIRED VARIABLES
d1_sort <- d1_sort[order(d1_sort[,"idno"], d1_sort[,"age_m1"]), ]  # SORT
library(tidyverse)
d1_sort1<- d1_sort
d1_sort1<- d1_sort1 %>% 
        group_by(idno) %>%
        count()
d1_sort1$cat <- ifelse(d1_sort1$n <= 4 , "less than 4",
                       ifelse(d1_sort1$n >= 5 , "5 and above",NA))

table(d1_sort1$cat)
d1_sort1<- filter(d1_sort1, cat %in%  c("5 and above"))

library(dplyr)
d1_sort2 <- merge(d1_sort1,height, by="idno") # merge by id # left joint

length(unique(d1_sort2$idno)) # 243 correct!
table(d1_sort2$cat)

head(d1_sort2)
height<- select(d1_sort2, idno,study,age_m1,age, weight,height,weand_cat,sex)




##########################################################

########### calculating residual standard deviation ######

##########################################################

head(height)
str(height)
height$age_m1<- as.character(height$age_m1)
height$age_m1<- as.numeric(height$age_m1)

random_id <- height[order(height[,"idno"], height[,"age_m1"]), ] # SORT


# CHECK THE DATA ON HEIGHT AND WEIGHT FOR RANDOM INDIVIDUAL
random_id1<-height[height$idno %in% sample(unique(height$idno),1),]
table(random_id1$idno)

summary(random_id1)
str(random_id1)

plot(age_m1, height)





###########################################################################

## to extract residual standard deviation for each children ###############

###########################################################################

# LINEAR REGRESSION
model1<- by(random_id, random_id$idno, 
            function(data)sigma(lm(height~age_m1,data=data)))

rse<- model1

# sum of all the rsd
linear<- sum(rse)
linear
nrow(na.omit(model1))


# calculating percent error through rsd
rsd <- data.frame(matrix(unlist(rse), nrow=length(rse), byrow=TRUE))
colnames(rsd)[1]<- 'x'
m1<- aggregate(x=random_id$height,by=list(random_id$idno),FUN=mean )
m1<- select(m1, x)
mod<- rsd*100/m1$x
table(is.na(mod$x))
rse_p1<- (sum(mod,na.rm = T)/24300)*100
rse_p1


# QUADRATIC POLYNOMIAL REGRESSION
model2<- by(random_id, random_id$idno, 
            function(data)sigma(lm(height~age_m1+I(age_m1^2),data=data)))

rse<- model2
quadratic<- sum(rse)
quadratic
nrow(na.omit(model2))

rsd <- data.frame(matrix(unlist(rse), nrow=length(rse), byrow=TRUE))
colnames(rsd)[1]<- 'x'
m1<- aggregate(x=random_id$height,by=list(random_id$idno),FUN=mean )
m1<- select(m1, x)
mod<- rsd*100/m1$x
rse_p1<- (sum(mod,na.rm = T)/24300)*100
rse_p1


# CUBIC POLYNOMIAL REGRESSION
model3<- by(random_id, random_id$idno, 
            function(data)sigma(lm(height~age_m1+I(age_m1^2)+
                                           +I(age_m1^2),data=data)))

rse<- model3
cubic<- sum(rse)
cubic
nrow(na.omit(model3))

rsd <- data.frame(matrix(unlist(rse), nrow=length(rse), byrow=TRUE))
colnames(rsd)[1]<- 'x'
m1<- aggregate(x=random_id$height,by=list(random_id$idno),FUN=mean )
m1<- select(m1, x)
mod<- rsd*100/m1$x
rse_p1<- (sum(mod,na.rm = T)/24300)*100
rse_p1


# LOG LINEAR REGRESSION
model4<- by(random_id, random_id$idno, 
            function(data)sigma(lm(log(height)~age_m1,data=data)))

rse<- model4
loglinear<- sum(rse)
loglinear
nrow(na.omit(model4))

rsd <- data.frame(matrix(unlist(rse), nrow=length(rse), byrow=TRUE))
colnames(rsd)[1]<- 'x'
m1<- aggregate(x=random_id$height,by=list(random_id$idno),FUN=mean )
m1<- select(m1, x)
mod<- rsd*100/m1$x
rse_p1<- (sum(mod,na.rm = T)/24300)*100
rse_p1


# LINEAR LOG REGRESSION
# TO PERFORM LOG 0 MUST BE EXCLUDED
random_id$age_cat <- ifelse(random_id$age_m1 <= 0 , "less than 0",
                            ifelse(random_id$age_m1 >= 1 , "others",NA))
table(random_id$age_cat)
##filter for multiple categories
height1 <- filter(random_id, age_cat %in% c("others"))

model5<- by(height1, height1$idno, 
            function(data)sigma(lm(height~log(age_m1),data=data)))

rse<- model5
#rse<- exp(rse) 
linearlog<- sum(rse)
linearlog
nrow(na.omit(model5))
#exp(linearlog) 

rsd <- data.frame(matrix(unlist(rse), nrow=length(rse), byrow=TRUE))
colnames(rsd)[1]<- 'x'
m1<- aggregate(x=height1$height,by=list(height1$idno),FUN=mean )
m1<- select(m1, x)
mod<- rsd*100/m1$x
rse_p1<- (sum(mod,na.rm = T)/24300)*100
rse_p1



# LOG LOG REGRESSION
model6<- by(height1, height1$idno, 
            function(data)sigma(lm(log(height)~log(age_m1),data=data)))

rse<- exp(rse)
rse<- model6
loglog<- sum(rse)
loglog
nrow(na.omit(model6))

rsd <- data.frame(matrix(unlist(rse), nrow=length(rse), byrow=TRUE))
colnames(rsd)[1]<- 'x'
m1<- aggregate(x=height1$height,by=list(height1$idno),FUN=mean )
m1<- select(m1, x)
mod<- rsd*100/m1$x
rse_p1<- (sum(mod,na.rm = T)/24300)*100
rse_p1


# SPLINE REGRESSION (PIESEWISE)
library(splines)
knots <- quantile(random_id$age_m1, p = c(0.25, 0.5, 0.75))

model7<- by(random_id, random_id$idno, 
            function(data)sigma(lm(height ~ bs(age_m1, knots = knots),
                                   data=data)))

rse<- model7
spline<- sum(rse,na.rm = T)
spline
nrow(na.omit(model7))

rsd <- data.frame(matrix(unlist(rse), nrow=length(rse), byrow=TRUE))
colnames(rsd)[1]<- 'x'
m1<- aggregate(x=random_id$height,by=list(random_id$idno),FUN=mean )
m1<- select(m1, x)
mod<- rsd*100/m1$x
rse_p1<- (sum(mod, na.rm = T)/24200)*100
rse_p1


# FRACTIONAL POLYNOMIAL
model8<- by(random_id, random_id$idno, 
            function(data)sigma(lm(log(height)~ I(age_m1^2)
                                   + I(age_m1^(1/2)),
                                   data=data)))

rse<- model8
fractional<- sum(rse,na.rm = T)
fractional
nrow(na.omit(model8))

rsd <- data.frame(matrix(unlist(rse), nrow=length(rse), byrow=TRUE))
colnames(rsd)[1]<- 'x'
m1<- aggregate(x=random_id$height,by=list(random_id$idno),FUN=mean )
m1<- select(m1, x)
mod<- rsd*100/m1$x
rse_p1<- (sum(mod,na.rm = T)/24300)*100
rse_p1


# FRACTIONAL POLYNOMIAL
model9<- by(height1, height1$idno, 
            function(data)sigma(lm(height ~ I(age_m1^(2))+log(age_m1)
                                   +I(age_m1^(2)),
                                   data=data)))

mod<- lm(weight ~ I(age_m1^2)+ I(age_m1^1/2),data=random_id)
summary(mod)

rse<- model9
fractional<- sum(rse,na.rm = T)
fractional
nrow(na.omit(model9))

rsd <- data.frame(matrix(unlist(rse), nrow=length(rse), byrow=TRUE))
colnames(rsd)[1]<- 'x'
m1<- aggregate(x=random_id$height,by=list(random_id$idno),FUN=mean )
m1<- select(m1, x)
mod<- rsd*100/m1$x
rse_p1<- (sum(mod,na.rm = T)/24300)*100
rse_p1



############################################################################

######## EXTRACT AVERAGE DISTANCE FOR AVG VS INDIVIDUAL PREDICTED VALUE ####

############################################################################

##################### spline regression ####################################

############################################################################
library(splines)
knots <- quantile(random_id$age_m1, p = c(0.25, 0.5, 0.75))



####################  spline regression model individual ########

####################### prediction value extraction ##############

###################################################################

# REMOVE NA'S FOR HEIGHT
random_id <- random_id[which(random_id$height != "NA"),]
model7<- by(random_id, random_id$idno, 
            function(data)fitted.values(lm(height ~ bs(age_m1, knots = knots),
                                           data=data)))
spline_fitted<- data.frame(t(sapply(model7, function(x) x[1:max(lengths(model7))])))

head(spline_fitted)


spline <- spline_fitted

head(spline)

sp<- melt(spline)
head(sp)
table(is.na(sp$value))

d <- sp[which(sp$value != "NA"),]
d1<- select(d,idno, value)

# select only id and age
df<- select(random_id, idno, age_m1)
df1<- select(df,idno,age_m1)


length(unique(d1$idno))
length(unique(df1$idno))

d1 <- d1[order(d1[,"idno"], d1[,"value"]), ] 
df1 <- df1[order(df1[,"idno"],df1[,"age_m1"]), ] 
head(df1)

d1_sum<- d1 %>% count(idno)
d1_sum$cat<- rep("predicted")

df1_sum<- df1 %>% count(idno)
df1_sum$cat<- rep("overall")

str(df1_sum)
str(d1_sum)
o <- rbind(df1_sum, d1_sum)
o1<- merge(df1_sum,d1_sum, by="idno")
str(o1)
# predicted minus overall
o1$dif<- (o1$n.x - o1$n.y)

sum(o1$dif)

str(df1)
str(d1)          

#write.csv(df1,"main/results/output/ps_output/new/observed_ido.csv")
#write.csv(d1,"main/results/output/ps_output/new/fitted.spline.csv")


mod<- lm(height ~ bs(age_m1, knots = knots),
         data=random_id)
random_id$fit<- predict(mod, random_id)

random_id2<- data.frame(unique(random_id$fit))
random_id2 <- random_id2[order(random_id2[,"unique.random_id.fit."]), ]


#write.csv(random_id2,"main/results/output/ps_output/new/spline_fitted.avg.csv")



## calculating distance for individual vs overall average ####

##############################################################
spline <- read.csv("main/results/output/ps_output/new/spline/fitted.spline1.csv")

head(spline)

spline$dis<- spline$avg_predict-spline$ind_predict


m1<- aggregate(x=spline$dis,by=list(spline$age_m1),FUN=sum )

#write.csv(m1,"main/results/output/ps_output/new/spline/spline_fit_bymonth.csv")


sum(m1$x)


########################################################



#########################################################

###################### linear regression ###############

#########################################################

model7<- by(random_id, random_id$idno, 
            function(data)fitted.values(lm(height ~ age_m1,
                                           data=data)))
spline_fitted<- data.frame(t(sapply(model7, function(x) x[1:max(lengths(model7))])))

head(spline_fitted)

spline <- spline_fitted
head(spline)

sp<- melt(spline)
head(sp)
table(is.na(sp$value))

d <- sp[which(sp$value != "NA"),]
d1<- select(d,idno, value)

# select only id and age
df<- select(random_id, idno, age_m1)
df1<- select(df,idno,age_m1)


length(unique(d1$idno))
length(unique(df1$idno))

d1 <- d1[order(d1[,"idno"], d1[,"value"]), ] 
df1 <- df1[order(df1[,"idno"],df1[,"age_m1"]), ] 
head(df1)

d1_sum<- d1 %>% count(idno)
d1_sum$cat<- rep("predicted")

df1_sum<- df1 %>% count(idno)
df1_sum$cat<- rep("overall")

str(df1_sum)
str(d1_sum)
o <- rbind(df1_sum, d1_sum)
o1<- merge(df1_sum,d1_sum, by="idno")
str(o1)
# predicted minus overall
o1$dif<- (o1$n.x - o1$n.y)

sum(o1$dif)

str(df1)
str(d1)          

write.csv(df1,"main/results/output/ps_output/new/observed_ido.csv")
write.csv(d1,"main/results/output/ps_output/new/fitted.spline.csv")


mod<- lm(height ~ age_m1,
         data=random_id)
random_id$fit<- predict(mod, random_id)

random_id2<- data.frame(unique(random_id$fit))
#random_id2 <- random_id2[order(random_id2[,"unique.random_id.fit."]), ]


write.csv(random_id2,"main/results/output/ps_output/new/spline_fitted.avg.csv")


## calculating distance for individual vs overall average ####

##############################################################
spline <- read.csv("main/results/output/ps_output/new/spline/fitted.spline2.csv")

head(spline)

spline$dis<- spline$avg_predict-spline$ind_predict

growth_vel <- filter(spline, age_m1 %in% c("8"))
sum(growth_vel$dis)

m1<- aggregate(x=spline$dis,by=list(spline$age_m1),FUN=sum )

write.csv(m1,"main/results/output/ps_output/new/spline/linear_fit_bymonth.csv")


sum(m1$x)




########################################################

#########################################################




#########################################################

###################### quadratic regression ###############

#########################################################
random_id <- random_id[which(random_id$height != "NA"),]
model7<- by(random_id, random_id$idno, 
            function(data)fitted.values(lm(height ~ I(age_m1)+I(age_m1^2),
                                           data=data)))
spline_fitted<- data.frame(t(sapply(model7, function(x) x[1:max(lengths(model7))])))

head(spline_fitted)

spline<- spline_fitted

head(spline)

sp<- melt(spline)
head(sp)
table(is.na(sp$value))

d <- sp[which(sp$value != "NA"),]
d1<- select(d,idno, value)

# select only id and age
df<- select(random_id, idno, age_m1)
df1<- select(df,idno,age_m1)


length(unique(d1$idno))
length(unique(df1$idno))

d1 <- d1[order(d1[,"idno"], d1[,"value"]), ] 
df1 <- df1[order(df1[,"idno"],df1[,"age_m1"]), ] 
head(df1)

d1_sum<- d1 %>% count(idno)
d1_sum$cat<- rep("predicted")

df1_sum<- df1 %>% count(idno)
df1_sum$cat<- rep("overall")

str(df1_sum)
str(d1_sum)
o <- rbind(df1_sum, d1_sum)
o1<- merge(df1_sum,d1_sum, by="idno")
str(o1)
# predicted minus overall
o1$dif<- (o1$n.x - o1$n.y)

sum(o1$dif)

str(df1)
str(d1)          

write.csv(df1,"main/results/output/ps_output/new/spline/observed_ido.csv")
write.csv(d1,"main/results/output/ps_output/new/spline/fitted.spline.csv")


mod<- lm(height ~ age_m1+I(age_m1^2),
         data=random_id)
summary(mod)
random_id$fit<- predict(mod, random_id)

random_id2<- data.frame(unique(random_id$fit))
#random_id2 <- random_id2[order(random_id2[,"unique.random_id.fit."]), ]


write.csv(random_id2,"main/results/output/ps_output/new/spline/spline_fitted.avg.csv")


## calculating distance for individual vs overall average ####

##############################################################
spline <- read.csv("main/results/output/ps_output/new/spline/fitted.spline3.csv")

head(spline)

spline$dis<- spline$avg_predict-spline$ind_predict

growth_vel <- filter(spline, age_m1 %in% c("8"))
sum(growth_vel$dis)

m1<- aggregate(x=spline$dis,by=list(spline$age_m1),FUN=sum )

write.csv(m1,"main/results/output/ps_output/new/spline/quadratic_fit_bymonth.csv")


sum(m1$x)



########################################################

#########################################################


#########################################################

###################### fractional regression ############

#########################################################
random_id <- random_id[which(random_id$height != "NA"),]
random_id$age_cat <- ifelse(random_id$age_m1 <= 0 , "less than 0",
                            ifelse(random_id$age_m1 >= 1 , "others",NA))
table(random_id$age_cat)
##filter for multiple categories
height1 <- filter(random_id, age_cat %in% c("others"))



model7<- by(height1, height1$idno, 
            function(data)fitted.values(lm(height ~ age_m1+log(age_m1)+I(age_m1^-2),
                                           data=data)))
spline_fitted<- data.frame(t(sapply(model7, function(x) x[1:max(lengths(model7))])))

head(spline_fitted)

spline <- spline_fitted
head(spline)
colnames(spline)[1]<- "idno"

sp<- melt(spline)
head(sp)
table(is.na(sp$value))

d <- sp[which(sp$value != "NA"),]
d1<- select(d,idno, value)

# select only id and age
df<- select(height1, idno, age_m1)
df1<- select(df,idno,age_m1)

#table(df1$age_m1)
#df1<- df1[!(df1$age_m1=="0"),] delete any value remove 0
#table(df1$age_m1)


length(unique(d1$idno))
length(unique(df1$idno))

d1 <- d1[order(d1[,"idno"], d1[,"value"]), ] 
df1 <- df1[order(df1[,"idno"],df1[,"age_m1"]), ] 
head(df1)

d1_sum<- d1 %>% count(idno)
d1_sum$cat<- rep("predicted")

df1_sum<- df1 %>% count(idno)
df1_sum$cat<- rep("overall")

str(df1_sum)
str(d1_sum)
o <- rbind(df1_sum, d1_sum)
o1<- merge(df1_sum,d1_sum, by="idno")
str(o1)



# predicted minus overall
o1$dif<- (o1$n.x - o1$n.y)

sum(o1$dif)

str(df1)
str(d1)          

write.csv(df1,"main/results/output/ps_output/new/spline/observed_ido.csv")
write.csv(d1,"main/results/output/ps_output/new/spline/fitted.spline.csv")


mod<- lm(height ~ age_m1+log(age_m1)+I(age_m1^-2),
         data=height1)
summary(mod)
height1$fit<- predict(mod, height1)

random_id2<- data.frame(unique(height1$fit))
#random_id2 <- random_id2[order(random_id2[,"unique.random_id.fit."]), ]


write.csv(random_id2,"main/results/output/ps_output/new/spline/spline_fitted.avg.csv")


## calculating distance for individual vs overall average ####

##############################################################
spline <- read.csv("main/results/output/ps_output/new/spline/fitted.spline4.csv")

head(spline)

spline$dis<- spline$avg_predict-spline$ind_predict

growth_vel <- filter(spline, age_m1 %in% c("8"))
sum(growth_vel$dis)

m1<- aggregate(x=spline$dis,by=list(spline$age_m1),FUN=sum )

write.csv(m1,"main/results/output/ps_output/new/spline/fractional_fit_bymonth.csv")


sum(m1$x)


########################################################

#########################################################






