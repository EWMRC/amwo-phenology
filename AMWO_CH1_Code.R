#load libraries
library(lme4)
library(AICcmodavg)
library(ggplot2)
library(wesanderson)
library(RColorBrewer)
library(tidyverse)
library(ggridges)



####################################################################################
##      Migration - Fall Initiation Analysis
####################################################################################

##upload migration initiation dataset and format for analysis
data19 <- read.csv("amwo.init.f19.csv")
data19$mig.init.between <- as.POSIXct(data19$mig.init.between, origin="1970-01-01 04:00:00")
data19$init.ord <- as.numeric(round(data19$mig.init.between - as.POSIXct("2019-10-01", origin="1970-01-01 04:00:00")))
data18 <- read.csv("amwo.init.f18.csv")
data18$mig.init.between <- as.POSIXct(data18$mig.init.between, origin="1970-01-01 04:00:00")
data18$init.ord <- as.numeric(round(data18$mig.init.between - as.POSIXct("2018-10-01", origin="1970-01-01 04:00:00")))
data17 <- read.csv("amwo.init.f17.csv")
data17$mig.init.between <- as.POSIXct(data17$mig.init.between, origin="1970-01-01 04:00:00")
data17$init.ord <- as.numeric(round(data17$mig.init.between - as.POSIXct("2017-10-01", origin="1970-01-01 04:00:00")))
data <- rbind(data18, data19, data17)
amwo.capt <- read.csv("G:/My Drive/PhD_AMWO/AMWo_Phenology/AMWO_Condition_Regression/amwo.capt_condition.csv")
#pre-processing data (indicator variables for age, sex)
data$sex <- as.numeric(if_else(data$sex == 'm', 0, 1))  #males=0; female=1
data$age <- as.numeric(if_else(data$age == 'juv', 0, 1))  #juv=o; ad=1
##creating condition indicie
amwo.capt$condition <- amwo.capt$resid
amwo.capt$ID <- amwo.capt$Movebank.ID
data <- left_join(data, amwo.capt[c(3,24:25)], by = "ID") 
#fall 2019 (remove duplicates)
data <- data[-c(40,54),]  # remove duplicate data from re-captured birds
#fall all initiation removal
data <- data[-c(35,39,87,103),] 
## done with data prep

##A priori Geospatial model... determining support for lat, long, st.start (Geospatial)
lat <- glm(init.ord ~ lat,data = data)
lon <- glm(init.ord ~ lon,data = data)
lat.plus.lon <- glm(init.ord ~ lat + lon,data = data)
lat.times.lon <- glm(init.ord ~ lat*lon,data = data)
state <- glm(init.ord ~ as.factor(st.start),data = data)
null <- glm(init.ord ~ 1,data = data)

cand.set <- list(lat, lon, lat.plus.lon, lat.times.lon, state, null)
modnames <- c("lat", "lon", "lat.plus.lon", "lat.times.lon", "state", "null")
aictab(cand.set=cand.set, modnames = modnames,  digits = 2, second.ord = T)
summary(lat.times.lon)$coefficients
summary(lat.plus.lon)
#lat.times.lon model has the greaters model support with 0.66 model weight

##A priori Demographic model... determing support for age, sex (Demographic)
age <- glm(init.ord ~ age + lat + lon,data = data)
sex <- glm(init.ord ~ sex + lat + lon,data = data)
age.plus.sex <- glm(init.ord ~ age + sex + lat + lon,data = data)
age.times.sex <- glm(init.ord ~ agesex + lat + lon,data = data)
null <- glm(init.ord ~ lat + lon,data = data)

cand.set <- list(age, sex, age.plus.sex, age.times.sex, null)
modnames <- c("age", "sex", "age.plus.sex", "age.times.sex", "null")
aictab(cand.set=cand.set, modnames = modnames,  digits = 2, second.ord = T)
summary(age)
summary(age.times.sex)$coefficients
#age model received the most support with 0.47 cumulative support

##A priori Biological model.... determing support for condition (Biological)
####  removing RI bird with NA for condition (subsetting data)
data2 <- subset(data,
                !(is.na(condition)))

null <- glm(init.ord ~ age + lat + lon,data = data2)
cond <- glm(init.ord ~ condition + age + lat + lon,data = data2)
cond.times.sex <- glm(init.ord ~ condition*sex + age + lat + lon,data = data2)
cond.times.age <- glm(init.ord ~ condition*age + age + lat + lon,data = data2)
cond.times.lat <- glm(init.ord ~ condition*lat + age + lat + lon,data = data2) 

cand.set <- list(null, cond, cond.times.sex, cond.times.age, cond.times.lat)
modnames <- c("null", "cond", "cond.times.sex", "cond.times.age", "cond.times.lat")
aictab(cand.set=cand.set, modnames = modnames,  digits = 2, second.ord = T)
summary(cond.times.age)$coefficients
summary(cond)$coefficients
summary(null)$coefficients
# null no non-significant parameters and had 0.22 cumulative weight

##Predicting initiation of fall Migration
### Model predictions
### using the inference model to predict distributions and calculate 95% CI
f.init.pred <- read.csv("fall_init_predict.csv")
prediction<- predict.glm(age, newdata=f.init.pred, se.fit=TRUE)
fall.pred<- cbind(f.init.pred, prediction)
colnames(fall.pred)<- c("st.start", "lon", "lat","age","init.pred","SE")
fall.pred$lower<- fall.pred$init.pred-(fall.pred$SE*1.96)
fall.pred$upper<- fall.pred$init.pred+(fall.pred$SE*1.96)
pred.adult <- subset(fall.pred, age=="1")
pred.juv <- subset(fall.pred, age=="0")



####################################################################################
##      Migration - Fall Terminiation Analysis
####################################################################################


##upload migration termination dataset and format for analysis
data19 <- read.csv("amwo.term.f19.csv")
data19$mig.term.between <- as.POSIXct(data19$mig.term.between, origin="1970-01-01 04:00:00")
data19$term.ord <- as.numeric(round(data19$mig.term.between - as.POSIXct("2019-10-01", origin="1970-01-01 04:00:00")))
data18 <- read.csv("amwo.term.f18.csv")
data18$mig.term.between <- as.POSIXct(data18$mig.term.between, origin="1970-01-01 04:00:00")
data18$term.ord <- as.numeric(round(data18$mig.term.between - as.POSIXct("2018-10-01", origin="1970-01-01 04:00:00")))
data17 <- read.csv("amwo.term.f17.csv")
data17$mig.term.between <- as.POSIXct(data17$mig.term.between, origin="1970-01-01 04:00:00")
data17$term.ord <- as.numeric(round(data17$mig.term.between - as.POSIXct("2017-10-01", origin="1970-01-01 04:00:00")))
data <- rbind(data18, data19, data17)
amwo.capt <- read.csv("G:/My Drive/PhD_AMWO/AMWo_Phenology/AMWO_Condition_Regression/amwo.capt_condition.csv")
#pre-processing data (indicator variables for age, sex) 
data$sex <- as.numeric(if_else(data$sex == 'm', 0, 1))  #males=0; female=1
data$age <- as.numeric(if_else(data$age == 'juv', 0, 1))  #ad=o; juv=1
##creating condition indicie
amwo.capt$condition <- amwo.capt$resid
amwo.capt$ID <- amwo.capt$Movebank.ID
data <- left_join(data, amwo.capt[c(3,24:25)], by = "ID") 

#importing initiation data to create column that has starting migration latitude
data192 <- read.csv("amwo.init.f19.csv")
data182 <- read.csv("amwo.init.f18.csv")
data172 <- read.csv("amwo.init.f17.csv")
data2 <- rbind(data182, data192, data172)
data2 <- data2 %>%
  rename(start.lat = lat,
         start.lon = lon)
data <- left_join(data, data2[c(2,6:7)], by = "ID") 
rm(data172, data182, data192, data2)
#fall 2019 (remove duplicates)
data <- data[-c(27:29, 33, 74:76, 81),]  # remove duplicate data from re-captured birds
## done with data prep


##A priori Geospatial model... determining support for lat, long, st.start (Geospatial)
      ## Note: lat long for both termination and initiation location are included in analysis
lat <- glm(term.ord ~ lat,data = data)
lat.start <- glm(term.ord ~ start.lat,data = data)
lon <- glm(term.ord ~ lon,data = data)
lon.start <- glm(term.ord ~ start.lon,data = data)
lat.plus.lon <- glm(term.ord ~ lat + lon,data = data)
lat.plus.lon.start <- glm(term.ord ~ start.lat + start.lon,data = data)
lat.times.lon <- glm(term.ord ~ lat*lon,data = data)
lat.times.lon.start <- glm(term.ord ~ start.lat*start.lon,data = data)
state.start <- glm(term.ord ~ as.factor(st.start),data = data)
state.stop <- glm(term.ord ~ as.factor(st.pr),data = data)
null <- glm(term.ord ~ 1,data = data)

cand.set <- list(lat, lat.start, lon, lon.start, lat.plus.lon, lat.plus.lon.start, lat.times.lon, lat.times.lon.start, state.start, state.stop, null)
modnames <- c("lat", "lat.start", "lon","start.lon", "lat.plus.lon","lat.plus.lon.start", "lat.times.lon","lat.times.lon.start", 
              "state.start", "state.stop", "null")
aictab(cand.set=cand.set, modnames = modnames,  digits = 2, second.ord = T)
summary(lat.times.lon.start)$coefficients
summary(lat.plus.lon.start)$coefficients
#lat.plus.lon.start model has no non-significant parameters and has 0.41 model weight; choosen over top model

##A priori Demographic model... determing support for age, sex (Demographic)
age <- glm(term.ord ~ age + start.lat + start.lon,data = data)
sex <- glm(term.ord ~ sex + start.lat + start.lon,data = data)
age.plus.sex <- glm(term.ord ~ age + sex + start.lat + start.lon,data = data)
age.times.sex <- glm(term.ord ~ agesex + start.lat + start.lon,data = data)
null <- glm(term.ord ~ start.lat + start.lon,data = data)

cand.set <- list(age, sex, age.plus.sex, age.times.sex, null)
modnames <- c("age", "sex", "age.plus.sex", "age.times.sex", "null")
aictab(cand.set=cand.set, modnames = modnames,  digits = 2, second.ord = T)
summary(age)$coefficients
#null model contained no non-significant parameters and had 0.30 support

##A priori Biological model.... determing support for condition (Biological)
####  removing RI bird with NA for condition (subsetting data)
data3 <- subset(data,
                !(is.na(condition)))

null <- glm(term.ord ~ start.lat + start.lon,data = data3)
cond <- glm(term.ord ~ condition + start.lat + start.lon,data = data3)
cond.times.sex <- glm(term.ord ~ condition*sex + start.lat + start.lon,data = data3)
cond.times.age <- glm(term.ord ~ condition*age + start.lat + start.lon,data = data3)
cond.times.lat <- glm(term.ord ~ condition*lat + start.lat + start.lon,data = data3)
cond.times.startlat <- glm(term.ord ~ condition*start.lat + start.lat + start.lon,data = data3)

cand.set <- list(null, cond, cond.times.sex, cond.times.age, cond.times.lat, cond.times.startlat)
modnames <- c("null", "cond", "cond.times.sex", "cond.times.age", "cond.times.lat", "cond.times.startlat")
aictab(cand.set=cand.set, modnames = modnames,  digits = 2, second.ord = T)
summary(cond)$coefficients
summary(null)$coefficients
## null model has the greatest model support with 0.27 sumulative weight

##No predication model for fall termination



####################################################################################
##      Migration - Spring Initiation Analysis
####################################################################################


##upload migration termination dataset and format for analysis
data19 <- read.csv("amwo.init.sp19.csv")
data19$mig.init.between <- as.POSIXct(data19$mig.init.between, origin="1970-01-01 04:00:00")
data19$init.ord <- as.numeric(round(data19$mig.init.between - as.POSIXct("2019-01-01", origin="1970-01-01 04:00:00")))
data20 <- read.csv("amwo.init.sp20.csv")
data20$mig.init.between <- as.POSIXct(data20$mig.init.between, origin="1970-01-01 04:00:00")
data20$init.ord <- as.numeric(round(data20$mig.init.between - as.POSIXct("2020-01-01", origin="1970-01-01 04:00:00")))
data <- rbind(data20, data19)
amwo.capt <- read.csv("amwo.capt_condition.csv")
amwo.capt <- read.csv("G:/My Drive/PhD_AMWO/AMWo_Phenology/AMWO_Condition_Regression/amwo.capt_condition.csv")
#pre-processing data (indicator variables for age, sex) 
data$sex <- as.numeric(if_else(data$sex == 'm', 0, 1))  #males=0; female=1
data$age <- as.numeric(if_else(data$age == 'juv', 0, 1))  #ad=o; juv=1
##creating condition indicie
amwo.capt$condition <- amwo.capt$resid
amwo.capt$ID <- amwo.capt$Movebank.ID
data <- left_join(data, amwo.capt[c(3,24:25)], by = "ID") 
#fall 2019 (remove duplicates)
data <- data[-c(57, 100, 111, 113, 114),]  # remove duplicate data from re-captured birds
#removing birds with no state listed
data <- subset(data,
               !(is.na(st.start)))
## done with data prep

##A priori Geospatial model... determining support for lat, long, st.start (Geospatial)
lat <- glm(init.ord ~ lat,data = data)
lon <- glm(init.ord ~ lon,data = data)
lat.plus.lon <- glm(init.ord ~ lat + lon,data = data)
lat.times.lon <- glm(init.ord ~ lat*lon,data = data)
state.start <- glm(init.ord ~ as.factor(st.start),data = data)
null <- glm(init.ord ~ 1,data = data)

cand.set <- list(lat, lon, lat.plus.lon, lat.times.lon, state.start, null)
modnames <- c("lat", "lon", "lat.plus.lon", "lat.times.lon", "state.start", "null")
aictab(cand.set=cand.set, modnames = modnames,  digits = 2, second.ord = T)
summary(lon)$coefficients
#lon model has the greaters model support with 0.55 model weight

##A priori Demographic model... determing support for age, sex (Demographic)
data <- data[-c(64, 108),]  # remove birds with no age sex data

age <- glm(init.ord ~ age + lon,data = data)
sex <- glm(init.ord ~ sex + lon,data = data)
age.plus.sex <- glm(init.ord ~ age + sex + lon,data = data)
age.times.sex <- glm(init.ord ~ agesex + lon,data = data)
null <- glm(init.ord ~ lon,data = data)

cand.set <- list(age, sex, age.plus.sex, age.times.sex, null)
modnames <- c("age", "sex", "age.plus.sex", "age.times.sex", "null")
aictab(cand.set=cand.set, modnames = modnames,  digits = 2, second.ord = T)
summary(sex)$coefficients
#sex model received the most support with 0.53 model weight

##A priori Biological model.... determing support for condition (Biological)
####  removing RI bird with NA for condition
data3 <- subset(data,
                !(is.na(condition)))
data3 <- data3[-c(30:33, 34:42, 50:58, 72, 77:86, 95),]  # remove AMWO marked in the fall, as conditon score is not longer acurate [NJ AWO 30:35 and 85:89]

null <- glm(init.ord ~ sex + lon,data = data3)
cond <- glm(init.ord ~ condition + sex + lon,data = data3)
cond.times.sex <- glm(init.ord ~ condition*sex + sex + lon,data = data3)
cond.times.age <- glm(init.ord ~ condition*age + sex + lon,data = data3)
cond.times.lat <- glm(init.ord ~ condition*lat + sex + lon,data = data3)
cond.times.lon <- glm(init.ord ~ condition*lon + sex + lon,data = data3)

cand.set <- list(null, cond, cond.times.sex, cond.times.age, cond.times.lat, cond.times.lon)
modnames <- c("null", "cond", "cond.times.sex", "cond.times.age", "cond.times.lat", "cond.times.lon")
aictab(cand.set=cand.set, modnames = modnames,  digits = 2, second.ord = T)
summary(cond.times.sex)$coefficients
summary(cond.times.sex2)
## cond.times.sex was the highest supported model with 0.90 cumulative weight


### predicting initiation date in spring
### using condition.sex as inference model
### condition used was mean=1.08 
sp.init.pred <- read.csv("sp.init.predict.csv")
prediction<- predict.glm(cond.times.sex, newdata=sp.init.pred, se.fit=TRUE)
sp.pred<- cbind(sp.init.pred, prediction)
colnames(sp.pred)<- c("st.start", "lon", "lat","condition","sex","init.pred","SE")
sp.pred$lower<- sp.pred$init.pred-(sp.pred$SE*1.96)
sp.pred$upper<- sp.pred$init.pred+(sp.pred$SE*1.96)
pred.male <- subset(sp.pred, sex=="0")
pred.female <- subset(sp.pred, sex=="1")



####################################################################################
##      Migration - Timing of Fall Migration Analysis
####################################################################################

##upload migration termination dataset and format for analysis
# initation data for 2019 and 2018; also included is calcuation for ordinal data
data19 <- read.csv("amwo.mig.f19.csv")
data19$date <- as.POSIXct(data19$date, origin="1970-01-01 04:00:00")
data19$date.ord <- as.numeric(round(data19$date - as.POSIXct("2019-10-01", origin="1970-01-01 04:00:00")))
data19$time <- as.POSIXct(data19$time, origin="1970-01-01 04:00:00")
data19$time.ord <- as.numeric(round(data19$time - as.POSIXct("2019-10-01", origin="1970-01-01 04:00:00")))
data18 <- read.csv("amwo.mig.f18.csv")
data18$date <- as.POSIXct(data18$date, origin="1970-01-01 04:00:00")
data18$date.ord <- as.numeric(round(data18$date - as.POSIXct("2018-10-01", origin="1970-01-01 04:00:00")))
data18$time <- as.POSIXct(data18$time, origin="1970-01-01 04:00:00")
data18$time.ord <- as.numeric(round(data18$time - as.POSIXct("2018-10-01", origin="1970-01-01 04:00:00")))
data17 <- read.csv("amwo.mig.f17.csv")
data17$date <- as.POSIXct(data17$date, origin="1970-01-01 04:00:00")
data17$date.ord <- as.numeric(round(data17$date - as.POSIXct("2017-10-01", origin="1970-01-01 04:00:00")))
data17$time <- as.POSIXct(data17$time, origin="1970-01-01 04:00:00")
data17$time.ord <- as.numeric(round(data17$time - as.POSIXct("2017-10-01", origin="1970-01-01 04:00:00")))
data <- rbind(data18, data19, data17)

data <- data[!(is.na(data$st.pr)),]
summary(data$st.pr)  ##verify no NA for st.pr exist

data$sex <- as.numeric(if_else(data$sex == 'm', 0, 1))  #males=0; female=1
data$age <- as.numeric(if_else(data$age == 'juv', 0, 1))  #ad=1; juv=0
## done with data prep

##A priori Geospatial model... determining support for lat, long, st.start (Geospatial)
lat <- lmer(date.ord~ lat + (1|ID),data = data, REML=FALSE)
lon <- lmer(date.ord ~ lon + (1|ID),data = data, REML=FALSE)
lat.plus.lon <- lmer(date.ord ~ lon + lat + (1|ID),data = data, REML=FALSE)
lat.times.lon <- lmer(date.ord ~ lon*lat + (1|ID),data = data, REML=FALSE)
state <- lmer(date.ord ~ as.factor(st.pr) + (1|ID),data = data, REML=FALSE)
null <- lmer(date.ord ~ 1 + (1|ID),data = data, REML=FALSE)

cand.set <- list(lat, lon, lat.plus.lon, lat.times.lon, state, null)
modnames <- c("lat", "lon", "lat.plus.lon", "lat.times.lon", "state", "null")
aictab(cand.set=cand.set, modnames = modnames,  digits = 2, second.ord = T)
summary(state)$coefficients
#state model has the greatest model support with 1.0 model weight

##A priori Demographic model... determing support for age, sex (Demographic)
age <- lmer(date.ord ~ age + as.factor(st.pr) + (1|ID),data = data, REML=FALSE)
sex <- lmer(date.ord ~ sex + as.factor(st.pr) + (1|ID),data = data, REML=FALSE)
age.plus.sex <-lmer(date.ord ~ age + sex + as.factor(st.pr) + (1|ID),data = data, REML=FALSE)
age.times.sex <- lmer(date.ord ~ agesex + as.factor(st.pr) + (1|ID),data = data, REML=FALSE)
null <- lmer(date.ord ~ as.factor(st.pr) + (1|ID),data = data, REML=FALSE)

cand.set <- list(age, sex, age.plus.sex, age.times.sex, null)
modnames <- c("age", "sex", "age.plus.sex", "age.times.sex", "null")
aictab(cand.set=cand.set, modnames = modnames,  digits = 2, second.ord = T)
summary(sex)$coefficients
summary(age.times.sex)
summary(age.plus.sex)$coefficients
#age.times.sex model received the most support with 1.0 cumulative support





###Also have code used to make migration predication rasters, but I did not think plotting code need be included
