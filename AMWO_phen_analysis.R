#calculate r-squared for models using the following....
library(rsq)
rsq(age.times.sex, adj = TRUE)




rm(list=ls())
setwd("G:/My Drive/PhD_AMWO/AMWo_Phenology/AMWO_phen_analysis/")

library(lme4)
library(AICcmodavg)
library(ggplot2)
library(wesanderson)
library(RColorBrewer)
library(tidyverse)
library(ggridges)

############################################################################
### Migration initiation analysis framework (combine 2017, 2018,  and 2019)

# initation data for 2019 and 2018; also included is calcuation for ordinal data
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

### date/time imported as factors, need to convert to POSIXct
  #data$mig.init.between <- as.POSIXct(data$mig.init.between, origin="1970-01-01 04:00:00")
  #data$mig.term.between <- as.POSIXct(data$mig.term.between, origin="1970-01-01 04:00:00")

#pre-processing data (indicator variables for age, sex) and convert date to 'ordinal date'?
data$sex <- as.numeric(if_else(data$sex == 'm', 0, 1))  #males=0; female=1
data$age <- as.numeric(if_else(data$age == 'juv', 0, 1))  #juv=o; ad=1
#data$m.year <- as.numeric(if_else(data$m.year == '2018', 0, 1))  #2018=o; 2019=1

  # data$init.ord <- as.numeric(round(data$mig.init.between - as.POSIXct("2019-10-01", origin="1970-01-01 04:00:00")))

##creating condition indicie
amwo.capt$condition <- amwo.capt$resid
#amwo.capt$condition[is.na(amwo.capt$condition)] <- 0  ## subset and remove legacy birds and condition=NA
amwo.capt$ID <- amwo.capt$Movebank.ID

data <- left_join(data, amwo.capt[c(3,24:25)], by = "ID") 

#data <- data %>%
#  group_by(ID) %>%
#  mutate(st.start = first(st.pr))
  


#fall 2019 (remove duplicates)
data <- data[-c(40,54),]  # remove duplicate data from re-captured birds
#fall all initiation removal
data <- data[-c(35,39,87,103),] 

##support for state vs lat long models (geospatial)
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

##support for age vs sex covarite models (demographic)
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

####  removing RI bird with NA for condition
#data$condition <- data$condition[!is.na(amwo.capt$condition)] ## subset and remove legacy birds and condition=NA
data2 <- subset(data,
                  !(is.na(condition)))

##a priori model set (biological)
   ## using lat.plus.lon and age fixed effect
   ## using lat.plus.lon as interaction effect
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


#summary of models and data
summary(mod6)
coef(mod5, complete = T)
summary(cond)$coefficients
summary(full)$coefficients
mod5$coefficients
summary(null)

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

##prducing figures that shown when woodcock initiated migration and when
  ## AMWO are predicted to initiate migration from different geopolitical boundaries

##need to remove single observation sites (KY and OH)
plot.data <- data %>%
  filter(st.start != 'KY',
         st.start != 'OH',
         st.start != 'VT')

nb.cols <- nlevels(as.factor(plot.data$st.start))
## stretch the 'Dark2' pallette into that number of colors
mycolors <- colorRampPalette(brewer.pal(8, "Blues"))(nb.cols)


#migration initation by state (data only, no model predictions)
state.plot3.1 <- ggplot () +
  geom_boxplot(data=plot.data, aes(x = reorder(as.factor(st.start),lat), # '-' before lat
                              y = as.Date(init.ord, origin = "2019-10-01"), 
                              fill = reorder(as.factor(st.start),lat)), outlier.shape=NA) + # '-' before lat
  geom_point(data=plot.data, aes(x = reorder(as.factor(st.start),lat), # '-' before lat
                                  y = as.Date(init.ord, origin = "2019-10-01")),
             position=position_jitter(width=0.15, height=1),
             size=1) +
  #scale_fill_manual(values=mycolors) +
  scale_fill_manual(values = c("steelblue4", "dodgerblue3", "deepskyblue3", "steelblue2",
                    "cadetblue3", "lightskyblue", "cadetblue2", "paleturquoise2", "lightblue1")) +
  scale_x_discrete(expand=c(0.1, 0)) +
  coord_flip() +
  theme_classic() +
  theme(legend.position = "none") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5), 
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12)) +
  labs(y = "", x = "State/Province") +
  ylim(as.Date(5, origin = "2019-10-01"),as.Date(75, origin = "2019-10-01"))
state.plot3.1

#migration initation by state (predictions only)
state.plot3.2 <- ggplot () +
  geom_pointrange(data=pred.adult, aes(x = reorder(as.factor(st.start),lat), # '-' before lat
                                       y = as.Date(init.pred, origin = "2019-10-01"), 
                                       ymin = as.Date(lower, origin = "2019-10-01"),
                                       ymax = as.Date(upper, origin = "2019-10-01")),
                  shape=15,   ##boxes for adults
                  position=position_nudge(x=0.15, y=0),
                  size=0.5) +
  geom_pointrange(data=pred.juv, aes(x = reorder(as.factor(st.start),lat), # '-' before lat
                                     y = as.Date(init.pred, origin = "2019-10-01"),
                                     ymin = as.Date(lower, origin = "2019-10-01"),
                                     ymax = as.Date(upper, origin = "2019-10-01")),
                  shape=18, ##diamonds are for juv
                  position=position_nudge(x=-0.15, y=0),
                  size=0.5) +
  #scale_fill_manual(values=mycolors) +
  scale_fill_manual(values = c("steelblue4", "dodgerblue3", "deepskyblue3", "steelblue2",
                               "cadetblue3", "lightskyblue", "cadetblue2", "paleturquoise2", "lightblue1")) +
  scale_x_discrete(expand=c(0.1, 0), limits=c("VA", "WV", "RI", "PA", "NY", "ME", "NS", "QUE", "ONT")) +
  coord_flip() +
  theme_classic() +
  theme(legend.position = "none") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5), 
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12)) +
  labs(y = "Fall initiation", x = "State/Province") +
  ylim(as.Date(5, origin = "2019-10-01"),as.Date(75, origin = "2019-10-01"))
state.plot3.2

library(cowplot)
state.plot3 <- plot_grid(
  state.plot3.1, state.plot3.2,
  labels = "AUTO", ncol = 1
)
state.plot3

ggsave("Fall_init_prediction.jpeg", device="jpeg",
       scale = 1, width = 3.25, height = 6, units = c("in"),
       dpi = 1200, limitsize = TRUE)

#scale_fill_manual(values=mycolors) +
#scale_fill_manual(values = c("midnightblue", "navyblue", "blue4", "mediumblue", "royalblue4",
#                  "dodgerblue4", "dodgerblue3", "deepskyblue2", "lightskyblue",
#                  "cadetblue3", "cadetblue2", "darkslategray1", "lightblue1")) +
scale_fill_manual(values = c("royalblue4",
                  "dodgerblue4", "dodgerblue3", "deepskyblue2", "lightskyblue",
                  "cadetblue3", "cadetblue2", "darkslategray1", "lightblue1")) +
  #scale_fill_manual(values = c("lightblue1", "powderblue", "lightblue", "blue4", "aliceblue", "lightskyblue2", "lightcyan", 
  #                             "skyblue2", "skyblue3", "paleturquoise1", "steelblue4", "blue4")) +




##not organized by latitude
state.plot4 <- ggplot () +
  geom_boxplot(data=data, aes(x = st.start, 
                              y = as.Date(init.ord, origin = "2019-10-01"), 
                              fill = st.start)) +
  scale_fill_manual(values = c("lightblue1", "lightblue", "aliceblue", "lightskyblue2", "lightcyan", 
                               "skyblue2", "skyblue3", "paleturquoise1", "steelblue4", "blue4", "navyblue", "blue4")) +
  coord_flip() +
  theme(legend.position = "none") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5), 
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold")) +
  labs(y = "Date", x = "State/Province")
state.plot4


## creating plots from the data
age.plot <- ggplot () +
  geom_boxplot(data=data, aes(x = as.factor(age), y = as.Date(init.ord, origin = "2019-10-01"))) +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5), 
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold")) +
  labs(y = "Date", x = "Age") +
  scale_x_discrete(labels=c("0" = "Adult", "1" = "Young"))
age.plot

age.plot <- ggplot () +
  geom_boxplot(data=data, aes(x = as.factor(age), y = as.Date(init.ord, origin = "2019-10-01"), 
                              fill = as.factor(age))) +
  scale_fill_brewer(palette = "Paired") +
  theme(legend.position = "none") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5), 
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold")) +
  labs(y = "Date", x = "Age") +
  scale_x_discrete(labels=c("0" = "Adult", "1" = "Young"))
age.plot

year.plot <- ggplot () +
  geom_boxplot(data=data, aes(x = as.factor(m.year), y = as.Date(init.ord, origin = "2019-10-01"))) +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5), 
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold")) +
  labs(y = "Date", x = "Year") +
  scale_x_discrete(labels=c("0" = "2018", "1" = "2019"))
year.plot

condition.plot <- ggplot () +
  geom_point(data=data, aes(x = condition, y = init.ord, colour = factor(m.state))) +
  scale_fill_manual(values = c("lightblue1", "lightblue", "aliceblue", "lightskyblue2", "lightcyan", 
                               "skyblue2", "skyblue3", "paleturquoise1", "steelblue4")) +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5), 
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold")) +
  labs(x = "Condition", y = "Migration Initiation (date)") +
  geom_smooth(data = data, aes(x = condition, y = init.ord), colour = "black", method = "lm", se = FALSE)
condition.plot

summary(data$st.start)

###########################################################################################
###########################################################################################
### Migration termination fall 
rm(list=ls())
setwd("G:/My Drive/PhD_AMWO/AMWo_Phenology/AMWO_phen_analysis/")

library(lme4)
library(AICcmodavg)
library(ggplot2)
library(wesanderson)
library(RColorBrewer)
library(tidyverse)
library(ggridges)

############################################################################
### Migration initiation analysis framework (combine 2017, 2018,  and 2019)

# initation data for 2019 and 2018; also included is calcuation for ordinal data
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

### date/time imported as factors, need to convert to POSIXct
#data$mig.init.between <- as.POSIXct(data$mig.init.between, origin="1970-01-01 04:00:00")
#data$mig.term.between <- as.POSIXct(data$mig.term.between, origin="1970-01-01 04:00:00")

#pre-processing data (indicator variables for age, sex) and convert date to 'ordinal date'?
data$sex <- as.numeric(if_else(data$sex == 'm', 0, 1))  #males=0; female=1
data$age <- as.numeric(if_else(data$age == 'juv', 0, 1))  #ad=o; juv=1
#data$m.year <- as.numeric(if_else(data$m.year == '2018', 0, 1))  #2018=o; 2019=1

# data$init.ord <- as.numeric(round(data$mig.init.between - as.POSIXct("2019-10-01", origin="1970-01-01 04:00:00")))

##creating condition indicie
amwo.capt$condition <- amwo.capt$resid
#amwo.capt$condition[is.na(amwo.capt$condition)] <- 0
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

#data <- data %>%
#  group_by(ID) %>%
#  mutate(st.start = first(st.pr))

#fall 2019 (remove duplicates)
data <- data[-c(27:29, 33, 74:76, 81),]  # remove duplicate data from re-captured birds

##support for state vs lat long models (geospatial)
   ## lat long for both termination and initiation location
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

##support for age vs sex covarite models (demographic)
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

####  removing RI bird with NA for condition
#data$condition <- data$condition[!is.na(amwo.capt$condition)] ## subset and remove legacy birds and condition=NA
data3 <- subset(data,
                !(is.na(condition)))

##a priori model set (biological)
   ## using lat.plus.long.start as fixed effect
   ## using lat plus lon so we can test for interactions between condition and lat
null <- glm(term.ord ~ start.lat + start.lon,data = data3)
cond <- glm(term.ord ~ condition + start.lat + start.lon,data = data3)
cond.times.sex <- glm(term.ord ~ condition*sex + start.lat + start.lon,data = data3)
cond.times.age <- glm(term.ord ~ condition*age + start.lat + start.lon,data = data3)
cond.times.lat <- glm(term.ord ~ condition*lat + start.lat + start.lon,data = data3)
cond.times.startlat <- glm(term.ord ~ condition*start.lat + start.lat + start.lon,data = data3)
## all three condition*

cand.set <- list(null, cond, cond.times.sex, cond.times.age, cond.times.lat, cond.times.startlat)
modnames <- c("null", "cond", "cond.times.sex", "cond.times.age", "cond.times.lat", "cond.times.startlat")
aictab(cand.set=cand.set, modnames = modnames,  digits = 2, second.ord = T)
summary(cond)$coefficients
summary(null)$coefficients
## null model has the greatest model support with 0.27 sumulative weight


#summary of models and data
summary(mod6)
coef(mod5, complete = T)
summary(mod3)$coefficients
summary(full)$coefficients
mod5$coefficients
summary(null)
summary(f19.dur)

plot.data2 <- data %>%
  filter(st.pr != 'PA',
         st.pr != 'TX')

nb.cols <- nlevels(as.factor(plot.data2$st.start))
## stretch the 'Dark2' pallette into that number of colors
mycolors <- colorRampPalette(brewer.pal(8, "Blues"))(nb.cols)

# migration termination plot (data only, no model predictions)
state.plot6 <- ggplot () +
  geom_boxplot(data=plot.data2, aes(x = reorder(as.factor(st.pr),lat), 
                              y = as.Date(term.ord, origin = "2019-10-01"), 
                              fill = reorder(as.factor(st.pr),lat)), outlier.shape=NA) +
  geom_point(data=plot.data2, aes(x = reorder(as.factor(st.pr),lat), # '-' before lat
                                 y = as.Date(term.ord, origin = "2019-10-01")),
             position=position_jitter(width=0.15, height=0.1),
             size=1) +
  scale_fill_manual(values = c("steelblue4", "dodgerblue3", "deepskyblue3", "steelblue2",
                               "cadetblue3", "lightskyblue", "cadetblue2", "lightskyblue2", "paleturquoise2", 
                               "lightblue1", "paleturquoise1", "lightcyan1")) +
  coord_flip() +
  theme_classic() +
  theme(legend.position = "none") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5), 
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12)) +
  labs(y = "Fall termination", x = "State")
state.plot6

ggsave("Fall_term_Final1.jpeg", device="jpeg",
       scale = 1, width = 3.25, height = 3, units = c("in"),
       dpi = 1200, limitsize = TRUE)


plot.data6 <- data %>%
  filter(st.start != 'VT',
         st.start != 'KY')

### Model predictions
### using the inference model to predict distributions and calculate 95% CI
f.term.pred <- read.csv("fall_term_predict.csv")
prediction<- predict.glm(null, newdata=f.term.pred, se.fit=TRUE)
fall.pred<- cbind(f.term.pred, prediction)
colnames(fall.pred)<- c("st.start", "start.lon", "start.lat","term.pred","SE")
fall.pred$lower<- fall.pred$term.pred-(fall.pred$SE*1.96)
fall.pred$upper<- fall.pred$term.pred+(fall.pred$SE*1.96)

# migration termination plot 2 (data only, no model predictions)
state.plot7.1 <- ggplot () +
  geom_boxplot(data=plot.data6, aes(x = reorder(as.factor(st.start), lat), 
                                    y = as.Date(term.ord, origin = "2019-10-01"),
                                    fill = reorder(as.factor(st.start),lat)), outlier.shape=NA) +
  geom_point(data=plot.data6, aes(x = reorder(as.factor(st.start),lat), # '-' before lat
                                  y = as.Date(term.ord, origin = "2019-10-01")),
             position=position_jitter(width=0.15, height=0.1),
             size=1) +
  scale_fill_manual(values = c("steelblue4", "dodgerblue3", "deepskyblue3", "steelblue2",
                               "cadetblue3", "lightskyblue", "cadetblue2", "lightskyblue2", "paleturquoise2", 
                               "lightblue1", "paleturquoise1", "lightcyan1")) +
  scale_x_discrete(limits=c("VA", "WV" ,"PA", "RI", "NY", "NS", "ME", "QUE", "ONT")) +
  coord_flip() +
  theme_classic() +
  theme(legend.position = "none") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5), 
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12)) +
  labs(y = "", x = "State/Province") +
  ylim(as.Date(35, origin = "2019-10-01"),as.Date(95, origin = "2019-10-01"))
state.plot7.1

state.plot7.2 <- ggplot () +
  geom_pointrange(data=fall.pred, aes(x = reorder(as.factor(st.start),start.lat), # '-' before lat
                                       y = as.Date(term.pred, origin = "2019-10-01"), 
                                       ymin = as.Date(lower, origin = "2019-10-01"),
                                       ymax = as.Date(upper, origin = "2019-10-01")),
                  size=0.5) +
  scale_fill_manual(values = c("steelblue4", "dodgerblue3", "deepskyblue3", "steelblue2",
                               "cadetblue3", "lightskyblue", "cadetblue2", "lightskyblue2", "paleturquoise2", 
                               "lightblue1", "paleturquoise1", "lightcyan1")) +
  coord_flip() +
  theme_classic() +
  theme(legend.position = "none") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5), 
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12)) +
  labs(y = "Fall termination", x = "State/Province") +
  ylim(as.Date(35, origin = "2019-10-01"),as.Date(95, origin = "2019-10-01"))
state.plot7.2

library(cowplot)
state.plot7 <- plot_grid(
  state.plot7.1, state.plot7.2,
  labels = "AUTO", ncol = 1
)
state.plot7

ggsave("Fall_term2_prediction.jpeg", device="jpeg",
       scale = 1, width = 3.25, height = 6, units = c("in"),
       dpi = 1200, limitsize = TRUE)


summary(data$st.pr)

###########################################################################################
###########################################################################################
### Spring migration initaition
rm(list=ls())
setwd("G:/My Drive/PhD_AMWO/AMWo_Phenology/AMWO_phen_analysis/")

library(lme4)
library(AICcmodavg)
library(ggplot2)
library(wesanderson)
library(RColorBrewer)
library(tidyverse)
library(ggridges)

############################################################################
### Migration initiation analysis framework (combine 2019 and 2020)

# initation data for 2019 and 2018; also included is calcuation for ordinal data
data19 <- read.csv("amwo.init.sp19.csv")
data19$mig.init.between <- as.POSIXct(data19$mig.init.between, origin="1970-01-01 04:00:00")
data19$init.ord <- as.numeric(round(data19$mig.init.between - as.POSIXct("2019-01-01", origin="1970-01-01 04:00:00")))
data20 <- read.csv("amwo.init.sp20.csv")
data20$mig.init.between <- as.POSIXct(data20$mig.init.between, origin="1970-01-01 04:00:00")
data20$init.ord <- as.numeric(round(data20$mig.init.between - as.POSIXct("2020-01-01", origin="1970-01-01 04:00:00")))
data <- rbind(data20, data19)

amwo.capt <- read.csv("amwo.capt_condition.csv")
amwo.capt <- read.csv("G:/My Drive/PhD_AMWO/AMWo_Phenology/AMWO_Condition_Regression/amwo.capt_condition.csv")

### date/time imported as factors, need to convert to POSIXct
#data$mig.init.between <- as.POSIXct(data$mig.init.between, origin="1970-01-01 04:00:00")
#data$mig.term.between <- as.POSIXct(data$mig.term.between, origin="1970-01-01 04:00:00")

#pre-processing data (indicator variables for age, sex) and convert date to 'ordinal date'?
data$sex <- as.numeric(if_else(data$sex == 'm', 0, 1))  #males=0; female=1
data$age <- as.numeric(if_else(data$age == 'juv', 0, 1))  #ad=o; juv=1
#data$m.year <- as.numeric(if_else(data$m.year == '2018', 0, 1))  #2018=o; 2019=1

# data$init.ord <- as.numeric(round(data$mig.init.between - as.POSIXct("2019-10-01", origin="1970-01-01 04:00:00")))

##creating condition indicie
amwo.capt$condition <- amwo.capt$resid
#amwo.capt$condition[is.na(amwo.capt$condition)] <- 0
amwo.capt$ID <- amwo.capt$Movebank.ID

data <- left_join(data, amwo.capt[c(3,24:25)], by = "ID") 

#data <- data %>%
#  group_by(ID) %>%
#  mutate(st.start = first(st.pr))

#fall 2019 (remove duplicates)
data <- data[-c(57, 100, 111, 113, 114),]  # remove duplicate data from re-captured birds

#removing birds with no state listed
data <- subset(data,
                !(is.na(st.start)))

##support for state vs lat long models (geospatial)
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

data <- data[-c(64, 108),]  # remove birds with no age sex data
##support for age vs sex covarite models (demographic)
age <- glm(init.ord ~ age + lon,data = data)
sex <- glm(init.ord ~ sex + lon,data = data)
age.plus.sex <- glm(init.ord ~ age + sex + lon,data = data)
age.times.sex <- glm(init.ord ~ agesex + lon,data = data)
null <- glm(init.ord ~ lon,data = data)

cand.set <- list(age, sex, age.plus.sex, age.times.sex, null)
modnames <- c("age", "sex", "age.plus.sex", "age.times.sex", "null")
aictab(cand.set=cand.set, modnames = modnames,  digits = 2, second.ord = T)
summary(sex)$coefficients
#sex model received the most support with 0.53 cumulative support

####  removing RI bird with NA for condition
#data$condition <- data$condition[!is.na(amwo.capt$condition)] ## subset and remove legacy birds and condition=NA
data3 <- subset(data,
                !(is.na(condition)))

data3 <- data3[-c(30:33, 34:42, 50:58, 72, 77:86, 95),]  # remove AMWO marked in the fall, as conditon score is not longer acurate [NJ AWO 30:35 and 85:89]
##a priori model set (biological)
  ## using state and sex fixed effect
  ## using lat plus lon as so that lat or lon can be included for interactions wiht condition
null <- glm(init.ord ~ sex + lon,data = data3)
cond <- glm(init.ord ~ condition + sex + lon,data = data3)
cond.times.sex <- glm(init.ord ~ condition*sex + sex + lon,data = data3)
cond.times.age <- glm(init.ord ~ condition*age + sex + lon,data = data3)
cond.times.lat <- glm(init.ord ~ condition*lat + sex + lon,data = data3)
cond.times.lon <- glm(init.ord ~ condition*lon + sex + lon,data = data3)
## all three condition*

cand.set <- list(null, cond, cond.times.sex, cond.times.age, cond.times.lat, cond.times.lon)
modnames <- c("null", "cond", "cond.times.sex", "cond.times.age", "cond.times.lat", "cond.times.lon")
aictab(cand.set=cand.set, modnames = modnames,  digits = 2, second.ord = T)
summary(cond.times.sex)$coefficients
summary(cond.times.sex2)
## cond.times.sex was the highest supported model with 0.90 cumulative weight
max(data3$condition)

## exploring sex and condition
library(modelr)
library(car)
mod.reis <- as.data.frame(resid(cond.times.sex))
resid.data <- add_residuals(data3, cond.times.sex)
resid.data <- subset(resid.data, sex =='1')
plot(resid.data$init.ord, resid.data$resid)
abline(h=0)


#summary of models and data
summary(mod6)
coef(mod5, complete = T)
summary(mod3)$coefficients
summary(full)$coefficients
mod5$coefficients
summary(null)
summary(f19.dur)

### predicting initiation data
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



plot.data3 <- data %>%
  filter(st.start != 'AR',
         st.start != 'NJ')

nb.cols <- nlevels(as.factor(plot.data3$st.start))
## stretch the 'Dark2' pallette into that number of colors
mycolors <- colorRampPalette(brewer.pal(8, "Blues"))(nb.cols)

state.plot1.1 <- ggplot () +
  geom_boxplot(data=plot.data3, aes(x = reorder(as.factor(st.start),lon), 
                              y = as.Date(init.ord, origin = "2019-01-01"), 
                              fill = reorder(as.factor(st.start),lon)), outlier.shape=NA) +
  geom_point(data=plot.data3, aes(x = reorder(as.factor(st.start),lat), # '-' before lat
                                 y = as.Date(init.ord, origin = "2019-01-01")),
             position=position_jitter(width=0.15, height=1),
             size=1) +
  scale_fill_manual(values = c("steelblue4", "dodgerblue3", "deepskyblue3", "steelblue2",
                               "cadetblue3", "lightskyblue", "cadetblue2", "paleturquoise2", "lightblue1")) +
  coord_flip() +
  theme_classic() + 
  theme(legend.position = "none") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5), 
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12)) +
  labs(y = "", x = "State/Province") +
  ylim(as.Date(15, origin = "2019-01-01"),as.Date(100, origin = "2019-01-01"))
state.plot1.1

state.plot1.2 <- ggplot () +
  geom_pointrange(data=pred.male, aes(x = reorder(as.factor(st.start),lat), # '-' before lat
                                      y = as.Date(init.pred, origin = "2019-01-01"),
                                      ymin = as.Date(lower, origin = "2019-01-01"), 
                                      ymax = as.Date(upper, origin = "2019-01-01")),
                  position=position_nudge(x=0.15, y=0),
                  size=0.5,
                  shape=15) +   ##boxes for males) +
  geom_pointrange(data=pred.female, aes(x = reorder(as.factor(st.start),lat), # '-' before lat
                                        y = as.Date(init.pred, origin = "2019-01-01"), 
                                        ymin = as.Date(lower, origin = "2019-01-01"), 
                                        ymax = as.Date(upper, origin = "2019-01-01")),
                  position=position_nudge(x=-0.15, y=0),
                  size=0.5,
                  shape=18) +   ##diamonds for females) +
  scale_fill_manual(values = c("steelblue4", "dodgerblue3", "deepskyblue3", "steelblue2",
                               "cadetblue3", "lightskyblue", "cadetblue2", "paleturquoise2", "lightblue1")) +
  scale_x_discrete(limits=c("LA", "MS", "AL", "GA", "SC", "NC", "VA", "MD", "RI")) +
  coord_flip() +
  theme_classic() + 
  theme(legend.position = "none") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5), 
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12)) +
  labs(y = "Spring initiation", x = "State/Province") +
  ylim(as.Date(15, origin = "2019-01-01"),as.Date(100, origin = "2019-01-01"))
state.plot1.2

library(cowplot)
state.plot1 <- plot_grid(
  state.plot1.1, state.plot1.2,
  labels = "AUTO", ncol = 1
)
state.plot1

ggsave("Spring_init_prediction.jpeg", device="jpeg",
       scale = 1, width = 3.25, height = 6, units = c("in"),
       dpi = 1200, limitsize = TRUE)



#### predicting how condition influenced migration initiation date
  ##1) predict relationship, 2) plot relationship
### predicting initiation data
### using condition.sex as inference model
### condition used was mean=1.08 
sp.condit.pred <- read.csv("sp.condit.predict.csv")
prediction<- predict.glm(cond.times.sex, newdata=sp.condit.pred, se.fit=TRUE)
sp.pred<- cbind(sp.condit.pred, prediction)
colnames(sp.pred)<- c("lon","condition","sex","init.pred","SE")
sp.pred$lower<- sp.pred$init.pred-(sp.pred$SE*1.96)
sp.pred$upper<- sp.pred$init.pred+(sp.pred$SE*1.96)
pred.male <- subset(sp.pred, sex=="0")
pred.female <- subset(sp.pred, sex=="1")

data3.male <- subset(data3, sex == "0")
data3.female <- subset(data3, sex == "1")

condition.male <- ggplot () +
  geom_line(data=pred.male, aes(y=as.Date(init.pred, origin = "2019-01-01"), x=condition), size = 1) +
  geom_ribbon(data=pred.male, aes(x=condition, ymin=as.Date(lower, origin = "2019-01-01"),
                                  ymax=as.Date(upper, origin = "2019-01-01")), alpha=0.3) +
  geom_point(data=data3.male, aes(x = condition, y = as.Date(init.ord, origin = "2019-01-01")), shape = 16, size = 2, color = "black") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5), 
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12)) +
  labs(x = "", y = "Migration initiation (male)") +
  xlim(-26, 28) +
  theme_classic()
condition.male

condition.female <- ggplot () +
  geom_line(data=pred.female, aes(y=as.Date(init.pred, origin = "2019-01-01"), x=condition), size = 1) +
  geom_ribbon(data=pred.female, aes(x=condition, ymin=as.Date(lower, origin = "2019-01-01"), 
                                    ymax=as.Date(upper, origin = "2019-01-01")), alpha=0.3) +
  geom_point(data=data3.female, aes(x = condition, y = as.Date(init.ord, origin = "2019-01-01")), shape = 17, size = 2, color = "black") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5), 
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12)) +
  labs(x = "Body condition", y = "Migration initiation (female)") +
  xlim(-26, 28) +
  theme_classic()
condition.female

library(cowplot)
cond.plot <- plot_grid(
  condition.male, condition.female,
  labels = "AUTO", ncol = 1
)
cond.plot

ggsave("Spring_condition_prediction.jpeg", device="jpeg",
       scale = 1, width = 3.25, height = 6, units = c("in"),
       dpi = 1200, limitsize = TRUE)




#############################################################################################
#############################################################################################
### All migratory (state 2 lcoations) for plotting
   ### when birds are in each geopolitical boundary
rm(list=ls())
setwd("G:/My Drive/PhD_AMWO/AMWo_Phenology/AMWO_phen_analysis/")

library(lme4)
library(AICcmodavg)
library(ggplot2)
library(wesanderson)
library(RColorBrewer)
library(tidyverse)
library(ggridges)

data <- read.csv("amwo.dur.f19.csv")
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
    ## need to determin minimal sample size for some states

###### removed extraction code from "NY_migration_esraction (bottom)

data$sex <- as.numeric(if_else(data$sex == 'm', 0, 1))  #males=0; female=1
data$age <- as.numeric(if_else(data$age == 'juv', 0, 1))  #ad=1; juv=0

##support for state vs lat long models (geospatial)
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

##support for age vs sex covarite models (demographic)
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
#nage.times.sex model received the most support with 1.0 cumulative support
summary(as.factor(data$agesex))
summary(data$st.pr)


### predicting stopover data
### using condition.sex as inference model
### condition used was mean=1.08 
stopover.pred <- read.csv("stopover.predict.csv")
prediction<- predict.glm(age.times.sex, newdata=stopover.pred, se.fit=TRUE)
sp.pred<- cbind(sp.init.pred, prediction)
colnames(sp.pred)<- c("st.start", "lon", "lat","condition","sex","init.pred","SE")
sp.pred$lower<- sp.pred$init.pred-(sp.pred$SE*1.96)
sp.pred$upper<- sp.pred$init.pred+(sp.pred$SE*1.96)
pred.male <- subset(sp.pred, sex=="0")
pred.female <- subset(sp.pred, sex=="1")




nb.cols <- nlevels(as.factor(data$st.pr))
## stretch the 'Dark2' pallette into that number of colors
mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.cols)


mig.boxplot1 <- ggplot () +
  geom_boxplot(data=data, aes(x = reorder(as.factor(st.pr),lat), 
                              y = as.Date(time.ord, origin = "2019-10-01"), 
                              fill = reorder(as.factor(st.pr),lat)), outlier.shape=NA) +
  geom_point(data=data, aes(x = reorder(as.factor(st.pr),lat), # '-' before lat
                                  y = as.Date(time.ord, origin = "2019-10-01")),
             position=position_jitter(width=0.15, height=3),
             size=0.25) +
  scale_fill_manual(values=mycolors) +
  #scale_fill_manual(values = c("lightblue1", "lightblue", "aliceblue", "lightskyblue2", "lightcyan", 
  #                             "skyblue2", "skyblue3", "paleturquoise1", "steelblue4")) +
  coord_flip() +
  theme_classic() +
  theme(legend.position = "none") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5), 
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12)) +
  labs(y = "Stopover locations", x = "State/Province")
mig.boxplot1

mig.boxplot2 <- ggplot () +
  geom_boxplot(data=data, aes(x = reorder(as.factor(st.pr),lat), 
                              y = as.Date(time.ord, origin = "2019-10-01"), 
                              fill = reorder(as.factor(st.pr),lat)), outlier.shape=NA) +
  geom_point(data=data, aes(x = reorder(as.factor(st.pr),lat), # '-' before lat
                            y = as.Date(time.ord, origin = "2019-10-01")),
             position=position_jitter(width=0.15, height=3),
             size=0.25) +
  scale_fill_manual(values=mycolors) +
  #scale_fill_manual(values = c("lightblue1", "lightblue", "aliceblue", "lightskyblue2", "lightcyan", 
  #                             "skyblue2", "skyblue3", "paleturquoise1", "steelblue4")) +
  coord_flip() +
  theme_classic() +
  theme(legend.position = "none") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5), 
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12)) +
  labs(y = "Stopover timing", x = "State/Province")
mig.boxplot2




ggsave("Fall_stopover_Final.jpeg", device="jpeg",
       scale = 1, width = 6.5, height = 6, units = c("in"),
       dpi = 1200, limitsize = TRUE)



mig.ridges <- ggplot () +
  geom_density_ridges2(data=data, aes(y = reorder(as.factor(st.pr),lat), 
                                     x = as.Date(time.ord, origin = "2019-10-01"), 
                                     fill = st.pr)) +
  #facet_wrap(~st.pr) +
  #scale_fill_manual(values = c("lightblue1", "lightblue", "aliceblue", "lightskyblue2", "lightcyan", 
  #                             "skyblue2", "skyblue3", "paleturquoise1", "steelblue4")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5), 
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold")) +
  labs(y = "Date", x = "State/Province")
mig.ridges

summary(data$st.pr)

###########################################################################################
###########################################################################################
###  Plotting the relationship between lat long and migration initation
    ### Output is a raster surface predicting migration initiation across capture range
    ###   basically a nice visualization for the manuscript

rm(list=ls())
setwd("G:/My Drive/PhD_AMWO/AMWo_Phenology/AMWO_phen_analysis/")
library(lme4)
library(AICcmodavg)
library(ggplot2)
library(wesanderson)
library(RColorBrewer)
library(tidyverse)
library(ggridges)
# initation data for 2019 and 2018; also included is calcuation for ordinal data
data19 <- read.csv("amwo.init.f19.csv")
data19$mig.init.between <- as.POSIXct(data19$mig.init.between, origin="1970-01-01 04:00:00")
data19$init.ord <- as.numeric(round(data19$mig.init.between - as.POSIXct("2019-10-01", origin="1970-01-01 04:00:00")))
data18 <- read.csv("amwo.init.f18.csv")
data18$mig.init.between <- as.POSIXct(data18$mig.init.between, origin="1970-01-01 04:00:00")
data18$init.ord <- as.numeric(round(data18$mig.init.between - as.POSIXct("2018-10-01", origin="1970-01-01 04:00:00")))
data17 <- read.csv("amwo.init.f17.csv")
data17$mig.init.between <- as.POSIXct(data17$mig.init.between, origin="1970-01-01 04:00:00")
data17$init.ord <- as.numeric(round(data17$mig.init.between - as.POSIXct("2017-10-01", origin="1970-01-01 04:00:00")))
data <- rbind(data18, data19)
data <- rbind(data18, data19, data17)
amwo.capt <- read.csv("G:/My Drive/PhD_AMWO/AMWo_Phenology/AMWO_Condition_Regression/amwo.capt_condition.csv")
#pre-processing data (indicator variables for age, sex) and convert date to 'ordinal date'?
data$sex <- as.numeric(if_else(data$sex == 'm', 0, 1))  #males=0; female=1
data$age <- as.numeric(if_else(data$age == 'juv', 0, 1))  #ad=o; juv=1
#data$m.year <- as.numeric(if_else(data$m.year == '2018', 0, 1))  #2018=o; 2019=1

##creating condition indicie
amwo.capt$condition <- amwo.capt$resid
amwo.capt$condition[is.na(amwo.capt$condition)] <- 0
amwo.capt$ID <- amwo.capt$Movebank.ID

data <- left_join(data, amwo.capt[c(3,24:25)], by = "ID") 
#fall 2019 (remove duplicates)
data <- data[-c(40,54),]  # remove duplicate data from re-captured birds
#fall all initiation removal
data <- data[-c(35,39,87,103),] 

######
#lat lon plots
mod2 <- lmer(init.ord ~ agesex + lat + lon + (1|Site),data = data)
mod3 <- glm(init.ord ~ lon + lat + lon*lat ,data = data)

cond.times.age <- glm(init.ord ~ condition*age + age + lat + lon,data = data2)

#lat lon interacton plot
lat.nd<- as.data.frame(seq(min(data$lat, na.rm = TRUE),max(data$lat, na.rm = TRUE),.1))
lon.nd<- as.data.frame(seq(min(data$lon, na.rm = TRUE),max(data$lon, na.rm = TRUE),.1))

lat.nd.all<- rep(lat.nd[,1], each=nrow(lon.nd))
lon.nd.all<- rep(lon.nd[,1], each=nrow(lat.nd))

df.latlon<- data.frame(expand.grid(lon.nd[,1], lat.nd[,1]))
colnames(df.latlon)<- c("lon", "lat")

mod3.predict<- predict(mod3, newdata=df.latlon)
#view(mod3.predict)

predict.df<-data.frame(lon=df.latlon$lon, lat=df.latlon$lat, predict=mod3.predict)

## plot of dataframe (fall migration initiation)
ggplot(predict.df, aes(lon, lat, fill= predict)) + 
  geom_tile()



#### trying to combine the data.frame erik predicted with a spatial area of interest
###plotting US and CAN; clipping raster to 
library(raster)
library(sf)
library(rgeos)
##convertine data.frame of migration initations into a raster object
dfr <- rasterFromXYZ(predict.df, crs ="+proj=longlat +datum=NAD83 +no_defs" )
plot(dfr)

#reading in shape files
can <- st_read("G:/My Drive/PhD_AMWO/Collaborator_Reports/Report_2_2019/R_Code_From_Erin/Canada/Canada.shp")
us <- st_read("G:/My Drive/PhD_AMWO/Collaborator_Reports/Report_2_2019/R_Code_From_Erin/tl_2018_us_state/tl_2018_us_state.shp")
gl <- st_read("G:/My Drive/PhD_AMWO/Collaborator_Reports/Report_2_2019/R_Code_From_Erin/Great_Lakes-shp/Great_Lakes.shp")

#states within the range of AMWO
us$STUPS <- as.factor(us$STUPS)
us4 <- us %>%
  dplyr::filter(STUSPS %in% c("ME", "VT", "NH", "MA", "RI", "PA", "NY", "NJ", "DC", "CT", 
                              "MD", "VA", "WV", "DE", "NC", "SC", "GA", "AL", "FL", "TN",
                              "KY", "MS", "LA", "TX", "OK", "AR", "MO", "OH", "IN", "IL",
                              "IA", "MN", "MI", "WI", "KS", "ND", "SD", "NE")) 
plot(us4[6])

#provinces in the range of AMWO
can4 <- can %>%
  dplyr::filter(NAME %in% c("Ontario", "Quebec", "New Brunswick", "Nova Scotia", "Prince Edward Island",
                            "Newfoundland and Labrador", "Manitoba", "Saskatchewan"))
plot(can4[1])

#attempting to clip/crop the raster with the shapefile, but no luck...
clip <- gIntersection(us4[6], dfr, byid = FALSE)
clip <- crop(us4[6], dfr)
plot(clip)
plot(us4[6])

##plot of shapefile layers... however cannot get the raster to dispaly or interact with these shape files
ggplot() +
  geom_sf(data=dfr, fill="Dark2", color = "black") +
  geom_sf(data=us4, fill="#99FF66", color= "black") +  
  geom_sf(data=can4, fill="#99FF66", color= "black") +
  geom_sf(data=gl, fill="#0000FF", color= "black") +
  annotation_scale(location = "bl", width_hint = 0.75) +
  annotation_north_arrow(location = "tr", which_north = "true", 
                         pad_x = unit(0.25, "in"), pad_y = unit(0.25, "in"),
                         style = north_arrow_fancy_orienteering) +
  xlab("Longitude") +
  ylab("Latitude") +
  theme_bw() 

summary(us)
crs(dfr)
crs(us)
crs(can)
#





























#### I beleive all code below is garbage... avenues pursued but in the end, not fruitful





## migration efficiency 
###########################################################################################
###########################################################################################
### Fall migration duration
rm(list=ls())
setwd("G:/My Drive/PhD_AMWO/AMWo_Phenology/AMWO_phen_analysis/")

library(lme4)
library(AICcmodavg)
library(ggplot2)
library(wesanderson)
library(RColorBrewer)
library(tidyverse)
library(ggridges)

# duration data for 2019, 2018, and 2017; also included is calcuation for ordinal data
data19 <- read.csv("amwo.dur.f19.csv")
data18 <- read.csv("amwo.dur.f18.csv")
data17 <- read.csv("amwo.dur.f17.csv")
data <- rbind(data19, data18, data17)

amwo.capt <- read.csv("G:/My Drive/PhD_AMWO/AMWo_Phenology/AMWO_Condition_Regression/amwo.capt_condition.csv")

#pre-processing data (indicator variables for age, sex) and convert date to 'ordinal date'?
data$sex <- as.numeric(if_else(data$sex == 'm', 0, 1))  #males=0; female=1
data$age <- as.numeric(if_else(data$age == 'juv', 0, 1))  #ad=o; juv=1
#data$m.year <- as.numeric(if_else(data$m.year == '2018', 0, 1))  #2018=o; 2019=1

# data$init.ord <- as.numeric(round(data$mig.init.between - as.POSIXct("2019-10-01", origin="1970-01-01 04:00:00")))

##creating condition indicie
amwo.capt$condition <- amwo.capt$resid
#amwo.capt$condition[is.na(amwo.capt$condition)] <- 0
amwo.capt$ID <- amwo.capt$Movebank.ID

data <- left_join(data, amwo.capt[c(3,24:25)], by = "ID") 

##inporting migration data, to calculate migration length for each bird
data19 <- read.csv("amwo.mig.f19.csv")
data19$time <- as.POSIXct(data19$time, origin="1970-01-01 04:00:00")
data19$time.ord <- as.numeric(round(data19$time - as.POSIXct("2019-10-01", origin="1970-01-01 04:00:00")))
data18 <- read.csv("amwo.mig.f18.csv")
data18$time <- as.POSIXct(data18$time, origin="1970-01-01 04:00:00")
data18$time.ord <- as.numeric(round(data18$time - as.POSIXct("2018-10-01", origin="1970-01-01 04:00:00")))
data17 <- read.csv("amwo.mig.f17.csv")
data17$time <- as.POSIXct(data17$time, origin="1970-01-01 04:00:00")
data17$time.ord <- as.numeric(round(data17$time - as.POSIXct("2017-10-01", origin="1970-01-01 04:00:00")))
data.mig <- rbind(data18, data19, data17)
data.mig <- subset(data.mig,!(is.na(step)))

data.mig2 <- data.mig %>%
  group_by(ID) %>%
  dplyr::summarize(mig.dist=sum(step))

data <- left_join(data, data.mig2[c(1:2)], by = "ID") 

#fall 2019 (remove duplicates)
data <- data[-c(41, 46, 49, 55, 92, 96),]  # remove duplicate data from re-captured birds

## caulcualting the average distance per day as a measure of migratory efficiency
data <- data %>%
  mutate(dist.day=mig.dist/dur.between)

##support for state vs lat long models (geospatial)
lat <- glm(dur.between ~ lat,data = data)
lon <- glm(dur.between ~ lon,data = data)
lat.plus.lon <- glm(dur.between ~ lon + lat,data = data)
lat.times.lon <- glm(dur.between ~ lon*lat,data = data)
state.start <- glm(dur.between ~ as.factor(mig.start),data = data)
null <- glm(dur.between ~ 1,data = data)

##support for state vs lat long models (geospatial) dist per day
lat <- glm(dist.day ~ lat,data = data)
lon <- glm(dist.day ~ lon,data = data)
lat.plus.lon <- glm(dist.day ~ lon + lat,data = data)
lat.times.lon <- glm(dist.day ~ lon*lat,data = data)
state.start <- glm(dist.day ~ as.factor(mig.start),data = data)
null <- glm(dist.day ~ 1,data = data)

cand.set <- list(lat, lon, lat.plus.lon, lat.times.lon, state.start, null)
modnames <- c("lat", "lon", "lat.plus.lon", "lat.times.lon", "state.start", "null")
aictab(cand.set=cand.set, modnames = modnames,  digits = 2, second.ord = T)
summary(lat.times.lon)$coefficients
#lat.times.lon model has the greaters model support with 0.71 model weight

##support for age vs sex covarite models (demographic)
age <- glm(dur.between ~ age + lon*lat,data = data)
sex <- glm(dur.between ~ sex + lon*lat,data = data)
age.plus.sex <- glm(dur.between ~ age + sex + lon*lat,data = data)
age.times.sex <- glm(dur.between ~ agesex + lon*lat,data = data)
null <- glm(dur.between ~ lon*lat,data = data)

##support for age vs sex covarite models (demographic) dist per day
age <- glm(dist.day ~ age,data = data)
sex <- glm(dist.day ~ sex,data = data)
age.plus.sex <- glm(dist.day ~ age + sex,data = data)
age.times.sex <- glm(dist.day ~ agesex,data = data)
null <- glm(dist.day ~ 1,data = data)

cand.set <- list(age, sex, age.plus.sex, age.times.sex, null)
modnames <- c("age", "sex", "age.plus.sex", "age.times.sex", "null")
aictab(cand.set=cand.set, modnames = modnames,  digits = 2, second.ord = T)
summary(age.plus.sex)$coefficients
summary(age.times.sex)$coefficients
#null model received the most support with 0.53 cumulative support

####  removing RI bird with NA for condition
#data$condition <- data$condition[!is.na(amwo.capt$condition)] ## subset and remove legacy birds and condition=NA
data2 <- subset(data,
                !(is.na(condition)))

data2 <- data2[-c(40:42),]  # remove AMWO marked in the spring, as conditon score is not longer acurate
##a priori model set (biological)
## using state and sex fixed effect
## using lat plus lon as so that lat or lon can be included for interactions wiht condition
null <- glm(dur.between ~ lon + lat,data = data2)
cond <- glm(dur.between ~ condition + lon + lat,data = data2)
cond.times.sex <- glm(dur.between ~ condition*sex + lon + lat,data = data2)
cond.times.age <- glm(dur.between ~ condition*age + lon + lat,data = data2)
cond.times.lat <- glm(dur.between ~ condition*lat + lon + lat,data = data2)
cond.times.lon <- glm(dur.between ~ condition*lon + lon + lat,data = data2)
## all three condition*

cand.set <- list(null, cond, cond.times.sex, cond.times.age, cond.times.lat, cond.times.lon)
modnames <- c("null", "cond", "cond.times.sex", "cond.times.age", "cond.times.lat", "cond.times.lon")
aictab(cand.set=cand.set, modnames = modnames,  digits = 2, second.ord = T)
summary(cond.times.sex)$coefficients
summary(cond.times.sex2)
## null was the highest supported model with 0.48 cumulative weight






###########################################################################
#####
##### Determining migration efficiency
#####
###########################################################################
library(tidyverse)

rm(list=ls())
setwd("G:/My Drive/PhD_AMWO/AMWo_Phenology/AMWO_phen_analysis/loc_efficiency/")

## ids from fall migration (initiation)
data19 <- read.csv("amwo.init.f19.csv")
data18 <- read.csv("amwo.init.f18.csv")
data17 <- read.csv("amwo.init.f17.csv")


## ids from spring migration (initiation)
data19 <- read.csv("amwo.init.sp19.csv")
data20 <- read.csv("amwo.init.sp20.csv")


#import all locations from specific date range
f17 <- read.csv("f17_loc_efficiency.csv")
f18 <- read.csv("f18_loc_efficiency.csv")
f19 <- read.csv("f19_loc_efficiency.csv")
s19 <- read.csv("s19_loc_efficiency.csv")
s20 <- read.csv("s20_loc_efficiency.csv")


#joing data
data <- left_join( data20[c(2,12,5)],s20, by = c("ID"))

#summarize, then mutate to obtain efficiency (#days/#locs)
data <- data %>%
  group_by(ID) %>%
  summarise(
    min.date=min(as.Date(time.y)),
    max.date=max(as.Date(time.y)),
    n.locs=n())
data <- data %>%
  group_by(ID) %>%
  mutate(
    n.days=max.date-min.date,
    eff=n.days/n.locs
  )

write.csv(data, "spring_20.csv")


## inport all fall and all spring, average to find efficiency
s20 <- read.csv("spring_20.csv")
s19 <- read.csv("spring_19.csv")
f19 <- read.csv("fall_19.csv")
f18 <- read.csv("fall_18.csv")
f17 <- read.csv("fall_17.csv")
data2 <- rbind(s20,s19,f19,f18,f17)
## efficincy including all locs and all days (rather than averaging a culculation)
sum(data2$n.days)/sum(data2$n.locs)



