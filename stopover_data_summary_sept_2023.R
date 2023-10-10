# Fall Migration Stopover Data summary
# Fish et al. Phenology paper
# September 2023
# Last Edit - SJC Sept 2023

library(leaflet)
library(tidyverse)
library(sf)
library(terra)
library(rnaturalearth)
library(wesanderson)
library(RColorBrewer)


# ********************* 01 Read in original data *******************************
# Alex's code for importing migration timing data
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
#View(data)

# *********************** 02 Identify Stopovers ********************************

# This section uses code adapted from the migration strategy paper
# Only modified to match Alex's data structure and to better handle birds with only 1 migration point

birds <- data

# Make a category for stopped vs not
stop_dist <- 16.1 # 16 km = distance threshold for movement to a new site
birds$step[which(is.na(birds$step)==T)] <- 0 # NAs have 0 steop length

birds <- birds %>% mutate(stopover=if_else(step<stop_dist, 1, 0)) 

# then separate into brd list, this time complete with both distance and stop threshold info
ids <- unique(birds$ID) 
brd <- list()
for (i in 1:length(ids)){
  brd[[i]] <- birds[which(birds$ID==ids[i]),]
}

# categorizing loop (taken from bbpl code/adapted for amwo)
# I don't think it's worth trying to make this run-length encoding thing in tidyverse because it depends on separating 
# so many different objects.
for (q in 1:length(brd)){
  tmp <- rle(brd[[q]]$stopover) # run length encoding based on distance, clusters by whether or not the next point is >16km away
  label <- seq(1:length(tmp$lengths)) # makes a number to label each cluster with
  count <- rep(0, times=length(tmp$lengths)) # empty column to fill count
  for (i in 1:length(count)){
    count[i] <- sum(tmp$lengths[1:i]) # count (number of points so far in the data frame) is there to be able to label clusters
  }
  lt <- data.frame(lengths=tmp$lengths,values=tmp$values, count=count, label=label) # makes a new table for each bird with information for assigning labels
  #if (lt$values==0&lt$lengths>1){
  #  
  #}
  brd[[q]]$label <- c(1, rep(0, times=nrow(brd[[q]])-1)) # empty column in bird for cluster labels - starts with 1
  brd[[q]]$lengths <- c(1, rep(0, times=nrow(brd[[q]])-1))# empty column in bird data for lengths (amount of points in each cluster) -starts with 1
  if (nrow(lt)>1){
  for (j in 2:nrow(lt)){
    brd[[q]]$label[(lt$count[j-1]+1):lt$count[j]] <- lt$label[j] # give consecutive labels to each cluster
    brd[[q]]$lengths[(lt$count[j-1]+1):lt$count[j]] <- lt$lengths[j] # add the amount of points within the same cluster as each point to the data
  }}
  
  brd[[q]]$class <- rep(NA, nrow(brd[[q]])) # empty column for class
  
  # this bit applies a distance rule to the back end - might not be necessary for Alex's analysis
  if(nrow(lt)>1){
  for (h in 2:nrow(brd[[q]])){
    if (brd[[q]]$step[h-1]<16.1){
      brd[[q]]$label[h] <- brd[[q]]$label[h-1]
   }
   }
  }
  # add class for each point 
  # need  to do these in order!
  for (v in 1:nrow(brd[[q]])){
    brd[[q]]$class[which(brd[[q]]$stopover==1&brd[[q]]$lengths>1)] <- 'stopover' # if <16km and more than 1 point, stopover
    brd[[q]]$class[is.na(brd[[q]]$class)] <- 'stop' # all others (points over 16km apart) are stops
  }
  
 # print(q) # completes 118 total birds
}


# *********************** 03 Check/Reclassify ********************************

# Leaflet map with stopover labels
stopover_id_map<- function(x){
  test <- x
  ll <- unique(x$label) # list of unique classes
  test$label <- factor(test$label) # make a factor for mapping
  pal <- colorFactor(rainbow(length(ll)), test$label) # color scheme, one color for each unique stopover
  m <- leaflet(test) %>% 
    addTiles() %>% 
    setView( lng = -85, lat = 40, zoom = 4) %>% 
    addProviderTiles("Esri.WorldImagery") %>%
    addScaleBar(position='bottomright') %>% 
    addPolylines(lat=test$lat, lng=test$lon, col='white') %>%
    addControl(paste(test$ID[1]), position = "topleft")%>%
    addCircleMarkers(lat=test$lat, lng=test$lon, color=~pal(label), radius=2, label=paste(test$time, test$label), labelOptions=labelOptions(textsize="20px"))%>%
    addLegend('topright', pal=pal, values=~label)
  m
} 

stopover_id_map(brd[[1]])
lapply(brd[1:12], FUN=stopover_id_map) # these are split up just to make them easier to keep track of while checking
lapply(brd[13:24], FUN=stopover_id_map)
lapply(brd[25:36], FUN=stopover_id_map)
lapply(brd[37:48], FUN=stopover_id_map)
lapply(brd[49:60], FUN=stopover_id_map)
lapply(brd[61:72], FUN=stopover_id_map)
lapply(brd[73:84], FUN=stopover_id_map)
lapply(brd[85:96], FUN=stopover_id_map)
lapply(brd[98:109], FUN=stopover_id_map)
lapply(brd[110:118], FUN=stopover_id_map)

# example of recursive 
which(ids=='PA-2018-11') 
stopover_id_map(brd[[90]])

# list of birds to check (from Kylie):
# ME-2018-08 (2)
# NY-2018-04 (9)
# NY-2018-06 (11)
# RI-2018-11 (34) ***has 2 years of data, cut out 2019...
# NY-2019-22 (55) # *** this one I'm going to actually leave as is bc just one terminal point
# RI-2019-21 (88)
# RI-2019-23 (89) # ***this one also leaving because it's more of a moving path overlapping a stopover
# VA-2019-21 (96)

### ones that weren't super weird but I wasn't sure about ***these we will roll with because there isn't really a better way to categorize them
# NS-2019-03 (41) -- it looked like overlapping colors (dark blue   and purple) at first but I think they're actually all the same   (purple) because dark blue isn't in the label legend
# VA-2018-01 (38) -- ending steps - 
# NY-2019-31 (62) -- ending steps 

# QUE-2019-13 (81) -- ending steps zig zag
# RI-2019-24 (90) -- no overlap but the ending steps are odd


# Manual reclassifications for birds that have the "recursive" stopovers (only 6 birds)
# ME-2018-08
brd[[2]]$label[which(brd[[2]]$label==10)] <- 8
brd[[2]]$label[which(brd[[2]]$label==12)] <- 8

# NY-2018-04
brd[[9]]$label[which(brd[[9]]$label==5)] <- 4
brd[[9]]$label[which(brd[[9]]$label==6)] <- 4

#NY-2018-06
brd[[11]]$label[which(brd[[11]]$label==5)] <- 4
brd[[11]]$label[which(brd[[11]]$label==6)] <- 4

# RI-2018-11
brd[[34]] <- brd[[34]][which(year(brd[[34]]$time)==2018),]

# RI-2019-21
brd[[88]]$label[which(brd[[88]]$label==6)] <- 4

# VA-2019-21
brd[[96]]$label[which(brd[[96]]$label==9)] <- 8
brd[[96]]$label[which(brd[[96]]$label==10)] <- 8


# ******************* 04 Stopover Data Calculations*****************************

# Calculate stopover metrics
# Create Metric table to fill in 
# Loop to fill in first set of metrics 
all_stopovers <- list()
for (i in 1:length(brd)){ 
    # Data frame for stops first (one point per stop)
      tmp_stop <- brd[[i]][which(brd[[i]]$class=='stop'),]
      tbl_stop <- data.frame(
        id=tmp_stop$ID,
        label=66:(66+(nrow(tmp_stop)-1)),
        start_time = tmp_stop$time,
        end_time = tmp_stop$time + 43200, # 12 hours later in seconds
        admin_unit = tmp_stop$m.state,
        start_lat = tmp_stop$lat,
        start_lon = tmp_stop$lon,
        end_lat = tmp_stop$lat, # start/end the same for stop
        end_lon = tmp_stop$lon,
        year=tmp_stop$m.year,
        st.pr=tmp_stop$st.pr,
        duration=rep(12, nrow(tmp_stop))
      )
    # Data frame for stopovers (more than 1 point per stop)
    if (length(which(brd[[i]]$class=='stopover'))>0){
      tmp <- brd[[i]][which(brd[[i]]$class=='stopover'),]
      tbl_stopover <- matrix(NA, nrow=length(unique(tmp$label[which(tmp$class=='stopover')])), ncol=11) #add columns here
      tbl_stopover <- data.frame(tbl_stopover)
      names(tbl_stopover) <- c('id', 'label', 'start_time', 'end_time', 'admin_unit', 'start_lat', 'start_lon',
                               'end_lat', 'end_lon', 'year', 'duration') #add column labels here
      tbl_stopover$id <- rep(tmp$ID[1], times=nrow(tbl_stopover))
      tbl_stopover$label <- unique(tmp$label)
      for (p in 1:nrow(tbl_stopover)){ # anything else you want to calculate by individual stopover you put in here
        yam <- tmp[which(tmp$label==tbl_stopover$label[p]),]
          tbl_stopover$start_time[p] <- yam$time[1]
          tbl_stopover$end_time[p] <- yam$time[nrow(yam)]
          tbl_stopover$admin_unit[p] <- NA # ignore this until later
          tbl_stopover$start_lat[p] <- yam$lat[1]
          tbl_stopover$end_lat[p] <- yam$lat[nrow(yam)]
          tbl_stopover$start_lon[p] <- yam$lon[1]
          tbl_stopover$end_lon[p] <- yam$lon[nrow(yam)]
          tbl_stopover$year[p] <- yam$m.year[1]
          tbl_stopover$st.pr[p] <- yam$st.pr[1]
        if (nrow(yam)>1){
          tbl_stopover$duration[p] <- abs(difftime(yam$time[1], yam$time[nrow(yam)], tz="EST", units='hours'))}
        else if (nrow(yam)==1){ 
          tbl_stopover$duration[p] <- 12}
      }}
      tbl_stopover$start_time <- as_datetime(tbl_stopover$start_time) # for some reason this converts here but not above
      tbl_stopover$end_time <- as_datetime(tbl_stopover$end_time) 
      
      all_stopovers[[i]] <- bind_rows(tbl_stop, tbl_stopover)
      all_stopovers[[i]] <- all_stopovers[[i]] %>% arrange(end_time)
      
      #print(i) completes 118 total birds
  }

# Combine all into a data frame
# This is 1 row for each stopover (id is bird id, label is stopover label)
# duration is in hours
fall_stopover_data <- bind_rows(all_stopovers)
fall_stopover_data <- fall_stopover_data %>% rename(duration_hours = duration)


# ********************* 05 Reclassify Admin Units ******************************

nj <- st_read(dsn='nj_hunting_zones', layer='nj_hunting_zones_revised')
qc <- st_read(dsn='quebec_hunting_zones', layer='quebec_hunting_zones_5070')
ont <- st_read(dsn='ontario_hunting_zones', layer='ontario_hunting_zones')

par(mfrow=c(1,1))

qc <- qc %>%
  select(DISTRICT, geometry) %>%
  rename(dstrct_=DISTRICT)

plot(nj$geometry)
plot(qc$geometry)
plot(ont$geometry)

nj_crs <- crs(nj)
qc_crs <- crs(qc)
ont_crs <- crs(ont)

# get usa and canada from rnaturalearth package
usa<-vect(ne_states(country="United States of America",returnclass="sf"))
can<-vect(ne_states(country="Canada",returnclass="sf"))
amwoext=ext(-95,-60,26,50)
usa <- crop(usa,amwoext) 
can <- crop(can, amwoext)
plot(can)
plot(usa, add=T)

usa <-st_as_sf(usa) %>% # convert to sf, get rid of new jersey
  filter(name != 'New Jersey') %>%
  select('geometry', 'name') %>%
  rename(dstrct_=name)
can <- st_as_sf(can) %>% # convert to sf, get rid of quebec and ontario
  filter(name != "Québec", name !="Ontario" ) %>%
  select('geometry', 'name') %>%
  rename(dstrct_=name)

crs_usa <- crs(usa)

# transform data
nj <- st_transform(nj, crs_usa)
qc <- st_transform(qc, crs_usa) 
ont <- st_transform(ont, crs_usa)

zones <- bind_rows(usa,can, nj, ont, qc) 
zones <- zones %>% rename(admin_unit=dstrct_)

head(zones)
ggplot(zones)+
  geom_sf(aes(fill=admin_unit))


amwo_sf <- st_as_sf(fall_stopover_data, coords=c('start_lon', 'start_lat'), remove=F) 
st_crs(amwo_sf) <- 4326
amwo_sf <- st_transform(amwo_sf, crs_usa)
amwo_sf <- st_join(amwo_sf, zones, join=st_within)

ggplot()+
  #coord_sf()+
  geom_sf(data=zones, fill='white')+
  geom_sf(data=amwo_sf, mapping=aes(color=admin_unit.y))

# looks like it worked
head(amwo_sf)
fall_stopovers <- amwo_sf %>%
  select(!admin_unit.x) %>%
  rename(admin_unit=admin_unit.y) %>%
  st_drop_geometry()

#save(fall_stopovers, file='amwo_fall_stopover_phenology_with_zones.rds')

# ABCDEF are the quebec ones - to match fall_stopover df
zone_dates <- read.csv('hunting_zone_dates.csv')
head(zone_dates)

zone_dates$open <- as_date(zone_dates$open, format='%m/%d')
zone_dates$close <- as_date(zone_dates$close, format='%m/%d')
zone_dates$reopen <- as_date(zone_dates$reopen, format='%m/%d')
zone_dates$reclose <- as_date(zone_dates$reclose, format='%m/%d')
zone_dates <- zone_dates %>%
  mutate(yday_open=yday(open),
         yday_close=yday(close),
         yday_reopen=yday(reopen),
         yday_reclose=yday(reclose))

fall_stopovers$yearless_date <- format(fall_stopovers$start_time, format='%m/%d')
fall_stopovers$yearless_date <- as_date(fall_stopovers$yearless_date, format='%m/%d')
fall_stopovers <- fall_stopovers %>%
  mutate(yday_start=yday(fall_stopovers$start_time),
         yday_end=yday(fall_stopovers$end_time)) %>%
  left_join(., zone_dates, by='admin_unit')

head(fall_stopovers)



sub <- fall_stopovers[which(is.na(fall_stopovers$yday_open)==F),] # subset of only places we had hunting informationf for

save(fall_stopovers, file='fall_stopover_plot_data_100623.rds')

# ************************ 06 Add plot *****************************************

fall_stopovers <- readRDS('fall_stopover_plot_data_100623.rds')

# remake Alex's "mig.boxplot2" from AMWO_CH1_Code script

# admin_unit is administrative unit that I pulled
# st.pr is from alex's original (I think)
# it looks like they're pretty much the same
 
data <- fall_stopovers
#data <-sub

data <- data %>%
  mutate(index=as.numeric(as.factor(admin_unit)),
         admin_unit=as.factor(admin_unit),
         st.pr=as.factor(st.pr))
# january isn't cooperating so we're just changing the year (sorry Liam if you don't like this piece of code)
data$yearless_date[which(data$yearless_date<as.Date('0000-02-01'))] <- 
  data$yearless_date[which(data$yearless_date<as.Date('0000-02-01'))] %m+% years(1)


length(unique(data$st.pr)) # without extra NJ, ONT, and QC zones
length(unique(data$admin_unit)) # with extra NJ, ONT, and QC Zones

xlabs <- as.character(unique(reorder(data$admin_unit, data$start_lat))) # this doesn't put them in the right order

nb.cols <- nlevels(as.factor(data$admin_unit))
## stretch the 'Dark2' pallette into that number of colors
mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.cols)

# what it looks like without the boxes, and admin unit as a factor
# letters by themselves are QC
mig.boxplot2 <- ggplot () +
  geom_boxplot(data=data, aes(x = reorder(admin_unit,start_lat), 
                              #y = as.Date(yday_start, origin = "2019-10-01"),
                              y=as.Date(yearless_date),
                              fill = reorder(admin_unit,start_lat)), outlier.shape=NA) +
  geom_point(data=data, aes(x = reorder(admin_unit,start_lat), # '-' before lat
                            # y = as.Date(yday_start, origin = "2019-10-01")),
                            y=as.Date(yearless_date)),
             position=position_jitter(width=0.15, height=3),
             size=0.25) +
  scale_fill_manual(values=mycolors) +
  scale_y_date(limits=c(as.Date('0000-10-01'), as.Date('0001-01-15')))+ # this worked? or did it cut off jan?
  
  coord_flip() +
  theme_classic() +
  theme(legend.position = "none") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5), 
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12)) +
  labs(y = "Stopover timing", x = "State/Province")
mig.boxplot2

# Modified Plot (with hunting season boxes, where we have the data in table 5)
# Boxes don't want to work unless x is numeric....
mig.boxplot2 <- ggplot (data, group=admin_unit) +
  geom_rect(data=data, aes(xmin=reorder(as.numeric(admin_unit), start_lat),
                           xmax=reorder(as.numeric(admin_unit)+0.3, start_lat),
                          # ymin = as.Date(yday_open, origin = "2019-10-01"), # changing these doesn't do anything
                          # ymax= as.Date(yday_close, origin = "2019-10-01")
                           ymin = as.Date(open),
                           ymax= as.Date(close)
  ), fill='gray')+
  
  geom_rect(data=data, aes(xmin=reorder(as.numeric(admin_unit), start_lat),
                           xmax=reorder(as.numeric(admin_unit)+0.3, start_lat),
                         #  ymin = as.Date(yday_reopen, origin = "2019-10-01"), 
                         #  ymax= as.Date(yday_reclose, origin = "2019-10-01")
                           ymin = as.Date(reopen),
                           ymax= as.Date(reclose)
  ), fill='gray')+
  geom_boxplot(data=data, aes(x = reorder(as.numeric(admin_unit),start_lat), 
                              #y = as.Date(yday_start, origin = "2019-10-01"),
                              y=as.Date(yearless_date),
                              fill = reorder(as.numeric(admin_unit),start_lat)), outlier.shape=NA) +
  geom_point(data=data, aes(x = reorder(as.numeric(admin_unit),start_lat), # '-' before lat
                           # y = as.Date(yday_start, origin = "2019-10-01")),
                           y=as.Date(yearless_date)),
             position=position_jitter(width=0.15, height=3),
             size=0.25) +
  scale_fill_manual(values=mycolors) +
  scale_y_date(limits=c(as.Date('0000-10-01'), as.Date('0001-01-15')))+ # this worked? or did it cut off jan?
 
  coord_flip() +
  theme_classic() +
  theme(legend.position = "none") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5), 
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12)) +
  labs(y = "Stopover timing", x = "State/Province")
mig.boxplot2




# Alex's original plot
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

