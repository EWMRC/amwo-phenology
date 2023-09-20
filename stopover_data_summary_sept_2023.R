# Fall Migration Stopover Data summary
# Fish et al. Phenology paper
# September 2023
# Last Edit - SJC Sept 2023

library(leaflet)
library(tidyverse)

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
View(data)

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
# I don't think it's worth trying to make this run-length encoding thing in tidyverse becasue it depends on separating 
# so many different objects.
for (q in 1:length(brd)){
  tmp <- rle(brd[[q]]$stopover) # run length encoding based on distance, clusters by whether or not the next point is >16km away
  label <- seq(1:length(tmp$lengths)) # makes a number to label each cluster with
  count <- rep(0, times=length(tmp$lengths)) # empty column to fill count
  for (i in 1:length(count)){
    count[i] <- sum(tmp$lengths[1:i]) # count (number of points so far in the data frame) is there to be able to label clusters
  }
  lt <- data.frame(lengths=tmp$lengths,values=tmp$values, count=count, label=label) # makes a new table for each bird with information for assigning labels
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
  
  print(q) # completes 118 total birds
}

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
lapply(brd[1:12], FUN=stopover_id_map)
lapply(brd[13:24], FUN=stopover_id_map)
lapply(brd[25:36], FUN=stopover_id_map)
lapply(brd[37:48], FUN=stopover_id_map)
lapply(brd[49:60], FUN=stopover_id_map)
lapply(brd[61:72], FUN=stopover_id_map)
lapply(brd[73:84], FUN=stopover_id_map)
lapply(brd[85:96], FUN=stopover_id_map)
lapply(brd[98:109], FUN=stopover_id_map)
lapply(brd[110:118], FUN=stopover_id_map)

# example of recursive - NY-2018-06
which(ids=='NY-2018-06') #11
stopover_id_map(brd[[11]])


