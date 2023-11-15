# New stopopver analysis
# Fish et al. migration phenology and hunting
# Las edit SJC Nov 2023


fst <- readRDS('fall_stopover_plot_data_100623.rds')
zone_dates <- readRDS('amwo_zone_dates_revised_111523.rds')
brd <- readRDS('birds_list_111523.rds')

head(fst)
head(zone_dates)
head(brd[[1]])

# stopover-specific metrics, including whether or not in hunting season
fst <- left_join(fst, zone_dates, by='admin_unit') 
fst$hunting <- 0
fst$hunting[which(fst$yearless_date>fst$open&fst$yearless_date<fst$close)] <- 1
fst$hunting[which(fst$yearless_date>fst$reopen&fst$yearless_date<fst$reclose)] <- 1

head(fst)
hist(fst$hunting)

# bird data frame with migration-wide metrics
dt <- data.frame(id=rep(NA, length(brd)), mig_duration_days=NA, start_lat=NA, start_lon=NA, end_lat=NA,
                 end_lon=NA, n_stop=NA, stop_hunting_proportion=NA, total_stop_duration_hours=NA, mean_stop_duration_hours=NA)
nrow(dt)
length(brd)
for (i in 1:nrow(dt)){
  dt$id[i] <- brd[[i]]$ID[1]
  dt$mig_duration_days[i] <- as.numeric(brd[[i]]$date[nrow(brd[[i]])]-brd[[i]]$date[1])
  dt$start_lat[i] <- brd[[i]]$lat[1]
  dt$start_lon[i] <- brd[[i]]$lon[1]
  dt$end_lat[i] <- brd[[i]]$lat[nrow(brd[[i]])]
  dt$end_lon[i]  <- brd[[i]]$lon[nrow(brd[[i]])]
  dt$n_stop[i] <- length(which(fst$id==dt$id[i]))
  dt$stop_hunting_proportion[i] <- sum(fst$duration_hours[which(fst$id==dt$id[i]&fst$hunting==1)])/sum(fst$duration_hours[which(fst$id==dt$id[i])])
  dt$total_stop_duration_hours[i] <- sum(fst$duration_hours[which(fst$id==dt$id[i])])
  dt$mean_stop_duration_hours[i] <- mean(fst$duration_hours[which(fst$id==dt$id[i])])
  print(i) # completes 118 total
}

# get rid of birds that didn't terminate
term17 <- read.csv('amwo.term.f17.csv')
term18 <- read.csv('amwo.term.f18.csv')
term19 <- read.csv('amwo.term.f19.csv')
term <- rbind(term17, term18, term19)
dt <- dt[which(dt$id%in%term$ID==T),]

hist(dt$mean_stop_duration_hours)
hist(dt$n_stop)
hist(dt$stop_hunting_proportion)

# save data frame
save(fst, file='stopover_data_for_new_analysis_111523')
save(dt, file='migration_data_for_new_analysis_111523')
