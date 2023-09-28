library(tidyverse)
library(sf)
library(mapview)

# NJ
nj <- st_read("nj_hunting_zones/nj_hunting_zones.shp")
mapview(nj)
#uniting areas where hunting is closed with the larger two zone

north_nj_names <- c("North", "Closed to Migratory Bird Hunting (Shark River)", "Closed to Migratory Bird Hunting (Oceanport Creek)", "Closed to Migratory Bird Hunting (Parkers Creek Branch)", "Closed to Migratory Bird Hunting (Manasquan River)")

south_nj_names <- c("South", "Closed to Migratory Bird Hunting (Cox Hall Creek WMA)", "Closed to Migratory Bird Hunting (Barnegat Inlet)",
                    "Closed to Migratory Bird Hunting (Herring Island)")

nj <- nj %>% 
  mutate(district_name = case_when(ZONE_NAME %in% north_nj_names ~ "New Jersey (North Zone)",
                                   ZONE_NAME %in% south_nj_names ~ "New Jersey (South Zone)",
                                   .default = "I'm broken"))

nj_grouped <- nj %>% 
  group_by(district_name) %>% 
  summarise()

st_write(nj_grouped, "nj_hunting_zones/nj_hunting_zones_revised.shp")

qc <- st_read(dsn='quebec_hunting_zones', layer='quebec_hunting_zones')

qc <- qc %>% 
  mutate(district_name = case_when(DISTRICT == "A" ~ "Quebec - District A",
                                   DISTRICT == "B" ~ "Quebec - District B",
                                   DISTRICT %in% c("C", "D", "E", "F") ~ "Quebec - District C-F",
                                   DISTRICT == "G" ~ "Quebec - District G"))

qc_grouped <- qc %>% 
  group_by(district_name) %>% 
  summarise()

st_write(nj_grouped, "quebec_hunting_zones/quebec_hunting_zones_revised.shp")

ont <- st_read(dsn='ontario_hunting_zones', layer='Wildlife_Management_Unit')

ont
# central (Ontario - Southern District H)
# 60:67, 69b

central_names <- c("60", "61", "62", "63A", "63B", "64A", "64B", "65", "66A", "66B", "67",  "69B")

# south (Ontario - Southern District I)
# 68, 69a, 70:95

south_names <- c("68A", "68B", "69A-1", "69A-2", "69A-3", "70", "71", "72A", "72B", "73", "74A", "74B", "75", "76A", "76B", "76C", "76D", "76E", "77A", "77B", "77C", "78A", "78B", "78C", "78D", "78E", "79C", "79D", "80", "81A", "81B", "82A", "82B", "82C", "83A", "83B", "83C", "84", "85A", "85B", "85C", "86A", "86B", "87A", "87B", "87C", "87D", "87E", "88", "89A", "89B", "90A", "90B", "91A", "91B", "92A", "92B", "92C", "92D", "93A", "93B", "93C", "94A", "94B", "95")

# north ("ontario")
# everything else

ont <- ont %>% 
  mutate(district_name = case_when(OFFICIAL_N %in% central_names ~ "Ontario - Southern District H",
                                   OFFICIAL_N %in% south_names ~ "Ontario - Southern District I",
                                   .default = "Ontario"))
ont_grouped <- ont %>% 
  group_by(district_name) %>% 
  summarise()

st_write(ont_grouped, "ontario_hunting_zones/ontario_hunting_zones.shp")

