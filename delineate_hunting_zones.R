library(tidyverse)
library(sf)
library(mapview)

# NJ
nj <- st_read("nj_hunting_zones/nj_hunting_zones.shp")
mapview(nj)
#uniting areas where hunting is closed with the larger two zones

nj %>% 
  mutate()