
#
# Script to tidy/prepare the data
# For the Southern Rocky Mountain ecoregion (EPA Level III)
# 

library(tidyverse)
library(sf)

getwd()

proj <- st_crs("ESRI:102039") # projection

# Environments
datadir <- '/Users/max/Library/CloudStorage/OneDrive-Personal/mcook/data/'


###############
# Load the data

# SRME boundary
srme <- st_read(paste0(datadir,'boundaries/ecological/ecoregion/na_cec_eco_l3.gpkg')) %>%
 filter(NA_L3NAME == 'Southern Rockies') %>%
 st_union() %>% st_as_sf() %>%
 st_transform(proj)

# MTBS perimeter data (2017-2021)
mtbs <- st_read('data/spatial/raw/MTBS/mtbs_perims_west_2017to2021.shp') %>%
 # First create the date/year field
 mutate(Ig_Date = as.Date(Ig_Date),
        Ig_Year = lubridate::year(Ig_Date)) %>%
 # For now, only fires from 2019 onward
 filter(
  Ig_Year >= 2019,
  Incid_Type != 'Prescribed Fire'
 ) %>%
 st_transform(proj) %>% 
 st_intersection(srme) # SRME fires

# FRP data, within MTBS
frp <- st_read('data/spatial/raw/VIIRS/DL_FIRE_SV-C2_361441/fire_archive_SV-C2_361441.shp') %>%
 st_transform(proj) %>%
 st_intersection(st_buffer(mtbs, 1000)) %>% # buffer fire perimeters 1km
 # Only retain nominal confidence observations
 filter(CONFIDENCE != 'l')

gc() # garbage clean

# Write out the files
st_write(srme,'data/spatial/mod/boundaries/na_cec_eco_l3_srme.gpkg',delete_dsn=T)
st_write(mtbs,'data/spatial/mod/MTBS/mtbs_perims_srme_2019to2021.gpkg',delete_dsn=T)
st_write(mtbs,'data/gee/imports/mtbs_perims_srme_2019to2021.shp',delete_dsn=T) # for GEE
st_write(frp,'data/spatial/mod/VIIRS/fire_archive_SV-C2_srme_2019to2021.gpkg',delete_dsn=T)
st_write(frp,'data/gee/imports/fire_archive_SV-C2_srme_2019to2021.shp',delete_dsn=T) # for GEE

