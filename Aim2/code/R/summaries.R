
##
# Summarize attributes within fires and FRP observations
##

library(tidyverse)
library(sf)


# Bring in the data
srme <- st_read('data/spatial/mod/boundaries/na_cec_eco_l3_srme.gpkg')
mtbs <- st_read('data/spatial/mod/MTBS/mtbs_perims_srme_2019to2021.gpkg')
frp <- st_read('data/spatial/mod/VIIRS/fire_archive_SV-C2_srme_2019to2021_plots.gpkg') %>%
 filter(Incid_Type != "Prescribed Fire") # only keep wildfires


####
# Prep the MTBS perimeters
####
mtbs_aspen <- read_csv('data/tabular/mod/aspen/mtbs_aspen_histo.csv') %>%
 # Calculate % aspen within fire perimeters
 mutate(Aspen_LF = aspen1*100 / (BurnBndAc * 4046.86) * 100) %>%
 select(Event_ID, Incid_Name, Aspen_LF)
# Filter fires w/ > 5% aspen
mtbs <- mtbs %>%
 left_join(mtbs_aspen, by="Event_ID") %>%
 rename(Incid_Name = Incid_Name.x) %>%
 select(-Incid_Name.y)
rm(mtbs_aspen)

# Bring in the TreeMap 2016 aspen summary
treemap <- read_csv("data/tabular/mod/TreeMap/treemap16_MTBS_west_foresttype.csv") %>%
 filter(ForTypName == 'Aspen') %>%
 mutate(
  ForTypM2 = Count * 900,
  # Acres
  ForTypAc = ForTypM2 / 4047,
  # Binary aspen presence
  AspenNoAspen = if_else(ForTypName == "Aspen", 1, 0)
 ) 
head(treemap)

# Join to MTBS
# filter fires where either has > 5%
mtbs <- mtbs %>%
 left_join(treemap,by="Event_ID") %>%
 mutate(Aspen_TM16 = ForTypAc / BurnBndAc * 100) %>%
 filter(Aspen_TM16 >= 5 | Aspen_TM16 > 5)
glimpse(mtbs)

rm(treemap) # clean up

####
# Prep the FRP observations
####

# Load the FRP summaries
 frp_aspen <- read_csv('data/tabular/mod/aspen/vnp14img_aspen_histo.csv') %>%
 rename(uid = fid) %>%
 select(uid, aspen_1)
# Join aspen pct to FRP observations
frp <- frp %>% 
 left_join(frp_aspen, by="uid") %>%
 mutate(Aspen_S2 = aspen_1*100 / area * 100) %>%
 select(-aspen_1) %>%
 filter(Event_ID %in% mtbs$Event_ID)
rm(frp_aspen)
glimpse(frp)


# Join LANDFIRE EVT summaries
lfevt <- read_csv('data/tabular/mod/LFEVT16/vnp14img_lf_evt_histo_phys_.csv') %>%
 rename(
  Hardwood = PHYS_2,
  Conifer = PHYS_3,
  Conifer_Hardwood = PHYS_4,
  Shrubland = PHYS_5,
  Grassland = PHYS_6,
  Riparian = PHYS_7
 ) %>%
 select(uid, Conifer, Conifer_Hardwood, Hardwood, Shrubland, Grassland, Riparian) %>%
 pivot_longer(-uid, names_to = "Class") %>%
 rename(Pixels = value) %>%
 mutate(AreaM2 = Pixels * 900,
        AreaM2 = if_else(AreaM2 > 141376, 141376, AreaM2),
        AreaPct = AreaM2 / 141376 * 100) %>%
 pivot_wider(id_cols=uid, names_from=Class, values_from=AreaPct)
head(lfevt,20)

# Join to the FRP obs.
frp <- frp %>%
 left_join(lfevt, by="uid") %>%
 # Tidy the data frame
 select(uid,Event_ID,Incid_Name,BurnBndAc,Ig_Date,Ig_Year,
        LATITUDE, LONGITUDE, BRIGHTNESS, FRP, ACQ_DATE, ACQ_TIME,
        CONFIDENCE, DAYNIGHT, Aspen_S2, Conifer, Shrubland, Grassland,
        Riparian, Hardwood, Conifer_Hardwood) %>%
 mutate(ACQ_YEAR = lubridate::year(ACQ_DATE)) %>%
 filter(ACQ_YEAR >= 2019)
glimpse(frp)

rm(lfevt)

# Export the tidy'd data frame and others
st_write(frp, 'data/spatial/mod/VIIRS/fire_archive_SV-C2_srme_19to21_plots_w_attr.gpkg',delete_dsn = TRUE)
st_write(mtbs, 'data/spatial/mod/MTBS/mtbs_perims_srme_2019to2021_wAspen.gpkg',delete_dsn = TRUE)

# Export a simplified version for GEE
st_write(
 frp %>% select(uid,ACQ_DATE,geom),
 'data/gee/imports/fire_archive_SV-C2_srme_19to21_plots_w_attr.shp',
 delete_dsn = TRUE
)
