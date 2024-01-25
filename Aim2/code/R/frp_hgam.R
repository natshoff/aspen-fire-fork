
# Libraries

library(tidyverse)
library(sf)
library(mgcv)

# Data

# Fire Radiative Power (FRP) obs. w/ attributes
frp <- st_read('data/spatial/mod/vnp14img_west_spatial_w_attr.gpkg') %>%
 # tidy the frame
 mutate(daynight = as.factor(daynight),
        acq_year = lubridate::year(acq_date),
        aspen = if_else(aspen_pct > 0, 1, 0)) %>%
 # Convert back to centroid (points) and ensure correct project
 st_centroid() %>% 
 st_transform(proj)

