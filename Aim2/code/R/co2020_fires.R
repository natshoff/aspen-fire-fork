
# Summary of aspen/wildfire for the 2020 wildfire (for press)

library(tidyverse)
library(sf)

frp <- st_read('data/spatial/mod/VIIRS/fire_archive_SV-C2_srme_19to21_plots_w_attr.gpkg') %>%
 # tidy the frame
 mutate(DAYNIGHT = as.factor(DAYNIGHT)) %>%
 rename(Aspen = Aspen_S2,
        MixedConifer = Conifer_Hardwood) %>%
 filter(ACQ_YEAR <= 2020) %>%
 # Convert back to centroid (points) and ensure correct project
 st_centroid() %>% st_transform(proj)
