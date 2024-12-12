# Load the required libraries
library(tidyverse)
library(sf) # spatial
library(INLA) # for spatial Bayes model
library(ggcorrplot)

# Environment variables
maindir <- '/Users/max/Library/CloudStorage/OneDrive-Personal/mcook/aspen-fire/Aim2/'

#=========Prep the grid data=========#

# Format the species composition data frame

# load the spatial grid
fp <- paste0(maindir,'data/spatial/mod/VIIRS/viirs_snpp_jpss1_afd_latlon_aspenfires_pixar_gridstats.gpkg')
grid <- st_read(fp) %>%
 select(c(grid_index, geom))

# load the aggregated FRP grid with TreeMap and climate/topography
fp <- paste0(maindir,'data/tabular/mod/viirs_snpp_jpss1_gridstats_fortypcd_climtopo.csv')
grid_fortyp <-  read_csv(fp) %>%
 # select the required columns
 select(c(grid_index, Fire_ID, frp_csum, frp_max, 
          grid_x, grid_y, SpeciesName, spp_pct, forest_pct,
          erc, erc_dv, vpd, vpd_dv, elev, slope, chili, tpi))
glimpse(grid_fortyp)

# join to the spatial data
grid_fortyp_sp <- inner_join(grid, grid_fortyp, by="grid_index")
head(grid_fortyp_sp)

# tidy the columns
grid_ <- grid_fortyp_sp %>%
 # remove missing FRP, prep columns
 filter(frp_max > 0) %>% # make sure FRP is not 0
 mutate(Fire_ID = as.factor(Fire_ID)) %>%
 distinct(grid_index, Fire_ID, SpeciesName, .keep_all = TRUE) # remove duplicates
head(grid_) # check the results

# reshape the data frame for modeling
grid_w <- grid_ %>% 
 # tidy the species names
 mutate(SpeciesName = str_replace_all(SpeciesName, "-", "_"),
        SpeciesName = str_to_lower(SpeciesName),
        spp_pct = as.numeric(spp_pct)) %>%
 pivot_wider(
  names_from = SpeciesName, 
  values_from = spp_pct, 
  values_fill = 0) %>% # pivot wider
 filter(frp_max > 0) %>%
 mutate(log_frp_max = log(frp_max + 1))
head(grid_w)
nrow(grid_w)

# create a conifer column
grid_w <- grid_w %>%
 as_tibble() %>%
 mutate(conifer = rowSums(select(., douglas_fir, lodgepole, ponderosa, spruce_fir, piñon_juniper), na.rm = TRUE))

# retain grid cells with some aspen component
grid_aspen <- grid_w %>%
 filter(aspen > 0)
nrow(grid_aspen)/nrow(grid_w)*100
rm(grid_)
gc()

