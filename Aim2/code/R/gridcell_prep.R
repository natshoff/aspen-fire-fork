#####################################################################
# Data prep for forest composition and structure modeling in R-INLA #

library(tidyverse)
library(sf)
library(ggcorrplot)
library(lubridate)
library(geosphere)

maindir <- '/Users/max/Library/CloudStorage/OneDrive-Personal/mcook/aspen-fire/Aim2/'

# load and tidy the fire data
# calculate the fire perimeter and fire centroid
fires_fp <- paste0(maindir,"data/spatial/mod/srm_fire_census_2017_to_2023_ics_perims.gpkg")
fires <- st_read(fires_fp) %>%
 st_transform(st_crs(4326)) %>% # make sure its in lat/lon
 mutate(
  fire_perim = st_cast(geom, "MULTILINESTRING")
 ) %>%
 as_tibble() %>%
 rename(
  fire_ig_dt = DISCOVERY_DATE,
  fire_acres = ICS_ACRES,
  fire_aspenpct = Aspen_Pct
 ) %>%
 select(Fire_ID, fire_ig_dt, fire_acres, fire_aspenpct, 
        fire_perim, geom) %>%
 mutate(
  Fire_ID = as.numeric(Fire_ID),
  fire_acres = as.numeric(fire_acres),
  log_fire_size = log(fire_acres),
  fire_aspenpct = as.numeric(fire_aspenpct))


#=========Prep the grid data=========#

# Format the forest composition and structure data frame 
# (long-format) each gridcell has rows for all species present
# each gridcell is assigned a majority forest type (FORTYPCD)
# climate and topography are summarized at the gridcell level
fp <- paste0(maindir,'data/tabular/mod/gridstats_fortypnm_gp_tm_ct_frp-cbi.csv')
grid_tm <-  read_csv(fp)  %>% # read in the file
 select(-c('...1')) %>%
 # join in some of the fire information
 left_join(fires%>%select(-geom), by="Fire_ID") %>%
 # tidy and create attributes ...
 mutate(
  # tidy the temporal fields
  fire_ig_dt = as.Date(fire_ig_dt),  
  fire_year = year(fire_ig_dt),              
  fire_month = month(fire_ig_dt),
  first_obs_date = as.Date(first_obs_date), # the day of first observation
  first_obs_month = month(first_obs_date), # month of first observation
  first_obs_doy = as.numeric(lubridate::yday(first_obs_date)),
  # Create a "season" variable based on the month of first observation
  season = case_when(
   first_obs_month %in% c(3, 4, 5) ~ "spring",
   first_obs_month %in% c(6, 7, 8) ~ "summer",
   first_obs_month %in% c(9, 10, 11) ~ "fall"
  ),
  # Year/season interaction
  year_season = interaction(fire_year, season, drop = TRUE),
  # Factorize the temporal attributes
  season = as.factor(season),
  fire_year = as.factor(fire_year),
  year_season = as.factor(year_season),
  first_obs_date = as.factor(first_obs_date),
  # Format species names consistently
  fortypnm_gp = str_to_lower(fortypnm_gp),
  fortypnm_gp = str_replace_all(fortypnm_gp, "-", "_"),
  fortypnm_gp = str_replace_all(fortypnm_gp, " ", "_"),
  # Format species names consistently
  species_gp_n = str_to_lower(species_gp_n),
  species_gp_n = str_replace_all(species_gp_n, "-", "_"),
  species_gp_n = str_replace_all(species_gp_n, " ", "_"),
  # make sure factor variables are set for species names
  fortypnm_gp = factor(fortypnm_gp),
  species_gp_n = factor(species_gp_n),
  Fire_ID = factor(Fire_ID),
  # Proportion of contributing AFD during daytime observations
  day_prop = day_count / afd_count,
  # log-scaled response variables
  log_frp_max = log(frp_max),
  log_frp_csum = log(frp_csum),
  log_frp_max_day = log(frp_max_day),
  log_frp_max_night = log(frp_max_night),
  log_frp_csum_day = log(frp_csum_day),
  log_frp_csum_night = log(frp_csum_night),
  # scale the percentages
  forest_pct = forest_pct / 100,
  fire_aspenpct = fire_aspenpct / 100,
  # create a 'northness' variable from aspect
  # -1 is south-facing, 1 is north-facing
  aspect_rad = aspect * pi / 180, # convert to radians
  northness = cos(aspect_rad) # northness
 ) %>%
 # fire-level aspen presence
 group_by(Fire_ID) %>%
 mutate(fire_aspen = ifelse(any(fortypnm_gp == "quaking_aspen"), 1, 0)) %>%
 ungroup() %>%
 # gridcell aspen presence
 group_by(grid_idx) %>%
 mutate(grid_aspen = ifelse(any(species_gp_n == "quaking_aspen"), 1, 0)) %>%
 ungroup() %>%
 # rename a few of the columns
 rename(
  dia_live = tree_dia_live,
  ht_live = tree_ht_live,
  tm_canopy = canopypct_mean,
  tm_balive = balive_sum,
  lf_canopy = lf_forest_cc_mean,
  lf_height = lf_forest_ch_mean
 ) %>%
 # aspen-specific metrics
 group_by(grid_idx) %>%
 mutate(
  # sumary stats
  aspen_ba_live = sum(if_else(species_gp_n == "quaking_aspen", ba_live, 0.0), na.rm = TRUE),
  aspen_tpp_live = sum(if_else(species_gp_n == "quaking_aspen", tpp_live, 0.0), na.rm = TRUE),
  aspen_ht_live = mean(if_else(species_gp_n == "quaking_aspen", ht_live, NA_real_), na.rm = TRUE),
  aspen_dia_live = mean(if_else(species_gp_n == "quaking_aspen", dia_live, NA_real_), na.rm = TRUE),
  aspen_qmd_live = mean(if_else(species_gp_n == "quaking_aspen", qmd_live, NA_real_), na.rm = TRUE),
  # proportions 
  aspen_ba_pr = if_else(ba_live_total > 0, aspen_ba_live / ba_live_total, 0.0),
  aspen_tpp_pr = if_else(tpp_live_total > 0, aspen_tpp_live / tpp_live_total, 0.0)
 ) %>%
 ungroup() %>%
 # be sure there are no duplicate rows for grid/fire/species
 distinct(grid_idx, species_gp_n, .keep_all = TRUE) # remove duplicates
glimpse(grid_tm)


############################################
# calculate the distance to fire perimeter #
grid_sf <- grid_tm %>%
 distinct(grid_idx, fortypnm_gp, .keep_all=T) %>%
 select(grid_idx, x, y, fire_perim, Fire_ID) %>%
 st_as_sf(coords = c("x", "y"), crs = 4326) %>%
 # Compute minimum distance from gridcell center to fire perimeter
 rowwise() %>%
 mutate(
  dist_to_perim = as.numeric(min(st_distance(geometry, fire_perim)))
 ) %>%
 ungroup() %>%
 # Normalize distance within each fire
 group_by(Fire_ID) %>%
 mutate(dist_to_perim = dist_to_perim / max(dist_to_perim, na.rm = TRUE)) %>%
 ungroup() %>%
 # Convert back to tibble for easy manipulation
 st_drop_geometry() %>%
 select(-fire_perim) %>% # remove this large geometry field
 select(grid_idx, dist_to_perim)
glimpse(grid_sf)
# join back to the gridcell data
grid_tm <- grid_tm %>%
 select(-fire_perim) %>% # remove this large geometry field
 left_join(., grid_sf, by = "grid_idx")
head(grid_tm%>%select(grid_idx,dist_to_perim))
rm(grid_sf) # tidy up
gc()

###################################################
# check how many grids are majority forested, etc #
print(length(unique(grid_tm$grid_idx))) # unique gridcells
# percent of gridcells majority forested
# Check how many grids are > 50% forested
print(paste0(
 "Percent forested gridcells (>50% forest pixels): ",
 round(dim(grid_tm %>% filter(forest_pct > 0.50) %>% distinct(grid_idx))[1]/
        dim(grid_tm %>% distinct(grid_idx))[1], 3) * 100, "%"
))
# # (optional) retain predominantly forested gridcells ...
# grid_tm <- grid_tm %>%
#  filter(forest_pct >= 0.50)

##############################################################
# Count how many times aspen co-occurs with other forest types
aspen_coo <- grid_tm %>%
 group_by(fortypnm_gp) %>%
 summarise(
  n_total = n(),  # Total occurrences of the forest type
  n_aspen = sum(grid_aspen == 1, na.rm = TRUE),  # How many co-occur with aspen
  aspen_ba_pr_mean = mean(aspen_ba_pr, na.rm=TRUE),
  aspen_ba_pr_p90 = quantile(aspen_ba_pr, probs = 0.9),
  aspen_tpp_pr_mean = mean(aspen_tpp_pr, na.rm=TRUE),
  aspen_tpp_pr_p90 = quantile(aspen_tpp_pr, probs = 0.9),
  pct_aspen_coo = round((n_aspen / n_total) * 100, 1)  # Convert to percentage
 ) %>%
 arrange(desc(pct_aspen_coo))  # Sort by co-occurrence percentage
print(aspen_coo)
# save the file
write_csv(aspen_coo, paste0(maindir,"data/tabular/mod/results/aspen_co-occurrence.csv"))
rm(aspen_coo)

##############################
# Flag duplicate gridcells ######
# These are mostly re-burns ... #
# flag duplicates in the data (reburns)
grid_tm <- grid_tm %>%
 # flag re-burn gridcells
 group_by(grid_index) %>%
 mutate(
  # Get earliest fire date for this gridcell
  first_fire_dt = min(fire_ig_dt, na.rm = TRUE),
  # Flag reburn: if current fire date is after earliest
  reburn = as.factor(fire_ig_dt > first_fire_dt)
 ) %>%
 ungroup()
print(dim(grid_tm%>%filter(reburn == "TRUE")))
##################
# plot CBI re-burns
grid_tm %>%
 distinct(grid_idx, reburn, CBIbc_p90) %>%
 ggplot(aes(x = reburn, y = CBIbc_p90)) +
 geom_boxplot() +
 labs(x = "Reburned", y = "90th Percentile CBI") +
 theme_minimal()
# plot FRP re-burns
grid_tm %>%
 distinct(grid_idx, reburn, log_frp_csum) %>%
 ggplot(aes(x = reburn, y = log_frp_csum)) +
 geom_boxplot() +
 labs(x = "Reburned", y = "Cumulative FRP") +
 theme_minimal()

###############################
# check on the species counts #
grid_tm %>% 
 group_by(species_gp_n) %>%
 summarise(n = length(species_gp_n)) %>%
 arrange(desc(n))

# check how many grids and fires
length(unique(grid_tm$grid_idx))
length(unique(grid_tm$Fire_ID))


###############################
# Save a fire perimeter dataset
fires <- fires %>% 
 filter(Fire_ID %in% unique(da$Fire_ID)) %>%
 select(-fire_perim) %>%
 st_as_sf() %>%
 st_transform(st_crs(5070))
st_write(fires,paste0(maindir,"data/spatial/mod/srm_fire_census_model_data.gpkg"),
         append=F)
# tidy up
rm(fires)


###################################################
# save the model data frame (tabular and spatial) #
# tabular data (model data)
out_fp = paste0(maindir,"data/tabular/mod/model_data_cleaned.csv")
write_csv(grid_tm, out_fp)
# spatial data (keep one row per gridcell)
grid_tm.w <- grid_tm %>%
 select(
  Fire_ID, Fire_Name, grid_idx, fortypnm_gp, species_gp_n, ba_live_pr,
  log_frp_csum, log_frp_csum_day, log_frp_csum_night,
  log_frp_max, log_frp_max_day, log_frp_max_night,
  CBIbc_mean, CBIbc_p90, CBIbc_p95, CBIbc_p99,
  fortyp_pct, forest_pct, lf_canopy, lf_height, H_tpp, H_ba,
  first_obs_doy, erc_dv, vpd, vs, elev, slope, aspect, tpi,
  ba_live_total, tpp_live_total, qmd_live_mean,
  tree_ht_live_mean, tree_dia_live_mean, 
  aspen_ba_live, aspen_tpp_live, aspen_ht_live, aspen_dia_live,  aspen_qmd_live,
  aspen_ba_pr, aspen_tpp_pr, dist_to_perim, overlap, day_prop, tmid_n 
 ) %>%
 pivot_wider(
  names_from = species_gp_n,
  values_from = ba_live_pr,
  values_fill = 0
 ) %>%
 left_join(., grid_tm %>% distinct(grid_idx, geometry), by="grid_idx") %>%
 distinct(grid_idx, fortypnm_gp, .keep_all = T) %>%
 mutate(geometry = st_as_sfc(geometry)) %>%
 st_as_sf(.) %>% 
 st_set_crs(st_crs(5070))
glimpse(grid_tm.w)
# save it out.
out_fp = paste0(maindir,"data/spatial/mod/gridcell_model_data_cleaned.gpkg")
st_write(grid_tm.w, out_fp, append=F)
rm(grid_tm.w)


#===============Explore Distributions, etc.================#

####################################
# distribution of response variables
resp_plot <- grid_tm %>%
 # pivot longer to facet plot
 pivot_longer(cols = c(log_frp_csum,
                       CBIbc_p90),
              names_to = "variable",
              values_to = "value") %>%
 # Plot with facets
 ggplot(aes(x = value)) +
 geom_histogram(bins = 30, fill = "orange", alpha = 0.7) +
 facet_wrap(
  ~ variable,
  scales = "free",
  labeller = as_labeller(c(log_frp_csum = "log(Cumulative FRP)",
                           CBIbc_p90 = "90th percentile CBIbc"))) +
 labs(x = "value",
      y = "Frequency") +
 theme_classic()
resp_plot
# save the plot.
out_png <- paste0(maindir,'figures/INLA_ResponseDistribution_FRP-CBI.png')
ggsave(out_png, plot = resp_plot, dpi=500, bg = 'white')
rm(resp_plot)

#==========EXPLORE CORRELATIONS==========#

#######################################
# species-specific correlation matrix #
sp_cor <- grid_tm %>%
 select(grid_idx, species_gp_n, ba_live) %>%
 spread(species_gp_n, ba_live, fill = 0)  # Reshape to wide format
# compute the correlation matrix
sp_cormat <- cor(sp_cor[,-1], use = "complete.obs", method = "spearman")
ggcorrplot(sp_cormat, method = "circle", type = "lower",
           lab = TRUE, lab_size = 3, colors = c("blue", "white", "red"))
rm(sp_cor, sp_cormat)
# save the plot.
out_png <- paste0(maindir,'figures/INLA_CorrelationMatrix_SpeciesBA.png')
ggsave(out_png, dpi=500, bg = 'white')
rm(sp_cor,sp_cormat)

# ########################################
# # correlation matrix for fixed effects #
# cor_da <- grid_tm %>%
#  mutate(
#   log_fire_size = log(fire_acres)
#  ) %>%
#  select(
#   fortypnm_gp,
#   # species structure metrics
#   tpp_live, ba_live, qmd_live, ht_live, dia_live, hdr_live, ba_per,
#   tpp_live_total, ba_live_total, qmd_live_mean,
#   tpp_dead_total, ba_dead_total, # dead BA and TPP
#   tpp_live_pr, ba_live_pr, # proportions of TPP and BA
#   forest_pct, fortyp_pct, # gridcell forest percent and majority forest type percent
#   H_ba, H_tpp, # gridcell species diversity (abundance- and dominance-based)
#   erc, erc_dv, vpd, vpd_dv, # day-of-burn climate
#   fm1000, rmin, tmmx, vs, # day-of-burn climate
#   elev, slope, tpi, northness,  # topography
#   lf_canopy, lf_height, # gridcell canopy percent/height and BA sum
#   day_prop, overlap, # VIIRS detection characteristics
#   log_fire_size, dist_to_perim # fire size and distance to perimeter
#  ) %>%
#  pivot_wider(
#   names_from = fortypnm_gp,
#   values_from = fortyp_pct,
#   values_fill = 0) %>%
#  mutate(across(everything(), ~ scale(.) %>% as.numeric()))  # Standardize variables
# 
# # Compute correlation matrix
# cor_mat <- cor(cor_da, use = "complete.obs", method = "spearman")
# # Plot correlation matrix
# ggcorrplot(cor_mat, method = "circle",
#  type = "lower", lab = TRUE, lab_size = 3, tl.cex = 10,
#  colors = c("blue", "white", "red")
# )
# rm(cor_da, cor_mat) # tidy up
# # save the plot.
# out_png <- paste0(maindir,'figures/INLA_CorrelationMatrix_FixedEffects.png')
# ggsave(out_png, dpi=500, width=12, height=12, bg = 'white')
# 









