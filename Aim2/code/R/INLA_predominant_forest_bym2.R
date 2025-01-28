
# Load the required libraries
library(tidyverse)
library(sf) # spatial
library(INLA) # for spatial Bayes model
library(ggcorrplot)
library(lubridate)
library(ggridges)
library(reshape2)
library(spdep)
library(patchwork)
library(forcats)

# Environment variables
maindir <- '/Users/max/Library/CloudStorage/OneDrive-Personal/mcook/aspen-fire/Aim2/'

# load the fire data
fires <- paste0(maindir,"data/spatial/mod/srm_fire_census_2017_to_2023_ics_perims.gpkg")
fires <- st_read(fires) %>%
 as_tibble() %>%
 rename(
  fire_ig_dt = DISCOVERY_DATE,
  fire_acres = ICS_ACRES,
  fire_name = Fire_Name,
  fire_aspenpct = Aspen_Pct
 ) %>%
 select(Fire_ID, fire_name, fire_ig_dt, fire_acres, fire_aspenpct) %>%
 mutate(
  Fire_ID = as.numeric(Fire_ID),
  fire_acres = as.numeric(fire_acres),
  fire_aspenpct = as.numeric(fire_aspenpct))
head(fires)


#=========Prep the grid data=========#

# Format the species composition data frame 
# (long-format) each grid has rows for species co-occurring with the dominant type
# climate and topography are summarized at the grid level

# load the aggregated FRP grid with TreeMap and climate/topography
fp <- paste0(maindir,'data/tabular/mod/gridstats_fortypnm_gp_tm_ct_frp-cbi.csv')
grid_tm <-  read_csv(fp) %>% # read in the file
 # join in some of the fire information
 left_join(fires, by="Fire_ID") %>%
 # filter out grids with no daytime observations
 filter(
  day_count > 0,
  CBIbc_p90 > 0 # and CBI 0 -> these are often boundary grids
 ) %>% 
 # create a numeric fire ID
 mutate(
  # create a new unique identifier
  grid_idx = paste0(Fire_ID, grid_index),
  # tidy the temporal fields
  fire_ig_dt = as.Date(fire_ig_dt),  
  fire_year = year(fire_ig_dt),              
  fire_month = month(fire_ig_dt),
  # Create a "season" variable based on the month
  season = case_when(
   fire_month %in% c(3, 4, 5) ~ "spring",
   fire_month %in% c(6, 7, 8) ~ "summer",
   fire_month %in% c(9, 10, 11) ~ "fall"
  ),
  # Year/season interaction
  year_season = interaction(fire_year, season, drop = TRUE),
  # Factorize the temporal attributes
  season = as.factor(season),
  fire_year = as.factor(fire_year),
  year_season = as.factor(year_season),
  # Format species names consistently
  fortypnm_gp = str_to_lower(fortypnm_gp),
  fortypnm_gp = str_replace_all(fortypnm_gp, "-", "_"),
  fortypnm_gp = str_replace_all(fortypnm_gp, " ", "_"),
  # make sure factor variables are set for species names
  fortypnm_gp = as.factor(fortypnm_gp),
  species_gp_n = as.factor(species_gp_n),
  Fire_ID = as.factor(Fire_ID),
  grid_idx = as.numeric(as.factor(grid_idx)),
  # proportion of contributing AFD during daytime observations
  day_prop = day_count / afd_count,
  # log-scaled outcomes
  log_frp_max = log(frp_max + 1e-5),
  log_frp_csum = log(frp_csum + 1e-5),
  log_frp_max_day = log(frp_max_day + 1e-5),
  log_frp_max_night = log(frp_max_night + 1e-5),
  log_frp_csum_day = log(frp_csum_day + 1e-5),
  log_frp_csum_night = log(frp_csum_night + 1e-5),
  log_frp_p90 = log(frp_p90_day + 1e-5),
  log_frp_p95 = log(frp_p95_day + 1e-5),
  # scale the percentages
  fortyp_pct = fortyp_pct / 100,
  forest_pct = forest_pct / 100,
  # create a species presence flag
  presence = factor(if_else(ba_live > 0 | tpp_live > 0, 1, 0))) %>%
 # be sure there are no duplicate rows for grid/fire/species
 distinct(grid_idx, fortypnm_gp, .keep_all = TRUE) # remove duplicates
glimpse(grid_tm) # check the results

rm(fires)

# Check how many grids are > 50% forested
print("Number of forested grids (>50% forest pixels):")
dim(grid_tm %>% filter(forest_pct > 0.50))[1]/dim(grid_tm)[1]

# Remove duplicate "grid_index"
# These are reburns ...
duplicate_grids <- grid_tm %>%
 group_by(grid_index) %>%
 filter(n() > 1) %>%
 arrange(grid_index, Fire_ID) 
# save a spatial polygon 
grid <- st_read(paste0(maindir,"data/spatial/mod/VIIRS/viirs_snpp_jpss1_afd_latlon_fires_pixar_gridstats.gpkg"))
grid <- grid %>%
 filter(grid_index %in% unique(duplicate_grids$grid_index)) %>%
 st_as_sf(.) %>% st_transform(st_crs(grid))
st_write(grid,paste0(maindir,"data/spatial/mod/grid_model_data_duplicates.gpkg"), append=F)
rm(grid)
gc()

# filter these out
idx <- duplicate_grids %>% pull(grid_index)
grid_tm <- grid_tm %>% 
 filter(!grid_index %in% idx)
rm(idx, duplicate_grids)


#===============Explore Distributions, etc.================#

####################################
# distribution of response variables
grid_tm %>%
 # pivot longer to facet plot
 pivot_longer(cols = c(log_frp_max_day, log_frp_csum, 
                       CBIbc_p90, CBIbc_p99),
              names_to = "variable",
              values_to = "value") %>%
 # Plot with facets
 ggplot(aes(x = value)) +
 geom_histogram(bins = 30, fill = "orange", alpha = 0.7) +
 facet_wrap(
  ~ variable, 
  scales = "free",
  labeller = as_labeller(c(log_frp_max_day = "log(Daytime Max FRP)",
                           log_frp_csum = "log(Cumulative FRP)",
                           CBIbc_p90 = "90th Percentile CBIbc",
                           CBIbc_p99 = "99th Percentile CBIbc"))) +
 labs(x = "value",
      y = "Frequency") +
 theme_minimal()

# save the plot.
out_png <- paste0(maindir,'figures/INLA_ResponseDistribution.png')
ggsave(out_png, dpi=500, bg = 'white')


####################
# correlation matrix
# Select only numeric columns and convert factors to dummy variables
# this correlation matrix is for this simple model (i.e., not forest composition)
cor_da <- grid_tm %>%
 select(
  fortypnm_gp, # forest type (factor)
  forest_pct, fortyp_pct, # forest and forest type percent
  canopypct_mean, balive_sum, # grid-level mean canopy percent and balive
  erc, erc_dv, vpd, vpd_dv, # climate
  fm1000, fm1000_dv, rmin, rmin_dv, # climate
  tmmx, tmmx_dv, vs, vs_dv, #climate
  elev, slope, tpi, chili  # topography
 ) %>%
 pivot_wider(
  names_from = fortypnm_gp, 
  values_from = fortyp_pct, 
  values_fill = 0) %>%
 mutate(across(everything(), ~ scale(.) %>% as.numeric()))  # Standardize variables

# Compute correlation matrix
cor_mat <- cor(cor_da, use = "complete.obs", method = "spearman")

# Plot correlation matrix
ggcorrplot(
 cor_mat,
 method = "circle",  # Circle or square for visualization
 type = "lower",  # Lower triangle of the correlation matrix
 lab = TRUE,  # Show correlation values
 lab_size = 3,
 tl.cex = 10,  # Text label size
 colors = c("blue", "white", "red")  # Color gradient
)

rm(cor_da, cor_mat)
gc()

# save the plot.
out_png <- paste0(maindir,'figures/INLA_CorrelationMatrix_Fortyp.png')
ggsave(out_png, dpi=500, bg = 'white')



#===========MODEL SETUP==============#

# list of species names
spps <- c("quaking_aspen", "mixed_conifer", "lodgepole", "ponderosa", 
          "spruce_fir", "piñon_juniper", "oak_woodland")

# force aspen to be the baseline
grid_tm <- grid_tm %>%
 mutate(
  fortypnm_gp = fct_relevel(fortypnm_gp, spps)
 )

# check the factor levels
# make sure aspen is first
levels(grid_tm$fortypnm_gp)

# create a data frame for just the predominant forest type
# one row per grid_index with the dominant forest type
# filter to where that type is at least 50% of the grid forested area
da <- grid_tm %>%
 # keep only grids where predominant species cover > 50%
 filter(forest_pct > 0.50) %>%
 # select the columns we need for modeling
 select(grid_index, grid_idx, Fire_ID, # ID columns
        fire_year, season, year_season, # temporal effects
        log_frp_max_day, log_frp_csum, log_frp_csum_day, # FRP response variables
        log_frp_p90, log_frp_p95,
        CBIbc_p90, CBIbc_p95, CBIbc_p99, CBIbc_mean, # CBI response variables
        fortypnm_gp, fortyp_pct, forest_pct, # forest type and percent cover
        canopypct_mean, balive_sum, # canopy percent and total live basal area
        erc, erc_dv, vpd, vpd_dv, # climate
        fm1000, fm1000_dv, rmin, rmin_dv, # climate
        tmmx, tmmx_dv, vs, vs_dv, #climate
        elev, slope, tpi, chili, # topography
        presence, # species presence (binary factor)
        x, y # grid centroid coordinate for spatial fields model
 ) %>%
 # center and scale continuous predictor variables
 mutate(
  across(
   c(
    forest_pct, fortyp_pct, # forest and forest type percent
    canopypct_mean, balive_sum, # canopy percent and total live basal area
    erc, vpd, vpd_dv, erc_dv, # climate
    fm1000, fm1000_dv, rmin, rmin_dv, # climate
    tmmx, tmmx_dv, vs, vs_dv, #climate
    elev, slope, tpi, chili # topography
   ),
   ~ as.numeric(scale(.))
  )
 ) %>%
 # rename the response variables
 rename(
  frp_day = log_frp_max_day,
  frp_csum = log_frp_csum,
  frp_csum_day = log_frp_csum_day,
  frp_p90 = log_frp_p90,
  frp_p95 = log_frp_p95,
  cbi_p90 = CBIbc_p90,
  cbi_p95 = CBIbc_p95,
  cbi_p99 = CBIbc_p99,
  cbi_mn = CBIbc_mean,
  fire_id = Fire_ID,
  canopypct = canopypct_mean,
  balive = balive_sum
 ) %>%
 # keep just one row per grid cell by predominant type
 distinct(grid_idx, fortypnm_gp, .keep_all = TRUE) %>%
 arrange(grid_idx)

########################################################
# Check on the grid cell counts for daytime observations
# After filtering, make sure we have enough grids per fire
grid_counts <- da %>%
 distinct(fire_id, grid_idx) %>% # keep only distinct rows
 group_by(fire_id) %>%
 summarise(n_grids = n()) %>%
 ungroup() %>% as_tibble()
summary(grid_counts$n_grids)
qt <- quantile(grid_counts$n_grids, probs = seq(.1, .9, by = .1))
qt

# how many fires have >= 10 grids?
dim(grid_counts %>% filter(n_grids >= 10))[1]
# how many fires have >= 50 grids?
dim(grid_counts %>% filter(n_grids >= 50))[1]
# how many fires have >= 100 grids?
dim(grid_counts %>% filter(n_grids >= 100))[1]

# Identify fires with n_grids below the Nth percentile
idx <- grid_counts %>%
 filter(n_grids < 10) %>%
 pull(fire_id)
length(idx)

# filter the data frame to remove these fires
da <- da %>%
 filter(!fire_id %in% idx)

# check how many grids and fires
length(unique(da$grid_idx))
length(unique(da$fire_id))

# check how many grids of each forest type there are
da %>%
 group_by(fortypnm_gp) %>%
 summarize(n = n()) %>%
 ungroup()

# check on counts by season
da %>%
 group_by(year_season) %>%
 summarize(n = n()) %>%
 ungroup()

# Subset species with too few observations
drop_spps = c("rocky_mountain_juniper", "oak_woodland")
print("Dropping small classes:")
print(dim(da%>%filter(fortypnm_gp%in%drop_spps)))
da <- da %>%
 filter(!fortypnm_gp %in% drop_spps) %>%
 mutate(fortypnm_gp = droplevels(fortypnm_gp))

# Tidy up!
rm(grid_tm, grid_counts, idx)
gc()

# save a spatial polygon 
grid <- st_read(paste0(maindir,"data/spatial/mod/VIIRS/viirs_snpp_jpss1_afd_latlon_fires_pixar_gridstats.gpkg"))
grid <- da %>%
 left_join(grid %>% select(grid_index, geom), by="grid_index") %>%
 st_as_sf(.)
st_write(grid,paste0(maindir,"data/spatial/mod/grid_model_data.gpkg"), append=F)
rm(grid)
gc()


#===========MODEL FITTING==============#

set.seed(435)

########################
# FIRE RADIATIVE POWER #

#######################################
# 1. Baseline model (no latent effects)

# Set up the model formula
mf.frp <- frp_csum ~ 
 fortypnm_gp + # predominant forest type
 forest_pct + # percent of grid of predominant forest type
 erc + tmmx_dv + vs + # climate/weather
 slope + tpi + chili # topography
# fit the model                     
model_bl.frp <- inla(
 mf.frp, data = da, 
 family = "gaussian",
 control.predictor = list(compute = TRUE),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
 control.fixed = list(
  prec = list(prior = "pc.prec", param = c(1, 0.01)) 
 ) 
)
summary(model_bl.frp)

##############################################
# 2. Baseline model + fire-level random effect

# specify a pc-prior for the random effect
hyper_pr.fire <- list(
 prec = list(prior = "pc.prec", param = c(1, 0.01))
)
# update the model formula
mf.frp.re <- update(
 mf.frp, . ~ 1 + . + 
  f(fire_id, model = "iid", hyper = hyper_pr.fire) # fire-level random effect
)
# fit the model                     
model_bl.frp.re <- inla(
 mf.frp.re, data = da, 
 family = "gaussian",
 control.predictor = list(compute = TRUE),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
 control.fixed = list(
  prec = list(prior = "pc.prec", param = c(1, 0.01)) 
 )   
)
summary(model_bl.frp.re)

#######################################################################
# 3. Baseline model + fire-level random effect + temporal random effect

# specify a pc-prior for the random effect
hyper_pr.time <- list(
 prec = list(prior = "pc.prec", param = c(1, 0.1))
)
# update the model formula
mf.frp.re2 <- update(
 mf.frp.re, . ~ 1 + . + 
  f(year_season, model = "iid", hyper = hyper_pr.time) # temporal random effect
)
# fit the model                     
model_bl.frp.re2 <- inla(
 mf.frp.re2, data = da, 
 family = "gaussian",
 control.predictor = list(compute = TRUE),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
 control.fixed = list(
  prec = list(prior = "pc.prec", param = c(1, 0.1)) 
 )   
)
summary(model_bl.frp.re2)

###########################################
# 4. Baseline model + spatial random effect

############################################
# create the adjacency structure
# spatial points representing grid centroids
grid_sf <- da %>%
 arrange(grid_index) %>%
 distinct(grid_index, x, y, .keep_all = TRUE) %>%
 st_as_sf(., coords = c("x", "y"), crs = 4326) %>%
 mutate(grid_index = as.numeric(as.factor(grid_index)))
coords <- grid_sf %>% st_coordinates(.)
st_write(grid_sf, paste0(maindir,"data/spatial/mod/model_grid_centroids.gpkg"))

# k-nearest neighbor and symmetric matrix
nbs <- knearneigh(coords, k = 7, longlat = T)
nbs <- knn2nb(nbs, row.names = grid_sf$grid_index, sym = T) #force symmetry!!

# plot the adjacency structure for one fire
fire_sf <- grid_sf %>% filter(fire_id == 51) # williams fork
plot(st_geometry(fire_sf), col = "grey")
plot(nbs, coords, add = TRUE, col = "red", lwd = 0.5)
title("KNN adjacency structure for example fire")

# check for disconnected components:
comp <- spdep::n.comp.nb(nbs)
print(table(comp$comp.id))

# setup for INLA
nb2INLA("cl_graph",nbs)
am_adj <- paste(getwd(),"data/spatial/cl_graph",sep="")
H <- inla.read.graph(filename="cl_graph")

##########################
# specify a pc-prior for the random effect
hyper_pr.fire <- list(
 prec = list(prior = "pc.prec", param = c(1, 0.1))
)
# update the model formula
mf.frp.sp <- update(
 mf.frp.re, . ~ 1 + . + 
  f(grid_index, 
    model = "bym2", 
    graph = H, 
    scale.model = T,
    constr = T,
    hyper = list(
     phi = list(prior = "pc", param = c(0.8, 2/3)),  # Proportion of spatial effect
     prec = list(prior = "pc.prec", param = c(1, 0.1))  # Overall variance
    )
   )
)
# make sure da has the same grid_index
da <- da %>%
 arrange(grid_index) %>%
 mutate(grid_index = as.numeric(as.factor(grid_index)))

# fit the model
model_bl.frp.sp <- inla(
 mf.frp.sp, 
 data = da, 
 family = "gaussian",
 control.predictor = list(compute = TRUE),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
 control.fixed = list(
  prec = list(prior = "pc.prec", param = c(1, 0.5)) 
 )
)
summary(model_bl.frp.sp)



#=================MODEL COMPARISON=================#

# Create a model comparison table
model_comparison <- tibble(
 Model = c("Baseline", 
           "W/Fire Random Effect", 
           "W/Fire + Temporal Random Effect",
           "W/Fire + Spatial Effect"),
 Response = "FRP",
 DIC = c(
  model_bl.frp$dic$dic,
  model_bl.frp.re$dic$dic,
  model_bl.frp.re2$dic$dic,
  model_bl.frp.sp$dic$dic
 ),
 WAIC = c(
  model_bl.frp$waic$waic,
  model_bl.frp.re$waic$waic,
  model_bl.frp.re2$waic$waic,
  model_bl.frp.sp$waic$waic
 ),
 Marginal_LogLikelihood = c(
  model_bl.frp$mlik[1],
  model_bl.frp.re$mlik[1],
  model_bl.frp.re2$mlik[1],
  model_bl.frp.sp$mlik[1]
 ),
 Effective_Params = c(
  model_bl.frp$dic$p.eff,
  model_bl.frp.re$dic$p.eff,
  model_bl.frp.re2$dic$p.eff,
  model_bl.frp.sp$dic$p.eff
 ),
 Mean_CPO = c(
  mean(model_bl.frp$cpo$cpo, na.rm = TRUE),
  mean(model_bl.frp.re$cpo$cpo, na.rm = TRUE),
  mean(model_bl.frp.re2$cpo$cpo, na.rm = TRUE),
  mean(model_bl.frp.sp$cpo$cpo, na.rm = TRUE)
 ),
 RMSE = c(
  sqrt(mean((model_bl.frp$summary.fitted.values$mean - da$frp_csum)^2, na.rm = TRUE)),
  sqrt(mean((model_bl.frp.re$summary.fitted.values$mean - da$frp_csum)^2, na.rm = TRUE)),
  sqrt(mean((model_bl.frp.re2$summary.fitted.values$mean - da$frp_csum)^2, na.rm = TRUE)),
  sqrt(mean((model_bl.frp.sp$summary.fitted.values$mean - da$frp_csum)^2, na.rm = TRUE))
 )
) %>% arrange(RMSE)
# Print the comparison table
print(model_comparison)


# Tidy up !
rm(model_bl.frp, model_bl.frp.re, model_bl.frp.re2)
gc()



#######################################
#######################################

########################
# COMPOSITE BURN INDEX #

# examine the distributions of CBIbc variables
da <- da %>%
 mutate(
  log_cbi_mn = log(cbi_mn),
  sqrt_cbi_p90 = sqrt(cbi_p90),
  sqrt_cbi_p99 = sqrt(cbi_p99)
 ) 
da %>%
 # pivot longer to facet plot
 pivot_longer(cols = c(cbi_p90, sqrt_cbi_p90,
                       cbi_p99, sqrt_cbi_p99,
                       cbi_mn, log_cbi_mn),
              names_to = "variable",
              values_to = "value") %>%
 # Plot with facets
 ggplot(aes(x = value)) +
 geom_histogram(bins = 30, fill = "orange", alpha = 0.7) +
 facet_wrap(
  ~ variable, 
  scales = "free",
  labeller = as_labeller(c(cbi_mn = "Average CBIbc",
                           log_cbi_mn = "log(Average CBIbc)",
                           cbi_p90 = "90th Percentile CBIbc",
                           cbi_p99 = "95th Percentile CBIbc",
                           sqrt_cbi_p90 = "sqrt(90th Percentile CBIbc)",
                           sqrt_cbi_p99 = "sqrt(99th Percentile CBIbc)"))) +
 labs(x = "value",
      y = "Frequency") +
 theme_minimal()

# save the plot.
out_png <- paste0(maindir,'figures/INLA_ResponseDistribution_CBI.png')
ggsave(out_png, dpi=500, bg = 'white')


#######################################
# 1. Baseline model (no latent effects)

# Set up the model formula
mf.cbi <- sqrt_cbi_p90 ~ 
 fortypnm_gp + # predominant forest type
 forest_pct + # percent of grid of predominant forest type
 erc + # energy release component (ERC)
 slope + tpi + chili # topography
# fit the model                     
model_bl.cbi <- inla(
 mf.cbi, data = da, 
 family = "gamma", 
 control.predictor = list(compute = TRUE),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
 control.fixed = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))
)
summary(model_bl.cbi)

##############################################
# 2. Baseline model + fire-level random effect

# update the model formula
hyper_pr <- list(prec = list(prior = "pc.prec", param = c(1, 0.1)))
mf.cbi.re <- update(
 mf.cbi, . ~ 1 + . + 
  f(fire_id, model = "iid", hyper = hyper_pr) # fire-level random effect
)
# fit the model                     
model_bl.cbi.re <- inla(
 mf.cbi.re, data = da, 
 family = "gamma",
 control.predictor = list(compute = TRUE),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
 control.fixed = list(
  prec = list(prior = "pc.prec", param = c(1, 0.1)) 
 )   
)
summary(model_bl.cbi.re)

########################################################
# 3. Baseline model + fire-level + spatial random effect

##########################
# update the model formula
mf.cbi.sp <- update(
 mf.cbi.re, . ~ 1 + . + 
  f(grid_index, 
    model = "bym2", 
    graph = H, 
    hyper = list(
     phi = list(prior = "pc", param = c(0.5, 0.3)),  # Proportion of spatial effect
     prec = list(prior = "pc.prec", param = c(1, 0.5))  # Overall variance
    )
  )
)
# fit the model
model_bl.cbi.sp <- inla(
 mf.cbi.sp, 
 data = da, 
 family = "gamma",
 control.predictor = list(compute = TRUE),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
 control.fixed = list(
  prec = list(prior = "pc.prec", param = c(1, 0.1)) 
 )  
)
summary(model_bl.cbi.sp)


#=================MODEL COMPARISON=================#

# Create a model comparison table
model_comparison <- tibble(
 Model = c("Baseline", 
           "W/Fire Random Effect", 
           "W/Fire + Temporal Random Effect",
           "W/Fire + Spatial Effect"),
 Response = "CBI",
 DIC = c(
  model_bl.frp$dic$dic,
  model_bl.frp.re$dic$dic,
  model_bl.frp.re2$dic$dic,
  model_bl.cbi.sp$dic$dic
 ),
 WAIC = c(
  model_bl.frp$waic$waic,
  model_bl.frp.re$waic$waic,
  model_bl.frp.re2$waic$waic,
  model_bl.cbi.sp$waic$waic
 ),
 Marginal_LogLikelihood = c(
  model_bl.frp$mlik[1],
  model_bl.frp.re$mlik[1],
  model_bl.frp.re2$mlik[1],
  model_bl.cbi.sp$mlik[1]
 ),
 Effective_Params = c(
  model_bl.frp$dic$p.eff,
  model_bl.frp.re$dic$p.eff,
  model_bl.frp.re2$dic$p.eff,
  model_bl.cbi.sp$dic$p.eff
 ),
 Mean_CPO = c(
  mean(model_bl.frp$cpo$cpo, na.rm = TRUE),
  mean(model_bl.frp.re$cpo$cpo, na.rm = TRUE),
  mean(model_bl.frp.re2$cpo$cpo, na.rm = TRUE),
  mean(model_bl.cbi.sp$cpo$cpo, na.rm = TRUE)
 ),
 RMSE = c(
  sqrt(mean((model_bl.frp$summary.fitted.values$mean - da$sqrt_cbi_p90)^2, na.rm = TRUE)),
  sqrt(mean((model_bl.frp.re$summary.fitted.values$mean - da$sqrt_cbi_p90)^2, na.rm = TRUE)),
  sqrt(mean((model_bl.frp.re2$summary.fitted.values$mean - da$sqrt_cbi_p90)^2, na.rm = TRUE)),
  sqrt(mean((model_bl.cbi.sp$summary.fitted.values$mean - da$sqrt_cbi_p90)^2, na.rm = TRUE))
 )
) %>% arrange(RMSE)
# Print the comparison table
print(model_comparison)


# Tidy up !
rm(coords, fire_sf, H, hyper_pr, model_bl.cbi, 
   model_bl.cbi.re, model_bl.frp, model_bl.frp.re, nbs)
gc()



#===========Plotting Posterior Effects============#

# Extract fixed effects related to fortypnm_gp

#####
# FRP
frp.eff <- as.data.frame(model_bl.frp.sp$summary.fixed) %>%
 rownames_to_column(var = "parameter") %>%
 filter(grepl("fortypnm_gp", parameter)) %>%
 mutate(
  response = "FRP",
  forest_type = gsub("fortypnm_gp", "", parameter),
  effect = mean,
  lower = `0.025quant`,
  upper = `0.975quant`
 )

#####
# CBI
cbi.eff <- as.data.frame(model_bl.cbi.sp$summary.fixed) %>%
 rownames_to_column(var = "parameter") %>%
 filter(grepl("fortypnm_gp", parameter)) %>%
 mutate(
  response = "CBIbc",
  forest_type = gsub("fortypnm_gp", "", parameter),
  effect = mean,
  lower = `0.025quant`,
  upper = `0.975quant`
 )

# Combine effects into one data frame
effects <- bind_rows(frp.eff, cbi.eff)

rm(cbi.eff, frp.eff) 
gc() # clean up

# Plot posterior effects with credible intervals
cmap <- c("CBIbc" = "#800026", "FRP" = "#FEB24C") # color map for the response
p1 <- ggplot(effects, aes(x = forest_type, y = effect, color = response)) +
 # Plot mean effects
 geom_point(position = position_dodge(width = 0.5), size = 3) +
 # Add credible intervals
 geom_errorbar(
  aes(ymin = lower, ymax = upper),
  position = position_dodge(width = 0.5),
  width = 0.2
 ) +
 # Add dashed horizontal line at y = 0
 geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
 # Customize the plot labels
 labs(
  x = "Predominant Forest Type",
  y = "Effect relative to aspen",
  color = "Response"
 ) +
 # Apply custom colors
 scale_color_manual(
  values = cmap,
  labels = c(
   "CBIbc" = "90th Percentile CBIbc", 
   "FRP" = "Maximum Daytime FRP")
 ) +
 # Adjust theme
 theme_bw() +
 theme(
  axis.text.x = element_text(angle = 45, hjust = 1, size=10),
  axis.text.y = element_text(angle = 0, hjust = 0, size=10),
  axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
  axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
  legend.position = c(0.18, 0.16),
  legend.background = element_rect(
   fill = scales::alpha("white", 0.2), 
   color = NA, size = 0.5),
  legend.title = element_text(size = 11),
  legend.text = element_text(size = 10)
 )
p1

# save the plot.
out_png <- paste0(maindir,'figures/INLA_FORTYPNM_PosteriorEffects.png')
ggsave(out_png, dpi=500, bg = 'white')


############
# Ridge plot
# Extract fixed effect marginals for FRP and CBI
frp_marginals <- model_bl.frp.sp$marginals.fixed
cbi_marginals <- model_bl.cbi.sp$marginals.fixed

# Define a function to tidy marginals
tidy_marginals <- function(marginals, response) {
 tibble::tibble(
  parameter = names(marginals),
  data = purrr::map(marginals, ~ as.data.frame(.x))
 ) %>%
  unnest(data) %>%
  mutate(response = response)
}

# Tidy FRP and CBI marginals
tidy_frp <- tidy_marginals(frp_marginals, "FRP")
tidy_cbi <- tidy_marginals(cbi_marginals, "CBI")

# Combine the data
tidy_combined <- bind_rows(tidy_frp, tidy_cbi)

# pull out the climate + topo
climate_topo <- tidy_combined %>%
 filter(parameter %in% c("erc_dv", "slope", "tpi", "chili")) %>%
 mutate(parameter = factor(parameter, levels = c("erc_dv", "slope", "tpi", "chili")))

# Filter for forest type effects
tidy_combined <- tidy_combined %>%
 filter(str_detect(parameter, "fortypnm_gp")) %>%
 mutate(forest_type = str_remove(parameter, "fortypnm_gp"))

# Compute mean effect sizes for reordering
effect_means <- tidy_combined %>%
 group_by(forest_type, response) %>%
 summarize(mean_effect = mean(x), .groups = "drop") %>%
 filter(response == "FRP")

tidy_combined <- tidy_combined %>%
 left_join(effect_means %>% select(forest_type, mean_effect), by = "forest_type") %>%
 mutate(forest_type = fct_reorder(forest_type, mean_effect, .desc = TRUE))

# create the plot
p2 <- ggplot(tidy_combined, aes(x = x, y = forest_type, height = y, fill = response)) +
 geom_density_ridges(stat = "identity", scale = 1.5, alpha = 0.7) +
 geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
 labs(
  x = "Effect relative to aspen",
  y = "Predominant Forest Type",
  fill = "Response",
 ) +
 scale_fill_manual(values = c("FRP" = "#FEB24C", "CBI" = "#800026"),
                   labels = c(
                    "CBI" = "90th Percentile CBIbc", 
                    "FRP" = "Maximum Daytime FRP")) +
 # Apply custom colors
 scale_color_manual(
  values = cmap,
  labels = c(
   "CBIbc" = "90th Percentile CBIbc", 
   "FRP" = "Maximum Daytime FRP")
 ) +
 theme_classic() +
 theme(axis.text.y = element_text(angle = 0, hjust = 0, size=11),
       axis.text.x = element_text(angle = 0, hjust = 0, size=11),
       axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
       axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
       legend.position = c(0.18, 0.16),
       legend.background = element_rect(
        fill = scales::alpha("white", 0.4), 
        color = NA, size = 0.8),
       legend.title = element_text(size = 11),
       legend.text = element_text(size = 10))
p2

# save the plot.
out_png <- paste0(maindir,'figures/INLA_FORTYPNM_PosteriorEffects_Ridge.png')
ggsave(out_png, dpi=300, bg = 'white')


# Plot the climate + topo effects
p_climate_topo <- ggplot(climate_topo, aes(x = x, y = parameter, height = y, fill = response)) +
 geom_density_ridges(stat = "identity", scale = 1.5, alpha = 0.7) +
 geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
 labs(
  x = "Effect size",
  y = "Climate/Topography Predictor",
  fill = "Response"
 ) +
 scale_fill_manual(
  values = c("FRP" = "#FEB24C", "CBI" = "#800026"),
  labels = c("FRP" = "Maximum Daytime FRP", "CBIbc" = "90th Percentile CBIbc")
 ) +
 theme_classic() +
 theme(
  axis.text.y = element_text(size = 11),
  axis.text.x = element_text(size = 11),
  axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
  axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
  legend.position = "top",
  legend.title = element_text(size = 11),
  legend.text = element_text(size = 10)
 )
p_climate_topo

# Save the plot
out_png <- paste0(maindir, 'figures/INLA_ClimateTopo_PosteriorEffects_Ridge.png')
ggsave(out_png, plot = p_climate_topo, dpi = 300, bg = 'white')



######################################
# Visualize the spatial random effects

# Extract spatial random effect estimates
spatial_effects <- model_bl.frp.sp$summary.random$grid_index %>%
 dplyr::select(ID, mean, sd, '0.025quant', '0.5quant', '0.975quant') %>%
 rename(grid_index = ID)

# Join spatial effects back to spatial data
grid_sf.eff <- grid_sf %>%
 left_join(spatial_effects, by = "grid_index")
sum(is.na(grid_sf.eff$mean))  # Count missing values

# check the histogram of spatial effects across fires
hist(grid_sf.eff$mean, breaks = 50, col = "skyblue")

# Plot the spatial effects (posterior mean)
# Use the Williams Fork Fire as an example
fire_sf <- grid_sf.eff %>% filter(fire_id == 51) %>% st_transform(st_crs(5070))
ggplot(fire_sf) +
 geom_sf(aes(color = mean), size=2.5) +
 scale_color_viridis_c(option = "magma") +
 theme_minimal()

# Save the plot
out_png <- paste0(maindir, 'figures/INLA_WilliamsFork_SpatialEff.png')
ggsave(out_png, dpi = 300, bg = 'white')


###########################################################
# Extract structured and unstructured components separately
n <- nrow(grid_sf.eff)
st_effects <- model_bl.frp.sp$summary.random$grid_index$mean[1:n] * 
 model_bl.frp.sp$summary.hyperpar$mean["Phi for grid_index"]
un_effects <- model_bl.frp.sp$summary.random$grid_index$mean[1:n] * 
 (1 - model_bl.frp.sp$summary.hyperpar$mean["Phi for grid_index"])

# Add components to spatial data
grid_sf.eff <- grid_sf.eff %>%
 mutate(structured = st_effects, unstructured = un_effects)

# Grab our test fire again
fire_sf <- grid_sf.eff %>% filter(fire_id == 51)

# Plot structured vs unstructured effects
p1 <- ggplot(fire_sf) +
 geom_sf(aes(color = structured), color = NA) +
 scale_color_viridis_c(name = "Structured Component") +
 labs(title = "Structured Spatial Component") +
 theme_minimal()
p2 <- ggplot(fire_sf) +
 geom_sf(aes(color = unstructured), color = NA) +
 scale_color_viridis_c(name = "Unstructured Component") +
 labs(title = "Unstructured Spatial Component") +
 theme_minimal()

# Combine the plots
p1 + p2


#==============MODEL DIAGNOSTICS================#

# Extract fitted values (posterior mean estimates)
fit_res <- data.frame(
 Fitted = model_bl.frp.sp$summary.fitted.values$mean,
 Residuals = da$frp - fitted
)
hist(fit_res$Residuals, breaks = 50, col = "skyblue", main = "Histogram of Residuals")

# Create a dataframe for plotting
comparison_df <- data.frame(
 type = c(rep("Observed", length(da$frp)), rep("Fitted", length(model_bl.frp.sp$summary.fitted.values$mean))),
 value = c(da$frp, model_bl.frp.sp$summary.fitted.values$mean)
)
ggplot(comparison_df, aes(x = value, color = type, fill = type)) +
 geom_density(alpha = 0.4) +
 labs(title = "Density Plot of Observed vs. Fitted Values", 
      x = "Value", 
      y = "Density") +
 scale_fill_manual(values = c("Observed" = "steelblue", "Fitted" = "orange")) +
 scale_color_manual(values = c("Observed" = "steelblue", "Fitted" = "orange")) +
 theme_minimal()


# Plot residuals vs fitted values
ggplot(fit_res, aes(x = Fitted, y = Residuals)) +
 geom_point(alpha = 0.5) +
 geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
 labs(title = "Residuals vs Fitted Values",
      x = "Fitted Values",
      y = "Residuals") +
 theme_minimal()






