
# libraries
library(tidyverse)
library(sf) 
library(INLA) 
library(ggcorrplot)
library(lubridate)
library(ggridges)
library(reshape2)
library(spdep)
library(patchwork)
library(forcats)
library(gstat)


maindir <- '/Users/max/Library/CloudStorage/OneDrive-Personal/mcook/aspen-fire/Aim2/'

# load and tidy the fire data
fires <- paste0(maindir,"data/spatial/mod/srm_fire_census_2017_to_2023_ics_perims.gpkg")
fires <- st_read(fires) %>%
 as_tibble() %>%
 rename(
  fire_ig_dt = DISCOVERY_DATE,
  fire_acres = ICS_ACRES,
  fire_aspenpct = Aspen_Pct
 ) %>%
 select(Fire_ID, fire_ig_dt, fire_acres, fire_aspenpct) %>%
 mutate(
  Fire_ID = as.numeric(Fire_ID),
  fire_acres = as.numeric(fire_acres),
  fire_aspenpct = as.numeric(fire_aspenpct))


#=========Prep the grid data=========#

# Format the species composition data frame 
# (long-format) each grid has rows for species co-occurring with the dominant type
# climate and topography are summarized at the grid level
fp <- paste0(maindir,'data/tabular/mod/gridstats_fortypnm_gp_tm_ct_frp-cbi.csv')
grid_tm <-  read_csv(fp)  %>% # read in the file
 select(-'...1') %>%
 # join in some of the fire information
 left_join(fires, by="Fire_ID") %>%
 # remove missing FRP, prep columns
 # Ensure daytime FRP and CBIbc is greater than 0
 filter(
  frp_csum > 0, # make sure cumulative FRP > 0
  CBIbc_p90 > 0 # make sure CBI > 0
 ) %>% 
 # create a numeric fire ID
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
  # create a gridcell-level live/dead proportion
  ba_ld_pr = ba_dead_total / (ba_live_total + ba_dead_total),
  tpp_ld_pr = tpp_dead_total / (tpp_live_total + tpp_dead_total),
 ) %>%
 # filter noise in species data
 filter(
  # # retain at least 50% forested gridcells
  # (forest_pct > 0.50) &
  # remove noise from  small species proportions
  (ba_live_pr > 0.10) # at least 10% of BA or TPP
 ) %>%
 group_by(Fire_ID) %>%
 mutate(aspen = ifelse(any(fortypnm_gp == "quaking_aspen"), 1, 0)) %>%
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
 # be sure there are no duplicate rows for grid/fire/species
 distinct(grid_idx, species_gp_n, .keep_all = TRUE) # remove duplicates
# tidy up
rm(fires)

################################################
# load and filter the duplicate (reburn) grids #
fp = paste0(maindir,"data/spatial/mod/gridcell_duplicates.gpkg")
dup_grids <- st_read(fp)
dup_idx <- unique(dup_grids$grid_index)
grid_tm <- grid_tm %>%
 filter(!grid_index %in% dup_idx)
print(paste0(
 "Removed [",dim(dup_grids)[1],"] duplicate gridcells."
))
rm(dup_grids, dup_idx) # tidy up

#############################
# check how many grids, etc #
print(length(unique(grid_tm$grid_idx))) # unique gridcells
# percent of gridcells majority forested
# Check how many grids are > 50% forested
print(paste0(
 "Percent forested gridcells (>50% forest pixels): ",
 round(dim(grid_tm %>% filter(forest_pct > 0.50) %>% distinct(grid_idx))[1]/
  dim(grid_tm %>% distinct(grid_idx))[1], 3) * 100, "%"
))
# retain predominantly forested gridcells ...
grid_tm <- grid_tm %>%
 filter(forest_pct >= 0.50)

##########################################
# filter fires with not enough gridcells #
# check on the grid cell counts
gridcell_counts <- grid_tm %>%
 distinct(Fire_ID, grid_idx) %>% # keep only distinct rows
 group_by(Fire_ID) %>%
 summarise(n_gridcells = n())
# check the distribution
summary(gridcell_counts$n_gridcells)
# calculate the quantile distribution
(qt <- tibble::enframe(
 round(quantile(gridcell_counts$n_gridcells, probs = seq(.1, .9, by = .1))),
 name = "qt", value = "val"
))
(qt10 = qt[1,]$val) # 10%
# filter fires below the 10th percentile
grid_tm <- grid_tm %>%
 group_by(Fire_ID) %>%
 filter(n() >= qt10) %>%
 ungroup()
# tidy up!
rm(gridcell_counts,qt)

# check how many grids and fires
length(unique(grid_tm$grid_idx))
length(unique(grid_tm$Fire_ID))

###############################
# check on the species counts #
grid_tm %>% 
 group_by(species_gp_n) %>%
 summarise(n = length(species_gp_n)) %>%
 arrange(desc(n))

###################################################
# save a model data frame for visualizations, etc #
# first, pivot wide for species columns (live basal area)
grid_tm.w <- grid_tm %>%
 select(
  grid_idx, fortypnm_gp, species_gp_n, ba_live_pr,
  log_frp_csum, log_frp_csum_day, log_frp_csum_night,
  log_frp_max, log_frp_max_day, log_frp_max_night,
  CBIbc_mean, CBIbc_p90, CBIbc_p95, CBIbc_p99,
  first_obs_doy, erc_dv, vpd, vs, elev, slope, aspect, tpi,
  ba_live_total, tpp_live_total, qmd_live_mean,
  lf_canopy, lf_height, H_tpp, H_ba, overlap, day_prop,
  tree_ht_live_mean, tree_dia_live_mean, ba_ld_pr, tpp_ld_pr,
  fortyp_pct, forest_pct, tmid_n
  ) %>%  
 pivot_wider(
  names_from = species_gp_n,  
  values_from = ba_live_pr,   
  values_fill = 0             
 ) %>%
 left_join(., grid_tm %>% distinct(grid_idx, geometry), by="grid_idx") %>%
 distinct(grid_idx, fortypnm_gp, geometry, .keep_all = T) %>%
 mutate(geometry = st_as_sfc(geometry)) %>%
 st_as_sf(.) %>% st_set_crs(st_crs(5070))
glimpse(grid_tm.w)
# save it out.
out_fp = paste0(maindir,"data/spatial/mod/gridcell_model_data.gpkg")
grid_tm.w %>% st_write(., out_fp, append=F)
rm(grid_tm.w)



#==========EXPLORE THE DATA==========#

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

########################################
# correlation matrix for fixed effects #
cor_da <- grid_tm %>%
 select(
  # species structure metrics
  tpp_live, ba_live, qmd_live, ht_live, dia_live,   
  tpp_live_pr, ba_live_pr, # proportions of TPP and BA
  ba_live_total, tpp_live_total, qmd_live_mean, # gridcell structural summaries
  ba_ld_pr, tpp_ld_pr, # live/dead proportions
  ba_dead_total, tpp_dead_total, # dead BA and TPP
  forest_pct, fortyp_pct, # gridcell forest percent and majority forest type percent
  H_ba, H_tpp, # gridcell species diversity (abundance- and dominance-based)
  erc, erc_dv, vpd, vpd_dv, # day-of-burn climate
  fm1000, rmin, tmmx, vs, # day-of-burn climate
  elev, slope, tpi, chili, aspect,  # topography
  tm_canopy, tm_balive, lf_canopy, lf_height, # gridcell canopy percent/height and BA sum
  day_prop, overlap, afd_count # VIIRS detection characteristics
 ) %>%
 mutate_if(is.factor, as.numeric)  # Convert factors to numeric (if needed)
# Compute correlation matrix
cor_mat <- cor(cor_da, use = "complete.obs", method = "spearman")
# Plot correlation matrix
ggcorrplot(cor_mat, method = "circle",  
 type = "lower", lab = TRUE, lab_size = 3, tl.cex = 10,  
 colors = c("blue", "white", "red")  
)
rm(cor_da, cor_mat) # tidy up
# save the plot.
out_png <- paste0(maindir,'figures/INLA_CorrelationMatrix_FixedEffects.png')
ggsave(out_png, dpi=500, bg = 'white')



#===========MODEL SETUP==============#

#########################################
# force aspen to be the baseline factor #
grid_tm <- grid_tm %>%
 mutate(
  species_gp_n = fct_relevel(species_gp_n, "quaking_aspen"),
  fortypnm_gp = fct_relevel(fortypnm_gp, "quaking_aspen")
 )
# check the factor levels
# make sure aspen is first
levels(grid_tm$species_gp_n)
levels(grid_tm$fortypnm_gp)

# prep the model data frame
# center and scale fixed effects
da <- grid_tm %>%
 mutate(
  # center/scale metrics / fixed effects
  across(
   c(ba_live, tpp_live, qmd_live, ht_live, dia_live, # species structure metrics
     ba_live_pr, tpp_live_pr, # proportion TPP and BA
     ba_ld_pr, tpp_ld_pr, # proportion dead BA/TPP
     ba_live_total, tpp_live_total, qmd_live_mean, # gridcell structural summaries
     ba_dead_total, tpp_dead_total, # dead structure
     H_ba, H_tpp, # gridcell diversity metrics
     forest_pct, fortyp_pct, # gridcell forest percent and majority forest type percent
     lf_canopy, lf_height, # gridcell forest characteristics
     erc_dv, vpd, # climate
     elev, slope, aspect, tpi, # topography
     day_prop, overlap # VIIRS detection characteristics
    ), ~ as.numeric(scale(.))
  ),
  log_fire_size = log(fire_acres) # log-scale fire size
 ) %>%
 # calculate the predominant species by raw BA and TPP
 group_by(grid_idx) %>%
 mutate(
  dom_sp_ba = as.factor(species_gp_n[which.max(ba_live)]),
  dom_sp_tpp = as.factor(species_gp_n[which.max(tpp_live)]),
  dom_sp_qmd = as.factor(species_gp_n[which.max(qmd_live)]),
 ) %>%
 ungroup() %>%
 arrange(grid_idx) # arrange by grid index
rm(grid_tm) # tidy up


##################################################
# check on the dominance/abundance distributions #
# get value counts for dominance/abundance by species
df1 <- da %>%
 group_by(fortypnm_gp) %>%
 distinct(grid_idx) %>%
 summarize(n_fortyp = n())# fortypcd
df2 <- da %>%
 group_by(dom_sp_ba) %>%
 distinct(grid_idx) %>%
 summarize(n_ba = n())# dominance
df3 <- da %>%
 group_by(dom_sp_tpp) %>%
 distinct(grid_idx) %>%
 summarize(n_tpp = n()) # abundance
(sp_dom_summary <- full_join(
 df1, df2, by = c("fortypnm_gp" = "dom_sp_ba")) %>%
  full_join(df3, by = c("fortypnm_gp" = "dom_sp_tpp")) %>%
  arrange(desc(n_fortyp)))
rm(df1,df2,df3)

###############################
# check on the species counts #
da %>% 
 group_by(species_gp_n) %>%
 summarise(n = length(species_gp_n)) %>%
 arrange(desc(n))


#===========MODEL FITTING==============#

set.seed(456)

# make aspen presence a factor
da$aspen <- as.factor(da$aspen)
levels(da$aspen)

########################
# FIRE RADIATIVE POWER #

#################################################
# 1. Baseline model (no random or latent effects)

# setup the model formula
mf.frp <- log_frp_csum ~ 1 + 
 f(fortypnm_gp, model="rw2") + # latent effect for primary forest type
 # fortypnm_gp:fortyp_pct + # majority forest type
 species_gp_n:ba_live + # species live BA
 species_gp_n:tpp_live + # species live BA
 # species_gp_n:qmd_live + # species live QMD
 species_gp_n:ht_live + # species live tree height
 species_gp_n:dia_live + # species live tree diameter
 # fortyp_pct + # percent cover of the dominant forest type
 H_ba + # gridcell structural diversity (dominance-based)
 lf_canopy + # gridcell forest canopy percent
 ba_dead_total + # proportion live/dead basal area
 erc_dv + vpd + vs + # day-of-burn climate/weather
 elev + slope + aspect + tpi + # topography 
 aspen + # fire-level aspen presence
 overlap + # gridcell VIIRS overlap (cumulative)
 day_prop # gridcell proportion daytime detections 

# fit the model
ml.frp <- inla(
 mf.frp, data = da,
 family = "gaussian",
 control.predictor = list(compute=T),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
 control.family = list(
  # relax the variance assumptions
  hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.1)))
 ),
 # controls for posterior approximation and hyperparameter search
 control.inla = list(strategy = "laplace", int.strategy = "grid")
)
summary(ml.frp)
# check on predictive power of the random effects models
mean(ml.frp$cpo$cpo, na.rm = TRUE)


#####################################
# 2. Baseline + Fire-ID Random Effect

# update the model formula
mf.frp.re <- update(
 mf.frp, . ~ 1 + . + 
  f(Fire_ID, model = "iid", 
    hyper = list(
     prec = list(prior = "pc.prec", param = c(1, 0.01))
    )) 
)

# fit the model
ml.frp.re <- inla(
 mf.frp.re, data = da,
 family = "gaussian",
 control.predictor = list(compute=T),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
 control.family = list(
  # relax the variance assumptions
  hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.1)))
 ),
 # controls for posterior approximation and hyperparameter search
 control.inla = list(strategy = "laplace", int.strategy = "grid")
)
summary(ml.frp.re)
# check on predictive power of the random effects models
mean(ml.frp.re$cpo$cpo, na.rm = TRUE)


###################################################
# 3. Baseline + Temporal + Fire-level Random Effect

# update the model formula
mf.frp.re2 <- update(
 mf.frp.re, . ~ 1 + . + 
  f(first_obs_date, 
    model = "iid", 
    hyper = list(
     prec = list(prior = "pc.prec", param = c(1, 0.01))
    )) 
)

# fit the model
ml.frp.re2 <- inla(
 mf.frp.re2, data = da,
 family = "gaussian",
 control.predictor = list(compute=T),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
 control.family = list(
  # relax the variance assumptions
  hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.1)))
 ),
 # controls for posterior approximation and hyperparameter search
 control.inla = list(strategy = "laplace", int.strategy = "grid")
)
summary(ml.frp.re2)
# check on predictive power of the random effects models
mean(ml.frp.re2$cpo$cpo, na.rm = TRUE)



#########################
# 4. Spatial SPDE model #

# extracting spatial coordinates for grid centroids
grid_sf <- da %>%
 arrange(grid_index) %>%
 st_as_sf(., coords = c("x", "y"), crs = 4326)
# extract coordinates
coords <- grid_sf %>% st_coordinates(.)
# convert coordinates to a matrix for INLA
coords_mat <- as.matrix(coords) 


# ##########################################################
# # fit a semivariogram to look at spatial dependence in FRP
# ##########################################################
# 
# # reproject
# sp <-  grid_sf %>%
#  st_transform(st_crs(5070))
# 
# # Function to compute semivariogram for each fire
# get_fire_vario <- function(Fire_ID, data) {
#  fire_sp <- data %>% filter(Fire_ID == !!Fire_ID)
#  
#  if (nrow(fire_sp) > 10) {  # Only compute if there are enough points
#   vario <- variogram(log_frp_csum ~ 1, data = fire_sp)
#   vario <- vario %>% mutate(Fire_ID = Fire_ID)  # Ensure fire_id is added correctly
#   print(paste("Computed variogram for fire:", Fire_ID, "with", nrow(fire_sp), "points"))
#   return(vario)
#  } else {
#   print(paste("Skipping fire:", Fire_ID, "- Not enough points"))
#   return(NULL)
#  }
# }
# 
# # Apply to each fire separately
# fires <- unique(sp$Fire_ID)
# variograms <- lapply(fires, get_fire_vario, data = sp)
# variograms <- do.call(rbind, variograms)  # Combine results
# 
# # Check if fire_id is still in the dataset
# print(unique(variograms$Fire_ID))
# 
# # Plot variograms grouped by Fire ID
# ggplot(variograms, aes(x = dist, y = gamma, 
#                        group = Fire_ID, color = Fire_ID)) +
#  geom_line(alpha = 0.3, size=0.9) +  # Light transparency to see overlap
#  labs(x = "Distance (meters)", y = "Semivariance", 
#       title = "Semivariograms for Individual Fires") +
#  theme_minimal() +
#  theme(legend.position = "none")
# # save the plot.
# out_png <- paste0(maindir,'figures/INLA_Fire_SemiVariograms_FRP.png')
# ggsave(out_png, dpi=500, bg = 'white')
# 
# # Find the approximate range (where semivariance plateaus) for each fire
# range_est <- variograms %>%
#  group_by(Fire_ID) %>%
#  summarize(range_m = max(dist[gamma < max(gamma) * 0.9], na.rm = TRUE))  # 90% of max
# qt <- quantile(range_est$range_m, probs = seq(.1, .9, by = .1))
# qt
# 
# # Compute and print the mean spatial range for FRP
# mean_range_frp <- mean(range_est$range_m, na.rm = TRUE)
# print(paste("Mean spatial range for FRP:", round(mean_range_frp, 2), "meters"))
# 
# # Histogram of estimated spatial ranges
# ggplot(range_est, aes(x = range_m)) +
#  geom_histogram(bins = 20, fill = "blue", alpha = 0.6) +
#  labs(x = "Estimated Range (meters)", y = "Count",
#       title = "Distribution of Within-Fire Spatial Correlation Ranges") +
#  theme_minimal()
# 
# # save the plot.
# out_png <- paste0(maindir,'figures/INLA_Fire_EstimatedRange_FRP.png')
# ggsave(out_png, dpi=500, bg = 'white')
# 
# rm(sp, range_est, variograms)
# gc()

##########################################################
##########################################################
##########################################################

# Define the spatial mesh
mesh <- inla.mesh.2d(
 loc = coords_mat, # Locations (grid centroids)
 max.edge = c(1, 20), # Maximum edge lengths (inner and outer)
 cutoff = 0.01, # Minimum distance between points (0.01deg = ~1.1km)
 offset = c(0.5, 0.1) # Boundary buffer
)
# Plot mesh to check
plot(mesh, main = "SPDE Mesh for FRP Model")
points(coords, col = "red", pch = 20)

rm(coords)
gc()

# ###################################################
# # Save a spatial file of the mesh grid and vertices
# ###################################################
# 
# ####### Save the mesh as spatial objects
# # Convert INLA mesh vertices to an sf object
# mesh_v <- data.frame(mesh$loc) # Extract coordinates
# colnames(mesh_v) <- c("x", "y") # set column names
# mesh_sf <- st_as_sf(mesh_v, coords = c("x", "y"), crs = 4326)
# # Export vertices to a gpkg
# st_write(mesh_sf, paste0(maindir,"data/spatial/mod/INLA_mesh_vertices.gpkg"), 
#          layer = "mesh_vertices", driver = "GPKG", append=FALSE)
# # Function to convert triangles to polygons
# mesh_t <- mesh$graph$tv  # Triangle indices
# tri_coords <- mesh$loc  # Mesh coordinates
# # Create list of polygons
# poly_list <- lapply(1:nrow(mesh_t), function(i) {
#  coords <- tri_coords[mesh_t[i, ], ]
#  coords <- rbind(coords, coords[1, ])  # Close the polygon
#  st_polygon(list(coords))
# })
# # Convert to sf object
# mesh_poly_sf <- st_sfc(poly_list, crs = 4326) %>% st_sf(geometry = .)
# # Export mesh triangles to GPKG
# st_write(mesh_poly_sf, paste0(maindir,"data/spatial/mod/INLA_mesh_triangles.gpkg"), 
#          layer = "mesh_triangles", driver = "GPKG", append=FALSE)
# rm(mesh_v, mesh_sf, mesh_t, tri_coords, poly_list, mesh_poly_sf)
# gc()

########################################
########################################

######################
# Build the SPDE model
spde.ml <- inla.spde2.pcmatern(
 # Mesh and smoothness parameter
 mesh = mesh, 
 alpha = 2,
 # P(practic.range < 0.3) = 0.5
 prior.range = c(0.3, 0.5),
 # P(sigma > 1) = 0.01
 prior.sigma = c(1, 0.01) # variance
)

# Compute the projector matrix (A)
A <- inla.spde.make.A(
 mesh = mesh,
 loc = coords_mat
)

# Assign the spatial index
field.idx <- inla.spde.make.index(
 name = "mesh.idx",
 n.spde = spde.ml$n.spde
)
str(field.idx)

# Extract the predictor variables
X <- da %>% 
 select(Fire_ID, first_obs_date, grid_index, 
        fortypnm_gp, fortyp_pct, species_gp_n,  
        tpp_live, tpp_live_pr, ba_live, ba_live_pr,    
        lf_canopy, H_tpp, H_ba, 
        ht_live, dia_live, qmd_live, 
        ba_ld_pr, tpp_ld_pr, ba_dead_total,
        erc_dv, vpd, vs, 
        elev, slope, aspect, tpi,
        aspen, day_prop, overlap)
head(X)

# Create the INLA data stack
stack.frp <- inla.stack(
 data = list(log_frp_csum = da$log_frp_csum),
 A = list(A, 1),  
 tag = 'est',
 effects = list(
  c(field.idx),
  list(as.data.frame(X))
 )
)
dim(inla.stack.A(stack.frp))

rm(grid_sf, X)
gc()

##########################
# update the model formula
mf.frp.re.sp <- update(
 mf.frp.re, . ~ 1 + . + 
  f(mesh.idx, model = spde.ml) # spatial process model
)
# fit the model
ml.frp.re.sp <- inla(
 mf.frp.re.sp, 
 data = inla.stack.data(stack.frp), ,
 family = "gaussian",
 control.predictor = list(A = inla.stack.A(stack.frp), compute=T),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
 control.family = list(
  # relax the variance assumptions
  hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.1)))
 ),
 # controls for posterior approximation and hyperparameter search
 control.inla = list(strategy = "laplace", int.strategy = "grid")
)
summary(ml.frp.re.sp)
mean(ml.frp.re.sp$cpo$cpo, na.rm = TRUE)


#==========EXTRACTING THE SPATIAL EFFECT===========#

# Compute spatial effect for each grid location
# mean
spat.eff <- inla.spde.make.A(mesh, coords_mat) %*% 
 ml.frp.re.sp$summary.random$mesh.idx$mean  # FRP spatial field
# credible intervals
# Compute 2.5% (lower bound) credible interval
lower <- inla.spde.make.A(mesh, coords_mat) %*% 
 ml.frp.re.sp$summary.random$mesh.idx$`0.025quant`
# Compute 97.5% (upper bound) credible interval
upper <- inla.spde.make.A(mesh, coords_mat) %*% 
 ml.frp.re.sp$summary.random$mesh.idx$`0.975quant`
# Compute uncertainty as the credible interval width (upper - lower)
uncertainty <- upper - lower

# convert to a spatial dataframe
spat.eff.frp <- data.frame(
 x = coords_mat[, 1],  # Longitude
 y = coords_mat[, 2],  # Latitude
 spat_effect = as.vector(spat.eff),  # Convert matrix to vector
 lower = as.vector(lower),  # 2.5% quantile
 upper = as.vector(upper),  # 97.5% quantile
 uncertainty = as.vector(uncertainty),  # CI width
 response = "FRP"
) %>%
 distinct(x, y, .keep_all=T) %>%
 st_as_sf(., coords = c("x", "y"), crs = 4326)
glimpse(spat.eff.frp)
# write this a spatial file
st_write(spat.eff.frp, paste0(maindir,"data/spatial/mod/spatial_effect_FRP_spp.gpkg"), append=F)

# Tidy up!
rm(spat.eff, lower, upper, uncertainty)
gc()


#=================MODEL COMPARISON=================#

# Create a model comparison table
ml_comparison.frp <- tibble(
 Model = c("Baseline", 
           "W/Fire Random Effect", 
           "W/Fire + Temporal Effect",
           "W/Fire + Temporal + Spatial Effect"),
 Response = "FRP",
 DIC = c(
  ml.frp$dic$dic,
  ml.frp.re$dic$dic,
  ml.frp.re2$dic$dic,
  ml.frp.re.sp$dic$dic
 ),
 WAIC = c(
  ml.frp$waic$waic,
  ml.frp.re$waic$waic,
  ml.frp.re2$waic$waic,
  ml.frp.re.sp$waic$waic
 ),
 Marginal_LogLikelihood = c(
  ml.frp$mlik[1],
  ml.frp.re$mlik[1],
  ml.frp.re2$mlik[1],
  ml.frp.re.sp$mlik[1]
 ),
 Effective_Params = c(
  ml.frp$dic$p.eff,
  ml.frp.re$dic$p.eff,
  ml.frp.re2$dic$p.eff,
  ml.frp.re.sp$dic$p.eff
 ),
 Mean_CPO = c(
  mean(ml.frp$cpo$cpo, na.rm = TRUE),
  mean(ml.frp.re$cpo$cpo, na.rm = TRUE),
  mean(ml.frp.re2$cpo$cpo, na.rm = TRUE),
  mean(ml.frp.re.sp$cpo$cpo, na.rm = TRUE)
 )
) %>% arrange(DIC)
# Print the comparison table
print(ml_comparison.frp)


# tidy up !
rm(A, coords_mat, field.idx, mesh, 
   ml.frp, ml.frp.re, ml.frp.re2, spde.ml, stack.frp)
gc()



#===========POSTERIOR EFFECTS===========#


#########################################
# Plot all of the posterior fixed effects
# Extract fixed effect marginals
frp_marginals <- ml.frp.re.sp$marginals.fixed
# Tidy marginals for all fixed effects
tidy.effects.frp <- tibble::tibble(
 parameter = names(frp_marginals),
 data = purrr::map(frp_marginals, ~ as.data.frame(.x))
) %>%
 unnest(data) %>%
 # Exclude the intercept
 filter(!parameter %in% c(
  "(Intercept)","aspen","day_prop"
  )
 ) %>%  
 mutate(
  effect = case_when(
   str_detect(parameter, "ba_live") ~ "Live basal area",
   # str_detect(parameter, "ba_live_pr") ~ "Proportion of\nlive basal area",
   # str_detect(parameter, "tpp_live") ~ "Live trees/pixel",
   str_detect(parameter, "dia_live") ~ "Live tree diameter",
   str_detect(parameter, "ht_live") ~ "Live tree height",
   # str_detect(parameter, "qmd_live") ~ "Live QMD",
   # str_detect(parameter, "H_tpp") ~ "Shannon diversity \n(abundance-based)",
   TRUE ~ parameter  # Default for all other fixed effects
  ),
  # Extract species names from parameters
  species = case_when(
   str_detect(parameter, "species_gp_n") ~ str_extract(parameter, "(?<=species_gp_n)[a-zA-Z_ñ]+"),
   str_detect(parameter, "fortypnm_gp") ~ str_extract(parameter, "(?<=fortypnm_gp)[a-zA-Z_ñ]+"),
   TRUE ~ NA_character_  # For non-species effects
  )
 ) %>%
 # Calculate mean effect size for ordering
 group_by(effect) %>%
 mutate(mean_effect = mean(x, na.rm = TRUE)) %>%
 ungroup() 

# Order effects: Species-specific effects first, 
# global effects by mean effect size
tidy.effects.frp <- tidy.effects.frp %>%
 mutate(
  effect_order = case_when(
   !is.na(species) ~ 2,  # Species-specific effects
   TRUE ~ 2              # Global effects
  ),
  effect = factor(effect, levels = tidy.effects.frp %>%
                   arrange(effect_order, desc(mean_effect)) %>%
                   pull(effect) %>%
                   unique()),
  fill_species = ifelse(is.na(species), "Global effect", species)
 ) %>%
 mutate(fill_species = recode(
  fill_species,
  "quaking_aspen" = "Quaking aspen",
  "lodgepole_pine" = "Lodgepole pine",
  "douglas_fir" = "Douglas-fir",
  "white_fir" = "White fir",
  "gambel_oak" = "Gambel oak",
  "piñon_juniper" = "Piñon-juniper",
  "ponderosa_pine" = "Ponderosa pine",
  "spruce_fir" = "Spruce-fir"
 )) %>%
 mutate(
  fill_species = factor(
   fill_species,
   levels = c("Lodgepole pine", "Douglas-fir", "White fir",
              "Gambel oak", "Piñon-juniper", "Ponderosa pine",
              "Spruce-fir", "Quaking aspen")))

# check on the species name extraction
unique(tidy.effects.frp$fill_species)
spps_breaks <- unique(tidy.effects.frp$fill_species)


###################################
# Extract the forest type effects #
fortyp_effects.frp <- tidy.effects.frp %>%
 filter(str_detect(parameter, "fortypnm_gp"),
        !str_detect(parameter, ":fortyp_pct")) %>%
 group_by(fill_species) %>%
 mutate(mean_effect = mean(x, na.rm = TRUE)) %>%
 ungroup() %>%
 mutate(fill_species = fct_reorder(fill_species, -mean_effect))
glimpse(fortyp_effects.frp)

# make the ridge plot
frp.p1 <- ggplot(fortyp_effects.frp,
       aes(x = x, y = fill_species, height = y,
           fill = fill_species)) +
 geom_density_ridges(stat = "identity", scale = 1.5, alpha = 0.7) +
 geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
 labs(
  x = "Effect Size",
  y = "Forest Type",
 ) +
 scale_fill_manual(
  values = c(
   "Quaking aspen" = "#e6ab02",
   "Lodgepole pine" = "#d95f02",
   "Douglas-fir" = "#7570b3",
   "White fir" = "#a6cee3",
   "Gambel oak" = "#e7298a",
   "Piñon-juniper" = "#66a61e",
   "Ponderosa pine" = "#1b9e77",
   "Spruce-fir" = "#1f78b4"
  )
 ) +
 theme_minimal() +
 theme(
  legend.position = "none",
  axis.text.y = element_text(size = 10),
  axis.title.y = element_text(size = 12)
 )
frp.p1


#######################################################
# create the ridge plot for species structure metrics #
sp_effects.frp <- tidy.effects.frp %>%
 filter(str_detect(parameter, "species_gp_n"))

# plot it
frp.p2 <- sp_effects.frp %>%
 ggplot(., aes(x = x, y = effect, 
               height = y, fill = fill_species, 
               alpha = fill_species)) +
 geom_density_ridges(
  stat = "identity", scale = 1.5, show.legend=F) +
 geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
 labs(
  x = "",
  y = "Fixed Effect",
  fill = "Species"
 ) +
 # # add a subplot label
 # annotate("text", x = -0.13, y = 5, 
 #          label = "(A)", fontface = "bold", 
 #          size = 4, hjust = 0) +
 scale_fill_manual(
  values = c(
   "Quaking aspen" = "#e6ab02",
   "Lodgepole pine" = "#d95f02",
   "Douglas-fir" = "#7570b3",
   "White fir" = "#a6cee3",
   "Gambel oak" = "#e7298a",
   "Piñon-juniper" = "#66a61e",
   "Ponderosa pine" = "#1b9e77",
   "Spruce-fir" = "#1f78b4"
  ),
  # Exclude "Global Effect" from the legend
  breaks = spps_breaks,  
  guide = guide_legend(title = "Species")
 ) +
 scale_alpha_manual(
  values = c(
   "Quaking aspen" = 0.9,
   "Lodgepole pine" = 0.6,
   "Douglas-fir" = 0.6,
   "White fir" = 0.6,
   "Gambel oak" = 0.6,
   "Piñon-juniper" = 0.6,
   "Ponderosa pine" = 0.6,
   "Spruce-fir" = 0.6
  ),
  guide = "none"  # Hide alpha from legend
 ) +
 # xlim(-0.13, 0.15) +
 theme_minimal() +
 theme(
  axis.text.y = element_text(angle = 0, hjust = 1, size=9),
  axis.text.x = element_text(angle = 0, hjust = 0, size=9),
  axis.title.y = element_text(size = 10, margin = margin(r = 8)),
  axis.title.x = element_text(size = 10, margin = margin(t = 8)),
  legend.position = c(0.85, 0.70),
  legend.background = element_rect(
   fill = scales::alpha("white", 0.6), 
   color = NA, size = 1),
  legend.title = element_text(size = 9),
  legend.text = element_text(size = 8)
 )
frp.p2

# Save the plot
out_png <- paste0(maindir, 'figures/INLA_ForestComposition_FixedEffects_FRP_species.png')
ggsave(out_png, dpi = 500, bg = 'white')



# ###############################################
# # create the ridge plot for all fixed effects #
# frp.p3 <- sp_effects.frp %>%
#  ggplot(., aes(x = x, y = effect, height = y, fill = fill_species)) +
#  geom_density_ridges(stat = "identity", scale = 1.5, alpha = 0.7) +
#  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
#  labs(
#   x = "Effect Size",
#   y = "Fixed Effect",
#   fill = "Species"
#  ) +
#  scale_fill_manual(
#   values = c(
#    "Global Effect" = "gray",  # Neutral color for global effects
#    "quaking_aspen" = "#e6ab02",
#    "lodgepole_pine" = "#d95f02",
#    "douglas_fir" = "#7570b3",
#    "white_fir" = "#a6cee3",
#    "piñon_juniper" = "#e7298a",
#    "ponderosa_pine" = "#66a61e",
#    "spruce_fir" = "#1b9e77"
#   ),
#   # Exclude "Global Effect" from the legend
#   breaks = spps_breaks,  
#   guide = guide_legend(title = "Species")
#  ) +
#  theme_minimal() +
#  theme(
#   legend.title = element_text(size = 10),
#   legend.text = element_text(size = 8),
#   axis.text.y = element_text(size = 10),
#   axis.title.y = element_text(size = 12)
#  )
# frp.p1
# 
# # Save the plot
# out_png <- paste0(maindir, 'figures/INLA_ForestComposition_FixedEffects_full.png')
# ggsave(out_png, dpi = 500, bg = 'white')

# # Tidy up!
# rm(frp_marginals, tidy.effects)
# gc()




########################
# COMPOSITE BURN INDEX #


#######################################
# 1. Baseline model (no latent effects)

# setup the model formula
mf.cbi <- CBIbc_p90 ~ 1 + 
 fortypnm_gp:fortyp_pct + # majority forest type
 species_gp_n:ba_live + # species live BA
 # species_gp_n:tpp_live + # species live BA
 # species_gp_n:qmd_live + # species live QMD
 # species_gp_n:ht_live + # species live tree height
 # species_gp_n:dia_live + # species live tree diameter
 # fortyp_pct + # percent cover of the dominant forest type
 H_tpp + # gridcell structural diversity (dominance-based)
 lf_canopy + # gridcell forest canopy percent
 ba_dead_total + # proportion live/dead basal area
 erc_dv + vpd + vs + # day-of-burn climate/weather
 elev + slope + aspect + tpi + # topography 
 aspen # fire-level aspen presence
# fit the model
ml.cbi <- inla(
 mf.cbi, data = da,
 family = "gamma",
 control.predictor = list(compute=T),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
 control.fixed = list(
  prec = list(prior = "loggamma", param = c(0.5, 0.1))
 ),
 control.family = list(
  hyper = list(prec = list(prior = "loggamma", param = c(0.5, 0.1)))
 ),
 control.inla = list(strategy = "adaptive", int.strategy = "grid")
)
summary(ml.cbi)
# check on predictive power of the random effects models
mean(ml.cbi$cpo$cpo, na.rm = TRUE)


##############################################
# 2. Baseline model + temporal random effect

# update the model formula
mf.cbi.re <- update(
 mf.cbi, . ~ 1 + . + 
  f(first_obs_date, model = "iid", 
    hyper = list(prec = list(prior = "pc.prec", param = c(0.5, 0.01)))
  ) # temporal random effect
)
# fit the model                     
ml.cbi.re <- inla(
 mf.cbi.re, data = da, 
 family = "gamma",
 control.predictor = list(compute = TRUE),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
 control.fixed = list(
  prec = list(prior = "loggamma", param = c(0.5, 0.1))
 ),
 control.family = list(
  hyper = list(prec = list(prior = "loggamma", param = c(0.5, 0.1)))
 ),
 control.inla = list(strategy = "adaptive", int.strategy = "grid")
)
summary(ml.cbi.re)
# check on predictive power of the random effects models
mean(ml.cbi.re$cpo$cpo, na.rm = TRUE)


#######################################################################
# 3. Baseline model + temporal random effect + fire-level random effect 

# update the model formula
mf.cbi.re2 <- update(
 mf.cbi.re, . ~ 1 + . + 
  f(fire_id, model = "iid", 
    hyper = list(prec = list(prior = "pc.prec", param = c(0.5, 0.01)))
  ) # fire-level random effect
)
# fit the model                     
ml.cbi.re2 <- inla(
 mf.cbi.re2, data = da, 
 family = "gamma",
 control.predictor = list(compute = TRUE),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
 control.fixed = list(
  prec = list(prior = "loggamma", param = c(0.5, 0.1))
 ),
 control.family = list(
  hyper = list(prec = list(prior = "loggamma", param = c(0.5, 0.1)))
 ),
 control.inla = list(strategy = "adaptive", int.strategy = "grid")
)
summary(ml.cbi.re2)
# check on predictive power of the random effects models
mean(ml.cbi.re2$cpo$cpo, na.rm = TRUE)


#######################
# 4. Spatial SPDE model

# extracting spatial coordinates for grid centroids
grid_sf <- da %>%
 arrange(grid_index) %>%
 st_as_sf(., coords = c("x", "y"), crs = 4326)
# extract coordinates
coords <- grid_sf %>% st_coordinates(.)
# convert coordinates to a matrix for INLA
coords_mat <- as.matrix(coords) 

# Define the spatial mesh
mesh <- inla.mesh.2d(
 loc = coords_mat, # Locations (grid centroids)
 max.edge = c(1, 20), # Maximum edge lengths (inner and outer)
 cutoff = 0.01, # Minimum distance between points (0.01deg = ~1.1km)
 offset = c(0.5, 0.1) # Boundary buffer
)
# Plot mesh to check
plot(mesh)
points(coords, col = "red", pch = 20)

rm(coords)
gc()

######################
# Build the SPDE model
spde.ml <- inla.spde2.pcmatern(
 # Mesh and smoothness parameter
 mesh = mesh, 
 alpha = 2,
 # P(practic.range < 0.3) = 0.5
 prior.range = c(0.1, 0.5),
 # P(sigma > 1) = 0.01
 prior.sigma = c(1, 0.01) # variance
)

# Compute the projector matrix (A)
A <- inla.spde.make.A(
 mesh = mesh,
 loc = coords_mat
)

# Assign the spatial index
field.idx <- inla.spde.make.index(
 name = "mesh.idx",
 n.spde = spde.ml$n.spde
)
str(field.idx)

# Extract the predictor variables
X <- da %>% 
 select(fire_id, first_obs_date, grid_index, 
        fortypnm_gp, species_gp_n, total_tpp, forest_pct,      
        canopypct, H_tpp, H_ba, ba_live_pr, tpp_live_pr, qmd_live, 
        tree_ht_live, live_dead_pr_tpp, live_dead_pr_ba,
        erc_dv, vpd, vs, slope, tpi, chili, elev,
        aspen, day_prop)
head(X)

# Create the INLA data stack
stack.cbi <- inla.stack(
 data = list(CBIbc_p90 = da$CBIbc_p90),
 A = list(A, 1),  
 tag = 'est',
 effects = list(
  c(field.idx),
  list(as.data.frame(X))
 )
)
dim(inla.stack.A(stack.cbi))

rm(grid_sf, X)
gc()

##########################
# update the model formula
mf.cbi.re.sp <- update(
 mf.cbi.re2, . ~ 1 + . + 
  f(mesh.idx, model = spde.ml) # spatial process model
)
# fit the model
ml.cbi.re.sp <- inla(
 mf.cbi.re.sp, 
 data = inla.stack.data(stack.cbi), ,
 family = "gamma",
 control.predictor = list(A = inla.stack.A(stack.cbi), compute=T),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
 control.fixed = list(
  prec = list(prior = "loggamma", param = c(0.5, 0.1))
 ),
 control.family = list(
  hyper = list(prec = list(prior = "loggamma", param = c(0.5, 0.1)))
 ),
 control.inla = list(strategy = "adaptive", int.strategy = "grid")
)
summary(ml.cbi.re.sp)

# check on predictive power of the random effects models
mean(ml.cbi.re.sp$cpo$cpo, na.rm = TRUE)



#==========EXTRACTING THE SPATIAL EFFECT===========#

# Compute spatial effect for each grid location
# mean
spat.eff <- inla.spde.make.A(mesh, coords_mat) %*% 
 ml.cbi.re.sp$summary.random$mesh.idx$mean  # FRP spatial field
# credible intervals
# Compute 2.5% (lower bound) credible interval
lower <- inla.spde.make.A(mesh, coords_mat) %*% 
 ml.cbi.re.sp$summary.random$mesh.idx$`0.025quant`
# Compute 97.5% (upper bound) credible interval
upper <- inla.spde.make.A(mesh, coords_mat) %*% 
 ml.cbi.re.sp$summary.random$mesh.idx$`0.975quant`
# Compute uncertainty as the credible interval width (upper - lower)
uncertainty <- upper - lower

# convert to a spatial dataframe
spat.eff.cbi <- data.frame(
 x = coords_mat[, 1],  # Longitude
 y = coords_mat[, 2],  # Latitude
 spat_effect = as.vector(spat.eff),  # Convert matrix to vector
 lower = as.vector(lower),  # 2.5% quantile
 upper = as.vector(upper),  # 97.5% quantile
 uncertainty = as.vector(uncertainty),  # CI width
 response = "CBI"
) %>%
 distinct(x, y, .keep_all=T) %>%
 st_as_sf(., coords = c("x", "y"), crs = 4326)
glimpse(spat.eff.cbi)

# write this a spatial file
st_write(spat.eff.cbi, paste0(maindir,"data/spatial/mod/spatial_effect_CBI_spp.gpkg"), append=F)

# Tidy up!
rm(spat.eff, lower, upper, uncertainty)
gc()



#=================MODEL COMPARISON=================#

# Create a model comparison table
ml_comparison.cbi <- tibble(
 Model = c("Baseline", 
           "W/Fire Random Effect", 
           "W/Fire + Temporal Effect",
           "W/Fire + Temporal + Spatial Effect"),
 Response = "CBI",
 DIC = c(
  ml.cbi$dic$dic,
  ml.cbi.re$dic$dic,
  ml.cbi.re2$dic$dic,
  ml.cbi.re.sp$dic$dic
 ),
 WAIC = c(
  ml.cbi$waic$waic,
  ml.cbi.re$waic$waic,
  ml.cbi.re2$waic$waic,
  ml.cbi.re.sp$waic$waic
 ),
 Marginal_LogLikelihood = c(
  ml.cbi$mlik[1],
  ml.cbi.re$mlik[1],
  ml.cbi.re2$mlik[1],
  ml.cbi.re.sp$mlik[1]
 ),
 Effective_Params = c(
  ml.cbi$dic$p.eff,
  ml.cbi.re$dic$p.eff,
  ml.cbi.re2$dic$p.eff,
  ml.cbi.re.sp$dic$p.eff
 ),
 Mean_CPO = c(
  mean(ml.cbi$cpo$cpo, na.rm = TRUE),
  mean(ml.cbi.re$cpo$cpo, na.rm = TRUE),
  mean(ml.cbi.re2$cpo$cpo, na.rm = TRUE),
  mean(ml.cbi.re.sp$cpo$cpo, na.rm = TRUE)
 )
) %>% arrange(DIC)
# Print the comparison table
print(ml_comparison.cbi)


# tidy up !
rm(A, coords_mat, field.idx, mesh, 
   ml.cbi, ml.cbi.re, ml.cbi.re2, spde.ml, stack.cbi)
gc()



#===========POSTERIOR EFFECTS===========#


#########################################
# Plot all of the posterior fixed effects
# Extract fixed effect marginals
cbi_marginals <- ml.cbi.re.sp$marginals.fixed
# Tidy marginals for all fixed effects
tidy.effects.cbi <- tibble::tibble(
 parameter = names(cbi_marginals),
 data = purrr::map(cbi_marginals, ~ as.data.frame(.x))
) %>%
 unnest(data) %>%
 filter(!parameter %in% c("(Intercept)","aspen","day_prop")) %>%  # Exclude the intercept
 mutate(
  effect = case_when(
   str_detect(parameter, "ba_live_pr") ~ "Proportion of \nLive Basal Area",
   str_detect(parameter, "qmd_live") ~ "Quadratic\nMean Diameter",
   str_detect(parameter, "H_tpp") ~ "Shannon Index\n(Abundance-based)",
   str_detect(parameter, "tree_ht_live") ~ "Tree height (m)",
   TRUE ~ parameter  # Default for all other fixed effects
  ),
  # Extract species names from parameters
  species = case_when(
   str_detect(parameter, "species_gp_n") ~ str_extract(parameter, "(?<=species_gp_n)[a-zA-Z_ñ]+"),
   str_detect(parameter, "fortypnm_gp") ~ str_extract(parameter, "(?<=fortypnm_gp)[a-zA-Z_ñ]+"),
   TRUE ~ NA_character_  # For non-species effects
  )
 ) %>%
 # Calculate mean effect size for ordering
 group_by(effect) %>%
 mutate(mean_effect = mean(x, na.rm = TRUE)) %>%
 ungroup() 
tidy.effects.cbi <- tidy.effects.cbi %>%
 # Order effects: Species-specific effects first, global effects by mean effect size
 mutate(
  effect_order = case_when(
   !is.na(species) ~ 2,  # Species-specific effects
   TRUE ~ 2              # Global effects
  ),
  effect = factor(effect, levels = tidy.effects.cbi %>%
                   arrange(effect_order, desc(mean_effect)) %>%
                   pull(effect) %>%
                   unique()),
  fill_species = ifelse(is.na(species), "Global Effect", species)
 ) %>%
 mutate(fill_species = recode(
  fill_species,
  "quaking_aspen" = "Quaking aspen",
  "lodgepole" = "Lodgepole",
  "mixed_conifer" = "Mixed-conifer",
  "piñon_juniper" = "Piñon-juniper",
  "ponderosa" = "Ponderosa",
  "spruce_fir" = "Spruce-fir"
 )) %>%
 mutate(
  fill_species = factor(fill_species, 
                        levels = c("Lodgepole", "Mixed-conifer", 
                                   "Piñon-juniper", "Ponderosa", 
                                   "Spruce-fir", "Quaking aspen")))


# check on the species name extraction
unique(tidy.effects.cbi$fill_species)
spps_breaks <- unique(tidy.effects.cbi$fill_species)


#######################
# create the ridge plot
cbi.p1 <- tidy.effects.cbi %>%
 filter(!parameter %in% c("aspen1"),
        !str_starts(parameter, "fortyp")) %>%
 ggplot(., aes(x = x, y = effect, height = y, fill = fill_species)) +
 geom_density_ridges(stat = "identity", scale = 1.5, alpha = 0.7) +
 geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
 labs(
  x = "Effect Size",
  y = "Fixed Effect",
  fill = "Species"
 ) +
 scale_fill_manual(
  values = c(
   "Global Effect" = "gray",  # Neutral color for global effects
   "Quaking aspen" = "#1b9e77",
   "Lodgepole" = "#d95f02",
   "Mixed-conifer" = "#7570b3",
   "Piñon-juniper" = "#e7298a",
   "Ponderosa" = "#66a61e",
   "Spruce-fir" = "#e6ab02"
  ),
  # Exclude "Global Effect" from the legend
  breaks = spps_breaks,  
  guide = guide_legend(title = "Species")
 ) +
 theme_minimal() +
 theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 8),
  axis.text.y = element_text(size = 10),
  axis.title.y = element_text(size = 12)
 )
cbi.p1
# Save the plot
out_png <- paste0(maindir, 'figures/INLA_ForestComposition_FixedEffects_CBI_full.png')
ggsave(out_png, dpi = 500, bg = 'white')


####################################
# plot with just the species metrics
cbi.p2 <- tidy.effects.cbi %>%
 filter(fill_species != "Global Effect") %>%
 ggplot(., aes(x = x, y = effect, height = y, fill = fill_species, alpha = fill_species)) +
 geom_density_ridges(
  stat = "identity", scale = 1.5,
  # alpha = 0.7,
  show.legend=T
 ) +
 geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
 labs(
  x = "Effect Size",
  y = "Fixed Effect",
  fill = "Species"
 ) +
 # add a subplot label
 annotate("text", x = -0.13, y = 5, 
          label = "(B)", fontface = "bold", 
          size = 4, hjust = 0) +
 scale_fill_manual(
  values = c(
   "Quaking aspen" = "#1b9e77",
   "Lodgepole" = "#d95f02",
   "Mixed-conifer" = "#7570b3",
   "Piñon-juniper" = "#e7298a",
   "Ponderosa" = "#66a61e",
   "Spruce-fir" = "#e6ab02"
  ),
  # Exclude "Global Effect" from the legend
  breaks = spps_breaks,  
  guide = guide_legend(title = "Species", override.aes = list(size = 2))
 ) +
 scale_alpha_manual(
  values = c(
   "Quaking aspen" = 0.9, 
   "Lodgepole" = 0.7,  
   "Mixed-conifer" = 0.7,
   "Piñon-juniper" = 0.7,
   "Ponderosa" = 0.7,
   "Spruce-fir" = 0.7
  ),
  guide = "none"  # Hide alpha from legend
 ) +
 coord_cartesian(xlim=c(-0.13,0.15)) +
 theme_minimal() +
 theme(
  axis.text.y = element_text(angle = 0, hjust = 1, size=9),
  axis.text.x = element_text(angle = 0, hjust = 0, size=9),
  axis.title.y = element_text(size = 10, margin = margin(r = 8)),
  axis.title.x = element_text(size = 10, margin = margin(t = 10)),
  legend.position = c(0.91, 0.74),
  legend.background = element_rect(
   fill = scales::alpha("white", 0.7), 
   color = NA, size = 1),
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 9),
  legend.key.height = unit(0.4, "cm"), 
  legend.key.width = unit(0.5, "cm"),
  legend.spacing.y = unit(0.8, "cm"),    
  legend.spacing.x = unit(0.6, "cm")
 )
cbi.p2

# Save the plot
out_png <- paste0(maindir, 'figures/INLA_ForestComposition_FixedEffects_CBI_species.png')
ggsave(out_png, dpi = 500, bg = 'white')
# 
# # Tidy up!
# rm(cbi_marginals, tidy.effects)
# gc()


# Make a combined figure (FRP and CBI)
# Load the two ridge plots (assuming they are `frp_plot` and `cbi_plot`)
frp.cbi.p3 <- frp.p2 / cbi.p2 # "/" stacks them, "|" would place them side-by-side
frp.cbi.p3

# Save the plot
out_png <- paste0(maindir, 'figures/INLA_ForestComposition_FixedEffects_FRP-CBI_species.png')
ggsave(out_png, dpi = 500, width = 8, height = 6, bg = 'white')

