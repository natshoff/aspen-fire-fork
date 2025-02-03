
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
grid_tm <-  read_csv(fp)  %>% # read in the file
 # join in some of the fire information
 left_join(fires, by="Fire_ID") %>%
 # remove missing FRP, prep columns
 # Ensure daytime FRP and CBIbc is greater than 0
 filter(
  day_count > 0, # only work with daytime FRP
  CBIbc_p90 > 0 # and CBI 0 -> these are often boundary grids
 ) %>% 
 # create a numeric fire ID
 mutate(
  # create a new unique identifier
  grid_idx = paste0(Fire_ID, grid_index),
  grid_idx = as.numeric(as.factor(grid_idx)),
  # tidy the temporal fields
  fire_ig_dt = as.Date(fire_ig_dt),  
  fire_year = year(fire_ig_dt),              
  fire_month = month(fire_ig_dt),
  first_obs_date = as.Date(first_obs_date), # the day of first observation
  first_obs_doy = as.numeric(lubridate::yday(first_obs_date)),
  fire_doy = as.factor(interaction(Fire_ID, first_obs_doy)),
  first_obs_date = as.factor(first_obs_date),
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
  # Format species names consistently
  species_gp_n = str_to_lower(species_gp_n),
  species_gp_n = str_replace_all(species_gp_n, "-", "_"),
  species_gp_n = str_replace_all(species_gp_n, " ", "_"),
  # make sure factor variables are set for species names
  fortypnm_gp = factor(fortypnm_gp),
  species_gp_n = factor(species_gp_n),
  Fire_ID = factor(Fire_ID),
  # proportion of contributing AFD during daytime observations
  day_prop = day_count / afd_count,
  # log-scaled response variables
  log_frp_max = log(frp_max + 1e-5),
  log_frp_csum = log(frp_csum + 1e-5),
  log_frp_max_day = log(frp_max_day + 1e-5),
  log_frp_max_night = log(frp_max_night + 1e-5),
  log_frp_csum_day = log(frp_csum_day + 1e-5),
  log_frp_csum_night = log(frp_csum_night + 1e-5),
  log_frp_p90 = log(frp_p90_day + 1e-5),
  log_frp_p95 = log(frp_p95_day + 1e-5),
  # scale the percentages
  forest_pct = forest_pct / 100,
  fire_aspenpct = fire_aspenpct / 100
 ) %>%
 # filter noise in species data
 filter(
  # retain forested grids (>50% forest cover)
  forest_pct >= 0.50,
  # filter where dominance/abundance is >1%
  # removes noise from extremely small proportions
  ba_live_pr > 0.01 | tpp_ld_pr > 0.05
 ) %>% 
 # Group by grid and sum dead BA and TPP separately
 group_by(grid_idx) %>%
 mutate(
  total_tpp_live = sum(tpp_live, na.rm = TRUE),
  total_ba_live = sum(ba_live, na.rm = TRUE),
  total_tpp_dead = sum(tpp_dead_pr, na.rm = TRUE),
  total_ba_dead = sum(ba_dead_pr, na.rm = TRUE),
  mean_qmd_live = mean(qmd_live, na.rm = TRUE),
  mean_qmd_dead = mean(qmd_dead, na.rm = TRUE),
  mean_tree_ht_l = mean(tree_ht_live, na.rm = TRUE)
 ) %>%
 ungroup() %>%
 group_by(Fire_ID) %>%
 mutate(aspen = ifelse(any(fortypnm_gp == "quaking_aspen"), 1, 0)) %>%
 ungroup() %>%
 # be sure there are no duplicate rows for grid/fire/species
 distinct(grid_idx, species_gp_n, .keep_all = TRUE) # remove duplicates

# load and filter the duplicate (reburn) grids
dup_grids <- st_read(paste0(maindir,"data/spatial/mod/grid_model_data_duplicates.gpkg"))
dup_idx <- unique(dup_grids$grid_index)
grid_tm <- grid_tm %>%
 filter(!grid_index %in% dup_idx)
rm(dup_grids, dup_idx)
gc()

glimpse(grid_tm)

rm(fires)

###########################################
# calculate the shannon index (H) for grids
# use both BA and TPP to calculate H
# also calculate an aspen presence column
shannon <- grid_tm %>%
 group_by(grid_idx) %>%
 mutate(
  # Replace 0 or NA proportions with a small value to avoid log issues
  ba_live_pr = ifelse(is.na(ba_live_pr) | ba_live_pr == 0, 1e-6, ba_live_pr),
  tpp_live_pr = ifelse(is.na(tpp_live_pr) | tpp_live_pr == 0, 1e-6, tpp_live_pr),
 ) %>%
 # Calculate Shannon index components
 summarise(
  H_ba = -sum(ba_live_pr * log(ba_live_pr), na.rm = TRUE),  # Based on basal area proportions
  H_tpp = -sum(tpp_live_pr * log(tpp_live_pr), na.rm = TRUE),  # Based on trees per pixel proportions
  .groups = "drop"
 )


#####################################
# get the grid-level aspen proportion
aspen_pr <- grid_tm %>%
 group_by(grid_idx) %>%
 reframe(
  # Aspen-specific live BA and TPP
  aspen_ba_live = sum(ba_live[species_gp_n == "quaking_aspen"], na.rm = TRUE),
  aspen_tpp_live = sum(tpp_live[species_gp_n == "quaking_aspen"], na.rm = TRUE),
  # Calculate proportions (avoid division by zero)
  aspen_ba_pr = if_else(total_ba_live > 0, aspen_ba_live / total_ba_live, 0),
  aspen_tpp_pr = if_else(total_tpp_live > 0, aspen_tpp_live / total_tpp_live, 0)
 ) %>%
 select(grid_idx, aspen_ba_pr, aspen_tpp_pr)
head(aspen_pr)


###########################################################
# Identify top two species by live basal area for each grid
# calculate the proportion of grids with that spp_pair
# the observed co-occurrence distribution...
top_spps <- grid_tm %>%
 mutate(species_gp_n = as.character(species_gp_n)) %>%
 group_by(grid_idx) %>%
 arrange(desc(ba_live_pr)) %>% # Sort by abundance or dominance
 slice_head(n = 2) %>% # select the top 2 species
 summarise(
  spp_1 = first(species_gp_n), # Most dominant species
  spp_2 = nth(species_gp_n, 2, default = first(species_gp_n))  # Second most dominant species
 ) %>%
 mutate(
  spp_2 = ifelse(is.na(spp_2), spp_1, spp_2),  # Handle pure stands by filling spp_2 with spp_1
  spp_pair = if_else(
   spp_1 == spp_2,
   paste0(spp_1, ".pure"),  # Label pure stands
   paste(spp_1, spp_2, sep = ".")  # Label mixed stands
  )
 ) %>%
 select(-c(spp_1,spp_2))
head(top_spps)

# Summarize unique spp pairs
n_grids <- length(unique(grid_tm$grid_idx))
spp_coo <- top_spps %>%
 group_by(spp_pair) %>%
 summarize(
  n = n(),  # Proportion of grids with co-occurrence
  .groups = "drop"
 ) %>%
 mutate(
  wt = n / n_grids,
  # Apply a smoothed scaling to balance rare/common pairs
  wt_inv = 1 / wt, # Inverse frequency scaling
  wt_log = log(wt + 1e-5),  # Optionally invert for comparison
  wt_log_inv = 1 / log(wt + 1e-5),
  wt_sq = sqrt(wt), # square root
  wt_sq_inv = 1 / sqrt(wt), # inverse square root
  wt_norm = (wt_log - min(wt_log)) / (max(wt_log) - min(wt_log))
 ) %>%
 arrange(n)
head(spp_coo,15)

dim(grid_tm)[1]
summary(spp_coo$wt)

# handle extremely rare cases ...
# these influence the variance too much on the latent effect
# gather list of rare pairings
rare_tr <- 10  # Threshold for rare species pairs (1% of grids)
rare_pairs <- spp_coo %>%
 filter(n < rare_tr) %>%
 pull(spp_pair)
rare_pairs

# check the weights distribution
ggplot(spp_coo%>%filter(!spp_pair %in% rare_pairs),
       aes(x = wt, y = wt_sq)) +
 geom_point(color = "blue") +
 labs(title = "Weight Scaling for Species Pairs", 
      x = "Fractional Occurrence (fr)", 
      y = "Scaled Weight (wt_sc)")

# merge the co-occurrence observed to the grid-level spp_pairs
spp_pairs <- top_spps %>%
 left_join(spp_coo, by="spp_pair") %>%
 arrange(n)
head(spp_pairs,10)

# Tidy up !
rm(top_spps, spp_coo)
gc()

# drop oak woodlands to align with predominant forest type model
grid_tm <- grid_tm %>%
 filter(fortypnm_gp != "oak_woodland",
        species_gp_n != "gambel_oak") %>%
 mutate(
  fortypnm_gp = droplevels(fortypnm_gp),
  species_gp_n = droplevels(species_gp_n)
 )


########################################
# merge attributes back to the grid data
grid_tm <- grid_tm %>%
 left_join(shannon, by="grid_idx") %>%
 left_join(aspen_pr, by = "grid_idx") %>%
 left_join(spp_pairs, by = "grid_idx") %>%
 # be sure there are no duplicate rows for grid/fire/species
 distinct(grid_idx, species_gp_n, .keep_all = TRUE) # remove duplicates
glimpse(grid_tm) # check the results

# Check on the grid cell counts for daytime observations
grid_counts <- grid_tm %>%
 distinct(Fire_ID, grid_idx) %>% # keep only distinct rows
 group_by(Fire_ID) %>%
 summarise(n_grids = n())
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
 pull(Fire_ID)
length(idx)

# filter the data frame to remove these fires
# also remove our rare_species
grid_tm <- grid_tm %>%
 filter(
  !Fire_ID %in% idx,
  # !spp_pair %in% rare_pairs
 )

# Check rows where spp_pair is NA
nrow(grid_tm %>% filter(is.na(spp_pair))) # should be 0

# tidy up!
rm(shannon, grid_counts, aspen_pr, spp_pairs, rare_pairs, idx)
gc()


# get the species counts
grid_tm %>% 
 group_by(species_gp_n) %>%
 summarise(n = length(species_gp_n))

# check how many grids and fires
length(unique(grid_tm$grid_idx))
length(unique(grid_tm$Fire_ID))


#==========EXPLORE THE DATA==========#

######################################
# correlation matrix for fixed effects
# Select only numeric columns
cor_da <- grid_tm %>%
 select(
  tpp_live_pr, ba_live_pr, qmd_live_pr, # species-level proportion live metrics
  aspen_ba_pr, aspen_tpp_pr, forest_pct, # grid-level aspen proportion and forest percent
  ba_live, tpp_live, qmd_live, tree_ht_live,  # species-level live forest composition metrics
  total_tpp_live, total_ba_live, mean_qmd_live, # grid-level live biomass
  total_ba_dead, total_tpp_dead, mean_qmd_dead, # grid-level dead biomass
  mean_tree_ht_l, # grid-level mean tree height
  H_ba, H_tpp, # species diversity
  erc, erc_dv, vpd, vpd_dv, # climate
  fm1000, fm1000_dv, rmin, rmin_dv, # climate
  tmmx, tmmx_dv, vs, vs_dv, #climate
  elev, slope, tpi, chili,  # topography
  forest_pct, canopypct_mean, balive_sum # grid-level forest characteristics
 ) %>%
 mutate_if(is.factor, as.numeric)  # Convert factors to numeric (if needed)

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
out_png <- paste0(maindir,'figures/INLA_CorrelationMatrix_SppComp.png')
ggsave(out_png, dpi=500, bg = 'white')



#===========MODEL SETUP==============#

# list of species names
spps <- c("quaking_aspen", "mixed_conifer", "lodgepole", 
          "ponderosa", "spruce_fir", "piñon_juniper")

# force aspen to be the baseline
grid_tm <- grid_tm %>%
 mutate(
  species_gp_n = fct_relevel(species_gp_n, spps)
 )

# check the factor levels
# make sure aspen is first
levels(grid_tm$species_gp_n)

# prep the model data frame
# center and scale fixed effects
da <- grid_tm %>%
 mutate(
  # create a grid-level live/dead proportion
  live_dead_pr_ba = total_ba_dead / (total_ba_live + total_ba_dead),
  live_dead_pr_tpp = total_tpp_dead / (total_tpp_live + total_tpp_dead),
  total_tpp = total_tpp_live + total_tpp_dead,
  # center/scale metrics / fixed effects
  across(
   c(ba_live, ba_live_pr, ba_ld, ba_ld_pr, # live basal area (dominance)
     ba_dead, ba_dead_pr, # dead basal area
     tpp_live, tpp_live_pr, tpp_ld, tpp_ld_pr, total_tpp, # tree/pixel (abundance)
     tpp_dead, tpp_dead_pr, # dead TPP
     qmd_live, qmd_live_pr, qmd_ld, qmd_ld_pr, # quadratic mean diameter
     qmd_dead, qmd_dead_pr, # dead QMD
     total_ba_live, total_tpp_live, mean_qmd_live, # grid-level live biomass
     total_ba_dead, total_tpp_dead, mean_qmd_dead, # grid-level dead biomass
     live_dead_pr_ba, live_dead_pr_tpp, # proportion dead BA/TPP
     mean_tree_ht_l, tree_ht_live, tree_ht_dead, # species- and grid-level tree height
     H_ba, H_tpp, # grid-level diversity metrics
     aspen_ba_pr, aspen_tpp_pr, # grid-level aspen proportion
     forest_pct, fortyp_pct, canopypct_mean, balive_sum, # grid-level forest characteristics
     erc, erc_dv, vpd, vpd_dv, # climate
     fm1000, fm1000_dv, rmin, rmin_dv, # climate
     tmmx, tmmx_dv, vs, vs_dv, #climate
     elev, slope, tpi, chili, # topography
     day_prop # proportion daytime observations
    ), ~ as.numeric(scale(.))
  ),
  log_fire_size = log(fire_acres) # log-scale fire size
 ) %>%
 # rename the response variables
 rename(
  fire_id = Fire_ID,
  canopypct = canopypct_mean,
  balive = balive_sum
 ) %>%
 # keep just one row per grid cell by predominant type
 distinct(grid_idx, species_gp_n, .keep_all = TRUE) %>%
 arrange(grid_idx)


# Check rows where spp_pair is NA
if (nrow(da %>% filter(is.na(spp_pair))) == 0) {
 print("No NaN values for species pair")
} # should be 0

rm(grid_tm)
gc()

# Create a wide-format data frame
da.w <- da %>%
 select(grid_index, fire_id,
        species_gp_n, fortypnm_gp,
        forest_pct, canopypct_mean, balive_sum,
        frp_day, frp_csum, cbi_p90, cbi_mn,
        ba_live_pr, qmd_live, tpp_live,
        erc_dv, slope, tpi, chili, # topography
       ) %>%
 # Pivot to wide format with species as a prefix
 pivot_wider(
  # Keep the grid-level attributes
  id_cols = c(grid_index, fire_id,
              frp_day, frp_csum, cbi_p90,
              fortypnm_gp, forest_pct,
              canopypct_mean, balive_sum,
              erc_dv, slope, tpi, chili),
  names_from = species_gp_n,
  values_from = c(ba_live_pr, qmd_live, tpp_live),
  names_sep = "_"
 ) %>%
 mutate(
  # Replace NA with 0 for all BA and TPP columns
  across(starts_with("ba_"), ~replace_na(., 0)),
  across(starts_with("tpp_"), ~replace_na(., 0)),
  across(starts_with("qmd_"), ~replace_na(., 0)),
  # Ensure values are numeric
  across(c(starts_with("ba_"), starts_with("tpp_"), starts_with("qmd_")), ~as.numeric(.))
 )
# Check the transformed dataset
glimpse(da.w)


#===========MODEL FITTING==============#

set.seed(456)

#################################################
# 1. Baseline model (no random or latent effects)

da$aspen <- as.factor(da$aspen)
levels(da$aspen)

#######
# FRP #

# setup the model formula
mf.frp <- log_frp_max_day ~ 1 + 
 fortypnm_gp:H_tpp + # species diversity (abundance) by predominant forest type
 species_gp_n:ba_live_pr + # species proportion of live basal area
 species_gp_n:qmd_live + # species-level live quadratic mean diameter
 species_gp_n:tree_ht_live + # species-level live tree height
 canopypct + # grid-level canopy percent
 erc_dv + vpd + # day-of-burn climate/weather
 slope + tpi + chili + # grid-level topography 
 day_prop  # grid-level proportion of daytime observations
 # aspen # fire-level aspen presence
# fit the model
ml.frp <- inla(
 mf.frp, data = da,
 family = "gaussian",
 control.predictor = list(compute=T),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
 control.fixed = list(
  prec = list(prior = "pc.prec", param = c(1, 0.1))
 )
)
summary(ml.frp)
# check on predictive power of the random effects models
mean(ml.frp$cpo$cpo, na.rm = TRUE)


######################################
# 2. Baseline + Temporal Random Effect

# update the model formula
mf.frp.re <- update(
 mf.frp, . ~ 1 + . + 
  f(first_obs_date, model = "iid", 
    hyper = list(
     prec = list(prior = "pc.prec", param = c(1, 0.5))
    )) # temporal random effect
)
# fit the model
ml.frp.re <- inla(
 mf.frp.re, data = da,
 family = "gaussian",
 control.predictor = list(compute=T),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
 control.fixed = list(
  prec = list(prior = "pc.prec", param = c(1, 0.1))
 )
)
summary(ml.frp.re)


###################################################
# 3. Baseline + Temporal + Fire-level Random Effect

# update the model formula
mf.frp.re2 <- update(
 mf.frp.re, . ~ 1 + . + 
  f(fire_id, model = "iid", 
    hyper = list(
     prec = list(prior = "pc.prec", param = c(1, 0.5))
    )) # fire-level random effect
)
# fit the model
ml.frp.re2 <- inla(
 mf.frp.re2, data = da,
 family = "gaussian",
 control.predictor = list(compute=T),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
 control.fixed = list(
  prec = list(prior = "pc.prec", param = c(1, 0.5))
 )
)
summary(ml.frp.re2)

gc()



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
plot(mesh, main = "SPDE Mesh for FRP Model")
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
 prior.sigma = c(1, 0.1) # variance
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
stack.frp <- inla.stack(
 data = list(log_frp_max_day = da$log_frp_max_day),
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
 mf.frp.re2, . ~ 1 + . + 
  f(mesh.idx, model = spde.ml) # spatial process model
)
# fit the model
ml.frp.re.sp <- inla(
 mf.frp.re.sp, 
 data = inla.stack.data(stack.frp), ,
 family = "gaussian",
 control.predictor = list(A = inla.stack.A(stack.frp), compute=T),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
 control.fixed = list(
  prec = list(prior = "pc.prec", param = c(1, 0.5))
 )
)
summary(ml.frp.re.sp)


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
rm(A, coords, coords_mat, field.idx, mesh, 
   ml.frp, ml.frp.re, ml.frp.re2, spde.ml, stack.frp)
gc()



#===========POSTERIOR EFFECTS===========#


#########################################
# Plot all of the posterior fixed effects
# Extract fixed effect marginals
frp_marginals <- ml.frp.re.sp$marginals.fixed
# Tidy marginals for all fixed effects
tidy.effects <- tibble::tibble(
 parameter = names(frp_marginals),
 data = purrr::map(frp_marginals, ~ as.data.frame(.x))
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
tidy.effects <- tidy.effects %>%
 # Order effects: Species-specific effects first, global effects by mean effect size
 mutate(
  effect_order = case_when(
   !is.na(species) ~ 2,  # Species-specific effects
   TRUE ~ 2              # Global effects
  ),
  effect = factor(effect, levels = tidy.effects %>%
                   arrange(effect_order, desc(mean_effect)) %>%
                   pull(effect) %>%
                   unique()),
  fill_species = ifelse(is.na(species), "Global Effect", species)
 )

# check on the species name extraction
unique(tidy.effects$fill_species)
spps_breaks <- unique(tidy.effects$fill_species)


#######################
# create the ridge plot
ggplot(tidy.effects, aes(x = x, y = effect, height = y, fill = fill_species)) +
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
   "quaking_aspen" = "#1b9e77",
   "lodgepole" = "#d95f02",
   "mixed_conifer" = "#7570b3",
   "piñon_juniper" = "#e7298a",
   "ponderosa" = "#66a61e",
   "spruce_fir" = "#e6ab02"
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

# Save the plot
out_png <- paste0(maindir, 'figures/INLA_ForestComposition_FixedEffects_full.png')
ggsave(out_png, dpi = 500, bg = 'white')

####################################
# plot with just the species metrics
tidy.effects %>%
 filter(fill_species != "Global Effect") %>%
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
   "quaking_aspen" = "#1b9e77",
   "lodgepole" = "#d95f02",
   "mixed_conifer" = "#7570b3",
   "piñon_juniper" = "#e7298a",
   "ponderosa" = "#66a61e",
   "spruce_fir" = "#e6ab02"
  ),
  # Exclude "Global Effect" from the legend
  breaks = spps_breaks,  
  guide = guide_legend(title = "Species")
 ) +
 theme_minimal() +
 theme(
  legend.title = element_text(size = 12),
  legend.text = element_text(size = 11),
  axis.text = element_text(size = 11),
  axis.title = element_text(size = 12)
 )

# Save the plot
out_png <- paste0(maindir, 'figures/INLA_ForestComposition_FixedEffects_FRP_species.png')
ggsave(out_png, dpi = 500, bg = 'white')

# Tidy up!
rm(frp_marginals, tidy.effects)
gc()




#####################
# Fit a model for CBI
# define a new model formula
mf.cbi.re.sp <- update(mf.frp.re.sp, cbi_p90 ~ 1 + .)

# fit the model
model_tm.cbi.re <- inla(
 mf.cbi.re.sp, data = da,
 family = "gaussian",
 control.predictor = list(compute=T),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
 control.fixed = list(prec = 0.01)  # Penalize high precision
)
summary(model_tm.cbi.re)


#########################################
# Plot all of the posterior fixed effects
# Extract fixed effect marginals
cbi_marginals <- model_tm.cbi.re$marginals.fixed
# Tidy marginals for all fixed effects
tidy.effects <- tibble::tibble(
 parameter = names(cbi_marginals),
 data = purrr::map(cbi_marginals, ~ as.data.frame(.x))
) %>%
 unnest(data) %>%
 filter(parameter != "(Intercept)") %>%  # Exclude the intercept
 mutate(
  effect = case_when(
   str_detect(parameter, "ba_live_pr") ~ "Proportion of\n Live BA",
   str_detect(parameter, "qmd_live") ~ "Mean. Live QMD",
   str_detect(parameter, "H_tpp") ~ "H (TPP)",
   str_detect(parameter, "forest_pct") ~ "Forest %",
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

tidy.effects <- tidy.effects %>%
 # Order effects: Species-specific effects first, global effects by mean effect size
 mutate(
  effect_order = case_when(
   !is.na(species) ~ 2,  # Species-specific effects
   TRUE ~ 2              # Global effects
  ),
  effect = factor(effect, levels = tidy.effects %>%
                   arrange(effect_order, desc(mean_effect)) %>%
                   pull(effect) %>%
                   unique()),
  fill_species = ifelse(is.na(species), "Global Effect", species)
 )

# check on the species name extraction
unique(tidy.effects$species)

# create the ridge plot
ggplot(tidy.effects, aes(x = x, y = effect, height = y, fill = fill_species)) +
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
   "aspen" = "#1b9e77",
   "lodgepole" = "#d95f02",
   "mixed_conifer" = "#7570b3",
   "piñon_juniper" = "#e7298a",
   "ponderosa" = "#66a61e",
   "spruce_fir" = "#e6ab02"
  ),
  # Exclude "Global Effect" from the legend
  breaks = c("aspen", "lodgepole", "mixed_conifer", "piñon_juniper", "ponderosa", "spruce_fir"),  
  guide = guide_legend(title = "Species")
 ) +
 theme_minimal() +
 theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 8),
  axis.text.y = element_text(size = 10),
  axis.title.y = element_text(size = 12)
 )

# Save the plot
out_png <- paste0(maindir, 'figures/INLA_Spp_COO_FixedEffects-CBI.png')
ggsave(out_png, dpi = 500, bg = 'white')

# Tidy up!
rm(cbi_marginals, tidy.effects)
gc()


#########################################################
# Extract posterior summaries for spp_pair latent effects

# create a lookup table for the species pair names
pair_lookup <- data.frame(
 spp_pair_id = seq_along(levels(da$spp_pair)),  # Numeric indices
 spp_pair_name = levels(da$spp_pair)  # Original species pair names
)

# Extract the posterior distribution of effects
spp_pair_effects <- as.data.frame(model_tm.cbi.re$summary.random$spp_pair)
spp_pair_effects$spp_pair_id <- as.numeric(rownames(spp_pair_effects))

# Add confidence intervals and significance
spp_pair_effects <- spp_pair_effects %>%
 mutate(
  lower_ci = `0.025quant`,
  upper_ci = `0.975quant`,
  significant = (lower_ci > 0 | upper_ci < 0),  # Significant if CI does not overlap 0
  spp_pair = rownames(spp_pair_effects)  # Extract spp_pair names
 ) %>%
 arrange(mean) %>%
 filter(significant) %>% # keep significant pairs
 left_join(pair_lookup, by="spp_pair_id") %>%
 mutate(spp_pair_name = factor(spp_pair_name, levels = spp_pair_name[order(mean)]))

# Create the plot
ggplot(spp_pair_effects, aes(x = mean, y = spp_pair_name, fill = mean > 0)) +
 geom_col(
  aes(x = mean), 
  position = "identity",
  width = 0.8,
  alpha = 0.7
 ) +
 geom_point(aes(x = mean), color = "black", size = 2) +
 geom_errorbarh(aes(xmin = `0.025quant`, xmax = `0.975quant`), height = 0.2, color = "black") +
 scale_fill_manual(values = c("TRUE" = "coral", "FALSE" = "skyblue"), guide = "none") +
 labs(
  x = "Effect Size (Mean ± 95% CI)",
  y = "Species Pair",
 ) +
 theme_minimal() +
 theme(
  legend.position = "none",
  axis.text.y = element_text(size = 8),
  axis.text.x = element_text(size = 10),
  axis.title = element_text(size = 12)
 )

# Save the plot
out_png <- paste0(maindir, 'figures/INLA_Spp_COO_LatentEffects-CBI.png')
ggsave(out_png, dpi = 500, bg = 'white')

rm(spp_pair_effects, pair_lookup)
gc()




