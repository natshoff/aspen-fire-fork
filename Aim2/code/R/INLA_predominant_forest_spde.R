
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
library(gstat)
library(terra)
library(viridis)

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
 # do some filtering of grids
 filter(
  day_count > 0, # only consider grids with some daytime observations
  frp_max > 0, # make sure some daytime FRP values
  CBIbc_p90 > 0, # and CBI 0 -> these are often boundary grids
  overlap >= 0.5 # (optional) filter only grids with at least 50% detection overlap
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
  # make sure factor variables are set for species names
  fortypnm_gp = as.factor(fortypnm_gp),
  Fire_ID = as.factor(Fire_ID),
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
  forest_pct = forest_pct / 100,
  fire_aspenpct = fire_aspenpct / 100,
  ) %>%
 group_by(Fire_ID) %>%
 mutate(aspen = ifelse(any(fortypnm_gp == "quaking_aspen"), 1, 0)) %>%
 ungroup() %>%
 # be sure there are no duplicate rows for grid/fire/species
 distinct(grid_idx, fortypnm_gp, .keep_all = TRUE) # remove duplicates
glimpse(grid_tm) # check the results

rm(fires) # clean up the fire perimeters

# Check how many grids are > 50% forested
print("Number of forested grids (>50% forest pixels):")
dim(grid_tm %>% filter(forest_pct > 0.50))[1]/dim(grid_tm)[1]

# Check the proportion of daytime observations
summary(grid_tm$day_prop)
dim(grid_tm %>% filter(day_prop > 0.50))[1]

# check the FRP stats
summary(grid_tm$log_frp_csum)
summary(grid_tm$log_frp_csum_day)
summary(grid_tm$log_frp_max)
summary(grid_tm$log_frp_max_day)


# Remove duplicate "grid_index"
# These are reburns ...
duplicate_grids <- grid_tm %>%
 group_by(grid_index) %>%
 filter(n() > 1) %>%
 arrange(grid_index, Fire_ID)
# save a spatial polygon of the duplicated grids
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

# calculate the dominant forest type for the fire
dom_fortyp <- grid_tm %>%
 group_by(Fire_ID) %>%
 count(fortypnm_gp) %>%
 slice_max(n, n = 1, with_ties = FALSE) %>%
 rename(dom_fortyp = fortypnm_gp)
# join back to grid_tm
grid_tm <- grid_tm %>%
 left_join(dom_fortyp, by = "Fire_ID") %>%
 mutate(fire_fortyp = as.factor(interaction(Fire_ID, dom_fortyp)))
head(grid_tm%>%select(Fire_ID,fortypnm_gp, dom_fortyp))
rm(dom_fortyp)



#===============Explore Distributions, etc.================#

####################################
# distribution of response variables
resp_plot <- grid_tm %>%
 # pivot longer to facet plot
 pivot_longer(cols = c(log_frp_max,
                       log_frp_max_day,
                       log_frp_csum, 
                       log_frp_csum_day,
                       CBIbc_p90, 
                       CBIbc_p99),
              names_to = "variable",
              values_to = "value") %>%
 # Plot with facets
 ggplot(aes(x = value)) +
 geom_histogram(bins = 30, fill = "orange", alpha = 0.7) +
 facet_wrap(
  ~ variable, 
  scales = "free",
  labeller = as_labeller(c(log_frp_max = "log(Max FRP)",
                           log_frp_max_day = "log(Daytime Max FRP)",
                           log_frp_csum = "log(Cumulative FRP)",
                           log_frp_csum_day = "log(Daytime Cumulative FRP)",
                           CBIbc_p90 = "90th Percentile CBIbc",
                           CBIbc_p99 = "99th Percentile CBIbc"))) +
 labs(x = "value",
      y = "Frequency") +
 theme_minimal()
resp_plot

# save the plot.
out_png <- paste0(maindir,'figures/INLA_ResponseDistribution.png')
ggsave(out_png, plot = resp_plot, dpi=500, bg = 'white')
rm(resp_plot)

####################
# correlation matrix
# Select only numeric columns and convert factors to dummy variables
# this correlation matrix is for this simple model (i.e., not forest composition)
cor_da <- grid_tm %>%
 select(
  fortypnm_gp, fire_acres, # forest type (factor)
  forest_pct, fortyp_pct, # forest and forest type percent
  canopypct_mean, balive_sum, # grid-level mean canopy percent and balive
  erc, erc_dv, vpd, vpd_dv, # climate
  fm1000, fm1000_dv, rmin, rmin_dv, # climate
  tmmx, tmmx_dv, vs, vs_dv, #climate
  elev, slope, tpi, chili,  # topography
  fire_aspenpct # fire-level aspen percent
 ) %>%
 pivot_wider(
  names_from = fortypnm_gp, 
  values_from = fortyp_pct, 
  values_fill = 0) %>%
 mutate(across(everything(), ~ scale(.) %>% as.numeric()))  # Standardize variables

# Compute correlation matrix
cor_mat <- cor(cor_da, use = "complete.obs", method = "spearman")

# Plot correlation matrix
cor_plot <- ggcorrplot(
 cor_mat,
 method = "circle",  # Circle or square for visualization
 type = "lower",  # Lower triangle of the correlation matrix
 lab = TRUE,  # Show correlation values
 lab_size = 3,
 tl.cex = 10,  # Text label size
 colors = c("blue", "white", "red")  # Color gradient
)
cor_plot

rm(cor_da, cor_mat)
gc()

# save the plot.
out_png <- paste0(maindir,'figures/INLA_CorrelationMatrix_Fortyp.png')
ggsave(out_png, plot = cor_plot, dpi=500, bg = 'white')
rm(cor_plot)


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

# check the distribution of dominant forest type percent
summary(grid_tm$fortyp_pct)
qt <- tibble::enframe(
 quantile(grid_tm$fortyp_pct, probs = seq(.1, .9, by = .1)),
 name = "qt", value = "val"
)
qt
qt[2,]$val

# create a data frame for just the predominant forest type
# one row per grid_index with the dominant forest type
# filter to where that type is at least 50% of the grid forested area
da <- grid_tm %>%
 # keep only grids where predominant species cover at least 50%
 # filter grids below the 20th percentile in forest type percent
 filter(
  (forest_pct > 0.50) & (fortyp_pct > qt[2,]$val),
 ) %>%
 # filter((forest_pct >= 0.50) |(fortyp_pct >= 0.50)) %>%
 # filter(fortyp_pct > qt[2,]$val) %>%
 # select the columns we need for modeling
 select(grid_index, grid_idx, Fire_ID, fire_acres, # ID columns
        day_prop, overlap, afd_count, # proportion daytime observations and percent detection overlap
        fire_year, season, year_season, first_obs_date, fire_doy, # temporal effects
        log_frp_max, log_frp_max_day, # FRP response variables
        log_frp_csum, log_frp_csum_day, # FRP response variables
        log_frp_p90, log_frp_p95, # FRP response variables
        CBIbc_p90, CBIbc_p95, CBIbc_p99, CBIbc_mean, # CBI response variables
        fortypnm_gp, fortyp_pct, forest_pct, # forest type and percent cover
        canopypct_mean, balive_sum, # canopy percent and total live basal area
        erc, erc_dv, vpd, vpd_dv, # climate
        fm1000, fm1000_dv, rmin, rmin_dv, # climate
        tmmx, tmmx_dv, vs, vs_dv, #climate
        elev, slope, tpi, chili, # topography
        x, y, # grid centroid coordinate for spatial fields model
        aspen, fire_aspenpct, # fire-level aspen percent/presence 
        dom_fortyp, fire_fortyp # fire-level dominant forest type
 ) %>%
 # center and scale continuous predictor variables
 mutate(
  across(
   c(
    forest_pct, fortyp_pct, day_prop, overlap, afd_count,
    canopypct_mean, balive_sum, # canopy percent and total live basal area
    erc, vpd, vpd_dv, erc_dv, # climate
    fm1000, fm1000_dv, rmin, rmin_dv, # climate
    tmmx, tmmx_dv, vs, vs_dv, #climate
    elev, slope, tpi, chili, # topography
    fire_aspenpct # fire-level aspen percent
   ),
   ~ as.numeric(scale(.))
  ),
  log_fire_size = log(fire_acres),
 ) %>%
 # rename the response variables
 rename(
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


#################################################
# Count the number of unique fortypnm_gp per fire
da %>%
 group_by(fire_id) %>%
 summarize(n_types = n_distinct(fortypnm_gp)) %>%
 summarize(
  mean_types = mean(n_types),
  median_types = median(n_types),
  min_types = min(n_types),
  max_types = max(n_types)
 )

# Compute the number of forest types per fire
fire_div <- da %>%
 group_by(fire_id) %>%
 summarize(n_types = n_distinct(fortypnm_gp),
           fire_size = first(log_fire_size),  # or another metric
           mean_FRP = mean(log_frp_csum))  # Avg fire radiative power

# Check if `fire_id` is absorbing forest type effects
ggplot(fire_div, aes(x = n_types, y = mean_FRP)) +
 geom_point(alpha = 0.5) +
 geom_smooth(method = "lm", se = FALSE, color = "red") +
 labs(x = "Number of Forest Types per Fire", y = "Mean FRP",
      title = "Does Fire-Level Forest Type Diversity Influence FRP?") +
 theme_minimal()

rm(fire_div)
gc()

#################################################
# Subset species with too few observations
print("Dropping small classes:")
print(dim(da%>%filter(fortypnm_gp == "oak_woodland")))[1]
da <- da %>%
 filter(!fortypnm_gp == "oak_woodland") %>%
 mutate(fortypnm_gp = droplevels(fortypnm_gp),
        fire_fortyp = droplevels(fire_fortyp))
levels(da$fortypnm_gp)

# save a spatial polygon 
grid <- st_read(paste0(maindir,"data/spatial/mod/VIIRS/viirs_snpp_jpss1_afd_latlon_fires_pixar_gridstats.gpkg"))
grid <- da %>%
 left_join(grid %>% select(grid_index, geom), by="grid_index") %>%
 st_as_sf(.)
st_write(grid,paste0(maindir,"data/spatial/mod/grid_model_data.gpkg"), append=F)

# Tidy up!
rm(grid, grid_tm, grid_counts, idx)
gc()



#===========MODEL FITTING==============#

set.seed(435)

########################
# FIRE RADIATIVE POWER #

#######################################
# 1. Baseline model (no latent effects)

da$aspen <- as.factor(da$aspen)
levels(da$aspen)

# Set up the model formula
mf.frp <- log_frp_csum_day ~ 1 +
 fortypnm_gp + # predominant forest type
 canopypct + # grid-level mean canopy percent
 erc_dv + vpd + vs + # climate/weather
 elev + slope + tpi + chili + # topography
 overlap + # grid-level VIIRS overlap (sum)
 day_prop + # proportion of daytime detections 
 afd_count + # number of contributing detections
 aspen # fire-level aspen presence
# fit the model                     
ml.frp <- inla(
 mf.frp, data = da, 
 family = "gaussian",
 control.predictor = list(compute = TRUE),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
 control.fixed = list(
  prec = list(prior = "pc.prec", param = c(1, 0.01))
 ),
 control.family = list(
  hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))
 ),
 control.inla = list(strategy = "adaptive", int.strategy = "grid")
)
summary(ml.frp)
# check on predictive power of the random effects models
mean(ml.frp$cpo$cpo, na.rm = TRUE)


############################################
# 2. Baseline model + temporal random effect
# temporal effect for the first burn day

# update the model formula
mf.frp.re <- update(
 mf.frp, . ~ . + 
  f(first_obs_date, model = "iid", 
    hyper = list(
     prec = list(prior = "pc.prec", param = c(0.5, 0.01))
   )) # temporal random effect
)
# fit the model                     
ml.frp.re <- inla(
 mf.frp.re, data = da, 
 family = "gaussian",
 control.predictor = list(compute = TRUE),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
 control.fixed = list(
  prec = list(prior = "pc.prec", param = c(1, 0.01))
 ),
 control.family = list(
  hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))
 ),
 control.inla = list(strategy = "adaptive", int.strategy = "grid")
)
summary(ml.frp.re)
# check on predictive power of the random effects models
mean(ml.frp.re$cpo$cpo, na.rm = TRUE)


##############################################
# 3. Baseline model + fire-level random effect

# update the model formula
mf.frp.re2 <- update(
 mf.frp.re, . ~ . + 
  f(fire_id, model = "iid", 
    hyper = list(
     prec = list(prior = "pc.prec", param = c(0.5, 0.01))
  )) # fire-level random effect
)
# fit the model                     
ml.frp.re2 <- inla(
 mf.frp.re2, data = da, 
 family = "gaussian",
 control.predictor = list(compute = TRUE),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
 control.fixed = list(
  prec = list(prior = "pc.prec", param = c(1, 0.01))
 ),
 control.family = list(
  hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))
 ),
 control.inla = list(strategy = "adaptive", int.strategy = "grid")
)
summary(ml.frp.re2)
# check on predictive power of the random effects models
mean(ml.frp.re2$cpo$cpo, na.rm = TRUE)



##########################################
# 4. Spatial model using the SPDE approach

# extracting spatial coordinates for grid centroids
grid_sf <- da %>%
 arrange(grid_index) %>%
 distinct(grid_index, x, y, .keep_all = TRUE) %>%
 st_as_sf(., coords = c("x", "y"), crs = 4326)
# extract coordinates
coords <- grid_sf %>% st_coordinates(.)
st_write(grid_sf, paste0(maindir,"data/spatial/mod/model_grid_centroids.gpkg"),
         append = F)

##########################################################
# fit a semivariogram to look at spatial dependence in FRP
##########################################################

# reproject
sp <-  grid_sf %>%
 st_transform(st_crs(5070))

# Function to compute semivariogram for each fire
get_fire_vario <- function(fire_id, data) {
 fire_sp <- data %>% filter(fire_id == !!fire_id)
 
 if (nrow(fire_sp) > 10) {  # Only compute if there are enough points
  vario <- variogram(log_frp_csum_day ~ 1, data = fire_sp)
  vario <- vario %>% mutate(fire_id = fire_id)  # Ensure fire_id is added correctly
  print(paste("Computed variogram for fire:", fire_id, "with", nrow(fire_sp), "points"))
  return(vario)
 } else {
  print(paste("Skipping fire:", fire_id, "- Not enough points"))
  return(NULL)
 }
}

# Apply to each fire separately
fires <- unique(sp$fire_id)
variograms <- lapply(fires, get_fire_vario, data = sp)
variograms <- do.call(rbind, variograms)  # Combine results

# Check if fire_id is still in the dataset
print(unique(variograms$fire_id))

# Plot variograms grouped by Fire ID
ggplot(variograms, aes(x = dist, y = gamma, 
                group = fire_id, color = fire_id)) +
 geom_line(alpha = 0.3, size=0.9) +  # Light transparency to see overlap
 labs(x = "Distance (meters)", y = "Semivariance", 
      title = "Semivariograms for Individual Fires") +
 theme_minimal() +
 theme(legend.position = "none")

# save the plot.
out_png <- paste0(maindir,'figures/INLA_Fire_SemiVariograms_FRP.png')
ggsave(out_png, dpi=500, bg = 'white')

# Find the approximate range (where semivariance plateaus) for each fire
range_est <- variograms %>%
 group_by(fire_id) %>%
 summarize(range_m = max(dist[gamma < max(gamma) * 0.9], na.rm = TRUE))  # 90% of max
qt <- quantile(range_est$range_m, probs = seq(.1, .9, by = .1))
qt
# Histogram of estimated spatial ranges
ggplot(range_est, aes(x = range_m)) +
 geom_histogram(bins = 20, fill = "blue", alpha = 0.6) +
 labs(x = "Estimated Range (meters)", y = "Count",
      title = "Distribution of Within-Fire Spatial Correlation Ranges") +
 theme_minimal()

# save the plot.
out_png <- paste0(maindir,'figures/INLA_Fire_EstimatedRange_FRP.png')
ggsave(out_png, dpi=500, bg = 'white')

rm(sp, range_est, variograms)
gc()

###################################################
###################################################

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


###################################################
# Save a spatial file of the mesh grid and vertices
###################################################

####### Save the mesh as spatial objects
# Convert INLA mesh vertices to an sf object
mesh_v <- data.frame(mesh$loc) # Extract coordinates
colnames(mesh_v) <- c("x", "y") # set column names
mesh_sf <- st_as_sf(mesh_v, coords = c("x", "y"), crs = 4326)
# Export vertices to a gpkg
st_write(mesh_sf, paste0(maindir,"data/spatial/mod/INLA_mesh_vertices.gpkg"), 
         layer = "mesh_vertices", driver = "GPKG", append=FALSE)
# Function to convert triangles to polygons
mesh_t <- mesh$graph$tv  # Triangle indices
tri_coords <- mesh$loc  # Mesh coordinates
# Create list of polygons
poly_list <- lapply(1:nrow(mesh_t), function(i) {
 coords <- tri_coords[mesh_t[i, ], ]
 coords <- rbind(coords, coords[1, ])  # Close the polygon
 st_polygon(list(coords))
})
# Convert to sf object
mesh_poly_sf <- st_sfc(poly_list, crs = 4326) %>% st_sf(geometry = .)
# Export mesh triangles to GPKG
st_write(mesh_poly_sf, paste0(maindir,"data/spatial/mod/INLA_mesh_triangles.gpkg"), 
         layer = "mesh_triangles", driver = "GPKG", append=FALSE)
rm(mesh_v, mesh_sf, mesh_t, tri_coords, poly_list, mesh_poly_sf)
gc()

########################################
########################################

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
# check on the A matrix
dim(A)  # Should be (n_data_locations, n_mesh_vertices)
sum(rowSums(A) == 0)  # Should be 0 (no unassigned rows)
table(rowSums(A > 0))  # Should not have zero rows
table(colSums(A) > 0)  # Ensure all mesh nodes are assigned

# Assign the spatial index
field.idx <- inla.spde.make.index(
 name = "mesh.idx",
 n.spde = spde.ml$n.spde
)
str(field.idx)

# Create the FRP INLA stack
# Isolate the predictor variables first
# Create a data frame for the predictors
# baseline model includes predominant forest type, climate, and topography
X <- da %>% 
 select(fire_id, first_obs_date, fire_doy, overlap, afd_count,
        forest_pct, fortyp_pct, fortypnm_gp, canopypct, 
        erc_dv, vpd, vs, slope, tpi, chili, elev,
        aspen, day_prop)
head(X)

# Create the INLA data stack
stack.frp <- inla.stack(
 data = list(log_frp_csum_day = da$log_frp_csum_day),
 A = list(A, 1),  
 tag = 'est',
 effects = list(
  c(field.idx),
  list(as.data.frame(X))
 )
)
dim(inla.stack.A(stack.frp))

##########################
# update the model formula
mf.frp.sp <- update(
 mf.frp.re2, . ~ . + 
  f(mesh.idx, 
    model = spde.ml
  )
)
# fit the model
ml.frp.re.sp <- inla(
 mf.frp.sp, 
 data = inla.stack.data(stack.frp), 
 family = "gaussian",
 control.predictor = list(A = inla.stack.A(stack.frp), compute = TRUE),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
 control.fixed = list(
  prec = list(prior = "pc.prec", param = c(1, 0.01))
 ),
 control.family = list(
  hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))
 ),
 control.inla = list(strategy = "adaptive", int.strategy = "grid")
)
summary(ml.frp.re.sp)
# check on predictive power of the random effects models
mean(ml.frp.re.sp$cpo$cpo, na.rm = TRUE)



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


#==========PLOTTING THE SPATIAL EFFECT===========#

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
 st_as_sf(., coords = c("x", "y"), crs = 4326)
glimpse(spat.eff.frp)
# write this a spatial file
st_write(spat.eff.frp, paste0(maindir,"data/spatial/mod/spatial_effect_FRP.gpkg"), append=F)

# Tidy up!
rm(spat.eff, lower, upper, uncertainty)
gc()


###########
# Tidy up !
rm(coords, coords_mat, mesh, spde.ml, grid_sf, field.idx, A, stack.frp,
   ml.frp, ml.frp.re, ml.frp.re2, X)
gc()



#######################################
#######################################

########################
# COMPOSITE BURN INDEX #

# examine the distributions of CBIbc variables
da <- da %>%
 mutate(
  log_cbi_mn = log(CBIbc_mean),
  sqrt_cbi_p90 = sqrt(CBIbc_p90),
  sqrt_cbi_p95 = sqrt(CBIbc_p95),
  sqrt_cbi_p99 = sqrt(CBIbc_p99)
 ) 
da %>%
 # pivot longer to facet plot
 pivot_longer(cols = c(CBIbc_mean, log_cbi_mn,
                       CBIbc_p90, sqrt_cbi_p90,
                       CBIbc_p95, sqrt_cbi_p95,
                       CBIbc_p99, sqrt_cbi_p99),
              names_to = "variable",
              values_to = "value") %>%
 # Plot with facets
 ggplot(aes(x = value)) +
 geom_histogram(bins = 30, fill = "orange", alpha = 0.7) +
 facet_wrap(
  ~ variable, 
  scales = "free",
  labeller = as_labeller(c(CBIbc_mean = "Average CBIbc",
                           log_cbi_mn = "log(Average CBIbc)",
                           CBIbc_p90 = "90th Percentile CBIbc",
                           CBIbc_p95 = "95th Percentile CBIbc",
                           CBIbc_p99 = "99th Percentile CBIbc",
                           sqrt_cbi_p90 = "sqrt(90th Percentile CBIbc)",
                           sqrt_cbi_p95 = "sqrt(95th Percentile CBIbc)",
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
mf.cbi <- CBIbc_p90 ~ 
 fortypnm_gp + # predominant forest type
 canopypct + # grid-level mean canopy percent
 erc_dv + vpd + vs + # climate/weather
 slope + tpi + chili + # topography
 aspen # fire-level aspen presence
# fit the model                     
ml.cbi <- inla(
 mf.cbi, data = da, 
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



########################################################
# 4. Baseline model + fire-level + spatial random effect

# extract spatial coordinates for grid centroids
grid_sf <- da %>%
 arrange(grid_index) %>%
 distinct(grid_index, x, y, .keep_all = TRUE) %>%
 st_as_sf(., coords = c("x", "y"), crs = 4326)
# extract coordinates
coords <- grid_sf %>% st_coordinates(.)


###################################################
# fit a semivariogram to look at spatial dependence
###################################################

# reproject
sp <-  grid_sf %>%
 st_transform(st_crs(5070))

# Function to compute semivariogram for each fire
get_fire_vario <- function(fire_id, data) {
 fire_sp <- data %>% filter(fire_id == !!fire_id)
 
 if (nrow(fire_sp) > 20) {  # Only compute if there are enough points
  vario <- variogram(CBIbc_p90 ~ 1, data = fire_sp)
  vario <- vario %>% mutate(fire_id = fire_id)  # Ensure fire_id is added correctly
  print(paste("Computed variogram for fire:", fire_id, "with", nrow(fire_sp), "points"))
  return(vario)
 } else {
  print(paste("Skipping fire:", fire_id, "- Not enough points"))
  return(NULL)
 }
}

# Apply to each fire separately
fires <- unique(sp$fire_id)
variograms <- lapply(fires, get_fire_vario, data = sp)
variograms <- do.call(rbind, variograms)  # Combine results

# Check if fire_id is still in the dataset
print(unique(variograms$fire_id))

# Plot variograms grouped by Fire ID
ggplot(variograms, aes(x = dist, y = gamma, 
                       group = fire_id, color = fire_id)) +
 geom_line(alpha = 0.3, size=0.9) +  # Light transparency to see overlap
 labs(x = "Distance (meters)", y = "Semivariance", 
      title = "Semivariograms for Individual Fires") +
 theme_minimal() +
 theme(legend.position = "none")

# save the plot.
out_png <- paste0(maindir,'figures/INLA_Fire_SemiVariograms_CBI.png')
ggsave(out_png, dpi=500, bg = 'white')

# Find the approximate range (where semivariance plateaus) for each fire
range_est <- variograms %>%
 group_by(fire_id) %>%
 summarize(range_m = max(dist[gamma < max(gamma) * 0.9], na.rm = TRUE))  # 90% of max
qt <- quantile(range_est$range_m, probs = seq(.1, .9, by = .1))
qt
# Histogram of estimated spatial ranges
ggplot(range_est, aes(x = range_m)) +
 geom_histogram(bins = 20, fill = "blue", alpha = 0.6) +
 labs(x = "Estimated Range (meters)", y = "Count",
      title = "Distribution of Within-Fire Spatial Correlation Ranges") +
 theme_minimal()

# save the plot.
out_png <- paste0(maindir,'figures/INLA_Fire_EstimatedRange_CBI.png')
ggsave(out_png, dpi=500, bg = 'white')

rm(sp, range_est, variograms)
gc()

###################################################
###################################################

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
plot(mesh, main = "SPDE Mesh for CBI Model")
points(coords, col = "red", pch = 20)


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
# check on the A matrix
dim(A)  # Should be (n_data_locations, n_mesh_vertices)
sum(rowSums(A) == 0)  # Should be 0 (no unassigned rows)
table(rowSums(A > 0))  # Should not have zero rows
table(colSums(A) > 0)  # Ensure all mesh nodes are assigned

# Assign the spatial index
field.idx <- inla.spde.make.index(
 name = "mesh.idx",
 n.spde = spde.ml$n.spde
)
str(field.idx)

# Create the FRP INLA stack
# Isolate the predictor variables first
# Create a data frame for the predictors
# baseline model includes predominant forest type, climate, and topography
X <- da %>% 
 select(fire_id, first_obs_date, forest_pct, fortyp_pct, 
        fortypnm_gp, canopypct, erc_dv, vpd, vs,
        slope, tpi, chili, aspen, day_prop)
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

##############################
# update the model formula
mf.cbi.sp <- update(
 mf.cbi.re2, . ~ . + 
  f(mesh.idx, 
    model = spde.ml
  )
)
# fit the model
ml.cbi.re.sp <- inla(
 mf.cbi.sp, 
 data = inla.stack.data(stack.cbi), 
 family = "gamma",
 control.predictor = list(A = inla.stack.A(stack.cbi), compute = TRUE),
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



#==========PLOTTING THE SPATIAL EFFECT===========#

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
 st_as_sf(., coords = c("x", "y"), crs = 4326)
glimpse(spat.eff.cbi)
# write this a spatial file
st_write(spat.eff.cbi, paste0(maindir,"data/spatial/mod/spatial_effect_CBI.gpkg"), append=F)

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
) %>% arrange(WAIC)
# Print the comparison table
print(ml_comparison.cbi)


# Tidy up !
rm(coords, coords_mat, ml.cbi, ml.cbi.re, ml.cbi.re2,
   stack.cbi, mesh, field.idx, A, spde.ml, grid_sf, X)
gc()



#===========Plotting Posterior Effects============#

# Extract fixed effects related to fortypnm_gp

#####
# FRP
frp.eff <- as.data.frame(ml.frp.re.sp$summary.fixed) %>%
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
cbi.eff <- as.data.frame(ml.cbi.re.sp$summary.fixed) %>%
 rownames_to_column(var = "parameter") %>%
 filter(grepl("fortypnm_gp", parameter)) %>%
 mutate(
  response = "CBI",
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
cmap <- c("CBI" = "#800026", "FRP" = "#FEB24C") # color map for the response
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
   "CBI" = "90th Percentile CBIbc", 
   "FRP" = "Cumulative Daytime FRP")
 ) +
 # Adjust theme
 theme_bw() +
 theme(
  axis.text.x = element_text(angle = 45, hjust = 1, size=10),
  axis.text.y = element_text(angle = 0, hjust = 0, size=10),
  legend.position.inside = c(0.18, 0.16),  # Updated syntax for legend position
  legend.background = element_rect(
   fill = scales::alpha("white", 0.2), 
   color = NA, linewidth = 0.5),  # Use `linewidth` instead of `size`
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
frp_marginals <- ml.frp.re.sp$marginals.fixed
cbi_marginals <- ml.cbi.re.sp$marginals.fixed

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

# Filter for forest type effects
fortyp_marginals <- tidy_combined %>%
 filter(str_detect(parameter, "fortypnm_gp")) %>%
 mutate(forest_type = str_remove(parameter, "fortypnm_gp"),
        forest_type = recode(
         forest_type,
         "mixed_conifer" = "Mixed Conifer",
         "lodgepole" = "Lodgepole",
         "ponderosa" = "Ponderosa",
         "spruce_fir" = "Spruce-Fir",
         "piñon_juniper" = "Piñon-Juniper"
        ))
 

# Compute mean effect sizes for reordering
effect_means <- fortyp_marginals %>%
 group_by(forest_type, response) %>%
 summarize(mean_effect = mean(x), .groups = "drop") %>%
 filter(response == "FRP")

fortyp_marginals <- fortyp_marginals %>%
 left_join(effect_means %>% select(forest_type, mean_effect), by = "forest_type") %>%
 mutate(forest_type = fct_reorder(forest_type, mean_effect, .desc = TRUE))

# create the plot
p2 <- ggplot(fortyp_marginals, aes(x = x, y = forest_type, height = y, fill = response)) +
 geom_density_ridges(stat = "identity", scale = 1.5, alpha = 0.7) +
 geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
 labs(
  x = "Effect relative to aspen",
  y = "Predominant Forest Type",
  fill = "Response",
 ) +
 scale_fill_manual(values = c("FRP" = "#FEB24C", "CBI" = "#800026"),
                   labels = c(
                    "FRP" = "Cumulative Daytime FRP",
                    "CBI" = "90th Percentile CBIbc")) +
 theme_classic() +
 theme(axis.text.y = element_text(angle = 0, hjust = 1, size=9),
       axis.text.x = element_text(angle = 0, hjust = 0, size=9),
       axis.title.y = element_text(size = 10, margin = margin(r = 12)),
       axis.title.x = element_text(size = 10, margin = margin(t = 12)),
       legend.position = c(0.20, 0.14),
       legend.background = element_rect(
        fill = scales::alpha("white", 0.4), 
        color = NA, size = 0.8),
       legend.title = element_text(size = 9),
       legend.text = element_text(size = 8))
p2

# save the plot.
out_png <- paste0(maindir,'figures/INLA_FORTYPNM_PosteriorEffects_Ridge_species.png')
ggsave(out_png, plot = p2, dpi = 300, width = 7, height = 4, bg = 'white')



#################################
# Plot the climate + topo effects
p3 <- tidy_combined %>%
 filter(!parameter %in% c("(Intercept)", "aspen1", "overlap", "day_prop", "afd_count"),
        !str_starts(parameter, "fortyp")) %>%
 ggplot(., aes(x = x, y = parameter, height = y, fill = response)) +
 geom_density_ridges(stat = "identity", scale = 1.5, alpha = 0.7) +
 geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
 labs(
  x = "Effect size",
  y = "Predictor",
  fill = "Response"
 ) +
 scale_fill_manual(
  values = c("FRP" = "#FEB24C", "CBI" = "#800026"),
  labels = c("FRP" = "Cumulative Daytime FRP", "CBI" = "90th Percentile CBIbc")
 ) +
 theme_classic() +
 theme(
  axis.text.y = element_text(size = 12),
  axis.text.x = element_text(size = 12),
  axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), size = 14),
  axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 14),
  legend.position = "top",
  legend.title = element_text(size = 14),
  legend.text = element_text(size = 12)
 )
p3

# Save the plot
out_png <- paste0(maindir, 'figures/INLA_FORTYP_PosteriorEffects_Ridge_other_vars.png')
ggsave(out_png, plot = p3, dpi = 500, bg = 'white')





# #==============MODEL DIAGNOSTICS================#
# 
# # Extract fitted values (posterior mean estimates)
# fit_res <- data.frame(
#  Fitted = model_bl.frp.sp$summary.fitted.values$mean,
#  Residuals = da$frp - fitted
# )
# hist(fit_res$Residuals, breaks = 50, col = "skyblue", main = "Histogram of Residuals")
# 
# # Create a dataframe for plotting
# comparison_df <- data.frame(
#  type = c(rep("Observed", length(da$frp)), rep("Fitted", length(model_bl.frp.sp$summary.fitted.values$mean))),
#  value = c(da$frp, model_bl.frp.sp$summary.fitted.values$mean)
# )
# ggplot(comparison_df, aes(x = value, color = type, fill = type)) +
#  geom_density(alpha = 0.4) +
#  labs(title = "Density Plot of Observed vs. Fitted Values", 
#       x = "Value", 
#       y = "Density") +
#  scale_fill_manual(values = c("Observed" = "steelblue", "Fitted" = "orange")) +
#  scale_color_manual(values = c("Observed" = "steelblue", "Fitted" = "orange")) +
#  theme_minimal()
# 
# 
# # Plot residuals vs fitted values
# ggplot(fit_res, aes(x = Fitted, y = Residuals)) +
#  geom_point(alpha = 0.5) +
#  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
#  labs(title = "Residuals vs Fitted Values",
#       x = "Fitted Values",
#       y = "Residuals") +
#  theme_minimal()






