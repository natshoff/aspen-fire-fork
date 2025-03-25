
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
library(terra)
library(viridis)

maindir <- '/Users/max/Library/CloudStorage/OneDrive-Personal/mcook/aspen-fire/Aim2/'


#=========Load the Prepped gridcell data=========#

# gridcell_prep.R
fp <- paste0(maindir,"data/tabular/mod/model_data_cleaned.csv")
grid_tm <- read_csv(fp) %>%
 # prep some of the attributes for modeling
 mutate(
  first_obs_date = as.factor(first_obs_date),
  log_fire_size = log(fire_acres) # log-scale fire size
 ) %>%
 # flag re-burn gridcells
 group_by(grid_index) %>%
 mutate(
  # Get earliest fire date for this gridcell
  first_fire_dt = min(fire_ig_dt, na.rm = TRUE),
  # Flag reburn: if current fire date is after earliest
  reburn = as.factor(fire_ig_dt > first_fire_dt)
 ) %>%
 filter(
  # remove re-burn gridcells
  reburn == "FALSE"
 ) %>%
 distinct(grid_idx, fortypnm_gp, .keep_all = T) %>%
 ungroup()
glimpse(grid_tm)

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


##########################################
# Subset species with too few observations

# check how many grids of each forest type there are
(fortyp_counts <- da %>%
 group_by(fortypnm_gp) %>%
 summarize(n = n()) %>%
 ungroup())


#===============Explore Distributions, etc.================#

# ####################################
# # distribution of response variables
# resp_plot <- grid_tm %>%
#  # pivot longer to facet plot
#  pivot_longer(cols = c(log_frp_max,
#                        log_frp_max_day,
#                        log_frp_csum, 
#                        log_frp_csum_day),
#               names_to = "variable",
#               values_to = "value") %>%
#  # Plot with facets
#  ggplot(aes(x = value)) +
#  geom_histogram(bins = 30, fill = "orange", alpha = 0.7) +
#  facet_wrap(
#   ~ variable, 
#   scales = "free",
#   labeller = as_labeller(c(log_frp_max = "log(Max FRP)",
#                            log_frp_max_day = "log(Daytime Max FRP)",
#                            log_frp_csum = "log(Cumulative FRP)",
#                            log_frp_csum_day = "log(Daytime Cumulative FRP)"))) +
#  labs(x = "value",
#       y = "Frequency") +
#  theme_minimal()
# resp_plot
# 
# # save the plot.
# out_png <- paste0(maindir,'figures/INLA_ResponseDistribution_FRP.png')
# ggsave(out_png, plot = resp_plot, dpi=500, bg = 'white')
# rm(resp_plot)


# ####################
# # correlation matrix
# # Select only numeric columns and convert factors to dummy variables
# # this correlation matrix is for this simple model (i.e., not forest composition)
# cor_da <- grid_tm %>%
#  select(
#   fortypnm_gp, fire_acres, # forest type (factor)
#   forest_pct, fortyp_pct, # forest and forest type percent
#   canopypct_mean, balive_sum, # grid-level mean canopy percent and balive
#   lf_forest_cc_mean, lf_forest_ch_mean, # LANDFIRE canopy cover and height
#   erc, erc_dv, vpd, vpd_dv, # climate
#   fm1000, fm1000_dv, rmin, rmin_dv, # climate
#   tmmx, tmmx_dv, vs, vs_dv, #climate
#   elev, slope, tpi, chili,  # topography
#   fire_aspenpct, aspen, # fire-level aspen percent
#   overlap, day_prop, afd_count # VIIRS detections information
#  ) %>%
#  pivot_wider(
#   names_from = fortypnm_gp,
#   values_from = fortyp_pct,
#   values_fill = 0) %>%
#  mutate(across(everything(), ~ scale(.) %>% as.numeric()))  # Standardize variables
# 
# # Compute correlation matrix
# cor_mat <- cor(cor_da, use = "complete.obs", method = "spearman")
# 
# # Plot correlation matrix
# cor_plot <- ggcorrplot(
#  cor_mat,
#  method = "circle",  # Circle or square for visualization
#  type = "lower",  # Lower triangle of the correlation matrix
#  lab = TRUE,  # Show correlation values
#  lab_size = 3,
#  tl.cex = 10,  # Text label size
#  colors = c("blue", "white", "red")  # Color gradient
# )
# cor_plot
# 
# rm(cor_da, cor_mat)
# gc()
# 
# # save the plot.
# out_png <- paste0(maindir,'figures/INLA_CorrelationMatrix_Fortyp.png')
# ggsave(out_png, plot = cor_plot, dpi=500, bg = 'white')
# rm(cor_plot)



#===========MODEL SETUP==============#

# list of species names
spps <- c("quaking_aspen", "douglas_fir", "white_fir", "gambel_oak",
          "lodgepole_pine", "ponderosa_pine", "spruce_fir", "piñon_juniper")

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
qt[1,]$val

# create a data frame for just the predominant forest type
# one row per grid_index with the dominant forest type
da <- grid_tm %>%
 # select the columns we need for modeling
 select(grid_index, grid_idx, Fire_ID, fire_acres, # ID columns
        day_prop, overlap, afd_count, # proportion daytime observations and percent detection overlap
        fire_year, season, year_season, first_obs_date, fire_doy, # temporal effects
        log_frp_max, log_frp_max_day, # FRP response variables
        log_frp_csum, log_frp_csum_day, # FRP response variables
        CBIbc_p90, CBIbc_p95, CBIbc_p99, CBIbc_mean, # CBI response variables
        fortypnm_gp, fortyp_pct, forest_pct, # forest type and percent cover
        # fire_dfortyp, fire_dfortyp_int, # fire-level dominant forest type
        canopypct_mean, balive_sum, # canopy percent and total live basal area
        ba_dead_total, tpp_dead_total, # dead metrics
        lf_forest_cc_mean, lf_forest_ch_mean, # LANDFIRE canopy cover and height
        erc, erc_dv, vpd, vpd_dv, # climate
        fm1000, fm1000_dv, rmin, rmin_dv, # climate
        tmmx, tmmx_dv, vs, vs_dv, #climate
        elev, slope, tpi, chili, aspect, # topography
        x, y, # grid centroid coordinate for spatial fields model
        H_tpp, H_ba, # diversity
        aspen, fire_aspenpct # fire-level aspen percent/presence 
 ) %>%
 # center and scale continuous predictor variables
 mutate(
  across(
   c(
    forest_pct, fortyp_pct, day_prop, overlap, afd_count,
    canopypct_mean, balive_sum, # canopy percent and total live basal area
    lf_forest_cc_mean, lf_forest_ch_mean, # LANDFIRE canopy cover and height
    erc, vpd, vpd_dv, erc_dv, # climate
    fm1000, fm1000_dv, rmin, rmin_dv, # climate
    tmmx, tmmx_dv, vs, vs_dv, #climate
    elev, slope, tpi, chili, aspect, # topography
    H_tpp, H_ba, # diversity
    fire_aspenpct # fire-level aspen percent
   ),
   ~ as.numeric(scale(.))
  ),
  log_fire_size = log(fire_acres),
 ) %>%
 # rename the response variables
 rename(
  fire_id = Fire_ID,
  tm_canopypct = canopypct_mean,
  tm_balive = balive_sum,
  lf_canopypct = lf_forest_cc_mean,
  lf_canopyht = lf_forest_ch_mean
 ) %>%
 # keep just one row per grid cell by predominant type
 distinct(grid_idx, fortypnm_gp, .keep_all = TRUE) %>%
 arrange(grid_idx)

rm(grid_tm,qt)
gc()


# # identify species to drop
# small_spps <- spp_counts %>% 
#  # (fewer than 100 observations)
#  filter(n < 100) %>% 
#  pull(fortypnm_gp) # get the names
# print(small_spps)
# # count how many gridcells we are dropping
# print(paste0("Dropping [",dim(da%>%filter(fortypnm_gp %in% small_spps))[1],"] gridcells"))
# # filter these out
# da <- da %>%
#  filter(!fortypnm_gp %in% small_spps) %>%
#  # reset the factor levels
#  mutate(fortypnm_gp = droplevels(fortypnm_gp),
#         fire_dfortyp = droplevels(fire_dfortyp))
# levels(da$fortypnm_gp) # check levels
# 
# # save a spatial polygon 
# grid <- st_read(paste0(maindir,"data/spatial/mod/VIIRS/viirs_snpp_jpss1_afd_latlon_fires_pixar_gridstats.gpkg"))
# grid <- da %>%
#  left_join(grid %>% select(grid_index, geom), by="grid_index") %>%
#  st_as_sf(.)
# st_write(grid,paste0(maindir,"data/spatial/mod/grid_model_data_predom.gpkg"), append=F)
# 
# # Tidy up!
# rm(grid, grid_counts, idx, spp_counts)
# gc()




#===========MODEL FITTING==============#

set.seed(435)

da$aspen <- as.factor(da$aspen)
levels(da$aspen)

########################
# FIRE RADIATIVE POWER #

#######################################
# 1. Baseline model (no latent effects)

# Set up the model formula
mf.frp <- log_frp_csum ~ 1 + # cumulative FRP
 fortypnm_gp + # predominant (majority) forest type
 fortypnm_gp:fortyp_pct + # interaction with forest type cover
 lf_canopypct + # gridcell mean canopy percent
 H_ba + # gridcell diversity (basal area)
 ba_dead_total + # proportion live/dead basal area
 erc_dv + vpd + vs + # climate/weather
 elev + slope + aspect + tpi + # topography
 overlap + # gridcell VIIRS overlap (sum)
 day_prop + # gridcell proportion of daytime detections 
 aspen # fire-level aspen presence

# fit the model                     
ml.frp <- inla(
 mf.frp, data = da, 
 family = "gaussian",
 control.predictor = list(compute = TRUE),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
 control.family = list(
  # relax the variance assumptions
  hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.1)))
 ),
 # controls for posterior approximation and hyperparameter search
 control.inla = list(strategy = "laplace", int.strategy = "grid")
)
# check the model summary
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
     prec = list(prior = "pc.prec", param = c(1, 0.1))
   )) # temporal random effect
)
# fit the model                     
ml.frp.re <- inla(
 mf.frp.re, data = da, 
 family = "gaussian",
 control.predictor = list(compute = TRUE),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
 control.family = list(
  # relax the variance assumptions
  hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.5)))
 ),
 # controls for posterior approximation and hyperparameter search
 control.inla = list(strategy = "laplace", int.strategy = "grid")
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
     prec = list(prior = "pc.prec", param = c(1, 0.1))
  )) # fire-level random effect
)
# fit the model                     
ml.frp.re2 <- inla(
 mf.frp.re2, data = da, 
 family = "gaussian",
 control.predictor = list(compute = TRUE),
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

# ##########################################################
# # fit a semivariogram to look at spatial dependence in FRP
# ##########################################################
# 
# # reproject
# sp <-  grid_sf %>%
#  st_transform(st_crs(5070))
# 
# # Function to compute semivariogram for each fire
# get_fire_vario <- function(fire_id, data) {
#  fire_sp <- data %>% filter(fire_id == !!fire_id)
#  
#  if (nrow(fire_sp) > 10) {  # Only compute if there are enough points
#   vario <- variogram(log_frp_max ~ 1, data = fire_sp)
#   vario <- vario %>% mutate(fire_id = fire_id)  # Ensure fire_id is added correctly
#   print(paste("Computed variogram for fire:", fire_id, "with", nrow(fire_sp), "points"))
#   return(vario)
#  } else {
#   print(paste("Skipping fire:", fire_id, "- Not enough points"))
#   return(NULL)
#  }
# }
# 
# # Apply to each fire separately
# fires <- unique(sp$fire_id)
# variograms <- lapply(fires, get_fire_vario, data = sp)
# variograms <- do.call(rbind, variograms)  # Combine results
# 
# # Check if fire_id is still in the dataset
# print(unique(variograms$fire_id))
# 
# # Plot variograms grouped by Fire ID
# ggplot(variograms, aes(x = dist, y = gamma, 
#                 group = fire_id, color = fire_id)) +
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
#  group_by(fire_id) %>%
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


###################################################
###################################################

# convert coordinates to a matrix for INLA
coords_mat <- as.matrix(coords) 

# Define the spatial mesh
mesh <- inla.mesh.2d(
 loc = coords_mat, # Locations (grid centroids)
 max.edge = c(1, 20), # Maximum edge lengths (inner and outer)
 cutoff = 0.008, # Minimum distance between points (0.01deg = ~1.1km)
 offset = c(0.5, 0.1) # Boundary buffer
)
# Plot mesh to check
plot(mesh, main = "SPDE Mesh for FRP Model")
points(coords, col = "red", pch = 20)


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
 prior.range = c(0.2, 0.5),
 # P(sigma > 1) = 0.01
 prior.sigma = c(1, 0.1) # variance
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
 select(fire_id, first_obs_date, forest_pct, fortyp_pct, fortypnm_gp, lf_canopypct, 
        erc_dv, vpd, vs, slope, aspect, tpi, chili, elev, H_tpp, H_ba,
        aspen, day_prop, overlap, afd_count, ba_dead_total)
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
 control.family = list(
  # relax the variance assumptions
  hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.5)))
 ),
 # controls for posterior approximation and hyperparameter search
 control.inla = list(strategy = "laplace", int.strategy = "grid")
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


#=================MODEL STATEMENTS=================#

ml.frp.re.sp$summary.fixed

# lodgepole
ml.frp.re.sp$summary.fixed["fortypnm_gplodgepole_pine", ]
exp(ml.frp.re.sp$summary.fixed["fortypnm_gplodgepole_pine", "mean"]) - 1
exp(ml.frp.re.sp$summary.fixed["fortypnm_gplodgepole_pine", c("0.025quant", "0.975quant")]) - 1
# spruce-fir
ml.frp.re.sp$summary.fixed["fortypnm_gpspruce_fir", ]
exp(ml.frp.re.sp$summary.fixed["fortypnm_gpspruce_fir", "mean"]) - 1
exp(ml.frp.re.sp$summary.fixed["fortypnm_gpspruce_fir", c("0.025quant", "0.975quant")]) - 1
# douglas-fir
ml.frp.re.sp$summary.fixed["fortypnm_gpdouglas_fir", ]
exp(ml.frp.re.sp$summary.fixed["fortypnm_gpdouglas_fir", "mean"]) - 1
exp(ml.frp.re.sp$summary.fixed["fortypnm_gpdouglas_fir", c("0.025quant", "0.975quant")]) - 1


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
  sqrt_cbi_mn = sqrt(CBIbc_mean),
  sqrt_cbi_p90 = sqrt(CBIbc_p90),
  sqrt_cbi_p95 = sqrt(CBIbc_p95),
  sqrt_cbi_p99 = sqrt(CBIbc_p99)
 ) 
da %>%
 # pivot longer to facet plot
 pivot_longer(cols = c(CBIbc_mean, log_cbi_mn, 
                       sqrt_cbi_mn,
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
                           sqrt_cbi_mn = "sqrt(Average CBIbc)",
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
 fortypnm_gp:fortyp_pct + # by percent cover
 lf_canopypct + # grid-level mean canopy percent
 erc_dv + vpd + vs + # climate/weather
 elev + slope + aspect + tpi + # topography
 aspen # fire-level aspen presence

# fit the model                     
ml.cbi <- inla(
 mf.cbi, data = da, 
 family = "gamma", 
 control.predictor = list(compute = TRUE),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
 control.family = list(
  # relax the variance assumptions
  hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.1)))
 ),
 # controls for posterior approximation and hyperparameter search
 control.inla = list(strategy = "laplace", int.strategy = "grid")
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
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.1)))
    ) # temporal random effect
)
# fit the model                     
ml.cbi.re <- inla(
 mf.cbi.re, data = da, 
 family = "gamma",
 control.predictor = list(compute = TRUE),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
 control.family = list(
  # relax the variance assumptions
  hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.1)))
 ),
 # controls for posterior approximation and hyperparameter search
 control.inla = list(strategy = "laplace", int.strategy = "grid")
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
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.1)))
    ) # fire-level random effect
)
# fit the model                     
ml.cbi.re2 <- inla(
 mf.cbi.re2, data = da, 
 family = "gamma",
 control.predictor = list(compute = TRUE),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
 control.family = list(
  # relax the variance assumptions
  hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.1)))
 ),
 # controls for posterior approximation and hyperparameter search
 control.inla = list(strategy = "laplace", int.strategy = "grid")
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
# Compute and print the mean spatial range for FRP
mean_range_frp <- mean(range_est$range_m, na.rm = TRUE)
print(paste("Mean spatial range for FRP:", round(mean_range_frp, 2), "meters"))

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
 prior.sigma = c(1, 0.1) # variance
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
        fortypnm_gp, lf_canopypct, erc_dv, vpd, vs, H_tpp, H_ba,
        elev, slope, aspect, tpi, chili, aspen, ba_dead_total)
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
 control.family = list(
  # relax the variance assumptions
  hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.1)))
 ),
 # controls for posterior approximation and hyperparameter search
 control.inla = list(strategy = "laplace", int.strategy = "grid")
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
 filter(str_detect(parameter, "fortypnm_gp"),
        !str_detect(parameter, ":fortyp_pct")) %>%
 mutate(
  response = "FRP",
  forest_type = gsub("fortypnm_gp", "", parameter),
  effect = mean,
  lower = `0.025quant`,
  upper = `0.975quant`
 )

#######
# CBI #
cbi.eff <- as.data.frame(ml.cbi.re.sp$summary.fixed) %>%
 rownames_to_column(var = "parameter") %>%
 filter(str_detect(parameter, "fortypnm_gp"),
        !str_detect(parameter, ":fortyp_pct")) %>%
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


###############################################
# Ridge plot of forest type relative to aspen #
# extract fixed effects for FRP and CBI
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
 filter(str_detect(parameter, "fortypnm_gp"),
        !str_detect(parameter, ":fortyp_pct")) %>%
 mutate(forest_type = str_remove(parameter, "fortypnm_gp"),
        forest_type = recode(
         forest_type,
         "douglas_fir" = "Douglas-fir",
         "lodgepole_pine" = "Lodgepole pine",
         "ponderosa_pine" = "Ponderosa pine",
         "spruce_fir" = "Spruce-fir",
         "piñon_juniper" = "Piñon-juniper",
         "gambel_oak" = "Gambel oak",
         "white_fir" = "White fir"
        )) %>%
 group_by(forest_type) %>%
 mutate(mean_effect = mean(x, na.rm = TRUE)) %>%
 ungroup() %>%
 mutate(forest_type = fct_reorder(forest_type, -mean_effect))

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
                    "FRP" = "Cumulative FRP",
                    "CBI" = expression("90"^"th" ~ "Percentile CBI"))) +
 coord_cartesian(xlim=c(-0.48,0.44)) +
 theme_classic() +
 theme(axis.text.y = element_text(angle = 0, hjust = 1, size=9),
       axis.text.x = element_text(angle = 0, hjust = 0, size=9),
       axis.title.y = element_text(size = 10, margin = margin(r = 12)),
       axis.title.x = element_text(size = 10, margin = margin(t = 12)),
       legend.position = c(0.20, 0.25),
       legend.background = element_rect(
        fill = scales::alpha("white", 0.4), 
        color = NA, size = 0.8),
       legend.title = element_text(size = 9),
       legend.text = element_text(size = 8))
p2

# save the plot.
out_png <- paste0(maindir,'figures/INLA_FORTYPNM_PosteriorEffects_Ridge_species.png')
ggsave(out_png, plot = p2, dpi = 500, width = 7, height = 4, bg = 'white')

##########################################
# version without gambel oak and white fir 
p2.1 <- fortyp_marginals %>%
 filter(!forest_type %in% c("White fir","Gambel oak")) %>%
 ggplot(., aes(x = x, y = forest_type,
               height = y, fill = response)) +
 geom_density_ridges(stat = "identity", scale = 1.5, alpha = 0.7) +
 geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
 labs(
  x = "Effect relative to aspen",
  y = "Predominant Forest Cover",
  fill = "Response",
 ) +
 scale_fill_manual(values = c("FRP" = "#FEB24C", "CBI" = "#800026"),
                   labels = c(
                    "FRP" = "Cumulative FRP",
                    "CBI" = expression("90"^"th" ~ "Percentile CBI"))) +
 coord_cartesian(xlim=c(-0.30,0.40)) +
 theme_classic() +
 theme(axis.text.y = element_text(angle = 0, hjust = 1, size=9),
       axis.text.x = element_text(angle = 0, hjust = 0, size=9),
       axis.title.y = element_text(size = 10, margin = margin(r = 12)),
       axis.title.x = element_text(size = 10, margin = margin(t = 12)),
       legend.position = c(0.20, 0.25),
       legend.background = element_rect(
        fill = scales::alpha("white", 0.4), 
        color = NA, size = 0.8),
       legend.title = element_text(size = 9),
       legend.text = element_text(size = 8))
p2.1

# save the plot.
out_png <- paste0(maindir,'figures/INLA_FORTYPNM_PosteriorEffects_Ridge_species_v2.png')
ggsave(out_png, plot = p2.1, dpi = 500, width = 7, height = 4, bg = 'white')


#########################################
# interaction with forest percent cover #
# Filter for forest type effects
fortyp_pct_marginals <- tidy_combined %>%
 filter(str_detect(parameter, ":fortyp_pct")) %>%
 mutate(forest_type = str_remove(parameter, ":fortyp_pct"),
        forest_type = recode(
         forest_type,
         "douglas_fir" = "Douglas-fir",
         "lodgepole_pine" = "Lodgepole pine",
         "ponderosa_pine" = "Ponderosa pine",
         "spruce_fir" = "Spruce-fir",
         "piñon_juniper" = "Piñon-juniper",
         "gambel_oak" = "Gambel oak",
         "white_fir" = "White fir"
        )) %>%
 group_by(forest_type) %>%
 mutate(mean_effect = mean(x, na.rm = TRUE)) %>%
 ungroup() %>%
 mutate(forest_type = fct_reorder(forest_type, -mean_effect))

# create the plot
p3 <- ggplot(fortyp_pct_marginals, 
             aes(x = x, y = forest_type, 
                 height = y, fill = response)) +
 geom_density_ridges(stat = "identity", scale = 1.5, alpha = 0.7) +
 geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
 labs(
  x = "Effect",
  y = "Predominant Forest Cover",
  fill = "Response",
 ) +
 scale_fill_manual(values = c("FRP" = "#FEB24C", "CBI" = "#800026"),
                   labels = c(
                    "FRP" = "Cumulative FRP",
                    "CBI" = expression("90"^"th" ~ "Percentile CBI"))) +
 coord_cartesian(xlim=c(-0.48,0.44)) +
 theme_classic() +
 theme(axis.text.y = element_text(angle = 0, hjust = 1, size=9),
       axis.text.x = element_text(angle = 0, hjust = 0, size=9),
       axis.title.y = element_text(size = 10, margin = margin(r = 12)),
       axis.title.x = element_text(size = 10, margin = margin(t = 12)),
       legend.position = c(0.20, 0.25),
       legend.background = element_rect(
        fill = scales::alpha("white", 0.4), 
        color = NA, size = 0.8),
       legend.title = element_text(size = 9),
       legend.text = element_text(size = 8))
p3

# save the plot.
out_png <- paste0(maindir,'figures/INLA_FORTYPNM_PosteriorEffects_Ridge_species_pct.png')
ggsave(out_png, plot = p3, dpi = 500, width = 7, height = 4, bg = 'white')




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
# p3

# Save the plot
out_png <- paste0(maindir, 'figures/INLA_FORTYPNM_PosteriorEffects_Ridge_full.png')
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






