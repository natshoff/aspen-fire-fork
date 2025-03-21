
# libraries
library(tidyverse)
library(sf) 
library(INLA) 
library(ggcorrplot)
library(ggridges)
library(reshape2)
library(spdep)
library(patchwork)
library(forcats)
library(gstat)

maindir <- '/Users/max/Library/CloudStorage/OneDrive-Personal/mcook/aspen-fire/Aim2/'


#=========Load the Prepped gridcell data=========#

# gridcell_prep.R
fp <- paste0(maindir,"data/tabular/mod/model_data_cleaned.csv")
grid_tm <-  read_csv(fp) %>%
 # prep some of the attributes for modeling
 mutate(
  first_obs_date = as.factor(first_obs_date),
  log_fire_size = log(fire_acres) # log-scale fire size
 )
glimpse(grid_tm)
 

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
# retain predominantly forested gridcells ...
grid_tm <- grid_tm %>%
 filter(forest_pct >= 0.50)


###########################################
# (Optional) filter our re-burn gridcells #
grid_tm <- grid_tm %>%
 filter(reburn == "FALSE")


##################################################
# Assess the proportion of live basal area and TPP
# filter potential noise ...
# check proportion of live BA and TPP
(qt.ba <- tibble::enframe(
 quantile(grid_tm$ba_live_pr, probs = seq(0, 1, by = .1)),
 name = "qt", value = "val"
))
(qt10.ba = qt.ba[2,]$val) # 10%
# plot the distribution and 10th percentile line
ggplot(grid_tm, aes(x = ba_live_pr)) +
 geom_histogram(
  bins = 30, 
  fill = "grey", 
  color="grey40", 
  alpha = 0.4
 ) +
 geom_vline(
  xintercept = qt10.ba, color = "darkred", linewidth=1, linetype = "dashed"
 ) +
 labs(x = "Proportion live basal area", y = "Count") +
 theme_classic()
# tpp
(qt.tpp <- tibble::enframe(
 quantile(grid_tm$tpp_live_pr, probs = seq(0, 1, by = .1)),
 name = "qt", value = "val"
))
(qt10.tpp = qt.tpp[2,]$val) # 10%
# plot the distribution and 10th percentile line
ggplot(grid_tm, aes(x = tpp_live_pr)) +
 geom_histogram(
  bins = 30, 
  fill = "grey", 
  color="grey40", 
  alpha = 0.4
 ) +
 geom_vline(
  xintercept = qt10.tpp, color = "darkred", linewidth=1, linetype = "dashed"
 ) +
 labs(x = "Proportion tree/pixel", y = "Count") +
 theme_classic()
# check percentage of rows below the threshold
dim(grid_tm %>% filter(
 ba_live_pr >= qt10.ba &
  tpp_live_pr >= qt10.tpp
))[1] / dim(grid_tm)[1]

# filter species contributing less than the 10th percentile
grid_tm <- grid_tm %>%
 # filter noise in species data
 filter(
  # remove noise from  small species proportions
  # either tpp or ba over their 10th percentile
  (ba_live_pr >= qt10.ba | tpp_live_pr >= qt10.tpp)
 )
rm(qt.ba, qt.tpp)

###############################
# check on the species counts #
grid_tm %>% 
 group_by(species_gp_n) %>%
 summarise(n = length(species_gp_n)) %>%
 arrange(desc(n))

##########################################
# filter fires with not enough gridcells #
# should have at least 10 (?)
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
# # filter fires below the 10th percentile
# grid_tm <- grid_tm %>%
#  group_by(Fire_ID) %>%
#  filter(n() >= qt10) %>%
#  ungroup()
# # tidy up!
# rm(gridcell_counts,qt)

# check how many grids and fires
length(unique(grid_tm$grid_idx))
length(unique(grid_tm$Fire_ID))



# #==========EXPLORE THE DATA==========#
# 
# #######################################
# # species-specific correlation matrix #
# sp_cor <- grid_tm %>%
#  select(grid_idx, species_gp_n, ba_live) %>%
#  spread(species_gp_n, ba_live, fill = 0)  # Reshape to wide format
# # compute the correlation matrix
# sp_cormat <- cor(sp_cor[,-1], use = "complete.obs", method = "spearman")
# ggcorrplot(sp_cormat, method = "circle", type = "lower",
#            lab = TRUE, lab_size = 3, colors = c("blue", "white", "red"))
# rm(sp_cor, sp_cormat)
# # save the plot.
# out_png <- paste0(maindir,'figures/INLA_CorrelationMatrix_SpeciesBA.png')
# ggsave(out_png, dpi=500, bg = 'white')
# 
# ########################################
# # correlation matrix for fixed effects #
# cor_da <- grid_tm %>%
#  select(
#   fortypnm_gp,
#   # species structure metrics
#   tpp_live, ba_live, qmd_live, ht_live, dia_live,
#   tpp_live_pr, ba_live_pr, # proportions of TPP and BA
#   ba_live_total, tpp_live_total, qmd_live_mean, # gridcell structural summaries
#   ba_ld_pr, tpp_ld_pr, # live/dead proportions
#   ba_dead_total, tpp_dead_total, # dead BA and TPP
#   forest_pct, fortyp_pct, # gridcell forest percent and majority forest type percent
#   H_ba, H_tpp, # gridcell species diversity (abundance- and dominance-based)
#   erc, erc_dv, vpd, vpd_dv, # day-of-burn climate
#   fm1000, rmin, tmmx, vs, # day-of-burn climate
#   elev, slope, tpi, northness,  # topography
#   tm_canopy, tm_balive, lf_canopy, lf_height, # gridcell canopy percent/height and BA sum
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



#===========MODEL SETUP==============#

# prep the model data frame
# center and scale fixed effects
da <- grid_tm %>%
 mutate(
  # force aspen to be the baseline factor #
  species_gp_n = fct_relevel(species_gp_n, "quaking_aspen"),
  fortypnm_gp = fct_relevel(fortypnm_gp, "quaking_aspen"),
  # center/scale metrics / fixed effects
  across(
   c(ba_live, tpp_live, qmd_live, ht_live, dia_live, # species structure metrics
     ba_live_pr, tpp_live_pr, # proportion TPP and BA
     ba_live_total, tpp_live_total, qmd_live_mean, # gridcell structural summaries
     ba_dead_total, tpp_dead_total, # dead structure
     H_ba, H_tpp, # gridcell diversity metrics
     forest_pct, fortyp_pct, # gridcell forest percent and majority forest type percent
     lf_canopy, lf_height, # gridcell forest characteristics
     erc_dv, vpd, vs, # fire weather
     elev, slope, northness, tpi, # topography
     day_prop, overlap, # VIIRS detection characteristics
     dist_to_perim # gridcell distance to perimeter
    ), ~ as.numeric(scale(.))
  )
 ) %>%
 # calculate the predominant species by raw BA and TPP
 group_by(grid_idx) %>%
 mutate(
  dom_sp_ba = as.factor(species_gp_n[which.max(ba_live)]),
  dom_sp_tpp = as.factor(species_gp_n[which.max(tpp_live)])
 ) %>%
 ungroup() %>%
 arrange(grid_idx) # arrange by grid index
rm(grid_tm) # tidy up
gc()

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
# save this table out
write_csv(sp_dom_summary,paste0(maindir,"data/tabular/mod/results/species_dominance_counts.csv"))
rm(df1,df2,df3,sp_dom_summary)

###############################
# check on the species counts #
da %>% 
 group_by(species_gp_n) %>%
 summarise(n = length(species_gp_n)) %>%
 arrange(desc(n))

###############################
# check on forest type counts #
fortyp_counts <- da %>%
 group_by(fortypnm_gp) %>%
 summarize(n = n()) %>%
 ungroup()
fortyp_counts

################################
# check the mean FRP for aspen #
# extract all aspen rows
aspen.df = da%>%filter(fortypnm_gp == "quaking_aspen")
(aspen_frp_mean = mean(aspen.df$log_frp_csum)) # mean cumulative frp
(aspen_cbi_mean = mean(aspen.df$CBIbc_p90))
rm(aspen.df)
# also grab the mean frp in general
(mean_frp = mean(da$log_frp_csum))
(mean_cbi = mean(da$CBIbc_p90))

##################################
# # (optional) remove gambel oak #
# da <- da %>%
#  filter(species_gp_n != "gambel_oak")


#===========MODEL FITTING==============#
set.seed(456)

########################
# FIRE RADIATIVE POWER #

#################################################
# 1. Baseline model (no random or latent effects)

# setup the model formula
mf.frp <- log_frp_csum ~ 1 + 
 # forest composition and structure with mediation by VPD
 vpd * (fortypnm_gp) + # majority forest type
 vpd * (fortypnm_gp:fortyp_pct) + # majority forest type by percent cover
 vpd * (species_gp_n:ba_live) + # species total live basal area
 vpd * (species_gp_n:ht_live) + # species average live tree height
 vpd * (species_gp_n:dia_live) + # species average live tree diameter
 # other fixed effects (global effects)
 # fortyp_pct + # majority forest type percent
 lf_canopy + # gridcell forest canopy cover
 H_ba + # gridcell structural diversity (dominance-based)
 tpp_dead_total + # proportion live/dead basal area
 tpp_live_total + # total trees/pixel
 erc_dv + vpd + vs + # day-of-burn fire weather
 elev + slope + northness + tpi + # topography 
 overlap + # gridcell VIIRS overlap (cumulative)
 day_prop +  # gridcell proportion daytime detections 
 log_fire_size + # log-scaled fire size
 reburn + # gridcell re-burn status
 dist_to_perim # gridcell distance to perimeter

# fit the model
ml.frp <- inla(
 mf.frp, data = da,
 family = "gaussian",
 control.predictor = list(compute=T),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
 # set priors for the fixed effects
 control.family = list(
  # relax the variance assumptions
  hyper = list(
   prec = list(prior = "pc.prec", param = c(1, 0.5)))
 ),
 # controls for posterior approximation and hyperparameter search
 control.inla = list(strategy = "laplace", int.strategy = "grid") 
)
summary(ml.frp) # model summary table
# check on predictive power of the random effects models
mean(ml.frp$cpo$cpo, na.rm = TRUE)


#####################################
# 2. Baseline + Fire-ID Random Effect

# update the model formula
mf.frp.re <- update(
 mf.frp, . ~ 1 + . + 
  f(Fire_ID, model = "iid", hyper = list(
     prec = list(prior = "pc.prec", param = c(1, 0.5))
   )) 
)

# fit the model
ml.frp.re <- inla(
 mf.frp.re, data = da,
 family = "gaussian",
 control.predictor = list(compute=T),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
 # set priors for the fixed effects
 control.family = list(
  # relax the variance assumptions on gaussian precision
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
  f(first_obs_date, model = "iid", 
    hyper = list(prec = list(prior = "pc.prec", param = c(2, 0.5)))
   ) 
)

# fit the model
ml.frp.re2 <- inla(
 mf.frp.re2, data = da,
 family = "gaussian",
 control.predictor = list(compute=T),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
 # set priors for the fixed effects
 control.family = list(
  # relax the variance assumptions on gaussian precision
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

# extract gridcell coordinates
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
#  if (nrow(fire_sp) >= 10) {  # Only compute if there are enough points
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
#  drop_na() %>%
#  group_by(Fire_ID) %>%
#  summarize(
#   range_m = ifelse(any(gamma < max(gamma) * 0.9, na.rm = TRUE),
#                    max(dist[gamma < max(gamma) * 0.9], na.rm = TRUE),
#                    NA_real_)  # Assign NA when no valid values exist
#  )
# qt <- quantile(range_est$range_m, probs = seq(.1, .9, by = .1), na.rm=TRUE)
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
# # Plot mesh to check
# plot(mesh, main = "SPDE Mesh for FRP Model")
# points(coords, col = "red", pch = 20)

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
 prior.range = c(0.1, 0.5),
 # P(sigma > 1) = 0.01
 prior.sigma = c(1, 0.5) # variance
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
 select(
  Fire_ID, first_obs_date, grid_index, 
  fortypnm_gp, species_gp_n, 
  forest_pct, fortyp_pct, lf_canopy, H_tpp, H_ba, 
  ba_live, tpp_live, ht_live, dia_live, qmd_live, 
  tpp_live_pr, ba_live_pr, aspen_ba_pr,
  ba_live_total, ba_dead_total, tpp_live_total, tpp_dead_total,
  erc_dv, erc, vpd, vpd_dv, vs, 
  elev, slope, northness, tpi,
  fire_aspen, grid_aspen,
  day_prop, overlap, reburn,
  log_fire_size, dist_to_perim
 )

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
 data = inla.stack.data(stack.frp),
 family = "gaussian",
 control.predictor = list(A = inla.stack.A(stack.frp), compute=T),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
 # set priors for the fixed effects
 control.family = list(
  # relax the variance assumptions on gaussian precision
  hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.1)))
 ),
 # controls for posterior approximation and hyperparameter search
 control.inla = list(strategy = "laplace", int.strategy = "grid")
)
# model summary
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


#=================MODEL STATEMENTS=================#

# compute the exponentiated effects
exp.frp <- ml.frp.re.sp$summary.fixed %>%
 rownames_to_column(var = "parameter") %>%
 mutate(
  exp_mean = exp(mean) - 1,  # Convert log(FRP) effect to % difference
  lower_ci = exp(`0.025quant`) - 1,  # 2.5% CI bound
  upper_ci = exp(`0.975quant`) - 1   # 97.5% CI bound
 )
# save this table
write_csv(exp.frp, paste0(maindir,"data/tabular/mod/results/INLA_exp_FRP.csv"))
# check results
exp.frp%>%select(parameter,exp_mean,lower_ci,upper_ci)


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
  "(Intercept)","aspen","day_prop","log_fire_size"
  )
 ) %>%  
 mutate(
  effect = case_when(
   str_detect(parameter, "ba_live") ~ "Live basal area",
   # str_detect(parameter, "ba_live_pr") ~ "Proportion of\nlive basal area",
   str_detect(parameter, "tpp_live") ~ "Live trees/pixel",
   str_detect(parameter, "tpp_dead") ~ "Dead trees/pixel",
   str_detect(parameter, "dia_live") ~ "Tree diameter",
   str_detect(parameter, "ht_live") ~ "Tree height",
   str_detect(parameter, "fortyp_pct") ~ "Percent cover",
   # str_detect(parameter, "qmd_live") ~ "Live QMD",
   # str_detect(parameter, "H_tpp") ~ "Shannon diversity \n(abundance-based)",
   str_detect(parameter, "vs") ~ "Wind speed",
   str_detect(parameter, "elev") ~ "Elevation",
   str_detect(parameter, "northness") ~ "Northness",
   str_detect(parameter, "slope") ~ "Slope",
   str_detect(parameter, "H_ba") ~ "Shannon diversity (H-BA)",
   str_detect(parameter, "tpi") ~ "Topographic position",
   str_detect(parameter, "lf_canopy") ~ "Forest canopy cover",
   str_detect(parameter, "erc_dv") ~ "Energy release component\n(15-year deviation)",
   str_detect(parameter, "erc") ~ "Energy release component",
   str_detect(parameter, "vpd") ~ "Vapor pressure deficit",
   str_detect(parameter, "log_fire_size") ~ "log(Fire size)",
   str_detect(parameter, "dist_to_perim") ~ "Distance to fire edge",
   str_detect(parameter, "reburn") ~ "Re-burn",
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
              "Spruce-fir", "Quaking aspen"))) %>%
 mutate(exp_effect = exp(x))

# check on the species name extraction
unique(tidy.effects.frp$fill_species)
spps_breaks <- unique(tidy.effects.frp$fill_species)

##########################
# create the color mapping
color_map <- c(
 "Global Effect" = "gray",
 "Quaking aspen" = "#e6ab02",
 "Lodgepole pine" = "#d95f02",
 "Douglas-fir" = "#7570b3",
 "White fir" = "#a6cee3",
 "Gambel oak" = "#e7298a",
 "Piñon-juniper" = "#66a61e",
 "Ponderosa pine" = "#1b9e77", #1b9e77,
 "Spruce-fir" = "#1f78b4"
)

alpha_map <- c(
 "Quaking aspen" = 0.95,
 "Lodgepole pine" = 0.2,
 "Douglas-fir" = 0.2,
 "White fir" = 0.2,
 "Gambel oak" = 0.2,
 "Piñon-juniper" = 0.2,
 "Ponderosa pine" = 0.2,
 "Spruce-fir" = 0.2
)

###################################
# Extract the forest type effects #
fortyp_effects.frp <- tidy.effects.frp %>%
 filter(
  str_detect(parameter, "fortypnm_gp"),
  !str_detect(parameter, ":fortyp_pct")
 ) %>%
 mutate(
  group = case_when(
   str_detect(parameter, "vpd:") ~ "VPD-mediated",
   TRUE ~ "Forest type"  # Default for all other fixed effects
 )) %>%
 group_by(fill_species) %>%
 mutate(mean_effect = mean(x, na.rm = TRUE)) %>%
 ungroup() %>%
 mutate(fill_species = fct_reorder(fill_species, -mean_effect))
glimpse(fortyp_effects.frp)

# make the ridge plot
frp.p1 <- ggplot(fortyp_effects.frp,
       aes(x = x, y = fill_species, height = y,
           fill = group, alpha=group)) +
 geom_density_ridges(stat = "identity", scale = 1.5, aes(alpha = group)) +
 geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
 labs(
  x = "Effect relative to aspen",
  y = "Forest Type"
 ) +
 scale_fill_manual(
  values = c(
   "VPD-mediated" = "#e6ab02",
   "Forest type" = "#e6ab02"
  )
 ) +
 scale_alpha_manual(
  values = c(
   "VPD-mediated" = 0.2,
   "Forest type" = 0.98
  )
 ) +
 theme_minimal() +
 theme(
  legend.position = "none",
  axis.text.y = element_text(size = 10),
  axis.title.y = element_text(size = 11)
 )
frp.p1


#######################################################
# create the ridge plot for species structure metrics #
# no VPD-mediation
# set the order for parameters
param_order <- c("Tree diameter","Tree height",
                 "Live basal area","Percent cover")
# filter to structure effects
sp_effects.frp <- tidy.effects.frp %>%
 filter(effect %in% param_order,
        !str_detect(parameter, "vpd"))

# plot it
frp.p2 <- sp_effects.frp %>%
 mutate(effect = factor(effect, levels = param_order)) %>%
 ggplot(., aes(x = x, y = effect, 
               height = y, fill = fill_species, 
               alpha = fill_species)) +
 geom_density_ridges(
  stat = "identity", scale = 1.5, show.legend=F,
  color = scales::alpha("black", 0.6)) +
 geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
 labs(
  x = "",
  y = "Fixed Effect",
  fill = "Species"
 ) +
 # add a subplot label
 annotate("text", x = -0.20, y = 5,
          label = expression(bold("(A)") ~ "FRPc"),
          size = 4, hjust = 0) +
 scale_fill_manual(
  values = color_map,
  # Exclude "Global Effect" from the legend
  breaks = spps_breaks,  
  guide = guide_legend(title = "Species")
 ) +
 scale_alpha_manual(
  values = alpha_map,
  guide = "none"  # Hide alpha from legend
 ) +
 xlim(-0.20, 0.22) +
 theme_classic() +
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


################################################
# version 2: removing Gambel oak and white fir #
frp.p2.1 <- sp_effects.frp %>%
 filter(!fill_species %in% c("Gambel oak","White fir")) %>%
 mutate(effect = factor(effect, levels = param_order)) %>%
 ggplot(., aes(x = x, y = effect, 
               height = y, fill = fill_species, 
               alpha = fill_species)) +
 geom_density_ridges(
  stat = "identity", scale = 1.5, show.legend=F,
  color = scales::alpha("black", 0.6)
 ) +
 geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
 labs(
  x = "",
  y = "Fixed Effect",
  fill = "Species"
 ) +
 # add a subplot label
 annotate("text", x = -0.15, y = 5,
          label = expression(bold("(A)") ~ "FRPc"),
          size = 4, hjust = 0) +
 scale_fill_manual(
  values = color_map,
  # Exclude "Global Effect" from the legend
  breaks = spps_breaks,  
  guide = guide_legend(title = "Species")
 ) +
 scale_alpha_manual(
  values = alpha_map,
  guide = "none"  # Hide alpha from legend
 ) +
 xlim(-0.15, 0.15) +
 theme_classic() +
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
frp.p2.1


#####################################################
# create the ridge plot for all other fixed effects #
fixed_effects.frp <- tidy.effects.frp %>%
 filter(!effect %in% c("aspen1","overlap"),
        !str_detect(parameter,"fortypnm"),
        !str_detect(parameter,"species"))

# plot the ridges
frp.p3 <- ggplot(fixed_effects.frp, 
                 aes(x = x, y = effect, height = y)) +
 geom_density_ridges(stat = "identity", scale = 1.5, alpha = 0.7) +
 geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
 labs(
  x = "Effect Size",
  y = "Fixed Effect",
 ) +
 theme_minimal() +
 theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 8),
  axis.text.y = element_text(size = 10),
  axis.title.y = element_text(size = 12)
 )
frp.p3




########################
# COMPOSITE BURN INDEX #

#######################################
# 1. Baseline model (no latent effects)

# setup the model formula
mf.cbi <- CBIbc_p90 ~ 1 + 
 # forest composition and structure with mediation by VPD
 vpd * (fortypnm_gp) + # majority forest type
 vpd * (fortypnm_gp:fortyp_pct) + # majority forest type:percent cover
 vpd * (species_gp_n:ba_live) + # species total live basal area
 vpd * (species_gp_n:ht_live) + # species average live tree height
 vpd * (species_gp_n:dia_live) + # species average live tree diameter
 # other fixed effects (global effects)
 # fortyp_pct + # majority forest type poercent cover
 lf_canopy + # gridcell forest canopy cover
 H_ba + # gridcell structural diversity (dominance-based)
 tpp_dead_total + # proportion live/dead basal area
 tpp_live_total + # total trees/pixel
 erc_dv + vpd + vs + # day-of-burn fire weather
 elev + slope + northness + tpi + # topography 
 overlap + # gridcell VIIRS overlap (cumulative)
 day_prop +  # gridcell proportion daytime detections 
 log_fire_size + # log-scaled fire size
 # reburn + # gridcell re-burn status
 dist_to_perim # gridcell distance to fire edge
 
# fit the model
ml.cbi <- inla(
 mf.cbi, data = da,
 family = "gamma",
 control.predictor = list(compute=T),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
 # set priors for the fixed effects
 control.family = list(
  # relax the variance assumptions
  hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.5)))
 ),
 # controls for posterior approximation and hyperparameter search
 control.inla = list(strategy = "laplace", int.strategy = "grid")
)
summary(ml.cbi)
# check on predictive power of the random effects models
mean(ml.cbi$cpo$cpo, na.rm = TRUE)


##############################################
# 2. Baseline model + fire random effect

# update the model formula
mf.cbi.re <- update(
 mf.cbi, . ~ 1 + . + 
  f(Fire_ID, model = "iid", hyper = list(
   prec = list(prior = "pc.prec", param = c(5, 0.5))
  )) 
)

# fit the model
ml.cbi.re <- inla(
 mf.cbi.re, data = da,
 family = "gamma",
 control.predictor = list(compute = TRUE),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
 # set priors for the fixed effects
 control.family = list(
  # relax the variance assumptions
  hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.5)))
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
  f(first_obs_date, model = "iid", 
    hyper = list(prec = list(prior = "pc.prec", param = c(2, 0.5)))
  ) 
)

# fit the model                     
ml.cbi.re2 <- inla(
 mf.cbi.re2, data = da, 
 family = "gamma",
 control.predictor = list(compute = TRUE),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
 # set priors for the fixed effects
 control.family = list(
  # relax the variance assumptions
  hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.5)))
 ),
 # controls for posterior approximation and hyperparameter search
 control.inla = list(strategy = "laplace", int.strategy = "grid")
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
 prior.sigma = c(2, 0.5) # variance
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
        fortypnm_gp, species_gp_n, forest_pct, fortyp_pct,     
        lf_canopy, H_tpp, H_ba, ba_live_pr, tpp_live_pr, qmd_live, 
        ba_live, tpp_live, ht_live, dia_live, 
        ba_dead_total, tpp_dead_total, tpp_live_total,
        erc_dv, vpd, vpd_dv, vs, slope, tpi, northness, elev,
        fire_aspen, grid_aspen, day_prop, overlap, 
        log_fire_size, dist_to_perim, reburn)
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
 # set priors for the fixed effects
 control.family = list(
  # relax the variance assumptions
  hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.5)))
 ),
 # controls for posterior approximation and hyperparameter search
 control.inla = list(strategy = "laplace", int.strategy = "grid")
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


#=================MODEL STATEMENTS=================#

fixed_effects <- ml.cbi.re.sp$summary.fixed

# Compute percentage change relative to aspen
forest_diff.cbi <- fixed_effects %>%
 rownames_to_column(var = "parameter") %>%
 mutate(
  exp_mean = exp(mean) - 1,  # Convert log(FRP) effect to % difference
  lower_ci = exp(`0.025quant`) - 1,  # 2.5% CI bound
  upper_ci = exp(`0.975quant`) - 1   # 97.5% CI bound
 )
# save this table
write_csv(forest_diff.cbi, paste0(maindir,"data/tabular/mod/results/forest_differences_CBI.csv"))
# check results
forest_diff.cbi%>%select(parameter,exp_mean,lower_ci,upper_ci)


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
 # Exclude the intercept
 filter(!parameter %in% c(
  "(Intercept)","aspen","day_prop","overlap","log_fire_size"
 )
 ) %>%  
 mutate(
  effect = case_when(
   str_detect(parameter, "ba_live") ~ "Live basal area",
   # str_detect(parameter, "ba_live_pr") ~ "Proportion of\nlive basal area",
   str_detect(parameter, "tpp_live") ~ "Live trees/pixel",
   str_detect(parameter, "tpp_dead") ~ "Dead trees/pixel",
   str_detect(parameter, "dia_live") ~ "Tree diameter",
   str_detect(parameter, "ht_live") ~ "Tree height",
   str_detect(parameter, "fortyp_pct") ~ "Percent cover",
   # str_detect(parameter, "qmd_live") ~ "Live QMD",
   # str_detect(parameter, "H_tpp") ~ "Shannon diversity \n(abundance-based)",
   str_detect(parameter, "vs") ~ "Wind speed",
   str_detect(parameter, "elev") ~ "Elevation",
   str_detect(parameter, "northness") ~ "Northness",
   str_detect(parameter, "slope") ~ "Slope",
   str_detect(parameter, "H_ba") ~ "Shannon diversity (H-BA)",
   str_detect(parameter, "tpi") ~ "Topographic position",
   str_detect(parameter, "lf_canopy") ~ "Forest canopy cover",
   str_detect(parameter, "erc_dv") ~ "Energy release component\n(15-year deviation)",
   str_detect(parameter, "vpd") ~ "Vapor pressure deficit",
   str_detect(parameter, "log_fire_size") ~ "log(Fire size)",
   str_detect(parameter, "dist_to_perim") ~ "Distance to fire edge",
   str_detect(parameter, "reburn") ~ "Re-burn",
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
tidy.effects.cbi <- tidy.effects.cbi %>%
 mutate(
  effect_order = case_when(
   !is.na(species) ~ 2,  # Species-specific effects
   TRUE ~ 2              # Global effects
  ),
  effect = factor(effect, levels = tidy.effects.cbi %>%
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
unique(tidy.effects.cbi$fill_species)
spps_breaks <- unique(tidy.effects.cbi$fill_species)

###################################
# Extract the forest type effects #
fortyp_effects.cbi <- tidy.effects.cbi %>%
 filter(str_detect(parameter, "fortypnm_gp"),
        !str_detect(parameter, ":fortyp_pct")) %>%
 mutate(
  group = case_when(
   str_detect(parameter, "vpd:") ~ "VPD-mediated",
   TRUE ~ "Forest type"  # Default for all other fixed effects
  )) %>%
 group_by(fill_species) %>%
 mutate(mean_effect = mean(x, na.rm = TRUE)) %>%
 ungroup() %>%
 mutate(fill_species = fct_reorder(fill_species, -mean_effect))
glimpse(fortyp_effects.cbi)

# make the ridge plot
cbi.p1 <- ggplot(fortyp_effects.cbi,
                 aes(x = x, y = fill_species, height = y,
                     fill = group, alpha=group)) +
 geom_density_ridges(stat = "identity", scale = 1.5, aes(alpha = group)) +
 geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
 labs(
  x = "Effect relative to aspen",
  y = "Forest Type"
 ) +
 scale_fill_manual(
  values = c(
   "VPD-mediated" = "#800026",
   "Forest type" = "#800026"
  )
 ) +
 scale_alpha_manual(
  values = c(
   "VPD-mediated" = 0.2,
   "Forest type" = 0.98
  )
 ) +
 theme_minimal() +
 theme(
  legend.position = "none",
  axis.text.y = element_text(size = 10),
  axis.title.y = element_text(size = 11)
 )
cbi.p1


#######################################################
# create the ridge plot for species structure metrics #
sp_effects.cbi <- tidy.effects.cbi %>%
 filter(effect %in% c("Live basal area","Tree height",
                      "Percent cover","Tree diameter"),
        !str_detect(parameter, "vpd"))

# set the order for parameters
param_order <- c("Tree diameter","Tree height",
                 "Live basal area","Percent cover")
# plot it
cbi.p2 <- sp_effects.cbi %>%
 mutate(effect = factor(effect, levels = param_order)) %>%
 ggplot(., aes(x = x, y = effect, 
               height = y, fill = fill_species, 
               alpha = fill_species)) +
 geom_density_ridges(
  stat = "identity", scale = 1.5, show.legend=F,
  color = scales::alpha("black", 0.6)) +
 geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
 labs(
  x = "",
  y = "Fixed Effect",
  fill = "Species"
 ) +
 # add a subplot label
 annotate("text", x = -0.22, y = 5,
          label = expression(bold("(B)") ~ "CBIbc"),
          size = 4, hjust = 0) +
 scale_fill_manual(
  values = color_map,
  # Exclude "Global Effect" from the legend
  breaks = spps_breaks,  
  guide = guide_legend(title = "Species")
 ) +
 scale_alpha_manual(
  values = alpha_map,
  guide = "none"  # Hide alpha from legend
 ) +
 xlim(-0.22, 0.20) +
 theme_classic() +
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
cbi.p2

# Save the plot
out_png <- paste0(maindir, 'figures/INLA_ForestComposition_FixedEffects_CBI_species.png')
ggsave(out_png, dpi = 500, bg = 'white')


################################################
# version 2: removing Gambel oak and white fir #
cbi.p2.1 <- sp_effects.cbi %>%
 filter(!fill_species %in% c("Gambel oak","White fir")) %>%
 mutate(effect = factor(effect, levels = param_order)) %>%
 ggplot(., aes(x = x, y = effect, 
               height = y, fill = fill_species, 
               alpha = fill_species)) +
 geom_density_ridges(
  stat = "identity", scale = 1.5, show.legend=T,
  color = scales::alpha("black", 0.6)
 ) +
 geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
 labs(
  x = "Effect size",
  y = "Fixed Effect",
  fill = "Species",
 ) +
 # add a subplot label
 annotate("text", x = -0.14, y = 5,
          label = expression(bold("(B)") ~ "CBIbc"),
          size = 4, hjust = 0) +
 scale_fill_manual(
  values = color_map,
  guide = "none"
 ) +
 scale_alpha_manual(
  values = alpha_map,
  guide = "none"
 ) +
 xlim(-0.14, 0.11) +
 theme_classic() +
 theme(
  axis.text.y = element_text(angle = 0, hjust = 1, size=9),
  axis.text.x = element_text(angle = 0, hjust = 0, size=9),
  axis.title.y = element_text(size = 10, margin = margin(r = 8)),
  axis.title.x = element_text(size = 10, margin = margin(t = 8)),
  legend.position = c(0.16, 0.30),
  legend.background = element_rect(
   fill = scales::alpha("white", 0.6), 
   color = NA, size = 1),
  legend.title = element_text(size = 9),
  legend.text = element_text(size = 8)
 )
cbi.p2.1


#####################################################
# create the ridge plot for all other fixed effects #
fixed_effects.cbi <- tidy.effects.cbi %>%
 filter(!effect %in% c("aspen1","overlap"),
        !str_detect(parameter,"fortypnm"),
        !str_detect(parameter,"species"),
        !fill_species %in% c("Gambel oak","White fir"))

# plot the ridges
cbi.p3 <- ggplot(fixed_effects.cbi, 
                 aes(x = x, y = effect, height = y)) +
 geom_density_ridges(stat = "identity", scale = 1.5, alpha = 0.7) +
 geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
 labs(
  x = "Effect Size",
  y = "Fixed Effect",
 ) +
 theme_minimal() +
 theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 8),
  axis.text.y = element_text(size = 10),
  axis.title.y = element_text(size = 12)
 )
cbi.p3

# Save the plot
out_png <- paste0(maindir, 'figures/INLA_ForestComposition_FixedEffects_other.png')
ggsave(out_png, dpi = 500, bg = 'white')


#######################################
# Make combined figures (FRP and CBI) #

# Create the combined forest type relative to aspen
fortyp_effects.frp$response <- "FRP"
fortyp_effects.cbi$response <- "CBI"
fortyp_effects <- bind_rows(fortyp_effects.frp, fortyp_effects.cbi) %>%
 # reorder the effects so main effect is drawn on top
 mutate(
  group = factor(group, levels = c("VPD-mediated", "Forest type")),
  fill_cat = factor(
   paste(response, group),
   levels = c(
    "CBI VPD-mediated", "FRP VPD-mediated", 
    "CBI Forest type", "FRP Forest type"
   ))
  )
# order the forest types from low->high on the main effect
forest_order <- fortyp_effects %>%
 filter(group == "Forest type") %>%
 group_by(fill_species) %>%
 summarize(mean_effect = mean(x, na.rm = TRUE)) %>%
 arrange(mean_effect)
# apply the ordering
fortyp_effects$fill_species <- factor(
 fortyp_effects$fill_species, 
 levels = rev(forest_order$fill_species)
)

# Plot posterior effects with credible intervals
cmap <- c("CBI" = "#800026", "FRP" = "#FEB24C") # color map for the response
p4 <- ggplot(fortyp_effects, aes(x = x, y = fill_species, height = y, 
                                 fill = fill_cat, 
                                 alpha = fill_cat,
                                 color = fill_cat)) +
 geom_density_ridges(stat = "identity", scale = 1.5) +
 geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
 labs(
  x = "Effect relative to aspen",
  y = "Predominant Forest Type",
  fill = "Effect Type",
 ) +
 scale_fill_manual(
  values = c(
   "FRP Forest type" = "#FEB24C",  
   "FRP VPD-mediated" = scales::alpha("#FEB24C", 0.3),
   "CBI Forest type" = "#800026",  
   "CBI VPD-mediated" = scales::alpha("#800026", 0.3)
  ),
  labels = c(
   "Cumulative FRP", "VPD-mediated (FRPc)",
   "90th percentile CBIbc", "VPD-mediated (CBIbc)"
  ),
  guide = guide_legend(
   override.aes = list(
    fill = c("#FEB24C", scales::alpha("#FEB24C", 0.3), 
             "#800026", scales::alpha("#800026", 0.3)),
    alpha = c(1, 0.3, 1, 0.3),  # Match transparency
    color = c(
     scales::alpha("black", 0.7),
     scales::alpha("#FEB24C", 0.3), 
     scales::alpha("black", 0.7), 
     scales::alpha("#800026", 0.3))
   ),
   order = 1  # Ensure this legend appears first
  )
 ) +
 scale_alpha_manual(
  values = c(
   "FRP Forest type" = 1,
   "FRP VPD-mediated" = 0.3,
   "CBI Forest type" = 1,
   "CBI VPD-mediated" = 0.3
  ),
  guide = "none"  
 ) +
 scale_color_manual(
  values = c(
   "FRP Forest type" = scales::alpha("black", 0.6),
   "FRP VPD-mediated" = scales::alpha("#FEB24C", 0.3),
   "CBI Forest type" = scales::alpha("black", 0.6),  
   "CBI VPD-mediated" = scales::alpha("#800026", 0.3)
  ),
  guide = "none"
 ) +
 coord_cartesian(xlim=c(-0.38,0.44)) +
 theme_classic() +
 theme(axis.text.y = element_text(angle = 0, hjust = 1, size=9),
       axis.text.x = element_text(angle = 0, hjust = 0, size=9),
       axis.title.y = element_text(size = 10, margin = margin(r = 12)),
       axis.title.x = element_text(size = 10, margin = margin(t = 12)),
       legend.position = c(0.18, 0.26),
       legend.background = element_rect(
        fill = scales::alpha("white", 0.4), 
        color = NA, size = 0.8),
       legend.title = element_text(size = 9),
       legend.text = element_text(size = 8))
p4
# Save the plot
out_png <- paste0(maindir, 'figures/INLA_ForestComposition_FORTYPCD_FRP-CBI_species_all.png')
ggsave(out_png, dpi = 500, width = 7, height = 4, bg = 'white')


#############################################
# version 2: without Gambel oak and white fir
p4.1 <- fortyp_effects %>%
 filter(!fill_species %in% c("Gambel oak")) %>%
 ggplot(., aes(x = x, y = fill_species, height = y, 
               fill = fill_cat, 
               alpha = fill_cat,
               color = fill_cat)) +
 geom_density_ridges(stat = "identity", scale = 1.5) +
 geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
 labs(
  x = "Effect relative to aspen",
  y = "Predominant Forest Type",
  fill = "Effect Type",
 ) +
 scale_fill_manual(
  values = c(
   "FRP Forest type" = "#FEB24C",  
   "FRP VPD-mediated" = scales::alpha("#FEB24C", 0.3),
   "CBI Forest type" = "#800026",  
   "CBI VPD-mediated" = scales::alpha("#800026", 0.3)
  ),
  labels = c(
   "Cumulative FRP", "VPD-mediated (FRPc)",
   "90th percentile CBIbc", "VPD-mediated (CBIbc)"
  ),
  guide = guide_legend(
   override.aes = list(
    fill = c("#FEB24C", scales::alpha("#FEB24C", 0.3), 
             "#800026", scales::alpha("#800026", 0.3)),
    alpha = c(1, 0.3, 1, 0.3),  # Match transparency
    color = c(
     scales::alpha("black", 0.7),
     scales::alpha("#FEB24C", 0.3), 
     scales::alpha("black", 0.7), 
     scales::alpha("#800026", 0.3))
   ),
   order = 1  # Ensure this legend appears first
  )
 ) +
 scale_alpha_manual(
  values = c(
   "FRP Forest type" = 1,
   "FRP VPD-mediated" = 0.3,
   "CBI Forest type" = 1,
   "CBI VPD-mediated" = 0.3
  ),
  guide = "none"  
 ) +
 scale_color_manual(
  values = c(
   "FRP Forest type" = scales::alpha("black", 0.6),
   "FRP VPD-mediated" = scales::alpha("#FEB24C", 0.3),
   "CBI Forest type" = scales::alpha("black", 0.6),  
   "CBI VPD-mediated" = scales::alpha("#800026", 0.3)
  ),
  guide = "none"
 ) +
 coord_cartesian(xlim=c(-0.38,0.44)) +
 theme_classic() +
 theme(axis.text.y = element_text(angle = 0, hjust = 1, size=9),
       axis.text.x = element_text(angle = 0, hjust = 0, size=9),
       axis.title.y = element_text(size = 10, margin = margin(r = 12)),
       axis.title.x = element_text(size = 10, margin = margin(t = 12)),
       legend.position = c(0.18, 0.26),
       legend.background = element_rect(
        fill = scales::alpha("white", 0.4), 
        color = NA, size = 0.8),
       legend.title = element_text(size = 9),
       legend.text = element_text(size = 8))
p4.1
# Save the plot
out_png <- paste0(maindir, 'figures/INLA_ForestComposition_FORTYPCD_FRP-CBI_species.png')
ggsave(out_png, dpi = 500, width = 8, height = 4, bg = 'white')


###################################
# Species composition / structure #
p5 <- frp.p2 / cbi.p2 # "/" stacks them, "|" would place them side-by-side
p5
# Save the plot
out_png <- paste0(maindir, 'figures/INLA_ForestComposition_FixedEffects_FRP-CBI_species.png')
ggsave(out_png, plot = p5, dpi = 500, width = 7, height = 6, bg = 'white')
# version 2
p5.1 <- frp.p2.1 / cbi.p2.1 # "/" stacks them, "|" would place them side-by-side
p5.1
# Save the plot
out_png <- paste0(maindir, 'figures/INLA_ForestComposition_FixedEffects_FRP-CBI_species_v2.png')
ggsave(out_png, plot = p5.1, dpi = 500, width = 8, height = 6, bg = 'white')


#######################
# All other effects ...
# Create the combined forest type relative to aspen
fixed_effects.frp$response <- "FRP"
fixed_effects.cbi$response <- "CBI"
fixed_effects <- bind_rows(fixed_effects.frp, fixed_effects.cbi)

# plot it
p6 <- fixed_effects %>%
 filter(effect != "Re-burn") %>%
 ggplot(., aes(x = x, y = effect, height = y, fill = response)) +
 geom_density_ridges(stat = "identity", scale = 1.5, alpha = 0.7) +
 geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
 labs(
  x = "Effect size",
  y = "Parameter",
  fill = "Response",
 ) +
 scale_fill_manual(values = c("FRP" = "#FEB24C", "CBI" = "#800026"),
                   labels = c(
                    "FRP" = "Cumulative FRP",
                    "CBI" = expression("90"^"th" ~ "Percentile CBI"))) +
 coord_cartesian(xlim=c(-0.16,0.20)) +
 theme_classic() +
 theme(axis.text.y = element_text(angle = 0, hjust = 1, size=9),
       axis.text.x = element_text(angle = 0, hjust = 0, size=9),
       axis.title.y = element_text(size = 10, margin = margin(r = 12)),
       axis.title.x = element_text(size = 10, margin = margin(t = 12)),
       legend.position = c(0.18, 0.20),
       legend.background = element_rect(
        fill = scales::alpha("white", 0.4), 
        color = NA, size = 0.8),
       legend.title = element_text(size = 9),
       legend.text = element_text(size = 8))
p6
# Save the plot
out_png <- paste0(maindir, 'figures/INLA_ForestComposition_FixedEffects_FRP-CBI_fixedeffects.png')
ggsave(out_png, dpi = 500, width = 8, height = 6, bg = 'white')

gc()


#====================RESULTS TABLES=======================#

# set the results directory
results_dir = paste0(maindir,"data/tabular/mod/results/")

ml_comparison.frp$response <- "FRP"
ml_comparison.cbi$response <- "CBI"
ml_comparison <- bind_rows(ml_comparison.frp, ml_comparison.cbi)
head(ml_comparison)

# save to CSV
write_csv(ml_comparison, paste0(results_dir,"model_comparison_FRP-CBI.csv"))


######################################################
# Extract fixed effects summary for both FRP and CBI #
summary.frp <- as.data.frame(ml.frp.re.sp$summary.fixed) %>%
 rownames_to_column("Parameter") %>%
 rename(
  Mean = mean,
  SD = sd,
  `2.5% CI` = `0.025quant`,
  `50% CI` = `0.5quant`,
  `97.5% CI` = `0.975quant`
 ) %>%
 mutate(Response = "FRP")
summary.cbi <- as.data.frame(ml.cbi.re.sp$summary.fixed) %>%
 rownames_to_column("Parameter") %>%
 rename(
  Mean = mean,
  SD = sd,
  `2.5% CI` = `0.025quant`,
  `50% CI` = `0.5quant`,
  `97.5% CI` = `0.975quant`
 ) %>%
 mutate(Response = "CBI")
# combine them
ml_summary.fixed <- bind_rows(summary.frp, summary.cbi)
head(ml_summary.fixed)

# save to CSV
write_csv(ml_summary.fixed, paste0(results_dir,"model_summaries_FRP-CBI.csv"))

gc()