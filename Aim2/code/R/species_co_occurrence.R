# Load the required libraries
library(tidyverse)
library(sf) # spatial
library(INLA) # for spatial Bayes model

# Environment variables
maindir <- '/Users/max/Library/CloudStorage/OneDrive-Personal/mcook/aspen-fire/Aim2/'

proj <- ""

#=========Prep the grid data=========#

# Format the species composition data frame

# load the aggregated FRP grid with TreeMap
fp <- paste0(maindir,'data/tabular/mod/viirs_snpp_jpss1_gridstats_fortypcd.csv')
grid_fortyp <-  read_csv(fp)
glimpse(grid_fortyp)


# tidy the columns
grid_ <- grid_fortyp %>%
 # select the required columns
 select(c(grid_index, Fire_ID, frp_csum, frp_max, 
          grid_x, grid_y, SpeciesName, spp_pct, forest_pct)) %>%
 filter(frp_max > 0) %>% # make sure FRP is not 0
 mutate(Fire_ID = as.factor(Fire_ID))
head(grid_) # check the results


# reshape the data frame for modeling
grid_w <- grid_ %>% 
 # tidy the species names
 mutate(SpeciesName = str_replace_all(SpeciesName, "-", "_"),
        SpeciesName = str_to_lower(SpeciesName)) %>%
 pivot_wider(
  names_from = SpeciesName, 
  values_from = spp_pct, 
  values_fill = 0) %>% # pivot wider
 filter(frp_max > 0) %>%
 mutate(log_frp_max = log(frp_max + 1))
head(grid_w)
nrow(grid_w)


# retain grid cells with some aspen component
grid_aspen <- grid_w %>%
 filter(aspen > 0)
nrow(grid_aspen)/nrow(grid_w)*100
rm(grid_)
gc()


#===========Explore Distributions====#

# distribution of raw frp_max and log-transformed

# distribution of forest type percent cover
grid_w %>%
 select(aspen, douglas_fir, spruce_fir, lodgepole, ponderosa, piñon_juniper) %>%
 pivot_longer(everything(), names_to = "species", values_to = "percent_cover") %>%
 ggplot(aes(x = percent_cover)) +
 geom_histogram() +
 facet_wrap(~ species, scales = "free") +
 theme_minimal()


#===========MODEL SETUP==============#

set.seed(456)

spp <- c("aspen", "douglas_fir", "lodgepole", "ponderosa", "spruce_fir", "piñon_juniper")

# 1. Baseline model without spatial component
mf1 <- log_frp_max ~ aspen + douglas_fir + spruce_fir + lodgepole + ponderosa + piñon_juniper # species composition
# fit the model                     
model_bl1 <- inla(
 mf1, data = grid_w,
 family = "gaussian",
 control.predictor = list(compute = TRUE),
 control.compute = list(dic = TRUE, waic = TRUE)
)
summary(model_bl1)

# 2. Baseline model without spatial component, adding a Fire_ID random effect
mf1 <- log_frp_max ~ aspen + douglas_fir + spruce_fir + lodgepole + ponderosa + piñon_juniper +  # species composition
                     f(Fire_ID, model = "iid")  # Random effect for fire-level variability
 # fit the model                   
model_bl1 <- inla(
 mf1, data = grid_w,
 family = "gaussian",
 control.predictor = list(compute = TRUE),
 control.compute = list(dic = TRUE, waic = TRUE)
)
summary(model_bl1)








# fires <- unique(grid_aspen$Fire_ID)[1:3]  # Get unique fire IDs
# # fire_data <- split(grid_aspen, grid_aspen$Fire_ID)  # Split data by Fire_ID
# fire_data <- grid_aspen %>% filter(Fire_ID %in% fires)
# 
# fire_meshes <- list()
# fire_spdes <- list()
# 
# for (fire_id in unique(fire_data$Fire_ID)) {
#  print(fire_id)
#  fire_grid <- grid_aspen %>% filter(Fire_ID == fire_id)
#  
#  # Create the mesh for this fire
#  loc <- as.matrix(fire_grid[, c("grid_x", "grid_y")])
#  mesh <- inla.mesh.2d(
#   loc = loc,
#   max.edge = c(30, 100),  # Adjust based on fire size
#   cutoff = 10           # Minimum distance between nodes
#  )
#  
#  # Create the SPDE model for the mesh
#  spde <- inla.spde2.matern(mesh)
#  
#  # Store the mesh and SPDE model
#  fire_meshes[[fire_id]] <- mesh
#  fire_spdes[[fire_id]] <- spde
# }
# 
# ##################
# # attached to fire
# fire_indices <- list()
# 
# for (fire_id in names(fire_spdes)) {
#  spde <- fire_spdes[[fire_id]]
#  n_spde <- spde$n.spde
#  fire_indices[[fire_id]] <- data.frame(
#   spatial_field = 1:n_spde,  # Spatial index
#   Fire_ID = fire_id          # Fire ID
#  )
# }
# 
# # Combine all spatial indices
# spatial_indices <- do.call(rbind, fire_indices)
# 
# ########################
# # Create the data stacks
# stack_list <- list()
# 
# for (fire_id in names(fire_meshes)) {
#  fire_grid <- grid_aspen %>% filter(Fire_ID == fire_id)
#  mesh <- fire_meshes[[fire_id]]
#  spde <- fire_spdes[[fire_id]]
#  
#  # Create A matrix
#  A <- inla.spde.make.A(
#   mesh = mesh,
#   loc = as.matrix(fire_grid[, c("grid_x", "grid_y")])
#  )
#  
#  # isolate effects
#  effects_df <- fire_grid %>%
#   select(-c(grid_index, frp_max, grid_x, grid_y, forest_pct))
#  
#  # Create the stack
#  stack <- inla.stack(
#   data = list(frp_max = fire_grid$frp_max),
#   A = list(A, 1),
#   effects = list(
#    spatial_field = spatial_indices[spatial_indices$Fire_ID == fire_id, "spatial_field"],
#    data.frame(effects_df)  # Covariates
#   )
#  )
#  
#  stack_list[[fire_id]] <- stack
# }
# 
# # Combine all stacks
# full_stack <- do.call(inla.stack, stack_list)
# 
# data_stack <- inla.stack.data(full_stack)
# data_stack$Fire_ID <- as.factor(data_stack$Fire_ID)
# str(data_stack)
# 
# any(is.na(data_stack))
# sapply(data_stack, function(x) any(is.nan(x)))
# 
# 
# # Run the model.
# model_formula <- frp_max ~ aspen + lodgepole + 
#  f(Fire_ID, model = "iid") +  # Random effect for fire-level variability
#  f(spatial_field, model = spde, group = Fire_ID)  # Fire-specific spatial fields
# 
# result <- inla(
#  formula = model_formula,
#  data = inla.stack.data(full_stack),
#  control.predictor = list(A = inla.stack.A(full_stack)),
#  control.compute = list(dic = TRUE, waic = TRUE),  # Enable model evaluation metrics
#  family = "gaussian"  # Use Gaussian for log-transformed FRP
# )
