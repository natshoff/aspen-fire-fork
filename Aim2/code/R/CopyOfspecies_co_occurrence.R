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
grid_fortyp <-  read_csv(paste0(maindir,'data/tabular/mod/viirs_snpp_jpss1_gridstats_fortypcd.csv'))
glimpse(grid_fortyp)

grid_ <- grid_fortyp %>%
 # select the required columns
 select(c(grid_index, Fire_ID, frp_csum, frp_max, 
          grid_x, grid_y, SpeciesName, spp_pct, forest_pct)) %>%
 filter(frp_max > 0) %>%
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
  values_fill = 0) # pivot wider
head(grid_w)
nrow(grid_w)

# retain grid cells with some aspen component
grid_aspen <- grid_w %>%
 filter(aspen > 0,
        frp_max > 0) # remove 0 FRP

nrow(grid_aspen)/nrow(grid_w)*100
rm(grid_, grid_w)
gc()



#===========MODEL SETUP==============#

spp <- c("aspen", "douglas_fir", "lodgepole", "ponderosa", "spruce_fir", "piñon_juniper")

fires <- unique(grid_aspen$Fire_ID)[1:3]  # Get unique fire IDs
# fire_data <- split(grid_aspen, grid_aspen$Fire_ID)  # Split data by Fire_ID
fire_data <- grid_aspen %>% filter(Fire_ID %in% fires)

fire_meshes <- list()
fire_spdes <- list()

for (fire_id in unique(fire_data$Fire_ID)) {
 print(fire_id)
 fire_grid <- grid_aspen %>% filter(Fire_ID == fire_id)
 
 # Create the mesh for this fire
 loc <- as.matrix(fire_grid[, c("grid_x", "grid_y")])
 mesh <- inla.mesh.2d(
  loc = loc,
  max.edge = c(30, 100),  # Adjust based on fire size
  cutoff = 10           # Minimum distance between nodes
 )
 
 # Create the SPDE model for the mesh
 spde <- inla.spde2.matern(mesh)
 
 # Store the mesh and SPDE model
 fire_meshes[[fire_id]] <- mesh
 fire_spdes[[fire_id]] <- spde
}

##################
# attached to fire
fire_indices <- list()

for (fire_id in names(fire_spdes)) {
 spde <- fire_spdes[[fire_id]]
 n_spde <- spde$n.spde
 fire_indices[[fire_id]] <- data.frame(
  spatial_field = 1:n_spde,  # Spatial index
  Fire_ID = fire_id          # Fire ID
 )
}

# Combine all spatial indices
spatial_indices <- do.call(rbind, fire_indices)

########################
# Create the data stacks
stack_list <- list()

for (fire_id in names(fire_meshes)) {
 fire_grid <- grid_aspen %>% filter(Fire_ID == fire_id)
 mesh <- fire_meshes[[fire_id]]
 spde <- fire_spdes[[fire_id]]
 
 # Create A matrix
 A <- inla.spde.make.A(
  mesh = mesh,
  loc = as.matrix(fire_grid[, c("grid_x", "grid_y")])
 )
 
 # isolate effects
 effects_df <- fire_grid %>%
  select(-c(grid_index, frp_max, grid_x, grid_y, forest_pct))
 
 # Create the stack
 stack <- inla.stack(
  data = list(frp_max = fire_grid$frp_max),
  A = list(A, 1),
  effects = list(
   spatial_field = spatial_indices[spatial_indices$Fire_ID == fire_id, "spatial_field"],
   data.frame(effects_df)  # Covariates
  )
 )
 
 stack_list[[fire_id]] <- stack
}

# Combine all stacks
full_stack <- do.call(inla.stack, stack_list)

data_stack <- inla.stack.data(full_stack)
data_stack$Fire_ID <- as.factor(data_stack$Fire_ID)
str(data_stack)

any(is.na(data_stack))
sapply(data_stack, function(x) any(is.nan(x)))


# Run the model.
model_formula <- frp_max ~ aspen + lodgepole + 
 f(Fire_ID, model = "iid") +  # Random effect for fire-level variability
 f(spatial_field, model = spde, group = Fire_ID)  # Fire-specific spatial fields

result <- inla(
 formula = model_formula,
 data = inla.stack.data(full_stack),
 control.predictor = list(A = inla.stack.A(full_stack)),
 control.compute = list(dic = TRUE, waic = TRUE),  # Enable model evaluation metrics
 family = "gaussian"  # Use Gaussian for log-transformed FRP
)





##########################
# Set up the model (hgam)
library(mgcv)

# Testing for one species interaction

df <- grid_aspen %>%
 filter(spruce_fir > 0) %>%
 mutate(aspen_sp = aspen + spruce_fir) %>%
 filter(aspen_sp >= 50)

model <- gam(
 frp_max ~ s(aspen) + s(spruce_fir) + te(aspen, spruce_fir) +
  s(grid_x, grid_y, bs = "gp") + s(Fire_ID, bs = "re"),
 data = df, 
 family = gaussian()
)

summary(model)

plot(model, pages = 1)

# Visualize interaction as a contour plot
vis.gam(model, view = c("aspen", "spruce_fir"), plot.type = "contour", main = "Aspen:Spruce-fir")
# 3D surface plot
vis.gam(model, view = c("aspen", "spruce_fir"), plot.type = "persp", main = "Aspen:Spruce-fir")

# Testing for one species interaction

df <- grid_aspen %>%
 filter(lodgepole > 0) %>%
 mutate(aspen_sp = aspen + lodgepole) %>%
 filter(aspen_sp >= 50)

model <- gam(
 frp_max ~ s(aspen) + s(lodgepole) + te(aspen, lodgepole) +
  s(grid_x, grid_y, bs = "gp") + s(Fire_ID, bs = "re"),
 data = df, 
 family = gaussian()
)

summary(model)

plot(model, pages = 1)

# Visualize interaction as a contour plot
vis.gam(model, view = c("aspen", "lodgepole"), plot.type = "contour", main = "Aspen:Lodgepole")
# 3D surface plot
vis.gam(model, view = c("aspen", "lodgepole"), plot.type = "persp", main = "Aspen:Lodgepole")

# ##########################################
# 
# # loop through species interactions
# spp <- c("aspen", "douglas_fir", "lodgepole", "ponderosa", "spruce_fir", "piñon_juniper")
# 
# models <- list()
# summaries <- list()
# 
# # Loop through species
# for (species in spp) {
#  print(species)
#  if (species == "aspen") next  # Skip if species is "aspen"
#  
#  # Filter and prepare the dataframe
#  df <- grid_aspen[grid_aspen[[species]] > 0 & 
#                    (grid_aspen[["aspen"]] + grid_aspen[[species]]) >= 50, ]
#  
#  # define the formula
#  formula <- as.formula(
#   paste0(
#    "frp_max ~ s(aspen) + s(", species, ") + te(aspen, ", species, ") +",
#    " s(grid_x, grid_y, bs = 'gp') + s(Fire_ID, bs = 're')"
#   )
#  )
#  
#  # Fit GAM
#  model <- gam(
#   formula = formula,
#   data = df, 
#   family = gaussian()
#  )
#  
#  # Store the model and summary
#  models[[species]] <- model
#  summaries[[species]] <- summary(model)
# }
# 
# summary(models[["spruce_fir"]])
# summary(models[["lodgepole"]])
# summary(models[["ponderosa"]])
# summary(models[["douglas_fir"]])
# summary(models[["piñon_juniper"]])
# 
# plot(models[["spruce_fir"]], pages=1)
# plot(models[["lodgepole"]])
# 
# vis.gam(models[["spruce_fir"]], view = c("aspen", "spruce_fir"), plot.type = "persp", main = "Aspen:Spruce-fir")
# vis.gam(models[["lodgepole"]], view = c("aspen", "lodgepole"), plot.type = "persp", main = "Aspen:Lodgepole")
# 
# # retrieve the model results
# results <- lapply(names(summaries), function(species) {
#  summary <- summaries[[species]]
#  data.frame(
#   species = species,
#   term = rownames(summary$s.table),  # Smooth terms
#   edf = summary$s.table[, "edf"],   # Effective degrees of freedom
#   F_value = summary$s.table[, "F"], # F-statistic
#   p_value = summary$s.table[, "p-value"] # P-value
#  )
# }) %>%
#  bind_rows()  # Combine into one dataframe
# 
# # View the combined results
# print(results)
# 
# # plot a specific model
# library(mgcViz)
# # Get mgcViz object for a specific model
# viz <- getViz(models[["spruce_fir"]])
# # Plot smooth for aspen
# plot(sm(viz, 1)) + l_fitLine() + l_ciLine() + theme_classic()
# # Plot smooth for spruce_fir
# plot(sm(viz, 2)) + l_fitLine() + l_ciLine() + theme_classic()
# # Plot interaction term (aspen:spruce_fir)
# plotSlice(sm(viz, 3)) + l_fitRaster() + l_fitContour() + theme_classic()

# #########################
# # Set up the model (brms)
# 
# spp <- c("aspen", "douglas_fir", "lodgepole", "ponderosa", "spruce_fir", "piñon_juniper")
# # gather the interaction terms (aspen:other species)
# sp_interactions <- paste("aspen:", spp[spp != "aspen"], sep = "")
# 
# # define the model formula
# model.fit <- brm(
#  log(frp_max) ~ s(aspen) + s(lodgepole) + aspen:lodgepole + 
#                 gp(grid_x, grid_y) + 
#                 (1 | Fire_ID),
#  data = grid_aspen, family = gaussian(),
#  prior = c(set_prior("normal(0, 1)", class = "b")),  # Example priors
#  chains = 4, cores = 4
# )
# 
# summary(model.fit)
# 
# # Fit the model
# model <- gam(formula=f, data=grid_aspen, family=gaussian())
# 
# # Model summary
# summary(model)