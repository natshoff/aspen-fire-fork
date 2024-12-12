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

# join to the spatial data
grid_fortyp_sp <- inner_join(grid, grid_fortyp, by="grid_index")
head(grid_fortyp_sp)

# tidy the columns
grid_ <- grid_fortyp_sp %>%
 as_tibble() %>% # operate without the geometry
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


rm(grid, grid_fortyp, grid_)
gc()


#====Explore Distributions, etc.====#

# distribution of raw frp_max and log-transformed
# Reshape the data to long format
dl <- grid_w %>%
 pivot_longer(cols = c(frp_max, log_frp_max),
              names_to = "variable",
              values_to = "value")

# Plot with facets
ggplot(dl, aes(x = value)) +
 geom_histogram(bins = 30, fill = "orange", alpha = 0.7) +
 facet_wrap(~ variable, scales = "free", 
            labeller = as_labeller(c(frp_max = "frp_max", log_frp_max = "log(frp_max)"))) +
 labs(x = "value",
      y = "Frequency") +
 theme_minimal()

# distribution of forest type percent cover
grid_w %>%
 select(aspen, douglas_fir, spruce_fir, lodgepole, ponderosa, piñon_juniper) %>%
 pivot_longer(everything(), names_to = "species", values_to = "percent_cover") %>%
 ggplot(aes(x = percent_cover)) +
 geom_histogram() +
 facet_wrap(~ species, scales = "free") +
 theme_minimal()

# check for NA values
any(is.na(grid_w))
# sapply(grid_w, length)


##############################
# correlation matrix on effects
effects_da <- grid_w %>%
 as_tibble() %>%
 select(aspen, douglas_fir, lodgepole, ponderosa, spruce_fir, piñon_juniper,
        vpd_dv, erc_dv, elev, slope, tpi, chili)

# create the matrix
cor_matrix <- cor(effects_da, use = "complete.obs")
# plot it with ggcorrplot
ggcorrplot(cor_matrix, method = "circle", type = "lower", lab = FALSE)



#===========MODEL SETUP==============#

set.seed(456)

spp <- c("aspen", "douglas_fir", "lodgepole", "ponderosa", "spruce_fir", "piñon_juniper")

# # scale the effects variables
# grid_w <- grid_w %>%
#  mutate(across(c(vpd, erc, elev, slope, tpi, chili,
#                  aspen, douglas_fir, spruce_fir, lodgepole, ponderosa, piñon_juniper), scale))

# scale just the climate/topography effects variables
grid_sc <- grid_w %>%
 mutate(across(c(vpd_dv, erc_dv, elev, slope, tpi, chili), scale),
        Fire_ID_nm = as.numeric(as.factor(Fire_ID))) %>%
 select(-c(Fire_ID, grid_index, frp_max, frp_csum, grid_x, grid_y, 
           forest_pct, vpd, erc, geom))

# define the model effects data
effects_da <- grid_sc %>%
 select(-c(log_frp_max))
colnames(effects_da)

dim(grid_sc)
dim(effects_da)

##############################################################
# 1. Baseline model without spatial component or random effect

# define the formula
mf <- 
 log_frp_max ~ aspen + douglas_fir + lodgepole + ponderosa + piñon_juniper + spruce_fir + # species composition
               vpd_dv + erc_dv + elev + slope + tpi + chili # climate & topography

# fit the model                     
model_bl1 <- inla(
 mf, data = grid_sc, # data with scaled climate/topo
 family = "gaussian",
 control.predictor = list(compute = TRUE),
 control.compute = list(dic = TRUE, waic = TRUE)
)
summary(model_bl1)


###########################################################
# 2. Adding between fires effect (random effect on fire ID)

# define the formula
mf2 <- 
 log_frp_max ~ aspen + douglas_fir + lodgepole + ponderosa + piñon_juniper + spruce_fir + # species composition
               vpd + erc + elev + slope + tpi + chili + # climate & topography
               f(Fire_ID_nm, model = "iid")  # Random effect for fire-level variability

# fit the model                     
model_bl2 <- inla(
 mf2, data = grid_sc,
 family = "gaussian",
 control.predictor = list(compute = TRUE),
 control.compute = list(dic = TRUE, waic = TRUE)
)

######################
# Compare DIC and WAIC
cat("Baseline Model: \n")
cat("DIC:", model_bl1$dic$dic, "\n")
cat("WAIC:", model_bl1$waic$waic, "\n\n")

cat("With Fire_ID Random Effect: \n")
cat("DIC:", model_bl2$dic$dic, "\n")
cat("WAIC:", model_bl2$waic$waic, "\n")

rm(model_bl1)
gc()

############################################
# 3. Adding "within-fire" spatial dependence

# Extract coordinates
coords <- as.matrix(grid_sc[, c("grid_x", "grid_y")])
# Create a shared spatial mesh
mesh <- inla.mesh.2d(
 loc = coords,
 max.edge = c(10, 100),  
 cutoff = 0.02 # Minimum distance between points
)
plot(mesh)

# define the stochastic partial difference equation (SPDE)
spde <- inla.spde2.pcmatern(
 mesh = mesh,
 alpha = 2,  # Smoothness parameter
 prior.range = c(10, 0.01),  # Prior for spatial range
 prior.sigma = c(5, 0.01)    # Prior for variance
)

# link the mesh to data (effects) (A-Matrix)
A <- inla.spde.make.A(
 mesh = mesh, # mesh grid
 loc = coords, # coordinate matrix
 group = grid_sc$Fire_ID_nm  # Group spatial effects by Fire_ID
)
dim(A) # this should match the number of rows in our data
nrow(grid_sc) # data rows

# create spatial field index
spatial_field_idx <- rep(1:spde$n.spde, times = length(unique(grid_sc$Fire_ID_nm)))
length(spatial_field_idx) == dim(A)[2]

# create an INLA stack
stack <- inla.stack(
 data = list(log_frp_max = grid_sc$log_frp_max),  # Response variable
 A = list(A, 1),  # Link spatial and fixed effects
 effects = list(
  spatial_field = 1:dim(A)[2],  # Correct spatial field index
  grid_sc %>% select(-c(log_frp_max))  # Use the entire dataset as fixed effects
 )
)

length(grid_sc$Fire_ID_nm) == nrow(grid_sc)  # Should return TRUE
dim(A)[1] == nrow(grid_sc)  # Should return TRUE

# define the formula.
mf3 <- log_frp_max ~ aspen + douglas_fir + lodgepole + ponderosa + piñon_juniper + spruce_fir +
                     vpd_dv + erc_dv + elev + slope + tpi + chili +
                     f(Fire_ID_nm, model = "iid") +
                     f(spatial_field, model = spde, group = group_fixed)

model_bl3 <- inla(
 formula = mf3,
 family = "gaussian",
 data = inla.stack.data(stack),
 control.predictor = list(A = inla.stack.A(stack), compute = TRUE),
 control.compute = list(dic = TRUE, waic = TRUE)
)

# Summarize the model results
summary(model_bl3)

group_fixed <- rep(grid_sc$Fire_ID_nm, each = 1)
length(group_fixed) == dim(A)[2]  # Should return TRUE

# Length of group matches columns in A
length(group_fixed) == dim(A)[2]

# Spatial field indices match columns in A
length(spatial_field_idx) == dim(A)[2]

# check
dim(inla.stack.A(stack))
dim(inla.stack.data(stack))

spde$n.spde * length(unique(grid_sc$Fire_ID_nm)) == dim(A)[2]  # Should return TRUE




length(grid_sc$Fire_ID_nm)
spde$n.spde
unique(grid_sc$Fire_ID_nm)



dim(effects_da)  # Should match nrow(grid_sc)
stopifnot(nrow(effects_da) == nrow(grid_sc))



spatial_field <- data.frame(spatial_field = 1:ncol(A))
head(spatial_field)
any(is.na(grid_w$Fire_ID_nm))
length(grid_w$Fire_ID_nm) == nrow(grid_w)  # Should return TRUE

# create an INLA  data stack
stack <- inla.stack(
 data = list(log_frp_max = grid_sc$log_frp_max),
 A = list(A, 1),
 effects = list(
  spatial_field = spatial_field,
  effects_df
 )
)



# Expand baseline model with spatial effects
mf2 <- log_frp_max ~ aspen + douglas_fir + spruce_fir + lodgepole + ponderosa + piñon_juniper +
 f(Fire_ID, model = "iid") + # Between-fire random effect
 f(spatial_field, model = spde, group = Fire_ID_nm) # Within-fire spatial effect

# fit the model.
model_bl2 <- inla(
 formula = mf2,
 data = inla.stack.data(stack),
 control.predictor = list(A = inla.stack.A(stack), compute = TRUE),
 control.compute = list(dic = TRUE, waic = TRUE),
 family = "gaussian"
)
summary(model_bl2)

 


# library(sp)
# grid.sp <- as(grid_aspen, "Spatial")
# idx.mapping <- as.vector(t(matrix(1:50, nrow = 10, ncol = 5)))
# grid.sp2 <- grid.sp[idx.mapping, ]



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
