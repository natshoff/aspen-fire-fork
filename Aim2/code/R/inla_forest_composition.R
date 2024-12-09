# Load the required libraries
library(tidyverse)
library(sf)

# Environment variabls
maindir <- '/Users/max/Library/CloudStorage/OneDrive-Personal/mcook/aspen-fire/Aim2/'

# load the aggregated FRP grid
grid <-  st_read(paste0(maindir,'data/tabular/mod/viirs_snpp_jpss1_gridstats_fortypcd.csv'))
glimpse(grid)

grid_ <- grid %>%
 select(c(grid_index, Fire_ID, frp_csum, frp_max, grid_x, grid_y, SpeciesName, spp_pct, forest_pct)) %>%
 mutate(frp_csum = as.numeric(frp_csum),
        frp_max = as.numeric(frp_max),
        grid_x = as.numeric(grid_x),
        grid_y = as.numeric(grid_y),
        forest_pct = as.numeric(forest_pct),
        spp_pct = as.numeric(spp_pct),
        Fire_ID = as.factor(Fire_ID))
head(grid_)

# reshape the dataframe for modeling
grid_w <- grid_ %>%
 mutate(SpeciesName = str_replace_all(SpeciesName, "-", "_"),
        SpeciesName = str_to_lower(SpeciesName)) %>%
 pivot_wider(names_from = SpeciesName, values_from = spp_pct, values_fill = 0)
head(grid_w)
nrow(grid_w)

# retain grid cells with some aspen component
grid_aspen <- grid_w %>%
 filter(aspen > 0)
nrow(grid_aspen)/nrow(grid_w)*100



##########################
# Set up the model (hgam)

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
vis.gam(model, view = c("aspen", "lodgepole"), plot.type = "contour", main = "Aspen:Spruce-fir")
# 3D surface plot
vis.gam(model, view = c("aspen", "lodgepole"), plot.type = "persp", main = "Aspen:Spruce-fir")

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