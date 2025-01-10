
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


#=========Prep the grid data=========#

# Format the species composition data frame 
# (long-format) each grid has rows for species co-occurring with the dominant type
# climate and topography are summarized at the grid level

# load the aggregated FRP grid with TreeMap and climate/topography
fp <- paste0(maindir,'data/tabular/mod/gridstats_fortypnm_gp_tm_ct_frp-cbi.csv')
grid_tm <-  read_csv(fp)  %>% # read in the file
 # get the acquisition year and month
 mutate(first_obs_date = as.Date(first_obs_date),  # Convert to Date
        year = year(first_obs_date),              # Extract year
        month = month(first_obs_date)) %>%
 # remove missing FRP, prep columns
 filter(
  frp_max_day > 0, # make sure daytime FRP is not 0
 ) %>% 
 # create a numeric fire ID
 mutate(
  # Set the random effects (IDs) as numeric factors
  Fire_ID = as.numeric(as.factor(Fire_ID)),
  grid_index = as.numeric(as.factor(grid_index)),
  # Format species names consistently
  fortypnm_gp = str_replace_all(fortypnm_gp, "-", "_"),
  fortypnm_gp = str_to_lower(fortypnm_gp),
  species_gp_n = str_replace_all(species_gp_n, "-", "_"),
  species_gp_n = str_to_lower(species_gp_n),
  # make sure factor variables are set for species names
  fortypnm_gp = as.factor(fortypnm_gp),
  species_gp_n = as.factor(species_gp_n),
  # proportion of contributing AFD during daytime observations
  day_prop = day_count / afd_count,
  # log-scaled outcomes
  log_frp_max = log(frp_max + 1e-5),
  log_frp_csum = log(frp_csum + 1e-5),
  log_frp_max_day = log(frp_max_day + 1e-5),
  log_frp_max_night = log(frp_max_night + 1e-5),
  # Create a "season" variable based on the month
  season = case_when(
   month %in% c(3, 4, 5) ~ "spring",
   month %in% c(6, 7, 8) ~ "summer",
   month %in% c(9, 10, 11) ~ "fall"
  ),
  # Year/season interaction
  year_season = interaction(year, season, drop = TRUE),
  # scale the percentages
  fortyp_pct = fortyp_pct / 100,
  forest_pct = forest_pct / 100,
  # center and scale continuous predictor variables
  across(
   c(ba_live, ba_dead, ba_ld, # live/dead basal area
     tpp_live, tpp_dead, tpp_ld, # live/dead trees/pixel
     qmd_live, qmd_dead, qmd_ld, # live/dead quadratic mean diameter
     tree_ht_live, tree_ht_dead, # live/dead tree height
     erc, vpd, vpd_dv, erc_dv, # climate
     elev, slope, tpi, chili # topography
    ),
   ~ as.numeric(scale(.))
  ),
  # center/scale abundance and dominance
  ba_ld_pr = ba_ld_pr - mean(ba_ld_pr, na.rm = TRUE),
  tpp_ld_pr = tpp_ld_pr - mean(tpp_ld_pr, na.rm = TRUE),
  qmd_ld_pr = qmd_ld_pr - mean(qmd_ld_pr, na.rm = TRUE),
  # create an aspen presence flag
  aspen_pres = if_else(species_gp_n == "aspen" & ba_live > 0, 1, 0)) %>%
 # be sure there are no duplicate rows for grid/fire/species
 distinct(grid_index, Fire_ID, species_gp_n, .keep_all = TRUE) # remove duplicates
glimpse(grid_tm) # check the results


########################################################
# Check on the grid cell counts for daytime observations
grid_counts <- grid_tm %>%
 distinct(Fire_ID, grid_index) %>% # keep only distinct rows
 group_by(Fire_ID) %>%
 summarise(n_grids = n())
summary(grid_counts)
quantile(grid_counts$n_grids, probs = seq(.1, .9, by = .1))

# Identify fires with n_grids below the 10th percentile
idx <- grid_counts %>%
 filter(n_grids < quantile(grid_counts$n_grids, probs = 0.2)) %>%
 pull(Fire_ID)
length(idx)
# filter the data frame to remove these fires
grid_tm <- grid_tm %>%
 filter(!Fire_ID %in% idx)

rm(grid_counts, idx)
gc()


#===============Explore Distributions, etc.================#

# list of species names
spps <- c("aspen", "mixed_conifer", "lodgepole", "ponderosa", "spruce_fir", "piñon_juniper")

####################################
# distribution of response variables
grid_tm %>%
 # Ensure daytime FRP and CBIbc is greater than 0
 filter(frp_max_day > 0,
        CBIbc_p90 > 0) %>%
 # pivot longer to facet plot
 pivot_longer(cols = c(log_frp_max_day, CBIbc_p90),
              names_to = "variable",
              values_to = "value") %>%
 # Plot with facets
 ggplot(aes(x = value)) +
 geom_histogram(bins = 30, fill = "orange", alpha = 0.7) +
 facet_wrap(~ variable, scales = "free", 
            labeller = as_labeller(c(log_frp_max_day = "log(Daytime Max FRP)", 
                                     CBIbc_p90 = "90th Percentile CBIbc"))) +
 labs(x = "value",
      y = "Frequency") +
 theme_minimal()

# save the plot.
out_png <- paste0(maindir,'figures/INLA_ResponseDistribution.png')
ggsave(out_png, dpi=500, bg = 'white')


####################
# correlation matrix
# Select only numeric columns and convert factors to dummy variables
cor_da <- grid_tm %>%
 select(
  tpp_ld_pr, ba_ld_pr, qmd_ld_pr, 
  ba_live, tpp_live, qmd_live, tree_ht_live,  # Forest composition metrics
  ba_dead, tpp_dead, qmd_dead, tree_ht_dead,
  erc, erc_dv, vpd, vpd_dv, elev, slope, tpi, chili  # Climate and topography metrics
 ) %>%
 mutate_if(is.factor, as.numeric)  # Convert factors to numeric (if needed)

# Compute correlation matrix
cor_mat <- cor(cor_da, use = "complete.obs")

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
out_png <- paste0(maindir,'figures/INLA_CorrelationMatrix.png')
ggsave(out_png, dpi=500, bg = 'white')


#===========MODEL SETUP==============#

set.seed(456)

# check the factor levels
# make sure aspen is first
levels(grid_tm$fortypnm_gp)

############################################
# 1. Baseline model: Predominant forest type
# ~ Effect of forest type on FRP and CBI relative to aspen 
# ~ FRP/CBIbc ~ dominant forest type + climate + topo
# ~ Fire_ID random effect
# ~ Spatial fields random effect

# create a data frame for just the predominant forest type
# one row per grid_index with the dominant forest type
# filter to where that type is at least 50% of the grid forested area
da <- grid_tm %>%
 select(grid_index, Fire_ID, year, year_season, # random effects 
        log_frp_max_day, CBIbc_p90, # response variables
        fortypnm_gp, fortyp_pct, forest_pct, # forest type and percent cover
        erc, erc_dv, vpd, vpd_dv, elev, slope, tpi, chili, # climate + topo
        aspen_pres, # whether aspen is present in the grid
        x, y # grid centroid coordinate for spatial fields model
        ) %>%
 filter(fortyp_pct > 0.50) %>% # only grids where predominant species cover > 50%
 # keep just one row per grid cell by predominant type
 distinct(grid_index, Fire_ID, fortypnm_gp, .keep_all = TRUE)

# Create a data frame for the predictors
# baseline model includes predominant forest type, climate, and topography
X <- da[,c(1:2,7,10:17)]
head(X)

########################################
# Create the spatial fields model (SPDE)
# Extract coordinates as a matrix
coords_mat <- as.matrix(select(da, x, y))
# Create a shared spatial mesh
mesh <- inla.mesh.2d(
 loc = coords_mat,                    # Locations (grid centroids)
 max.edge = c(10, 30),             # Maximum edge lengths (inner and outer)
 cutoff = 0.1,                      # Minimum distance between points
 offset = c(1, 3)               # Boundary buffer
)
plot(mesh)
points(coords_mat, col = "red", pch = 19, cex = 0.5)

# save the plot.
out_png <- paste0(maindir,'figures/INLA_MeshGrid.png')
ggsave(out_png, dpi=500, bg = 'white')

# Build the SPDE model
spde.bl <- inla.spde2.pcmatern(
 # Mesh and smoothness parameter
 mesh = mesh, 
 alpha = 2,
 # P(practic.range < 0.3) = 0.5
 prior.range = c(0.3, 0.5),
 # P(sigma > 1) = 0.01
 prior.sigma = c(10, 0.01)
)

# Compute the projector matrix (A)
A <- inla.spde.make.A(
 mesh, # the spatial mesh created above
 loc = coords_mat
)
dim(A) # should be equal to # of data locations X number of vertices
table(rowSums(A > 0)) 
table(rowSums(A))
table(colSums(A) > 0) # triangles with no point locations

# Assign the spatial index
field.idx <- inla.spde.make.index(
 "spatial_field", n.spde=mesh$n
)


#===========MODEL FITTING==============#

#####
# FRP

# Create the FRP INLA stack
stk.frp <- inla.stack(
 data = list(log_frp_max_day = da$log_frp_max_day),
 A = list(A, 1), 
 effects = list(c(field.idx),
                list(X)),
 tag = 'est')
dim(inla.stack.A(stk.frp))

# Set up the model formula
mf.frp <- log_frp_max_day ~ 
 fortypnm_gp + # predominant forest type
 erc + vpd + elev + slope + tpi + chili + # climate+topography
 f(Fire_ID, model="iid") + # fire-level random effect
 f(grid_index, model="iid") + # grid-level random effect
 f(spatial_field, model=spde.bl) # spatial model

# fit the model                     
model_bl.frp <- inla(
 mf.frp, data = inla.stack.data(stk.frp), 
 family = "gaussian",
 control.predictor = list(A = inla.stack.A(stk.frp)),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
)
summary(model_bl.frp)


#####
# CBI

# Create the CBI INLA stack
stk.cbi <- inla.stack(
 data = list(CBIbc_p90 = da$CBIbc_p90),
 A = list(A, 1), 
 effects = list(c(field.idx),
                list(X)),
 tag = 'est')

# adjust the priors to encourage smoother effect
spde.bl.alt <- inla.spde2.pcmatern(
 mesh = mesh, 
 alpha = 2,
 prior.range = c(1, 0.01), # Increase the practical range
 prior.sigma = c(1, 0.01) # Reduce the allowed standard deviation
)

# Set up the model formula
mf.cbi <- CBIbc_p90 ~ 
 fortypnm_gp + # dominant forest type
 erc + vpd + elev + slope + tpi + chili + # climate+topography
 f(Fire_ID, model="iid") + # fire-level random effect
 f(grid_index, model="iid") + # grid-level random effect
 f(spatial_field, model=spde.bl.alt) # spatial model

# fit the model                     
model_bl.cbi <- inla(
 mf.cbi, data = inla.stack.data(stk.cbi), 
 family = "gaussian",
 control.predictor = list(A = inla.stack.A(stk.cbi)),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
)
summary(model_bl.cbi)


#===========Plotting Posterior Effects============#

# Extract fixed effects related to fortypnm_gp

#####
# FRP
frp.eff <- as.data.frame(model_bl.frp$summary.fixed) %>%
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
cbi.eff <- as.data.frame(model_bl.cbi$summary.fixed) %>%
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
frp_marginals <- model_bl.frp$marginals.fixed
cbi_marginals <- model_bl.cbi$marginals.fixed

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
                    "CBI" = "CBI", 
                    "FRP" = "FRP")) +
 theme_classic() +
 theme(axis.text.y = element_text(angle = 0, hjust = 0, size=11),
       axis.text.x = element_text(angle = 0, hjust = 0, size=11),
       axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
       axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
       legend.position = c(0.89, 0.86),
       legend.background = element_rect(
        fill = scales::alpha("white", 0.4), 
        color = NA, size = 0.8),
       legend.title = element_text(size = 11),
       legend.text = element_text(size = 10))
p2

# save the plot.
out_png <- paste0(maindir,'figures/INLA_FORTYPNM_PosteriorEffects_Ridge.png')
ggsave(out_png, dpi=300, bg = 'white')


######################################
# Visualize the spatial random effects
spat.eff.frp <- inla.spde.make.A(mesh, coords_mat) %*% 
 model_bl.frp$summary.random$spatial_field$mean
spat.eff.cbi <- inla.spde.make.A(mesh, coords_mat) %*% 
 model_bl.cbi$summary.random$spatial_field$mean
# Add spatial effects to the data frame
spat.eff.df <- cbind(as.data.frame(coords_mat), spat_eff_frp = as.vector(spat.eff.frp))
spat.eff.df <- cbind(spat.eff.df, spat_eff_cbi = as.vector(spat.eff.cbi))
colnames(spat.eff.df) <- c("x", "y", "spat_eff_frp", "spat_eff_cbi")

# Convert to an sf object for mapping
spat.sf <- st_as_sf(spat.eff.df, coords = c("x", "y"), crs = 4326)  # Adjust CRS if needed

# Plot spatial effects
# FRP Map
frp_map <- ggplot(spat.sf) +
 geom_sf(aes(color = spat_eff_frp), size = 2) +
 scale_color_viridis_c(option = "plasma", name = "Spatial Effect (FRP)") +
 theme_minimal() +
 labs(title = "Spatial Random Effects (FRP)", x = "Longitude", y = "Latitude")

# CBI Map
cbi_map <- ggplot(spat.sf) +
 geom_sf(aes(color = spat_eff_cbi), size = 2) +
 scale_color_viridis_c(option = "plasma", name = "Spatial Effect (CBI)") +
 theme_minimal() +
 labs(title = "Spatial Random Effects (CBI)", x = "Longitude", y = "Latitude")

# Combine maps
frp_map + cbi_map

# save the plot.
out_png <- paste0(maindir,'figures/INLA_FORTYPNM_SpatialFieldsMap.png')
ggsave(out_png, dpi=500, bg = 'white')


########################################
# Density plot of spatial random effects
p.frp <- ggplot(spat.eff.df, aes(x = spat_eff_frp)) +
 geom_density(fill = "blue", alpha = 0.5) +
 geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
 theme_minimal() +
 labs(
  x = "Spatial Effect",
  y = "Density"
 )

p.cbi <- ggplot(spat.eff.df, aes(x = spat_eff_cbi)) +
 geom_density(fill = "blue", alpha = 0.5) +
 geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
 theme_minimal() +
 labs(
  x = "Spatial Effect",
  y = "Density"
 )

p.frp + p.cbi

# save the plot.
out_png <- paste0(maindir,'figures/INLA_FORTYPNM_SpatialFieldsDistribution.png')
ggsave(out_png, dpi=500, bg = 'white')

# Tidy up
rm(spat.eff.frp, spat.eff.cbi, spat.eff.df, spat.sf, p.frp, p.cbi)
gc()

############################################
# Contributions to explaining model variance

# Extract variance components for random effects
############################################
# Contributions to Explaining Model Variance

# Compute variance contributions for FRP
res.prec <- model_bl.frp$summary.hyperpar["Precision for the Gaussian observations", "mean"]
res.var <- 1 / res.prec

grid.prec <- model_bl.frp$summary.hyperpar["Precision for grid_index", "mean"]
grid.var <- 1 / grid.prec

fire.prec <- model_bl.frp$summary.hyperpar["Precision for Fire_ID", "mean"]
fire.var <- ifelse(is.na(fire.prec), 0, 1 / fire.prec) # Handle Fire_ID precision if zero or undefined

spat.sd <- model_bl.frp$summary.hyperpar["Stdev for spatial_field", "mean"]
spat.var <- spat.sd^2

fixed.eff <- model_bl.frp$summary.fitted.values$mean
fixed.var <- var(fixed.eff, na.rm = TRUE)

total.var <- res.var + fixed.var + grid.var + fire.var + spat.var

# Compute variance contributions for CBI
res.prec.cbi <- model_bl.cbi$summary.hyperpar["Precision for the Gaussian observations", "mean"]
res.var.cbi <- 1 / res.prec.cbi

grid.prec.cbi <- model_bl.cbi$summary.hyperpar["Precision for grid_index", "mean"]
grid.var.cbi <- 1 / grid.prec.cbi

fire.prec.cbi <- model_bl.cbi$summary.hyperpar["Precision for Fire_ID", "mean"]
fire.var.cbi <- ifelse(is.na(fire.prec.cbi), 0, 1 / fire.prec.cbi) # Handle Fire_ID precision if zero or undefined

spat.sd.cbi <- model_bl.cbi$summary.hyperpar["Stdev for spatial_field", "mean"]
spat.var.cbi <- spat.sd.cbi^2

fixed.eff.cbi <- model_bl.cbi$summary.fitted.values$mean
fixed.var.cbi <- var(fixed.eff.cbi, na.rm = TRUE)

total.var.cbi <- res.var.cbi + fixed.var.cbi + grid.var.cbi + fire.var.cbi + spat.var.cbi

# Combine FRP and CBI variance contributions into a single data frame
var.contr.full <- bind_rows(
 data.frame(
  Response = "FRP",
  Component = c("Residual", "Fixed Effects", "Grid-Level IID", "Fire-Level IID", "Spatial Field"),
  Variance = c(res.var, fixed.var, grid.var, fire.var, spat.var),
  Proportion = c(
   res.var / total.var,
   fixed.var / total.var,
   grid.var / total.var,
   fire.var / total.var,
   spat.var / total.var
  )
 ),
 data.frame(
  Response = "CBI",
  Component = c("Residual", "Fixed Effects", "Grid-Level IID", "Fire-Level IID", "Spatial Field"),
  Variance = c(res.var.cbi, fixed.var.cbi, grid.var.cbi, fire.var.cbi, spat.var.cbi),
  Proportion = c(
   res.var.cbi / total.var.cbi,
   fixed.var.cbi / total.var.cbi,
   grid.var.cbi / total.var.cbi,
   fire.var.cbi / total.var.cbi,
   spat.var.cbi / total.var.cbi
  )
 )
)

# Plot the variance contributions for both FRP and CBI
ggplot(var.contr.full, aes(x = Component, y = Proportion, fill = Response)) +
 geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
 scale_fill_manual(values = c("FRP" = "#FEB24C", "CBI" = "#800026")) +
 theme_light() +
 labs(
  title = "Variance Contributions by Model Component",
  x = "Model Component",
  y = "Proportion of Variance Explained",
  fill = "Response"
 )

# Save the plot
out_png <- paste0(maindir, 'figures/INLA_VarianceContributions_FRP_CBI.png')
ggsave(out_png, dpi = 500, bg = 'white')


###########
# tidy up !
rm(frp_marginals, cbi_marginals, tidy_marginals, tidy_frp, tidy_cbi, tidy_combined,
   model_bl.cbi, model_bl.frp, da, effects, res.var, res.prec, grid.var, grid.prec,
   spat.var, spat.sd, fixed.var, fixed.eff, total.var, var.contr)
gc() # garbage clean up












