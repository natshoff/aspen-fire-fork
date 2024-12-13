
# Load the required libraries
library(tidyverse)
library(sf) # spatial
library(INLA) # for spatial Bayes model
library(ggcorrplot)
library(lubridate)

# Environment variables
maindir <- '/Users/max/Library/CloudStorage/OneDrive-Personal/mcook/aspen-fire/Aim2/'


#=========Prep the grid data grid=========#

# Format the species composition data frame

# load the aggregated FRP grid with TreeMap and climate/topography
fp <- paste0(maindir,'data/tabular/mod/viirs_snpp_jpss1_gridstats_fortypcd_climtopo.csv')
grid_fortyp <-  read_csv(fp)  %>% # read in the file
 # get the acquisition year and month
 mutate(first_obs_date = as.Date(first_obs_date),  # Convert to Date
        year = year(first_obs_date),              # Extract year
        month = month(first_obs_date)) %>%
 # select the required columns
 select(c(grid_index, year, month, Fire_ID, first_obs_date, 
          frp_csum, frp_max, frp_max_day, frp_max_night,
          afd_count, day_count, night_count, overlap,
          grid_x, grid_y, SpeciesName, spp_pct, forest_pct,
          erc, erc_dv, vpd, vpd_dv, elev, slope, chili, tpi)) %>%
 # remove missing FRP, prep columns
 filter(frp_max > 0) %>% # make sure FRP is not 0
 # create a numeric fire ID
 mutate(Fire_ID = as.factor(Fire_ID),
        Fire_ID_nm = as.numeric(as.character(Fire_ID)),
        spp_pct = spp_pct / 100, # set cover 0-1
        # proportion of contributing AFD during daytime observations
        day_prop = day_count / afd_count) %>% 
 # be sure there are no duplicate rows
 distinct(grid_index, Fire_ID_nm, SpeciesName, .keep_all = TRUE) # remove duplicates
glimpse(grid_fortyp) # check the results


#####################################
# reshape the data frame for modeling
grid_fortyp_w <- grid_fortyp %>% 
 # tidy the species names
 mutate(SpeciesName = str_replace_all(SpeciesName, "-", "_"),
        SpeciesName = str_to_lower(SpeciesName),
        spp_pct = as.numeric(spp_pct)) %>%
 # pivot wider to get species as columns
 pivot_wider(
  names_from = SpeciesName, 
  values_from = spp_pct, 
  values_fill = 0) %>% # pivot wider
 # create the log-scaled outcome variable
 mutate(log_frp_max = log(frp_max + 1e-5),
        log_frp_csum = log(frp_csum + 1e-5),
        log_frp_max_day = log(frp_max_day + 1e-5),
        log_frp_max_night = log(frp_max_night + 1e-5),
        # create a conifer column
        conifer = rowSums(select(., douglas_fir, lodgepole, ponderosa, 
                                 spruce_fir, piñon_juniper), na.rm = TRUE),
        # create a "season" based on month
        season = case_when(
         month %in% c(3, 4, 5) ~ "spring",
         month %in% c(6, 7, 8) ~ "summer",
         month %in% c(9, 10, 11) ~ "fall",
        ),
        # year/season interaction
        year_season = interaction(year, season, drop = TRUE))

# double check some dimension
colnames(grid_fortyp_w)
nrow(grid_fortyp_w)


#############################################
# retain grid cells with some aspen component
grid_aspen_w <- grid_fortyp_w %>%
 filter(aspen > 0) # greater than 0 percent aspen cover
print("Percent of grids with some aspen component:")
nrow(grid_aspen_w)/nrow(grid_fortyp_w)*100

# tidy up
rm(grid_fortyp)
gc()



#===============Explore Distributions, etc.================#


#################################################
# distribution of raw frp_max and log-transformed
# Reshape the data to long format
grid_fortyp_w %>%
 # Ensure daytime FRP is greater than 0
 filter(frp_max_day > 0) %>%
 # pivot longer to facet plot
 pivot_longer(cols = c(frp_max, log_frp_max, log_frp_max_day),
              names_to = "variable",
              values_to = "value") %>%
 # Plot with facets
 ggplot(aes(x = value)) +
  geom_histogram(bins = 30, fill = "orange", alpha = 0.7) +
  facet_wrap(~ variable, scales = "free", 
             labeller = as_labeller(c(frp_max = "frp_max", 
                                      log_frp_max = "log(frp_max)",
                                      log_frp_max_day = "log(frp_max_day)"))) +
  labs(x = "value",
       y = "Frequency") +
  theme_minimal()


###########################################
# distribution of forest type percent cover
grid_fortyp_w %>%
 select(aspen, douglas_fir, spruce_fir, lodgepole, ponderosa, piñon_juniper) %>%
 pivot_longer(everything(), names_to = "species", values_to = "percent_cover") %>%
 ggplot(aes(x = percent_cover)) +
 geom_histogram() +
 facet_wrap(~ species, scales = "free") +
 theme_minimal()

# check for NA values
any(is.na(grid_fortyp_w)) # should be FALSE


###############################
# correlation matrix on effects
cor_matrix <- cor(
 grid_fortyp_w %>%
  as_tibble() %>%
  select(aspen, douglas_fir, lodgepole, ponderosa, spruce_fir, piñon_juniper,
         vpd_dv, erc_dv, elev, slope, tpi, chili), 
 use = "complete.obs")
# plot it with ggcorrplot
ggcorrplot(cor_matrix, method = "circle", type = "lower", 
           lab = TRUE, lab_size = 2.5)
# tidy up
rm(cor_matrix)
gc()



#===========MODEL SETUP==============#

set.seed(456)

#############################
# Set up the model data frame

# list of species names
spp <- c("aspen", "douglas_fir", "lodgepole", "ponderosa", "spruce_fir", "piñon_juniper")

# scale just the climate/topography effects variables
grid_sc <- grid_fortyp_w %>%
 # scale effects variables
 mutate(across(c(vpd_dv, erc_dv, elev, slope, tpi, chili), scale)) %>%
 # make sure the response variable > 0
 filter(frp_max_day > 0)

# define the model effects data
effects_da <- grid_sc %>%
 # remove all but the model effects
 select(-c(Fire_ID, grid_index, grid_x, grid_y, 
           frp_max, frp_csum, frp_max_night, frp_max_day,
           log_frp_csum, log_frp_max, log_frp_max_day, log_frp_max_night,
           first_obs_date, afd_count, day_count, night_count, overlap,
           forest_pct, vpd, erc, month)) %>%
 as_tibble()
colnames(effects_da)

# make sure the dimensions match
dim(grid_sc)
dim(effects_da)


####################################
# Generate species covariance matrix
sp_mat <- grid_sc %>%
 select(aspen, douglas_fir, lodgepole, ponderosa, spruce_fir, piñon_juniper)
# caompute the pearson correlation matrix
sp_cov <- cor(sp_mat, method = "pearson")  # Pearson correlation
# ensure positive semi-definitiveness
eigenvalues <- eigen(sp_cov)$values
print(eigenvalues)
# plot the matrix
ggcorrplot(sp_cov, method = "circle", lab = TRUE)


################################
# Calculate the "species effect"

sp_l <- grid_sc %>%
 select(grid_index, aspen, douglas_fir, lodgepole, ponderosa, spruce_fir, piñon_juniper) %>%
 pivot_longer(cols=-grid_index,
              names_to = "species", 
              values_to = "sp_cover") %>%
 mutate(sp_id = as.numeric(as.factor(species)),
        cidx = interaction(grid_index, sp_id, drop = TRUE))


# merge back the climate + topography
clim_topo <- grid_sc %>%
 mutate(across(c(vpd_dv, erc_dv, elev, slope, tpi, chili), as.numeric))

sp_l <- sp_l %>%
 left_join(
  clim_topo %>% 
   select(grid_index, Fire_ID_nm, year_season, year, grid_x, grid_y,
          vpd_dv, erc_dv, elev, slope, tpi, chili), 
  by = "grid_index")
glimpse(sp_l)

rm(clim_topo)
gc()


##############################################################
# 1. Baseline model without spatial component or random effect

# set the formula.
mf <- 
 log_frp_max_day ~ (aspen + douglas_fir + lodgepole + ponderosa + piñon_juniper + spruce_fir) + # species composition
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
 log_frp_max_day ~ (aspen + douglas_fir + lodgepole + ponderosa + piñon_juniper + spruce_fir) + # species composition
                   vpd_dv + erc_dv + slope + tpi + chili + # climate & topography
                   f(Fire_ID_nm, model = "iid")  # Random effect for fire-level variability

# fit the model                     
model_bl2 <- inla(
 mf2, data = grid_sc,
 family = "gaussian",
 control.predictor = list(compute = TRUE),
 control.compute = list(dic = TRUE, waic = TRUE)
)
summary(model_bl2)


########################
# Compare DIC and WAIC #

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
 cutoff = 0.01 # Minimum distance between points
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
 mesh = mesh,
 loc = coords
)
dim(A) # this should match the number of rows in our data
nrow(grid_sc) # data rows

# create an INLA stack
stack <- inla.stack(
 data = list(log_frp_max_day = grid_sc$log_frp_max_day),
 A = list(A, 1), # Link spatial field and fixed effects
 effects = list(
  spatial_field = 1:spde$n.spde,
  data.frame(effects_da)
 )
)

# define the formula.
mf3 <- log_frp_max_day ~ (aspen + douglas_fir + lodgepole + ponderosa + piñon_juniper + spruce_fir) +
                         vpd_dv + erc_dv + slope + tpi + chili +
                         f(Fire_ID_nm, model = "iid") +
                         f(spatial_field, model = spde)

model_bl3 <- inla(
 formula = mf3,
 family = "gaussian",
 data = inla.stack.data(stack),
 control.predictor = list(A = inla.stack.A(stack)),
 control.compute = list(dic = TRUE, waic = TRUE)
)

# Summarize the model results
summary(model_bl3)


##############################################
# 4. Adding temporal random effect (fire year)

# define the formula.
mf4 <- log_frp_max_day ~ aspen + douglas_fir + lodgepole + ponderosa + piñon_juniper + spruce_fir +
                         vpd_dv + erc_dv + slope + tpi + chili +
                         f(Fire_ID_nm, model = "iid") + # fire random effect
                         f(spatial_field, model = spde) + # spatial effect
                         f(year_season, model = "rw1")  # temporal random effect

model_bl4 <- inla(
 formula = mf4,
 family = "gaussian",
 data = inla.stack.data(stack),
 control.predictor = list(A = inla.stack.A(stack)),
 control.compute = list(dic = TRUE, waic = TRUE)
)

# Summarize the model results
summary(model_bl4)


######################
# Compare DIC and WAIC

cat("Spatial model: \n")
cat("DIC:", model_bl3$dic$dic, "\n")
cat("WAIC:", model_bl3$waic$waic, "\n\n")

cat("Spatial-temporal model (year): \n")
cat("DIC:", model_bl4$dic$dic, "\n")
cat("WAIC:", model_bl4$waic$waic, "\n")

# rm(model_bl1)
# gc()


#===========Plotting==============#

# Extract fixed effects
fixed_effects <- as.data.frame(model_bl4$summary.fixed)
fixed_effects$Variable <- rownames(fixed_effects)

# Exponentiate coefficients for interpretation
fixed_effects$mean_exp <- exp(fixed_effects$mean)
fixed_effects$lower_exp <- exp(fixed_effects$`0.025quant`)
fixed_effects$upper_exp <- exp(fixed_effects$`0.975quant`)

# Plot exponentiated coefficients (excluding intercept for clarity)
ggplot(fixed_effects %>% filter(Variable != "(Intercept)"), 
       aes(x = Variable, y = mean_exp, ymin = lower_exp, ymax = upper_exp)) +
 geom_pointrange() +
 coord_flip() +
 labs(y = "Effect on FRP", x = "Variable") +
 theme_minimal()


# Extract spatial field posterior mean
spatial_effects <- inla.spde.make.index("spatial_field", n.spde = spde$n.spde)
spatial_field_mean <- model_bl4$summary.random$spatial_field$mean

# Add to mesh
mesh_df <- data.frame(
 x = mesh$loc[, 1],
 y = mesh$loc[, 2],
 spatial_effect = spatial_field_mean
)
ggplot(mesh_df, aes(x = x, y = y, color = spatial_effect)) +
 geom_point(size = 1) +
 scale_color_viridis_c() +
 coord_equal() +
 labs(title = "Spatial Field (Posterior Mean)", x = "Longitude", y = "Latitude") +
 theme_minimal()


# Extract year random effect
year_effects <- model_bl4$summary.random$year_season

ggplot(year_effects, aes(x = ID, y = mean, ymin = `0.025quant`, ymax = `0.975quant`)) +
 geom_pointrange() +
 geom_line() +
 labs(title = "Year-Season Effect (Posterior Mean)", x = "Year-Season", y = "Effect Size") +
 theme_minimal() +
 theme(axis.text.x = element_text(angle = 45, hjust = 1))



#===========Model Setup (Interactions)==============#


