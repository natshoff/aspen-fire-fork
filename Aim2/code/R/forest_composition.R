
# Load the required libraries
library(tidyverse)
library(sf) # spatial
library(INLA) # for spatial Bayes model
library(ggcorrplot)
library(lubridate)

# Environment variables
maindir <- '/Users/max/Library/CloudStorage/OneDrive-Personal/mcook/aspen-fire/Aim2/'


#=========Prep the grid data=========#

# Format the species composition data frame 
# (long-format) each grid has rows for species occurrence
# climate and topography are summarized at the grid level

# load the aggregated FRP grid with TreeMap and climate/topography
fp <- paste0(maindir,'gridstats_fortypnm_gp_tm_ct_frp-cbi.csv')
grid_tm <-  read_csv(fp)  %>% # read in the file
 # get the acquisition year and month
 mutate(first_obs_date = as.Date(first_obs_date),  # Convert to Date
        year = year(first_obs_date),              # Extract year
        month = month(first_obs_date)) %>%
 # remove missing FRP, prep columns
 filter(frp_max_day > 0) %>% # make sure daytime FRP is not 0
 # create a numeric fire ID
 mutate(Fire_ID = as.factor(Fire_ID),
        Fire_ID_nm = as.numeric(as.character(Fire_ID)),
        # Format species names consistently
        fortypnm_gp = str_replace_all(fortypnm_gp, "-", "_"),
        fortypnm_gp = str_to_lower(fortypnm_gp),
        species_gp_n = str_replace_all(species_gp_n, "-", "_"),
        species_gp_n = str_to_lower(species_gp_n),
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
        year_season = interaction(year, season, drop = TRUE)) %>%
 # calculate the percent conifer
 group_by(grid_index) %>%
 mutate(
  conifer = max(
   fortyp_pct[fortypnm_gp %in% c("mixed_conifer", "lodgepole", "ponderosa", 
                                 "spruce_fir", "piñon_juniper")], na.rm = TRUE)
  ) %>% ungroup() %>%
 # be sure there are no duplicate rows
 distinct(grid_index, Fire_ID_nm, species_gp_n, .keep_all = TRUE) # remove duplicates
glimpse(grid_tm) # check the results


df <- grid_tm %>% slice(1:100)
write.csv(df, paste0(maindir,'data/tabular/mod/TEMP_100rows_gridstats_example.csv'))


#############################
# create a wider format table 
# (one row for each grid, species percent as columns)
grid_fortyp_w <- grid_fortyp %>% 
 # pivot wider to get species as columns
 pivot_wider(
  names_from = species_gp_n, 
  values_from = spp_pct, 
  values_fill = 0)
# double check some dimension
colnames(grid_fortyp_w)
nrow(grid_fortyp_w)



#===============Explore Distributions, etc.================#

# list of species names
spp <- c("aspen", "douglas_fir", "lodgepole", "ponderosa", "spruce_fir", "piñon_juniper")

###########################################
# distribution of forest type percent cover
grid_fortyp_w %>%
 select(aspen, douglas_fir, spruce_fir, lodgepole, ponderosa, piñon_juniper) %>%
 pivot_longer(everything(), names_to = "species", values_to = "percent_cover") %>%
 ggplot(aes(x = percent_cover)) +
 geom_histogram() +
 facet_wrap(~ species, scales = "free") +
 theme_minimal()

##############################
# Distribution of forest types
grid_fortyp_w %>%
 select(all_of(spp)) %>%
 summarize(across(everything(), ~ sum(. > 0))) %>%
 pivot_longer(cols = everything(), names_to = "species", values_to = "count") %>%
 ggplot(aes(x = reorder(species, -count), y = count, fill = species)) +
  geom_bar(stat = "identity") +
  labs(
   x = "Species",
   y = "Number of Grid Cells"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

#################################################
# distribution of raw frp_max and log-transformed
# Reshape the data to long format
grid_fortyp_w %>%
 # Ensure daytime FRP is greater than 0
 filter(frp_max_day > 0) %>%
 # pivot longer to facet plot
 pivot_longer(cols = c(frp_max, log_frp_max_day, log_frp_csum),
              names_to = "variable",
              values_to = "value") %>%
 # Plot with facets
 ggplot(aes(x = value)) +
  geom_histogram(bins = 30, fill = "orange", alpha = 0.7) +
  facet_wrap(~ variable, scales = "free", 
             labeller = as_labeller(c(frp_max = "frp_max", 
                                      log_frp_max_day = "log(frp_csum)",
                                      log_frp_csum = "log(frp_max_day)"))) +
  labs(x = "value",
       y = "Frequency") +
  theme_minimal()

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

##########################################
# Generate diversity indices for each grid
# 'shannon_H' = dominance-adjusted Shannon index
spp_effect <- grid_fortyp_w %>%
 rowwise() %>%
 mutate(
  shannon = -sum(c_across(all_of(spp)) * log(c_across(all_of(spp))), na.rm = TRUE),
  simpson = 1 - sum(c_across(all_of(spp))^2, na.rm = TRUE),
  shannon_H = shannon * simpson
 ) %>%
 ungroup() %>%
 mutate(across(c(vpd_dv, erc_dv, elev, slope, tpi, chili), scale))
head(spp_effect%>%select(grid_index,shannon,simpson,shannon_H))

####################################
# Generate species covariance matrix
sp_mat <- grid_fortyp_w %>%
 select(all_of(spp))
# compute the pearson correlation matrix
sp_cov <- cor(sp_mat, method = "spearman")  # Spearman rank correlation
# ensure positive semi-definite
eigenvalues <- eigen(sp_cov)$values
print(eigenvalues)
# plot the matrix
ggcorrplot(sp_cov, method = "circle", lab = TRUE)

# tidy up
rm(sp_mat, eigenvalues, grid_fortyp_w)
gc()

###############################################
# Update the covariance matrix for grid-species
# Determine the number of grid cells
n_grids <- length(unique(grid_fortyp$grid_index))
# Construct a block-diagonal covariance matrix
sp_cov_grid <- Matrix::bdiag(replicate(n_grids, sp_cov, simplify = FALSE))  # Block-diagonal matrix
dim(sp_cov_grid)  # Should be (n_grids * 6) x (n_grids * 6)
Matrix::isSymmetric(sp_cov_grid)  # Should return TRUE

# tidy up.
rm(sp_cov)
gc()

#############################################
# Calculate the Latent Field "species effect"
lat_spp_effect <- grid_fortyp %>%
 # filter to remove 0 spp cover
 filter(spp_pct > 0) %>%
 mutate(
  spp_id = as.numeric(as.factor(SpeciesName)),  # Assign numeric IDs to species
  cidx = as.numeric(as.factor(interaction(grid_index, spp_id, drop = TRUE))),  # Grid × species interaction index
  # Scale climate/topography variables
  across(c(vpd_dv, erc_dv, elev, slope, tpi, chili), scale)
 )
glimpse(lat_spp_effect)
length(unique(lat_spp_effect$cidx))  # Should match the number of unique grid × species combinations

# tidy up.
rm(grid_fortyp)
gc()





#===========MODEL FITTING==============#

##############################################################
# 1. Baseline model without spatial component or random effect
# Species cover, climate, topography as fixed effects
# Dominance-adjusted Shannon diversity index by species contributions

mf <- log_frp_max_day ~ shannon_H * (aspen + lodgepole + ponderosa + spruce_fir + piñon_juniper + douglas_fir) + # species effects
                        vpd_dv + erc_dv + elev + slope + tpi + chili + # climate+topography
 f(as.factor(Fire_ID), model="iid") # random effect for fire ID
                    
# fit the model                     
model_bl1 <- inla(
 mf, data = spp_effect, 
 family = "gaussian",
 control.predictor = list(compute = TRUE),
 control.compute = list(dic = TRUE, waic = TRUE)
)

summary(model_bl1)


##############################################################
# 1. Baseline model without spatial component or random effect

# set the formula.
mf <- log_frp_max_day ~ vpd_dv + erc_dv + elev + slope + tpi + chili +  # Climate/topography
                        # Latent fields for species composition
                        f(cidx, spp_pct, model = "generic0", Cmatrix = sp_cov_grid, 
                          hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))
# fit the model                     
model_bl1 <- inla(
 mf, data = spp_effect, 
 family = "gaussian",
 control.predictor = list(compute = TRUE),
 control.compute = list(dic = TRUE, waic = TRUE)
)
summary(model_bl1)


###########################################################
# 2. Adding between fires effect (random effect on fire ID)

# set the formula.
mf2 <- log_frp_max_day ~ vpd_dv + erc_dv + elev + slope + tpi + chili +  # Climate/topography
                         f(cidx, spp_pct, model = "generic0", Cmatrix = sp_cov_grid, 
                           hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
                         f(Fire_ID_nm, model = "iid") # Fire ID random effect

# fit the model                     
model_bl2 <- inla(
 mf2, data = spp_effect, 
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

print("Keeping better model")
if (model_bl1$waic$waic > model_bl2$waic$waic) {
 rm(model_bl1) # clean up
 gc()
} else {
 rm(model_bl2) # clean up
 gc()
}


############################################
# 3. Adding "within-fire" spatial dependence

# Extract coordinates from wide data frame
coords <- spp_effect %>% distinct(grid_index, grid_x, grid_y)
coords_mat <- as.matrix(coords[, c("grid_x", "grid_y")])
# Create a shared spatial mesh
mesh <- inla.mesh.2d(
 loc = coords_mat,
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

# create the A-matrix (linking mesh to coords)
A <- inla.spde.make.A(
 mesh = mesh,
 loc = coords_mat
)

# Create the INLA stack
stack <- inla.stack(
 data = list(log_frp_max_day = spp_effect$log_frp_max_day),  # Response variable
 A = list(A, diag(nrow(spp_effect))),  # Sparse spatial field and identity matrix
 effects = list(
  spatial_field = 1:spde$n.spde,  # Spatial random field
  spp_effect %>% select(-log_frp_max_day)  # All covariates except response
 )
)

#####################
# define the formula.
mf3 <- log_frp_max_day ~ vpd_dv + erc_dv + elev + slope + tpi + chili +  # Climate/topography
                         f(cidx, spp_pct, model = "generic0", Cmatrix = sp_cov_grid, 
                           hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
                         f(Fire_ID_nm, model = "iid") + # Fire ID random effect
                         f(spatial_field, model = spde) # spatial effect

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


