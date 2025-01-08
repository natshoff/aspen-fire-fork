
# Load the required libraries
library(tidyverse)
library(sf) # spatial
library(INLA) # for spatial Bayes model
library(ggcorrplot)
library(lubridate)
library(ggridges)
library(reshape2)

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
  frp_max_day > 1, # make sure daytime FRP is not 0
 ) %>% 
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
        year_season = interaction(year, season, drop = TRUE),
        # make sure factor variables are set
        fortypnm_gp = as.factor(fortypnm_gp),
        species_gp_n = as.factor(species_gp_n),
        Fire_ID = as.factor(Fire_ID),
        # create an aspen presence flag
        aspen_pres = if_else(species_gp_n == "aspen" & balive > 0, 1, 0),
        # scale continuous predictors
        across(
         c(balive, badead, tpa_live, tpa_dead,
           vpd_dv, erc_dv, elev, slope, tpi, chili),
         ~ as.numeric(scale(.))
        ),
        # center/scale abundance and dominance
        sp_abundance_ld_c = sp_abundance_ld - mean(sp_abundance_ld, na.rm = TRUE),
        sp_dominance_ld_c = sp_dominance_ld - mean(sp_dominance_ld, na.rm = TRUE)) %>%
 # be sure there are no duplicate rows
 distinct(grid_index, Fire_ID, species_gp_n, .keep_all = TRUE) # remove duplicates
glimpse(grid_tm) # check the results



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

####################
# correlation matrix
# Select only numeric columns and convert factors to dummy variables
cor_da <- grid_tm %>%
 select(
  sp_abundance_ld, sp_dominance_ld, balive, tpa_live, tree_ht_live,  # Forest composition metrics
  erc_dv, vpd_dv, elev, slope, tpi, chili  # Climate and topography metrics
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


#===========MODEL SETUP==============#

set.seed(456)

# check the factor levels
# make sure aspen is first
levels(grid_tm$fortypnm_gp)


#===========MODEL FITTING==============#

############################################
# 1. Baseline model: Predominant forest type
# ~ Effect of forest type on FRP and CBI relative to aspen 
# ~ FRP/CBIbc ~ dominant forest type + climate + topo
# ~ Fire_ID random effect
# ~ Spatial fields random effect

# create a data frame for just the predominant forest type
da <- grid_tm %>%
 select(grid_index, Fire_ID, year, # random effects 
        log_frp_max_day, CBIbc_p90, # response variables
        fortypnm_gp, fortyp_pct, forest_pct, # forest type and percent cover
        erc_dv, vpd_dv, elev, slope, tpi, chili, # climate + topo
        x, y # grid centroid coordinate for spatial fields model
        ) %>%
 filter(fortyp_pct > 50) %>% # only grids where predominant species cover > 50%
 # keep just one row per grid cell
 distinct(grid_index, Fire_ID, fortypnm_gp, .keep_all = TRUE)
 
# Create the spatial fields model (SPDE)
# Extract coordinates from wide data frame
coords <- da %>% distinct(grid_index, x, y)
coords_mat <- as.matrix(coords[, c("x", "y")])
# Create a shared spatial mesh
mesh <- inla.mesh.2d(
 loc = coords_mat,
 max.edge = c(1, 5),  
 cutoff = 0.01 # Minimum distance between points
)
plot(mesh)

# define the stochastic partial difference equation (SPDE)
spde <- inla.spde2.pcmatern(
 mesh = mesh,
 alpha = 2,  # Smoothness parameter
 prior.range = c(10, 0.01),  # Prior for spatial range
 prior.sigma = c(1, 0.01)    # Prior for variance
)

# create the A-matrix (linking mesh to coords)
A <- inla.spde.make.A(
 mesh = mesh,
 loc = coords_mat,
 group = as.numeric(as.factor(da$Fire_ID))
)

# Define the spatial index
spatial_index <- inla.spde.make.index(
 name = "spatial",
 n.spde = spde$n.spde,
 group = as.numeric(as.factor(da$Fire_ID)),
 n.group = length(unique(da$Fire_ID))
)

# Rename the group variable explicitly
names(spatial_index) <- c("spatial", "spatial.group", "spatial.repl")

# Create the INLA stack
stack.frp <- inla.stack(
 data = list(log_frp_max_day = da$log_frp_max_day),  # Response variable
 A = list(A, 1),  # Sparse spatial field and identity matrix
 effects = list(
  spatial_index,  # Spatial field
  da %>% 
   select(fortypnm_gp, erc_dv, vpd_dv, elev, slope, tpi, chili, Fire_ID) %>%  # Covariates
   as.data.frame()
 )
)

# FRP

# Set up the model formula
mf.frp <- log_frp_max_day ~ 
 fortypnm_gp + # dominant forest type
 erc_dv + vpd_dv + elev + slope + tpi + chili + # climate+topography
 f(Fire_ID, model = "iid") + # Fire ID random effect
 f(spatial, model = spde, group = spatial.group) # spatial effect

# fit the model                     
model_bl.frp <- inla(
 mf.frp, data = inla.stack.data(stack.frp), 
 family = "gaussian",
 control.predictor = list(A = inla.stack.A(stack.frp), compute = TRUE),
 control.compute = list(dic = TRUE, waic = TRUE)
)
summary(model_bl.frp)


# CBI

# Re-create the INLA stack for CBI
stack.frp <- inla.stack(
 data = list(CBIbc_p90 = da$CBIbc_p90),  # Response variable
 A = list(A, 1),  # Sparse spatial field and identity matrix
 effects = list(
  spatial_index,  # Spatial field
  da %>% 
   select(fortypnm_gp, erc_dv, vpd_dv, elev, slope, tpi, chili, Fire_ID) %>%  # Covariates
   as.data.frame()
 )
)

# Set up the model formula
mf.cbi <- CBIbc_p90 ~ 
 fortypnm_gp + # dominant forest type
 erc_dv + vpd_dv + elev + slope + tpi + chili + # climate+topography
 f(Fire_ID, model = "iid") # Fire ID random effect

# fit the model                     
model_bl.cbi <- inla(
 mf.cbi, data = da, 
 family = "gaussian",
 control.predictor = list(compute = TRUE),
 control.compute = list(dic = TRUE, waic = TRUE)
)
summary(model_bl.cbi)


#===========Plotting Effects============#

# Extract fixed effects related to fortypnm_gp

# fixed effects for FRP
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

# fixed effects for CBIbc
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
out_png <- paste0(maindir,'figures/INLA_PredominantFORTYP_PosteriorEffects.png')
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
 theme(axis.text.y = element_text(angle = 0, hjust = 0, size=10),
       axis.text.x = element_text(angle = 0, hjust = 0, size=10),
       axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
       axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
       legend.position = c(0.10, 0.86),
       legend.background = element_rect(
        fill = scales::alpha("white", 0.4), 
        color = NA, size = 0.8),
       legend.title = element_text(size = 11),
       legend.text = element_text(size = 10))
p2

# save the plot.
out_png <- paste0(maindir,'figures/INLA_PredominantFORTYP_PosteriorEffects_Ridge.png')
ggsave(out_png, dpi=500, bg = 'white')

rm(frp_marginals, cbi_marginals, tidy_marginals, tidy_frp, tidy_cbi, tidy_combined,
   model_bl.cbi, model_bl.frp, da, effects)
gc() # clean up


##############################################################
# 2. Species dominance (basal area) interactions
# ~ FRP/CBIbc ~ (species_balive)^2 + climate + topo

# Prepare the grid data for the model
da <- grid_tm %>%
 # Pivot the data to wide format: each species is a column
 pivot_wider(
  id_cols = c(grid_index, Fire_ID, log_frp_max_day, CBIbc_p90, fortypnm_gp,
              erc_dv, vpd_dv, elev, slope, tpi, chili),  # Retain relevant variables
  names_from = species_gp_n,  # Use species name as column names
  values_from = balive,  # Use BALIVE as values for each species
  names_prefix = "balive_"  # Prefix columns with "balive_"
 ) %>%
 # Replace NA with 0 (indicates species absence in grid cell)
 mutate(across(starts_with("balive_"), ~ replace_na(., 0))) %>%
 # keep just one row per grid cell
 distinct(grid_index, Fire_ID, fortypnm_gp, .keep_all = TRUE)
glimpse(da)


#########################
# setup the model formula
mf.frp <- log_frp_max_day ~ 
 (balive_aspen + balive_lodgepole + balive_mixed_conifer + 
  balive_ponderosa + balive_spruce_fir + balive_piñon_juniper)^2 +  # Main effects and interactions
 erc_dv + vpd_dv + elev + slope + tpi + chili +  # Climate and topography predictors
 f(Fire_ID, model = "iid")  # Random effect for Fire ID

# fit the model                     
model_bl_tm.frp <- inla(
 mf.frp, data = da, 
 family = "gaussian",
 control.predictor = list(compute = TRUE, link = 1),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE)
)
summary(model_bl_tm.frp)


###########################
# Extract interaction terms
effects <- as.data.frame(model_bl_tm.frp$summary.fixed) %>%
 rownames_to_column(var = "parameter") %>%
 filter(str_detect(parameter, ":")) %>%
 separate(parameter, into = c("species1", "species2"), sep = ":") %>%
 mutate(across(everything(), ~ str_replace(., "balive_", "")))

# Prepare data for heatmap
heatmap_da <- effects %>%
 select(species1, species2, mean) %>%
 pivot_wider(names_from = species2, values_from = mean) %>%
 column_to_rownames(var = "species1")

# Convert to long format for ggplot
heatmap_l <- melt(heatmap_da, varnames = c("Species1", "Species2"), value.name = "Effect")

# Plot heatmap
ggplot(heatmap_l, aes(x = Species2, y = Species1, fill = Effect)) +
 geom_tile(color = "white") +
 scale_fill_gradient2(
  low = "blue", mid = "white", high = "red", midpoint = 0,
  limits = c(min(heatmap_l$Effect, na.rm = TRUE), max(heatmap_l$Effect, na.rm = TRUE))
 ) +
 labs(
  title = "Pairwise Interaction Effects on FRP",
  x = "Co-occurring Species",
  y = "Dominant Species",
  fill = "Effect"
 ) +
 theme_minimal() +
 theme(axis.text.x = element_text(angle = 45, hjust = 1))








##########################################################
# Extract fixed effect marginals for the interaction terms

# mf.frp <- log_frp_max_day ~ 
#  fortypnm_gp + # predominant forest type
#  (species_gp_n * balive) + # species dominance (balive)
#  erc_dv + vpd_dv + elev + slope + tpi + chili + # climate + topo
#  f(Fire_ID, model = "iid") # random effect for Fire ID

marginals <- model_bl_tm.frp$marginals.fixed

# Define a function to tidy marginals for interaction terms
tidy_int <- function(marginals) {
 tibble::tibble(
  parameter = names(marginals),
  data = purrr::map(marginals, ~ as.data.frame(.x))
 ) %>%
  unnest(data) %>%
  filter(str_detect(parameter, "species_gp_n") & str_detect(parameter, ":balive")) %>%
  mutate(
   species = str_remove(parameter, ":balive"),
   species = str_remove(species, "species_gp_n"),
   species = str_replace(species, "_", " "),
   effect_type = "Interaction"
  )
}

# Tidy the marginals for the interaction terms
tidy_int_df <- tidy_int(marginals)

# Create ridge plot for interaction terms
ggplot(tidy_int_df, aes(x = x, y = species, height = y, fill = effect_type)) +
 geom_density_ridges(stat = "identity", scale = 1.5, alpha = 0.7) +
 geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
 labs(
  x = "Effect Magnitude",
  y = "Co-occurring Species",
  fill = "Effect Type"
 ) +
 scale_fill_manual(values = c("Interaction" = "#FEB24C")) +
 theme_classic() +
 theme(
  axis.text.y = element_text(size = 10),
  axis.title.y = element_text(size = 11, margin = margin(r = 10)),
  axis.title.x = element_text(size = 11, margin = margin(t = 10)),
  legend.position = "top"
 )

# save the plot.
out_png <- paste0(maindir,'figures/INLA_Species-BALIVE_PosteriorEffects_Ridge.png')
ggsave(out_png, dpi=500, bg = 'white')



##############################################################
# 3. Aspen-specific model
# ~ Effect of aspen basal area on FRP and CBI for predominant forest types

da.aspen <- grid_tm %>%
 # Ensure that only relevant species are included
 filter(
  species_gp_n %in% c("aspen", "mixed_conifer", "lodgepole", 
                      "ponderosa", "spruce_fir", "piñon_juniper")) %>%
 filter(
  # filter to aspen grids
  aspen_pres == 1,
  # filter to at least 1% dominance or abundance
  sp_dominance_ld >= 0.01 | sp_abundance_ld >= 0.01
 ) %>%
 # Pivot the data
 pivot_wider(
  id_cols = c(grid_index, Fire_ID, log_frp_max_day, CBIbc_p90, fortypnm_gp,
              erc_dv, vpd_dv, elev, slope, tpi, chili),  # Keep these variables as-is
  names_from = species_gp_n,  # Pivot by species name
  values_from = balive,  # live basal area
  names_prefix = "balive_"  
 ) %>%
 # Replace NA with 0 for `balive` columns (no presence of species in the grid)
 mutate(across(starts_with("balive_"), ~replace_na(., 0)))
glimpse(da.aspen)


##########################
# Set up the model formula
mf.frp.aspen <- log_frp_max_day ~ 
 fortypnm_gp * balive_aspen + # dominant forest type * aspen live basal area
 erc_dv + vpd_dv + elev + slope + tpi + chili + # climate+topography
 f(Fire_ID, model = "iid") # Fire ID random effect

# fit the model                     
model_bl.frp.aspen <- inla(
 mf.frp.aspen, data = da.aspen, 
 family = "gaussian",
 control.predictor = list(compute = TRUE),
 control.compute = list(dic = TRUE, waic = TRUE)
)
summary(model_bl.frp.aspen)









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
mf4 <- log_frp_max_day ~ aspen + mixed_conifer + lodgepole + ponderosa + piñon_juniper + spruce_fir +
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


