
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
  # create an aspen presence flag based on abundance
  aspen_flag = if_else(species_gp_n == "aspen" & tpp_live > 0, 1, 0)
 ) %>%
 # be sure there are no duplicate rows for grid/fire/species
 distinct(grid_index, Fire_ID, species_gp_n, .keep_all = TRUE) # remove duplicates

# calculate the shannon index (H) for grids
# use both BA and TPP to calculate H
# also calculate an aspen presence column
shannon <- grid_tm %>%
 group_by(grid_index) %>%
 mutate(
  # Replace 0 or NA proportions with a small value to avoid log issues
  ba_live_pr = ifelse(is.na(ba_live_pr) | ba_live_pr == 0, 1e-6, ba_live_pr),
  tpp_live_pr = ifelse(is.na(tpp_live_pr) | tpp_live_pr == 0, 1e-6, tpp_live_pr),
 ) %>%
 # Calculate Shannon index components
 summarise(
  H_ba = -sum(ba_live_pr * log(ba_live_pr), na.rm = TRUE),  # Based on basal area proportions
  H_tpp = -sum(tpp_live_pr * log(tpp_live_pr), na.rm = TRUE),  # Based on trees per pixel proportions
  aspen = max(aspen_flag),
  .groups = "drop"
 )

# get the proportion aspen
aspen_pr <- grid_tm %>%
 group_by(grid_index) %>%
 summarise(
  # Total live BA and TPP in the grid
  total_ba_live = sum(ba_live, na.rm = TRUE),
  total_tpp_live = sum(tpp_live, na.rm = TRUE),
  # Aspen-specific live BA and TPP
  aspen_ba_live = sum(ba_live[species_gp_n == "aspen"], na.rm = TRUE),
  aspen_tpp_live = sum(tpp_live[species_gp_n == "aspen"], na.rm = TRUE),
  # Calculate proportions (avoid division by zero)
  aspen_ba_pr = if_else(total_ba_live > 0, aspen_ba_live / total_ba_live, 0),
  aspen_tpp_pr = if_else(total_tpp_live > 0, aspen_tpp_live / total_tpp_live, 0)
 ) %>%
 select(grid_index, aspen_ba_pr, aspen_tpp_pr) 

# merge back to the grid data
# center and scale the fixed effects
grid_tm <- grid_tm %>%
 left_join(shannon, by="grid_index") %>%
 left_join(aspen_pr, by = "grid_index") %>%
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

# tidy up!
rm(shannon, grid_counts, aspen_pr, idx)
gc()


#===========EXPLORE THE DATA==============#

######################################
# correlation matrix for fixed effects
# Select only numeric columns
cor_da <- grid_tm %>%
 select(
  tpp_live_pr, ba_live_pr, qmd_live_pr, # proportion live metrics
  ba_dead_pr, tpp_dead_pr, qmd_dead_pr, # proportion dead metrics
  aspen_ba_pr, aspen_tpp_pr, forest_pct,
  ba_live, tpp_live, qmd_live, tree_ht_live,  # Live forest composition metrics
  H_ba, H_tpp, # species diversity
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
levels(grid_tm$species_gp_n)

# filter where dominance | abundance is >= 5% of grid
# center and scale fixed effects
da <- grid_tm %>%
 filter(
  ba_live_pr >= 0.05
 ) %>%
 mutate(
  # center/scale metrics / fixed effects
  across(
   c(ba_live, ba_live_pr, ba_ld, ba_ld_pr,
     tpp_live, tpp_live_pr, tpp_ld, tpp_ld_pr,
     qmd_live, qmd_live_pr, qmd_ld, qmd_ld_pr,
     aspen_ba_pr, aspen_tpp_pr, H_ba, H_tpp, tree_ht_live, forest_pct,
     erc, erc_dv, vpd, vpd_dv, elev, slope, tpi, chili),
   ~ as.numeric(scale(.))
  ))


#===========MODEL FITTING==============#

###################################################################
# Model 2. Effect of species metrics and composition on FRP and CBI
# ~ FRP/CBIbc ~ species * (ba_live + qmd_live + shannon) + climate + topo

#####
# FRP

# setup the model formula
mf.frp <- log_frp_max_day ~ 1 + 
 fortypnm_gp:H_tpp + # species diversity by predominant forest type
 species_gp_n:ba_live_pr + # species TPP proportion
 species_gp_n:qmd_live + # species QMD
 forest_pct + # proportion of basal area that is aspen
 erc + vpd_dv + elev + slope + tpi + chili # climate + topography

# fit the model
model_tm.frp1 <- inla(
 mf.frp, data = da,
 family = "gaussian",
 control.predictor = list(compute=T),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
)
summary(model_tm.frp1)

################################
# Add random effects and compare

# update the model formula
mf.frp <- update(mf.frp, . ~ . + f(Fire_ID, model = "iid"))

# fit the model
model_tm.frp2 <- inla(
 mf.frp, data = da,
 family = "gaussian",
 control.predictor = list(compute=T),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
)
summary(model_tm.frp2)

# Compare the two models using WAIC
########################
# Compare DIC and WAIC #

cat("Baseline Model: \n")
cat("DIC:", model_tm.frp1$dic$dic, "\n")
cat("WAIC:", model_tm.frp1$waic$waic, "\n\n")

cat("With Fire_ID Random Effect: \n")
cat("DIC:", model_tm.frp2$dic$dic, "\n")
cat("WAIC:", model_tm.frp2$waic$waic, "\n")

print("Keeping better model")
if (model_tm.frp1$waic$waic > model_tm.frp2$waic$waic) {
 model_tm.frp = model_tm.frp2
} else {
 model_tm.frp = model_tm.frp1
}

rm(model_tm.frp1, model_tm.frp2) # clean up
gc()


#===========POSTERIOR EFFECTS=============#

########################
# Plotting fixed effects
# Extract fixed effect marginals
frp_marginals <- model_tm.frp$marginals.fixed
# Tidy up marginals for key predictors
tidy_marginals <- function(marginals, predictor) {
 tibble::tibble(
  parameter = names(marginals),
  data = purrr::map(marginals, ~ as.data.frame(.x))
 ) %>%
  unnest(data) %>%
  filter(str_detect(parameter, predictor))
}
# Tidy up species-specific effects
tidy_tpp <- tidy_marginals(frp_marginals, "ba_live_pr")
tidy_qmd <- tidy_marginals(frp_marginals, "qmd_live")
# Combine and add species labels
tidy_combined <- bind_rows(
 tidy_tpp %>% mutate(effect = "ba_live_pr"),
 tidy_qmd %>% mutate(effect = "qmd_live")
) %>%
 # Add species or forest type names
 mutate(
  species = case_when(
   str_detect(parameter, "species_gp_n") ~ str_extract(parameter, "species_gp_n\\w+"),
  )
 )

# Plot the ridge plot
ggplot(tidy_combined, aes(x = x, y = effect, height = y, fill = species)) +
 geom_density_ridges(stat = "identity", scale = 1.5, alpha = 0.7) +
 geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
 labs(
  x = "Effect Size",
  y = "Fixed Effect",
  fill = "Species"
 ) +
 theme_minimal()

rm(frp_marginals, tidy_marginals, tidy_tpp, tidy_qmd, tidy_combined)
gc()







#####
# CBI

# setup the model formula
mf.cbi <- CBIbc_p90 ~
 species_gp_n * (ba_live + qmd_live) + # Species and their structure/diversity
 erc + vpd + elev + slope + tpi + chili + # Climate and topography predictors
 f(Fire_ID, model = "iid") + # Fire-level random effect
 f(grid_index, model = "iid") # Grid-level random effect

# fit the model
model_bl_tm.cbi <- inla(
 mf.cbi, data = grid_tm,
 family = "gaussian",
 control.predictor = list(compute=T),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
)
summary(model_bl_tm.cbi)



# # center/scale abundance and dominance proportions
# ba_live_pr = ba_live_pr - mean(ba_live_pr, na.rm = TRUE),
# tpp_live_pr = tpp_live_pr - mean(tpp_live_pr, na.rm = TRUE),
# qmd_live_pr = qmd_live_pr - mean(qmd_live_pr, na.rm = TRUE),
# ba_dead_pr = ba_dead_pr - mean(ba_dead_pr, na.rm = TRUE),
# tpp_dead_pr = tpp_dead_pr - mean(tpp_dead_pr, na.rm = TRUE),
# qmd_dead_pr = qmd_dead_pr - mean(qmd_dead_pr, na.rm = TRUE)
# )
# # Prep a wide-format data frame
# grid_tm_w <- grid_tm %>%
#  # Pivot the data
#  pivot_wider(
#   id_cols = c(grid_index, Fire_ID, log_frp_max_day, CBIbc_p90, 
#               fortypnm_gp, H_ba, H_tpp,
#               erc, vpd, elev, slope, tpi, chili),  # Keep these variables as-is
#   names_from = species_gp_n,  # Pivot by species name
#   values_from = ba_live,  # live basal area
#   names_prefix = "balive_"
#  ) %>%
#  # Replace NA with 0 for `balive` columns (no presence of species in the grid)
#  mutate(across(starts_with("balive_"), ~replace_na(., 0)))
# glimpse(grid_tm_w)




# #################################################################
# # Expand the grid matrix for full presence/absence of all species
# # Get all unique species
# spps <- unique(grid_tm$species_gp_n)
# # Generate all combinations of grid_index and species
# complete_spps <- expand.grid(
#  grid_index = unique(grid_tm$grid_index),
#  species_gp_n = spps
# ) %>% as_tibble()
# 
# # isolate the grid-level variables for after
# grid.s <- grid_tm %>%
#  select(grid_index, Fire_ID, log_frp_max_day, CBIbc_p90, 
#         erc, vpd, elev, slope, tpi, chili, x, y) %>%
#  distinct()
# 
# # Create the model data frame and scale forest metrics
# grid_tm_c <- complete_spps %>%
#  left_join(grid_tm, by = c("grid_index", "species_gp_n")) %>%
#  mutate(
#   # Fill in missing values with 0
#   across(
#    c(ba_live, qmd_live, tpp_live, H_ba, H_tpp,
#      ba_live_pr, tpp_live_pr, qmd_live_pr,
#      ba_dead_pr, tpp_dead_pr, qmd_dead_pr
#    ), ~ replace_na(.x, 0))) %>%
#  arrange(grid_index) %>%
#  select(grid_index, species_gp_n, 
#         ba_live, tpp_live, qmd_live, H_ba, H_tpp,
#         ba_live_pr, tpp_live_pr, qmd_live_pr,
#         ba_dead_pr, tpp_dead_pr, qmd_dead_pr
#  ) %>%
#  mutate(
#   # create the species presence flag based on BA and TPP
#   presence = ifelse(ba_live == 0 & tpp_live == 0, 0, 1),
#   # center/scale metrics / fixed effects
#   across(
#    c(ba_live, tpp_live, qmd_live, H_ba, H_tpp),
#    ~ as.numeric(scale(.x, center = TRUE, scale = TRUE))
#   ),
#   # center/scale abundance and dominance proportions
#   ba_live_pr = ba_live_pr - mean(ba_live_pr, na.rm = TRUE),
#   tpp_live_pr = tpp_live_pr - mean(tpp_live_pr, na.rm = TRUE),
#   qmd_live_pr = qmd_live_pr - mean(qmd_live_pr, na.rm = TRUE),
#   ba_dead_pr = ba_dead_pr - mean(ba_dead_pr, na.rm = TRUE),
#   tpp_dead_pr = tpp_dead_pr - mean(tpp_dead_pr, na.rm = TRUE),
#   qmd_dead_pr = qmd_dead_pr - mean(qmd_dead_pr, na.rm = TRUE)
#  ) %>%
#  # merge back to get the response variables and grid-level variables
#  left_join(grid.s, by="grid_index")
# glimpse(grid_tm_c)
# 
# # Tidy up !
# rm(complete_spps, grid.s, grid_tm)
# gc()

# ########################################
# # Create the spatial fields model (SPDE)
# 
# # arrange by grid_index
# grid_tm_c <- grid_tm_c %>%
#  arrange(grid_index)
# 
# # Extract coordinates as a matrix
# # extract only unique rows
# coords_mat <- select(grid_tm_c, x, y) %>% as.matrix()
# 
# # Create a shared spatial mesh
# mesh <- inla.mesh.2d(
#  loc = coords_mat,                    # Locations (grid centroids)
#  max.edge = c(10, 30),             # Maximum edge lengths (inner and outer)
#  cutoff = 0.1,                      # Minimum distance between points
#  offset = c(1, 3)               # Boundary buffer
# )
# 
# # Compute the projector matrix using unique coordinates
# n.rep <- length(unique(grid_tm_c$species_gp_n)) # six species per grid
# m <- length(unique(grid_tm_c$grid_index))
# A <- inla.spde.make.A(
#  mesh = mesh, 
#  loc = coords_mat,
#  index = rep(1:m, times = n.rep),
#  repl = rep(1:n.rep, each = m)
# )
# 
# # Build the SPDE model
# spde.ml <- inla.spde2.pcmatern(
#  # Mesh and smoothness parameter
#  mesh = mesh, 
#  alpha = 2,
#  # P(practic.range < 0.3) = 0.5
#  prior.range = c(0.3, 0.5),
#  # P(sigma > 1) = 0.01
#  prior.sigma = c(10, 0.01)
# )
# 
# # Create the spatial index with replicates
# field.idx <- inla.spde.make.index(
#  name = "spatial_field",
#  mesh = mesh,
#  n.spde = spde.ml$n.spde,
#  n.repl = n.rep
# )
# 
# # Validate the length of the spatial field index
# n_spde <- spde.ml$n.spde
# length(field.idx$spatial_field) == n_spde * n.rep
# 
# # tidy up !
# rm(coords_mat, mesh)
# gc()

# # Define the response variable
# response <- list(log_frp_max_day = grid_tm_c$log_frp_max_day)
# # Define the fixed effects (covariates)
# X <- grid_tm_c[, c("species_gp_n", "Fire_ID", "grid_index",
#                    "ba_live", "qmd_live", "H_tpp",
#                    "erc", "vpd", "elev", "slope", "tpi", "chili")]
# 
# # Create the INLA stack
# stk.frp <- inla.stack(
#  data = response,
#  A = list(A, 1),
#  effects = list(field.idx, X),
#  tag = 'est'
# )


# #########################

# 
# 
# ###########################
# # Extract interaction terms
# effects <- as.data.frame(model_bl_tm.frp$summary.fixed) %>%
#  rownames_to_column(var = "parameter") %>%
#  filter(str_detect(parameter, ":")) %>%
#  separate(parameter, into = c("species1", "species2"), sep = ":") %>%
#  mutate(across(everything(), ~ str_replace(., "balive_", "")))
# 
# # Prepare data for heatmap
# heatmap_da <- effects %>%
#  select(species1, species2, mean) %>%
#  pivot_wider(names_from = species2, values_from = mean) %>%
#  column_to_rownames(var = "species1")
# 
# # Convert to long format for ggplot
# heatmap_l <- melt(heatmap_da, varnames = c("Species1", "Species2"), value.name = "Effect")
# 
# # Plot heatmap
# ggplot(heatmap_l, aes(x = Species2, y = Species1, fill = Effect)) +
#  geom_tile(color = "white") +
#  scale_fill_gradient2(
#   low = "blue", mid = "white", high = "red", midpoint = 0,
#   limits = c(min(heatmap_l$Effect, na.rm = TRUE), max(heatmap_l$Effect, na.rm = TRUE))
#  ) +
#  labs(
#   title = "Pairwise Interaction Effects on FRP",
#   x = "Co-occurring Species",
#   y = "Dominant Species",
#   fill = "Effect"
#  ) +
#  theme_minimal() +
#  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# 
# ##########################################################
# 
# 
# 
# 
# 
# 
# ##########################################################
# # Extract fixed effect marginals for the interaction terms
# 
# # mf.frp <- log_frp_max_day ~ 
# #  fortypnm_gp + # predominant forest type
# #  (species_gp_n * balive) + # species dominance (balive)
# #  erc_dv + vpd_dv + elev + slope + tpi + chili + # climate + topo
# #  f(Fire_ID, model = "iid") # random effect for Fire ID
# 
# marginals <- model_bl_tm.frp$marginals.fixed
# 
# # Define a function to tidy marginals for interaction terms
# tidy_int <- function(marginals) {
#  tibble::tibble(
#   parameter = names(marginals),
#   data = purrr::map(marginals, ~ as.data.frame(.x))
#  ) %>%
#   unnest(data) %>%
#   filter(str_detect(parameter, "species_gp_n") & str_detect(parameter, ":balive")) %>%
#   mutate(
#    species = str_remove(parameter, ":balive"),
#    species = str_remove(species, "species_gp_n"),
#    species = str_replace(species, "_", " "),
#    effect_type = "Interaction"
#   )
# }
# 
# # Tidy the marginals for the interaction terms
# tidy_int_df <- tidy_int(marginals)
# 
# # Create ridge plot for interaction terms
# ggplot(tidy_int_df, aes(x = x, y = species, height = y, fill = effect_type)) +
#  geom_density_ridges(stat = "identity", scale = 1.5, alpha = 0.7) +
#  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
#  labs(
#   x = "Effect Magnitude",
#   y = "Co-occurring Species",
#   fill = "Effect Type"
#  ) +
#  scale_fill_manual(values = c("Interaction" = "#FEB24C")) +
#  theme_classic() +
#  theme(
#   axis.text.y = element_text(size = 10),
#   axis.title.y = element_text(size = 11, margin = margin(r = 10)),
#   axis.title.x = element_text(size = 11, margin = margin(t = 10)),
#   legend.position = "top"
#  )
# 
# # save the plot.
# out_png <- paste0(maindir,'figures/INLA_Species-BALIVE_PosteriorEffects_Ridge.png')
# ggsave(out_png, dpi=500, bg = 'white')
# 
# 
# 
# ##############################################################
# # 3. Aspen-specific model
# # ~ Effect of aspen basal area on FRP and CBI for predominant forest types
# 
# da.aspen <- grid_tm %>%
#  # Ensure that only relevant species are included
#  filter(
#   species_gp_n %in% c("aspen", "mixed_conifer", "lodgepole", 
#                       "ponderosa", "spruce_fir", "piñon_juniper")) %>%
#  filter(
#   # filter to aspen grids
#   aspen_pres == 1,
#   # filter to at least 1% dominance or abundance
#   sp_dominance_ld >= 0.01 | sp_abundance_ld >= 0.01
#  ) %>%
#  # Pivot the data
#  pivot_wider(
#   id_cols = c(grid_index, Fire_ID, log_frp_max_day, CBIbc_p90, fortypnm_gp,
#               erc_dv, vpd_dv, elev, slope, tpi, chili),  # Keep these variables as-is
#   names_from = species_gp_n,  # Pivot by species name
#   values_from = balive,  # live basal area
#   names_prefix = "balive_"  
#  ) %>%
#  # Replace NA with 0 for `balive` columns (no presence of species in the grid)
#  mutate(across(starts_with("balive_"), ~replace_na(., 0)))
# glimpse(da.aspen)
# 
# 
# ##########################
# # Set up the model formula
# mf.frp.aspen <- log_frp_max_day ~ 
#  fortypnm_gp * balive_aspen + # dominant forest type * aspen live basal area
#  erc_dv + vpd_dv + elev + slope + tpi + chili + # climate+topography
#  f(Fire_ID, model = "iid") # Fire ID random effect
# 
# # fit the model                     
# model_bl.frp.aspen <- inla(
#  mf.frp.aspen, data = da.aspen, 
#  family = "gaussian",
#  control.predictor = list(compute = TRUE),
#  control.compute = list(dic = TRUE, waic = TRUE)
# )
# summary(model_bl.frp.aspen)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ########################
# # Compare DIC and WAIC #
# 
# cat("Baseline Model: \n")
# cat("DIC:", model_bl1$dic$dic, "\n")
# cat("WAIC:", model_bl1$waic$waic, "\n\n")
# 
# cat("With Fire_ID Random Effect: \n")
# cat("DIC:", model_bl2$dic$dic, "\n")
# cat("WAIC:", model_bl2$waic$waic, "\n")
# 
# print("Keeping better model")
# if (model_bl1$waic$waic > model_bl2$waic$waic) {
#  rm(model_bl1) # clean up
#  gc()
# } else {
#  rm(model_bl2) # clean up
#  gc()
# }
# 
# ##############################################################
# # 1. Baseline model without spatial component or random effect
# 
# 
# # set the formula.
# mf <- log_frp_max_day ~ vpd_dv + erc_dv + elev + slope + tpi + chili +  # Climate/topography
#  # Latent fields for species composition
#  f(cidx, spp_pct, model = "generic0", Cmatrix = sp_cov_grid, 
#    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))
# # fit the model                     
# model_bl1 <- inla(
#  mf, data = spp_effect, 
#  family = "gaussian",
#  control.predictor = list(compute = TRUE),
#  control.compute = list(dic = TRUE, waic = TRUE)
# )
# summary(model_bl1)
# 
# 
# ###########################################################
# # 2. Adding between fires effect (random effect on fire ID)
# 
# # set the formula.
# mf2 <- log_frp_max_day ~ vpd_dv + erc_dv + elev + slope + tpi + chili +  # Climate/topography
#  f(cidx, spp_pct, model = "generic0", Cmatrix = sp_cov_grid, 
#    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
#  f(Fire_ID_nm, model = "iid") # Fire ID random effect
# 
# # fit the model                     
# model_bl2 <- inla(
#  mf2, data = spp_effect, 
#  family = "gaussian",
#  control.predictor = list(compute = TRUE),
#  control.compute = list(dic = TRUE, waic = TRUE)
# )
# summary(model_bl2)
# 
# 
# ########################
# # Compare DIC and WAIC #
# 
# cat("Baseline Model: \n")
# cat("DIC:", model_bl1$dic$dic, "\n")
# cat("WAIC:", model_bl1$waic$waic, "\n\n")
# 
# cat("With Fire_ID Random Effect: \n")
# cat("DIC:", model_bl2$dic$dic, "\n")
# cat("WAIC:", model_bl2$waic$waic, "\n")
# 
# print("Keeping better model")
# if (model_bl1$waic$waic > model_bl2$waic$waic) {
#  rm(model_bl1) # clean up
#  gc()
# } else {
#  rm(model_bl2) # clean up
#  gc()
# }
# 
# 
# ############################################
# # 3. Adding "within-fire" spatial dependence
# 
# # Extract coordinates from wide data frame
# coords <- spp_effect %>% distinct(grid_index, grid_x, grid_y)
# coords_mat <- as.matrix(coords[, c("grid_x", "grid_y")])
# # Create a shared spatial mesh
# mesh <- inla.mesh.2d(
#  loc = coords_mat,
#  max.edge = c(10, 100),  
#  cutoff = 0.01 # Minimum distance between points
# )
# plot(mesh)
# 
# # define the stochastic partial difference equation (SPDE)
# spde <- inla.spde2.pcmatern(
#  mesh = mesh,
#  alpha = 2,  # Smoothness parameter
#  prior.range = c(10, 0.01),  # Prior for spatial range
#  prior.sigma = c(5, 0.01)    # Prior for variance
# )
# 
# # create the A-matrix (linking mesh to coords)
# A <- inla.spde.make.A(
#  mesh = mesh,
#  loc = coords_mat
# )
# 
# # Create the INLA stack
# stack <- inla.stack(
#  data = list(log_frp_max_day = spp_effect$log_frp_max_day),  # Response variable
#  A = list(A, diag(nrow(spp_effect))),  # Sparse spatial field and identity matrix
#  effects = list(
#   spatial_field = 1:spde$n.spde,  # Spatial random field
#   spp_effect %>% select(-log_frp_max_day)  # All covariates except response
#  )
# )
# 
# #####################
# # define the formula.
# mf3 <- log_frp_max_day ~ vpd_dv + erc_dv + elev + slope + tpi + chili +  # Climate/topography
#  f(cidx, spp_pct, model = "generic0", Cmatrix = sp_cov_grid, 
#    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
#  f(Fire_ID_nm, model = "iid") + # Fire ID random effect
#  f(spatial_field, model = spde) # spatial effect
# 
# model_bl3 <- inla(
#  formula = mf3,
#  family = "gaussian",
#  data = inla.stack.data(stack),
#  control.predictor = list(A = inla.stack.A(stack)),
#  control.compute = list(dic = TRUE, waic = TRUE)
# )
# 
# # Summarize the model results
# summary(model_bl3)
# 
# 
# ##############################################
# # 4. Adding temporal random effect (fire year)
# 
# # define the formula.
# mf4 <- log_frp_max_day ~ aspen + mixed_conifer + lodgepole + ponderosa + piñon_juniper + spruce_fir +
#  vpd_dv + erc_dv + slope + tpi + chili +
#  f(Fire_ID_nm, model = "iid") + # fire random effect
#  f(spatial_field, model = spde) + # spatial effect
#  f(year_season, model = "rw1")  # temporal random effect
# 
# model_bl4 <- inla(
#  formula = mf4,
#  family = "gaussian",
#  data = inla.stack.data(stack),
#  control.predictor = list(A = inla.stack.A(stack)),
#  control.compute = list(dic = TRUE, waic = TRUE)
# )
# 
# # Summarize the model results
# summary(model_bl4)
# 
# 
# ######################
# # Compare DIC and WAIC
# 
# cat("Spatial model: \n")
# cat("DIC:", model_bl3$dic$dic, "\n")
# cat("WAIC:", model_bl3$waic$waic, "\n\n")
# 
# cat("Spatial-temporal model (year): \n")
# cat("DIC:", model_bl4$dic$dic, "\n")
# cat("WAIC:", model_bl4$waic$waic, "\n")
# 
# # rm(model_bl1)
# # gc()
# 
# 
# #===========Plotting==============#
# 
# # Extract fixed effects
# fixed_effects <- as.data.frame(model_bl4$summary.fixed)
# fixed_effects$Variable <- rownames(fixed_effects)
# 
# # Exponentiate coefficients for interpretation
# fixed_effects$mean_exp <- exp(fixed_effects$mean)
# fixed_effects$lower_exp <- exp(fixed_effects$`0.025quant`)
# fixed_effects$upper_exp <- exp(fixed_effects$`0.975quant`)
# 
# # Plot exponentiated coefficients (excluding intercept for clarity)
# ggplot(fixed_effects %>% filter(Variable != "(Intercept)"), 
#        aes(x = Variable, y = mean_exp, ymin = lower_exp, ymax = upper_exp)) +
#  geom_pointrange() +
#  coord_flip() +
#  labs(y = "Effect on FRP", x = "Variable") +
#  theme_minimal()
# 
# 
# # Extract spatial field posterior mean
# spatial_effects <- inla.spde.make.index("spatial_field", n.spde = spde$n.spde)
# spatial_field_mean <- model_bl4$summary.random$spatial_field$mean
# 
# # Add to mesh
# mesh_df <- data.frame(
#  x = mesh$loc[, 1],
#  y = mesh$loc[, 2],
#  spatial_effect = spatial_field_mean
# )
# ggplot(mesh_df, aes(x = x, y = y, color = spatial_effect)) +
#  geom_point(size = 1) +
#  scale_color_viridis_c() +
#  coord_equal() +
#  labs(title = "Spatial Field (Posterior Mean)", x = "Longitude", y = "Latitude") +
#  theme_minimal()
# 
# 
# # Extract year random effect
# year_effects <- model_bl4$summary.random$year_season
# 
# ggplot(year_effects, aes(x = ID, y = mean, ymin = `0.025quant`, ymax = `0.975quant`)) +
#  geom_pointrange() +
#  geom_line() +
#  labs(title = "Year-Season Effect (Posterior Mean)", x = "Year-Season", y = "Effect Size") +
#  theme_minimal() +
#  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# 
# 
# #===========Model Setup (Interactions)==============#