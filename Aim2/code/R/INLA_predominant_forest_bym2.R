
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
grid_tm <-  read_csv(fp) %>% # read in the file
 # get the acquisition year and month
 mutate(day_max_frp = as.Date(day_max_frp),  # Convert to Date
        fire_year = year(day_max_frp),              # Extract year
        fire_month = month(day_max_frp),
        # create a new unique identifier
        grid_idx = paste0(Fire_ID, grid_index)) %>%
 # remove missing FRP, prep columns
 filter(
  frp_max_day > 0, # only work with daytime FRP
  CBIbc_p90 > 0 # and CBI 0 -> these are often boundary grids
 ) %>% 
 # create a numeric fire ID
 mutate(
  # Format species names consistently
  fortypnm_gp = str_to_lower(fortypnm_gp),
  fortypnm_gp = str_replace_all(fortypnm_gp, "-", "_"),
  fortypnm_gp = str_replace_all(fortypnm_gp, " ", "_"),
  # make sure factor variables are set for species names
  fortypnm_gp = factor(fortypnm_gp),
  species_gp_n = factor(species_gp_n),
  Fire_ID = factor(Fire_ID),
  grid_idx = as.numeric(as.factor(grid_idx)),
  # proportion of contributing AFD during daytime observations
  day_prop = day_count / afd_count,
  # log-scaled outcomes
  log_frp_max = log(frp_max + 1e-5),
  log_frp_csum = log(frp_csum + 1e-5),
  log_frp_max_day = log(frp_max_day + 1e-5),
  log_frp_max_night = log(frp_max_night + 1e-5),
  log_frp_csum_day = log(frp_csum_day + 1e-5),
  log_frp_csum_night = log(frp_csum_night + 1e-5),
  # Create a "season" variable based on the month
  season = case_when(
   fire_month %in% c(3, 4, 5) ~ "spring",
   fire_month %in% c(6, 7, 8) ~ "summer",
   fire_month %in% c(9, 10, 11) ~ "fall"
  ),
  # Year/season interaction
  year_season = interaction(fire_year, season, drop = TRUE),
  # Factorize the temporal attributes
  season = factor(season),
  fire_year = factor(fire_year),
  year_season = factor(year_season),
  # scale the percentages
  fortyp_pct = fortyp_pct / 100,
  forest_pct = forest_pct / 100,
  # create a species presence flag
  presence = factor(if_else(ba_live > 0 | tpp_live > 0, 1, 0))) %>%
 # be sure there are no duplicate rows for grid/fire/species
 distinct(grid_idx, fortypnm_gp, .keep_all = TRUE) # remove duplicates
glimpse(grid_tm) # check the results

# Check how many grids are > 50% forested
print("Number of forested grids (>50% forest pixels):")
dim(grid_tm %>% filter(forest_pct > 0.50))[1]/dim(grid_tm)[1]


#===============Explore Distributions, etc.================#

####################################
# distribution of response variables
grid_tm %>%
 # pivot longer to facet plot
 pivot_longer(cols = c(log_frp_max_day, CBIbc_p90),
              names_to = "variable",
              values_to = "value") %>%
 # Plot with facets
 ggplot(aes(x = value)) +
 geom_histogram(bins = 30, fill = "orange", alpha = 0.7) +
 facet_wrap(
  ~ variable, 
  scales = "free",
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
# this correlation matrix is for this simple model (i.e., not forest composition)
cor_da <- grid_tm %>%
 select(
  fortypnm_gp, # forest type (factor)
  forest_pct, fortyp_pct, # forest and forest type percent
  canopypct_mean, balive_sum, # grid-level mean canopy percent and balive
  erc, erc_dv, vpd, vpd_dv, # climate
  fm1000, fm1000_dv, rmin, rmin_dv, # climate
  tmmx, tmmx_dv, vs, vs_dv, #climate
  elev, slope, tpi, chili  # topography
 ) %>%
 pivot_wider(
  names_from = fortypnm_gp, 
  values_from = fortyp_pct, 
  values_fill = 0) %>%
 mutate(across(everything(), ~ scale(.) %>% as.numeric()))  # Standardize variables

# Compute correlation matrix
cor_mat <- cor(cor_da, use = "complete.obs", method = "spearman")

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
out_png <- paste0(maindir,'figures/INLA_CorrelationMatrix_Fortyp.png')
ggsave(out_png, dpi=500, bg = 'white')



#===========MODEL SETUP==============#

# list of species names
spps <- c("quaking_aspen", "mixed_conifer", "lodgepole", 
          "ponderosa", "spruce_fir", "piñon_juniper",
          "rocky_mountain_juniper", "oak_woodland")

# force aspen to be the baseline
grid_tm <- grid_tm %>%
 mutate(
  fortypnm_gp = fct_relevel(fortypnm_gp, spps)
 )

# check the factor levels
# make sure aspen is first
levels(grid_tm$fortypnm_gp)

# create a data frame for just the predominant forest type
# one row per grid_index with the dominant forest type
# filter to where that type is at least 50% of the grid forested area
da <- grid_tm %>%
 # keep only grids where predominant species cover > 50%
 filter(fortyp_pct > 0.50) %>%
 # select the columns we need for modeling
 select(grid_index, grid_idx, Fire_ID, # ID columns
        fire_year, season, year_season, # temporal effects
        log_frp_max_day, log_frp_csum, log_frp_csum_day, # FRP response variables
        CBIbc_p90, CBIbc_p95, CBIbc_p99, CBIbc_mean, # CBI response variables
        fortypnm_gp, fortyp_pct, forest_pct, # forest type and percent cover
        canopypct_mean, balive_sum, # canopy percent and total live basal area
        erc, erc_dv, vpd, vpd_dv, # climate
        fm1000, fm1000_dv, rmin, rmin_dv, # climate
        tmmx, tmmx_dv, vs, vs_dv, #climate
        elev, slope, tpi, chili, # topography
        presence, # species presence (binary factor)
        x, y # grid centroid coordinate for spatial fields model
 ) %>%
 # center and scale continuous predictor variables
 mutate(
  across(
   c(
    forest_pct, fortyp_pct, # forest and forest type percent
    canopypct_mean, balive_sum, # canopy percent and total live basal area
    erc, vpd, vpd_dv, erc_dv, # climate
    fm1000, fm1000_dv, rmin, rmin_dv, # climate
    tmmx, tmmx_dv, vs, vs_dv, #climate
    elev, slope, tpi, chili # topography
   ),
   ~ as.numeric(scale(.))
  )
 ) %>%
 # rename the response variables
 rename(
  frp_day = log_frp_max_day,
  frp_csum = log_frp_csum,
  frp_csum_day = log_frp_csum_day,
  cbi_p90 = CBIbc_p90,
  cbi_p95 = CBIbc_p95,
  cbi_p99 = CBIbc_p99,
  cbi_mn = CBIbc_mean,
  fire_id = Fire_ID
 ) %>%
 # keep just one row per grid cell by predominant type
 distinct(grid_idx, fortypnm_gp, .keep_all = TRUE) %>%
 arrange(grid_idx)


########################################################
# Check on the grid cell counts for daytime observations
# After filtering, make sure we have enough grids per fire
grid_counts <- da %>%
 distinct(fire_id, grid_idx) %>% # keep only distinct rows
 group_by(fire_id) %>%
 summarise(n_grids = n())
summary(grid_counts$n_grids)
quantile(grid_counts$n_grids, probs = seq(.1, .9, by = .1))

# Identify fires with n_grids below the 20th percentile
idx <- grid_counts %>%
 filter(n_grids < quantile(grid_counts$n_grids, probs = 0.2)) %>%
 pull(fire_id)
length(idx)

# filter the data frame to remove these fires
da <- da %>%
 filter(!fire_id %in% idx)
# tidy up
rm(grid_counts, idx)
gc()

# check how many grids and fires
length(unique(da$grid_idx))
length(unique(da$fire_id))

# check how many grids of each forest type there are
da %>%
 group_by(fortypnm_gp) %>%
 summarize(n = n()) %>%
 ungroup()

# check on counts by season
da %>%
 group_by(season) %>%
 summarize(n = n()) %>%
 ungroup()

# Subset species with too few observations
drop_spps = c("rocky_mountain_juniper", "oak_woodland")
print("Dropping small classes:")
print(dim(da%>%filter(fortypnm_gp%in%drop_spps)))
da <- da %>%
 filter(!fortypnm_gp %in% drop_spps) %>%
 mutate(fortypnm_gp = droplevels(fortypnm_gp))

# Tidy up!
rm(grid_tm)
gc()



#===========MODEL FITTING==============#

set.seed(435)

#######
# FRP #

#######################################
# 1. Baseline model (no latent effects)

# Set up the model formula
mf.frp <- frp_day ~ 
 fortypnm_gp + # predominant forest type
 forest_pct + # percent of grid of predominant forest type
 erc_dv + # energy release component (ERC) and wind speed
 slope + tpi + chili # topography

# fit the model                     
model_bl.frp <- inla(
 mf.frp, data = da, 
 family = "gaussian",
 control.predictor = list(compute = TRUE),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
 control.fixed = list(
  prec = list(prior = "pc.prec", param = c(1, 0.01)) 
 ) 
)
summary(model_bl.frp)

 
##############################################
# 2. Baseline model + fire-level random effect

# specify a pc-prior for the random effect
hyper_pr <- list(
 prec = list(prior = "pc.prec", param = c(1, 0.1))
)
# update the model formula
mf.frp.re <- update(
 mf.frp, . ~ 1 + . + 
  f(fire_id, model = "iid", hyper = hyper_pr) # fire-level random effect
)

# fit the model                     
model_bl.frp.re <- inla(
 mf.frp.re, data = da, 
 family = "gaussian",
 control.predictor = list(compute = TRUE),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
 control.fixed = list(
  prec = list(prior = "pc.prec", param = c(1, 0.1)) 
 )   
)
summary(model_bl.frp.re)


#######################################################################
# 3. Baseline model + fire-level random effect + temporal random effect

# specify a pc-prior for the random effect
hyper_pr <- list(
 prec = list(prior = "pc.prec", param = c(1, 0.01))
)
# update the model formula
mf.frp.re2 <- update(
 mf.frp.re, . ~ 1 + . + 
  f(year_season, model = "iid", hyper = hyper_pr) # temporal random effect
)

# fit the model                     
model_bl.frp.re2 <- inla(
 mf.frp.re2, data = da, 
 family = "gaussian",
 control.predictor = list(compute = TRUE),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
 control.fixed = list(
  prec = list(prior = "pc.prec", param = c(1, 0.01)) 
 )   
)
summary(model_bl.frp.re2)


##################
# Model comparison

model_comparison <- tibble(
 Model = c("Baseline", "Baseline + Fire RE", "Baseline + Fire + Temporal RE"),
 DIC = c(
  model_bl.frp$dic$dic,
  model_bl.frp.re$dic$dic,
  model_bl.frp.re2$dic$dic
 ),
 WAIC = c(
  model_bl.frp$waic$waic,
  model_bl.frp.re$waic$waic,
  model_bl.frp.re2$waic$waic
 ),
 Marginal_LogLikelihood = c(
  model_bl.frp$mlik[1],
  model_bl.frp.re$mlik[1],
  model_bl.frp.re2$mlik[1]
 ),
 Effective_Params = c(
  model_bl.frp$dic$p.eff,
  model_bl.frp.re$dic$p.eff,
  model_bl.frp.re2$dic$p.eff
 )
) %>% arrange(WAIC)
 
# Print the comparison table
print(model_comparison)

# # Tidy up !
# rm(model_bl.frp, model_bl.frp.re, model_bl.frp.re2)
# gc()


###########################################
# 4. Baseline model + spatial random effect

# spatial points representing grid centroids
da$grid_idx <- as.numeric(as.factor(da$grid_idx))
grid_sf <- da %>%
 distinct(grid_idx, x, y, .keep_all = TRUE) %>%
 st_as_sf(., coords = c("x", "y"), crs = 4326)
coords <- grid_sf %>% st_coordinates(.)

# k-nearest neighbor and symmetric matrix
nbs <- knearneigh(coords, k = 3, longlat = TRUE)
nbs <- knn2nb(nbs, row.names = grid_sf$grid_idx, sym = TRUE) #force symmetry!!
# plot
fire_sf <- grid_sf %>% filter(fire_id == 51)
plot(st_geometry(fire_sf), col = "grey")
plot(nbs, coords, add = TRUE, col = "red", lwd = 0.5)
title("KNN Adjacency Structure")

# setup for INLA
nb2INLA("cl_graph",nbs)
am_adj <-paste(getwd(),"data/spatial/cl_graph",sep="")
H <- inla.read.graph(filename="cl_graph")


##########################
# update the model formula
mf.frp.sp <- update(
 mf.frp.re, . ~ 1 + . + 
  f(grid_idx, 
    model = "bym2", 
    graph = H, 
    hyper = list(
     phi = list(prior = "pc", param = c(0.5, 2/3)),  # Proportion of spatial effect
     prec = list(prior = "pc.prec", param = c(1, 0.1))  # Overall variance
    )
   )
)

# fit the model
model_bl.frp.sp <- inla(
 mf.frp.sp, 
 data = da, 
 family = "gaussian",
 control.predictor = list(compute = TRUE),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
 control.fixed = list(
  prec = list(prior = "pc.prec", param = c(1, 0.5)) 
 )  
)

summary(model_bl.frp.sp)



#######################################
#######
# CBIbc

# handle 0s in the data
# we looked at these ... many are boundary grids
da.cbi <- da %>%
 filter(cbi_p90 > 0)

# examine the distributions of CBIbc variables
da.cbi %>%
 mutate(log_cbi_mn = log(cbi_mn)) %>%
 # pivot longer to facet plot
 pivot_longer(cols = c(cbi_p90, cbi_p95, cbi_mn, log_cbi_mn),
              names_to = "variable",
              values_to = "value") %>%
 # Plot with facets
 ggplot(aes(x = value)) +
 geom_histogram(bins = 30, fill = "orange", alpha = 0.7) +
 facet_wrap(
  ~ variable, 
  scales = "free",
  labeller = as_labeller(c(cbi_mn = "Average CBIbc",
                           log_cbi_mn = "log(Average CBIbc)",
                           cbi_p90 = "90th Percentile CBIbc",
                           cbi_p95 = "95th Percentile CBIbc"))) +
 labs(x = "value",
      y = "Frequency") +
 theme_minimal()

# save the plot.
out_png <- paste0(maindir,'figures/INLA_ResponseDistribution.png')
ggsave(out_png, dpi=500, bg = 'white')



#######################################
# 1. Baseline model (no latent effects)

# Set up the model formula
mf.cbi <- cbi_p90 ~ 
 fortypnm_gp + # predominant forest type
 forest_pct + # percent of grid of predominant forest type
 erc_dv + # energy release component (ERC)
 slope + tpi + chili # topography

# fit the model                     
model_bl.cbi <- inla(
 mf.cbi, data = da.cbi, 
 family = "gamma", 
 control.predictor = list(compute = TRUE),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
 control.fixed = list(prec = list(prior = "pc.prec", param = c(1, 0.1)))
)
summary(model_bl.cbi)


##############################################
# 2. Baseline model + fire-level random effect

# update the model formula
hyper_pr <- list(prec = list(prior = "pc.prec", param = c(1, 0.5)))
mf.cbi.re <- update(
 mf.cbi, . ~ 1 + . + 
  f(fire_id, model = "iid", hyper = hyper_pr) # fire-level random effect
)

# fit the model                     
model_bl.cbi.re <- inla(
 mf.cbi.re, data = da.cbi, 
 family = "gamma",
 control.predictor = list(compute = TRUE),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
 control.fixed = list(
  prec = list(prior = "pc.prec", param = c(1, 0.1)) 
 )   
)
summary(model_bl.cbi.re)


###########################################
# 3. Baseline model + spatial random effect

##########################
# update the model formula
mf.cbi.sp <- update(
 mf.cbi, . ~ 1 + . + 
  f(grid_index, 
    model = "bym2", 
    graph = H, 
    hyper = list(
     phi = list(prior = "pc", param = c(0.5, 0.5)),  # Proportion of spatial effect
     prec = list(prior = "pc.prec", param = c(1, 0.01))  # Overall variance
    ), 
    constr = T
  )
)

# fit the model
model_bl.cbi.sp <- inla(
 mf.cbi.sp, 
 data = da.cbi, 
 family = "gamma",
 control.predictor = list(compute = TRUE),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
 control.fixed = list(
  prec = list(prior = "pc.prec", param = c(1, 0.1)) 
 )  
)

summary(model_bl.cbi.sp)


###########################################
# 3. Baseline model + spatial random effect

da.cbi <- da.cbi %>%
 mutate(fire_id = as.numeric(as.factor(fire_id)))
##########################
# update the model formula
mf.cbi.sp.re <- update(
 mf.cbi, . ~ 1 + . + 
  f(fire_id)
)

# fit the model
model_bl.cbi.sp.re <- inla(
 mf.cbi.sp.re, 
 data = da.cbi, 
 family = "gamma",
 control.predictor = list(compute = TRUE),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
 control.fixed = list(
  prec = list(prior = "pc.prec", param = c(1, 0.1)) 
 )  
)

summary(model_bl.cbi.sp.re)


# Tidy up !
rm(coords, fire_sf, H, hyper_pr, model_bl.cbi, 
   model_bl.cbi.re, model_bl.frp, model_bl.frp.re, nbs)
gc()



#===========Plotting Posterior Effects============#

# Extract fixed effects related to fortypnm_gp

#####
# FRP
frp.eff <- as.data.frame(model_bl.frp.sp$summary.fixed) %>%
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
cbi.eff <- as.data.frame(model_bl.cbi.sp$summary.fixed) %>%
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
frp_marginals <- model_bl.frp.sp$marginals.fixed
cbi_marginals <- model_bl.cbi.sp$marginals.fixed

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

# pull out the climate + topo
climate_topo <- tidy_combined %>%
 filter(parameter %in% c("erc_dv", "slope", "tpi", "chili")) %>%
 mutate(parameter = factor(parameter, levels = c("erc_dv", "slope", "tpi", "chili")))

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
                    "CBI" = "90th Percentile CBIbc", 
                    "FRP" = "Maximum Daytime FRP")) +
 # Apply custom colors
 scale_color_manual(
  values = cmap,
  labels = c(
   "CBIbc" = "90th Percentile CBIbc", 
   "FRP" = "Maximum Daytime FRP")
 ) +
 theme_classic() +
 theme(axis.text.y = element_text(angle = 0, hjust = 0, size=11),
       axis.text.x = element_text(angle = 0, hjust = 0, size=11),
       axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
       axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
       legend.position = c(0.18, 0.16),
       legend.background = element_rect(
        fill = scales::alpha("white", 0.4), 
        color = NA, size = 0.8),
       legend.title = element_text(size = 11),
       legend.text = element_text(size = 10))
p2

# save the plot.
out_png <- paste0(maindir,'figures/INLA_FORTYPNM_PosteriorEffects_Ridge.png')
ggsave(out_png, dpi=300, bg = 'white')


# Plot the climate + topo effects
p_climate_topo <- ggplot(climate_topo, aes(x = x, y = parameter, height = y, fill = response)) +
 geom_density_ridges(stat = "identity", scale = 1.5, alpha = 0.7) +
 geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
 labs(
  x = "Effect size",
  y = "Climate/Topography Predictor",
  fill = "Response"
 ) +
 scale_fill_manual(
  values = c("FRP" = "#FEB24C", "CBI" = "#800026"),
  labels = c("FRP" = "Maximum Daytime FRP", "CBIbc" = "90th Percentile CBIbc")
 ) +
 theme_classic() +
 theme(
  axis.text.y = element_text(size = 11),
  axis.text.x = element_text(size = 11),
  axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
  axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
  legend.position = "top",
  legend.title = element_text(size = 11),
  legend.text = element_text(size = 10)
 )
p_climate_topo

# Save the plot
out_png <- paste0(maindir, 'figures/INLA_ClimateTopo_PosteriorEffects_Ridge.png')
ggsave(out_png, plot = p_climate_topo, dpi = 300, bg = 'white')



######################################
# Visualize the spatial random effects

# Extract spatial random effect estimates
spatial_effects <- model_bl.frp.sp$summary.random$grid_index %>%
 dplyr::select(ID, mean, sd, '0.025quant', '0.5quant', '0.975quant') %>%
 rename(grid_index = ID)

# Join spatial effects back to spatial data
grid_sf.eff <- grid_sf %>%
 left_join(spatial_effects, by = "grid_index")
sum(is.na(grid_sf.eff$mean))  # Count missing values

# check the histogram of spatial effects across fires
hist(grid_sf.eff$mean, breaks = 50, col = "skyblue")

# Plot the spatial effects (posterior mean)
# Use the Williams Fork Fire as an example
fire_sf <- grid_sf.eff %>% filter(fire_id == 51) %>% st_transform(st_crs(5070))
ggplot(fire_sf) +
 geom_sf(aes(color = mean), size=2.5) +
 scale_color_viridis_c(option = "magma") +
 theme_minimal()

# Save the plot
out_png <- paste0(maindir, 'figures/INLA_WilliamsFork_SpatialEff.png')
ggsave(out_png, dpi = 300, bg = 'white')


###########################################################
# Extract structured and unstructured components separately
n <- nrow(grid_sf.eff)
st_effects <- model_bl.frp.sp$summary.random$grid_index$mean[1:n] * 
 model_bl.frp.sp$summary.hyperpar$mean["Phi for grid_index"]
un_effects <- model_bl.frp.sp$summary.random$grid_index$mean[1:n] * 
 (1 - model_bl.frp.sp$summary.hyperpar$mean["Phi for grid_index"])

# Add components to spatial data
grid_sf.eff <- grid_sf.eff %>%
 mutate(structured = st_effects, unstructured = un_effects)

# Grab our test fire again
fire_sf <- grid_sf.eff %>% filter(fire_id == 51)

# Plot structured vs unstructured effects
p1 <- ggplot(fire_sf) +
 geom_sf(aes(color = structured), color = NA) +
 scale_color_viridis_c(name = "Structured Component") +
 labs(title = "Structured Spatial Component") +
 theme_minimal()
p2 <- ggplot(fire_sf) +
 geom_sf(aes(color = unstructured), color = NA) +
 scale_color_viridis_c(name = "Unstructured Component") +
 labs(title = "Unstructured Spatial Component") +
 theme_minimal()

# Combine the plots
p1 + p2


#==============MODEL DIAGNOSTICS================#

# Extract fitted values (posterior mean estimates)
fit_res <- data.frame(
 Fitted = model_bl.frp.sp$summary.fitted.values$mean,
 Residuals = da$frp - fitted
)
hist(fit_res$Residuals, breaks = 50, col = "skyblue", main = "Histogram of Residuals")

# Create a dataframe for plotting
comparison_df <- data.frame(
 type = c(rep("Observed", length(da$frp)), rep("Fitted", length(model_bl.frp.sp$summary.fitted.values$mean))),
 value = c(da$frp, model_bl.frp.sp$summary.fitted.values$mean)
)
ggplot(comparison_df, aes(x = value, color = type, fill = type)) +
 geom_density(alpha = 0.4) +
 labs(title = "Density Plot of Observed vs. Fitted Values", 
      x = "Value", 
      y = "Density") +
 scale_fill_manual(values = c("Observed" = "steelblue", "Fitted" = "orange")) +
 scale_color_manual(values = c("Observed" = "steelblue", "Fitted" = "orange")) +
 theme_minimal()


# Plot residuals vs fitted values
ggplot(fit_res, aes(x = Fitted, y = Residuals)) +
 geom_point(alpha = 0.5) +
 geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
 labs(title = "Residuals vs Fitted Values",
      x = "Fitted Values",
      y = "Residuals") +
 theme_minimal()






