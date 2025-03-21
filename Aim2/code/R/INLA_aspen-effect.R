
# Load the required libraries
library(tidyverse)
library(sf) 
library(INLA) 
library(ggridges)
library(reshape2)
library(patchwork)
library(forcats)

# Environment variables
maindir <- '/Users/max/Library/CloudStorage/OneDrive-Personal/mcook/aspen-fire/Aim2/'

# load the cleaned model data frame
grid_tm <- paste0(maindir,"data/tabular/mod/model_data_cleaned.csv")
grid_tm <- read_csv(grid_tm) %>%
 # filter to fires with some aspen
 filter(
  fire_aspen == 1
 ) %>%
 # make sure factors are set
 mutate(
  first_obs_date = factor(first_obs_date)
 ) %>%
 group_by(grid_idx) %>%
 mutate(
  
  # Compute aspen-specific metrics only for quaking aspen
  aspen_ba_live = sum(if_else(species_gp_n == "quaking_aspen", ba_live, 0.0), na.rm = TRUE),
  aspen_tpp_live = sum(if_else(species_gp_n == "quaking_aspen", tpp_live, 0.0), na.rm = TRUE),
  aspen_ht_live = mean(if_else(species_gp_n == "quaking_aspen", ht_live, NA_real_), na.rm = TRUE),
  aspen_dia_live = mean(if_else(species_gp_n == "quaking_aspen", dia_live, NA_real_), na.rm = TRUE),
  
  # Compute proportions (avoiding division by zero)
  aspen_ba_pr = if_else(ba_live_total > 0, aspen_ba_live / ba_live_total, 0.0),
  aspen_tpp_pr = if_else(tpp_live_total > 0, aspen_tpp_live / tpp_live_total, 0.0)
 ) %>%
 ungroup() %>%
 # be sure there are no duplicate rows for grid/fire/species
 distinct(grid_idx, species_gp_n, .keep_all = TRUE) # remove duplicates

###############################
# check on the species counts #
grid_tm %>% 
 group_by(species_gp_n) %>%
 summarise(n = length(species_gp_n)) %>%
 arrange(desc(n))

##################################
# check on the aspen proportions #
summary(grid_tm$aspen_ba_pr)
summary(grid_tm$aspen_tpp_pr)

################################
# check how many grids and fires
length(unique(grid_tm$grid_idx))
length(unique(grid_tm$Fire_ID))

##############################################################
# Count how many times aspen co-occurs with other forest types
aspen_coo <- grid_tm %>%
 group_by(fortypnm_gp) %>%
 summarise(
  n_total = n(),  # Total occurrences of the forest type
  n_aspen = sum(grid_aspen == 1, na.rm = TRUE),  # How many co-occur with aspen
  aspen_ba_pr_mean = mean(aspen_ba_pr, na.rm=TRUE),
  aspen_ba_pr_med = median(aspen_ba_pr),
  aspen_ba_pr_p90 = quantile(aspen_ba_pr, probs = 0.9),
  pct_aspen_coo = round((n_aspen / n_total) * 100, 1)  # Convert to percentage
 ) %>%
 arrange(desc(pct_aspen_coo))  # Sort by co-occurrence percentage
print(aspen_coo)
rm(aspen_coo)


#===========MODEL SETUP==============#

# list of species names
spps <- c("quaking_aspen", "douglas_fir", "lodgepole_pine", 
          "ponderosa_pine", "spruce_fir", "piñon_juniper",
          "white_fir", "gambel_oak")
# force aspen to be the baseline
grid_tm <- grid_tm %>%
 mutate(
  species_gp_n = fct_relevel(species_gp_n, spps)
 )
# check the factor levels
# make sure aspen is first
levels(grid_tm$species_gp_n)

# prep the model data frame
# center and scale fixed effects
da <- grid_tm %>%
 mutate(
  # center/scale metrics / fixed effects
  across(
   c(aspen_ba_live, aspen_tpp_live, aspen_ht_live, # gridcell aspen metrics
     aspen_dia_live, aspen_ba_pr, aspen_tpp_pr, # gridcell aspen metrics
     H_ba, H_tpp, # grid-level diversity metrics
     ba_live, dia_live, ba_live, # structure metrics
     forest_pct, fortyp_pct, lf_canopy, fire_aspenpct, # gridcell forest characteristics
     ba_dead_total, tpp_live_total, tpp_dead_total, # gridcell forest characteristics
     erc, erc_dv, vpd, vs, #climate
     elev, slope, northness, tpi, # topography
     day_prop, overlap, # VIIRS detection characteristics
   ), ~ as.numeric(scale(., center=T, scale=T))
  ),
  log_fire_size = log(fire_acres), # log-scale fire size
  grid_aspen = as.factor(grid_aspen), # factor for grid aspen presence
  fire_aspen = as.factor(fire_aspen) # fire-level aspen presence
 ) %>%
 # keep just one row per grid cell by predominant type
 distinct(grid_idx, fortypnm_gp, .keep_all = TRUE) %>%
 arrange(grid_idx)

gc()

#########################################
# force aspen to be the baseline factor #
da <- da %>%
 mutate(
  species_gp_n = fct_relevel(species_gp_n, "quaking_aspen"),
  fortypnm_gp = fct_relevel(fortypnm_gp, "quaking_aspen")
 )
# check the factor levels
# make sure aspen is first
levels(da$species_gp_n)
levels(da$fortypnm_gp)

################################
# check the mean FRP for aspen #
# extract all aspen rows
aspen = da%>%filter(fortypnm_gp == "quaking_aspen")
aspen_mean_frpc = mean(aspen$log_frp_csum) # mean cumulative frp
aspen_mean_cbibc = mean(aspen$CBIbc_p90) # mean cumulative frp
rm(aspen)

#===========MODEL FITTING==============#

set.seed(456)

########################
# FIRE RADIATIVE POWER #

#################################################
# 1. Baseline model (no random or latent effects)

# setup the model formula
mf.frp <- log_frp_csum ~ 1 +
 vpd * (fortypnm_gp * aspen_ba_pr) + 
 fortyp_pct + # forest type percent cover
 H_ba + # gridcell structural diversity (dominance-based)
 lf_canopy + # gridcell forest canopy percent
 tpp_dead_total + # proportion live/dead basal area
 tpp_live_total + # total trees/pixel
 erc_dv + # day-of-burn climate/weather
 elev + slope + northness + tpi + # topography 
 overlap + # gridcell VIIRS overlap (cumulative)
 day_prop + # gridcell proportion daytime detections 
 log_fire_size + # log-scaled fire size
 fire_aspenpct + # fire-level aspen percent
 grid_aspen + # gridcell aspen presence
 # random effects
 f(Fire_ID, model = "iid", 
   hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.5)))
 ) + 
 f(first_obs_date, model = "iid",
   hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.5)))
 )

# fit the model
ml.frp <- inla(
 mf.frp, data = da,
 family = "gaussian",
 control.predictor = list(compute=T),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
 # set priors for the fixed effects
 control.fixed = list(
  mean.intercept = aspen_mean_frpc, # prior on intercept
  prec.intercept = 0.5,  # prior on intercept precision
  mean = 0, # mean on the fixed effects
  prec = 0.1  # shrinkage prior for all fixed effects
 ),
 control.family = list(
  # relax the variance assumptions
  hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.5)))
 ),
 # controls for posterior approximation and hyperparameter search
 control.inla = list(strategy = "laplace", int.strategy = "grid")
)
summary(ml.frp)
# check on predictive power of the random effects models
mean(ml.frp$cpo$cpo, na.rm = TRUE)


######################
# Spatial SPDE model #

# extracting spatial coordinates for grid centroids
grid_sf <- da %>%
 arrange(grid_index) %>%
 st_as_sf(., coords = c("x", "y"), crs = 4326)
# extract coordinates
coords <- grid_sf %>% st_coordinates(.)
# convert coordinates to a matrix for INLA
coords_mat <- as.matrix(coords) 

# Define the spatial mesh
mesh <- inla.mesh.2d(
 loc = coords_mat, # Locations (grid centroids)
 max.edge = c(1, 20), # Maximum edge lengths (inner and outer)
 cutoff = 0.01, # Minimum distance between points (0.01deg = ~1.1km)
 offset = c(0.5, 0.1) # Boundary buffer
)

# Build the SPDE model
spde.ml <- inla.spde2.pcmatern(
 # Mesh and smoothness parameter
 mesh = mesh, 
 alpha = 2,
 # P(practic.range < 0.3) = 0.5
 prior.range = c(0.1, 0.5),
 # P(sigma > 1) = 0.01
 prior.sigma = c(1, 0.1) # variance
)

# Compute the projector matrix (A)
A <- inla.spde.make.A(
 mesh = mesh,
 loc = coords_mat
)

# Assign the spatial index
field.idx <- inla.spde.make.index(
 name = "mesh.idx",
 n.spde = spde.ml$n.spde
)
str(field.idx)

# Extract the predictor variables
X <- da %>% 
 select(Fire_ID, first_obs_date, grid_index, fire_aspenpct,
        aspen_ba_live, aspen_tpp_live, aspen_ht_live, 
        aspen_dia_live, aspen_ba_pr, aspen_tpp_pr, 
        fortypnm_gp, fortyp_pct, species_gp_n,  
        tpp_live, tpp_live_pr, ba_live, ba_live_pr,    
        lf_canopy, H_tpp, H_ba, ba_dead_total,
        erc_dv, vpd, vs, elev, slope, northness, tpi,
        tpp_live_total, tpp_dead_total,
        grid_aspen, day_prop, overlap, log_fire_size)

# Create the INLA data stack
stack.frp <- inla.stack(
 data = list(log_frp_csum = da$log_frp_csum),
 A = list(A, 1),  
 tag = 'est',
 effects = list(
  c(field.idx),
  list(as.data.frame(X))
 )
)
dim(inla.stack.A(stack.frp))
rm(grid_sf, X)
gc()

##########################
# update the model formula
mf.frp.sp <- update(
 mf.frp, . ~ 1 + . + 
  f(mesh.idx, model = spde.ml) # spatial process model
)
# fit the model
ml.frp.re.sp <- inla(
 mf.frp.sp, 
 data = inla.stack.data(stack.frp), ,
 family = "gaussian",
 control.predictor = list(A = inla.stack.A(stack.frp), compute=T),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
 # set priors for the fixed effects
 control.fixed = list(
  mean.intercept = 0, # prior on intercept
  prec.intercept = 0.5,  # prior on intercept precision
  mean = 0, # mean on the fixed effects
  prec = 0.5  # shrinkage prior for all fixed effects
 ),
 control.family = list(
  # relax the variance assumptions
  hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.5)))
 ),
 # controls for posterior approximation and hyperparameter search
 control.inla = list(strategy = "laplace", int.strategy = "grid")
)
summary(ml.frp.re.sp)
mean(ml.frp.re.sp$cpo$cpo, na.rm = TRUE)


#===========POSTERIOR EFFECTS===========#

#########################################
# Plot all of the posterior fixed effects
# Extract fixed effect marginals
frp_marginals <- ml.frp.re.sp$marginals.fixed
# Tidy marginals for all fixed effects
tidy.effects.frp <- tibble::tibble(
 parameter = names(frp_marginals),
 data = purrr::map(frp_marginals, ~ as.data.frame(.x))
) %>%
 unnest(data) %>%
 filter(!parameter %in% c("(Intercept)","aspen1","day_prop","overlap","log_fire_size")) %>%  
 mutate(
  # Extract species names from parameters
  fortyp = case_when(
   str_detect(parameter, "fortypnm_gp") ~ str_extract(parameter, "(?<=fortypnm_gp)[a-zA-Z_ñ]+"),
   TRUE ~ NA_character_  # For non-species effects
  )
 ) %>%
 drop_na() %>%
 # tidy the parameter
 mutate(
  effect = case_when(
   str_detect(parameter, "vpd:") ~ "VPD-mediated",  
   str_detect(parameter, ":aspen_ba_pr") & !str_detect(parameter, "vpd:")  ~ "Aspen proportion",
   TRUE ~ "Gloabl effect"  # For non-species effects
  ),
  # handle forest type names
  fortyp = recode(
   fortyp,
   "quaking_aspen" = "Quaking aspen",
   "lodgepole_pine" = "Lodgepole pine",
   "douglas_fir" = "Douglas-fir",
   "white_fir" = "White fir",
   "gambel_oak" = "Gambel oak",
   "piñon_juniper" = "Piñon-juniper",
   "ponderosa_pine" = "Ponderosa pine",
   "spruce_fir" = "Spruce-fir"
  )
 ) %>%
 # Calculate mean effect size for ordering
 group_by(effect) %>%
 mutate(mean_effect = mean(x, na.rm = TRUE)) %>%
 ungroup() 
glimpse(tidy.effects.frp)
 
# reorder the factors
tidy.effects.frp <- tidy.effects.frp %>%
 mutate(
  effect_order = case_when(
   !is.na(fortyp) ~ 2,  # Species-specific effects
   TRUE ~ 2              # Global effects
  ),
  effect = factor(effect, levels = tidy.effects.frp %>%
                   arrange(effect_order, desc(mean_effect)) %>%
                   pull(effect) %>%
                   unique()),
  fill_species = ifelse(is.na(fortyp), "Global effect", fortyp)
 ) %>%
 mutate(fill_species = recode(
  fill_species,
  "quaking_aspen" = "Quaking aspen",
  "lodgepole_pine" = "Lodgepole pine",
  "douglas_fir" = "Douglas-fir",
  "white_fir" = "White fir",
  "gambel_oak" = "Gambel oak",
  "piñon_juniper" = "Piñon-juniper",
  "ponderosa_pine" = "Ponderosa pine",
  "spruce_fir" = "Spruce-fir"
 )) %>%
 mutate(
  fill_species = factor(
   fill_species,
   levels = c("Lodgepole pine", "Douglas-fir", "White fir",
              "Gambel oak", "Piñon-juniper", "Ponderosa pine",
              "Spruce-fir", "Quaking aspen")))

# color map
color_map = c(
 "Aspen proportion" = "#e6ab02",  # Keep orange for direct effect
 "VPD-mediated" = "gray80"  # Keep gray for mediation effect
)

# make the ridge plot
frp.p1 <- ggplot(tidy.effects.frp, 
                 aes(x = x, y = fill_species, height = y, 
                     fill = factor(effect, levels = c("VPD-mediated","Aspen proportion")),
                     alpha = factor(effect, levels = c("VPD-mediated","Aspen proportion")))) +
 geom_density_ridges(
  stat = "identity", scale = 1.5, show.legend=T
 ) +
 geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
 labs(
  x = "",
  y = "Predominant Forest Type",
  fill = "Model effect"
 ) +
 scale_fill_manual(
  values = color_map,
  guide = guide_legend(title = "Model effect")
 ) +
 scale_alpha_manual(
  values = c(
   "Aspen proportion" = 0.95,
   "VPD-mediated" = 0.6
  ),
  guide = "none"  # Hide alpha from legend
 ) +
 coord_cartesian(xlim=c(-0.26,0.51)) +
 theme_classic() +
 theme(
  axis.text.y = element_text(angle = 0, hjust = 1, size=9),
  axis.text.x = element_text(angle = 0, hjust = 0, size=9),
  axis.title.y = element_text(size = 10, margin = margin(r = 12)),
  axis.title.x = element_blank(),
  legend.position = c(0.85, 0.80),
  legend.background = element_rect(
   fill = scales::alpha("white", 0.6), 
   color = NA, linewidth = 1),
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 9)
 )
frp.p1

#######################################
# Version 2: no Gambel oak or White fir
frp.p1.1 <- tidy.effects.frp %>%
 filter(!fill_species %in% c("Gambel oak", "White fir", "Piñon-juniper")) %>%
 ggplot(., aes(x = x, y = fill_species, height = y,
               fill = factor(effect, levels = c("VPD-mediated","Aspen proportion")),
               alpha = factor(effect, levels = c("VPD-mediated","Aspen proportion")))) +
 geom_density_ridges(
  stat = "identity", scale = 1.5, show.legend=T) +
 geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
 labs(
  x = "",
  y = "Predominant Forest Type",
  fill = "Model effect"
 ) +
 scale_fill_manual(
  values = color_map,
  guide = guide_legend(title = "Model effect")
 ) +
 scale_alpha_manual(
  values = c(
   "Aspen proportion" = 0.95,
   "VPD-mediated" = 0.6
  ),
  guide = "none"  # Hide alpha from legend
 ) +
 coord_cartesian(xlim=c(-0.268,0.51)) +
 theme_minimal() +
 theme(
  axis.text.y = element_text(angle = 0, hjust = 1, size=9),
  axis.text.x = element_text(angle = 0, hjust = 0, size=9),
  axis.title.y = element_text(size = 10, margin = margin(r = 12)),
  axis.title.x = element_blank(),
  legend.position = c(0.85, 0.80),
  legend.background = element_rect(
   fill = scales::alpha("white", 0.6), 
   color = NA, size = 1),
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 9)
 )
frp.p1.1





#==========#

########################
# COMPOSITE BURN INDEX #

#################################################
# 1. Baseline model (no random or latent effects)

# setup the model formula
mf.cbi <- CBIbc_p90 ~ 1 + 
 grid_aspen + # grid-level aspen presence
 # vpd * aspen_qmd_live + # grid-level aspen QMD
 vpd * (fortypnm_gp:aspen_ba_pr) + # maj. forest type by aspen live basal area * VPD
 H_tpp + # grid-level diversity metric
 canopypct + # grid-level canopy percent
 erc_dv + vpd + vs + # day-of-burn climate/weather
 elev + slope + tpi + chili # grid-level topography
# fit the model
ml.cbi <- inla(
 mf.cbi, data = da,
 family = "gamma",
 control.predictor = list(compute=T),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
 control.fixed = list(
  prec = list(prior = "loggamma", param = c(0.5, 0.1))
 ),
 control.family = list(
  hyper = list(prec = list(prior = "loggamma", param = c(0.5, 0.1)))
 ),
 control.inla = list(strategy = "adaptive", int.strategy = "grid")
)
summary(ml.cbi)
# check on predictive power of the random effects models
mean(ml.cbi$cpo$cpo, na.rm = TRUE)


######################################
# 2. Baseline + Temporal Random Effect

# update the model formula
mf.cbi.re <- update(
 mf.cbi, . ~ 1 + . + 
  f(first_obs_date, model = "iid", 
    hyper = list(
     prec = list(prior = "pc.prec", param = c(0.5, 0.01))
    )) # temporal random effect
)
# fit the model
ml.cbi.re <- inla(
 mf.cbi.re, data = da,
 family = "gamma",
 control.predictor = list(compute=T),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
 control.fixed = list(
  prec = list(prior = "loggamma", param = c(0.5, 0.1))
 ),
 control.family = list(
  hyper = list(prec = list(prior = "loggamma", param = c(0.5, 0.1)))
 ),
 control.inla = list(strategy = "adaptive", int.strategy = "grid")
)
summary(ml.cbi.re)
# check on predictive power of the random effects models
mean(ml.cbi.re$cpo$cpo, na.rm = TRUE)


###################################################
# 3. Baseline + Temporal + Fire-level Random Effect

# update the model formula
mf.cbi.re2 <- update(
 mf.cbi.re, . ~ 1 + . + 
  f(fire_id, model = "iid", 
    hyper = list(
     prec = list(prior = "pc.prec", param = c(0.5, 0.01))
    )) # fire-level random effect
)
# fit the model
ml.cbi.re2 <- inla(
 mf.cbi.re2, data = da,
 family = "gaussian",
 control.predictor = list(compute=T),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
 control.fixed = list(
  prec = list(prior = "loggamma", param = c(0.5, 0.1))
 ),
 control.family = list(
  hyper = list(prec = list(prior = "loggamma", param = c(0.5, 0.1)))
 ),
 control.inla = list(strategy = "adaptive", int.strategy = "grid")
)
summary(ml.cbi.re2)
# check on predictive power of the random effects models
mean(ml.cbi.re2$cpo$cpo, na.rm = TRUE)



#######################
# 4. Spatial SPDE model

# extracting spatial coordinates for grid centroids
grid_sf <- da %>%
 arrange(grid_index) %>%
 st_as_sf(., coords = c("x", "y"), crs = 4326)
# extract coordinates
coords <- grid_sf %>% st_coordinates(.)
# convert coordinates to a matrix for INLA
coords_mat <- as.matrix(coords) 

# Define the spatial mesh
mesh <- inla.mesh.2d(
 loc = coords_mat, # Locations (grid centroids)
 max.edge = c(1, 20), # Maximum edge lengths (inner and outer)
 cutoff = 0.01, # Minimum distance between points (0.01deg = ~1.1km)
 offset = c(0.5, 0.1) # Boundary buffer
)
# Plot mesh to check
plot(mesh, main = "SPDE Mesh for FRP Model")
points(coords, col = "red", pch = 20)

rm(coords)
gc()

######################
# Build the SPDE model
spde.ml <- inla.spde2.pcmatern(
 # Mesh and smoothness parameter
 mesh = mesh, 
 alpha = 2,
 # P(practic.range < 0.3) = 0.5
 prior.range = c(0.1, 0.5),
 # P(sigma > 1) = 0.01
 prior.sigma = c(0.5, 0.01) # variance
)

# Compute the projector matrix (A)
A <- inla.spde.make.A(
 mesh = mesh,
 loc = coords_mat
)

# Assign the spatial index
field.idx <- inla.spde.make.index(
 name = "mesh.idx",
 n.spde = spde.ml$n.spde
)
str(field.idx)

# Extract the predictor variables
X <- da %>% 
 select(fire_id, first_obs_date, grid_index, grid_aspen,
        fortypnm_gp, total_tpp, forest_pct, canopypct, 
        H_tpp, H_ba, aspen_ba_pr, aspen_qmd_live, aspen_tpp_pr, 
        erc_dv, vpd, vs, slope, tpi, chili, elev,
        day_prop, afd_count, overlap)
head(X)

# Create the INLA data stack
stack.cbi <- inla.stack(
 data = list(CBIbc_p90 = da$CBIbc_p90),
 A = list(A, 1),  
 tag = 'est',
 effects = list(
  c(field.idx),
  list(as.data.frame(X))
 )
)
dim(inla.stack.A(stack.cbi))

rm(grid_sf, X)
gc()

##########################
# update the model formula
mf.cbi.re.sp <- update(
 mf.cbi.re2, . ~ 1 + . + 
  f(mesh.idx, model = spde.ml) # spatial process model
)
# fit the model
ml.cbi.re.sp <- inla(
 mf.cbi.re.sp, 
 data = inla.stack.data(stack.cbi), ,
 family = "gaussian",
 control.predictor = list(A = inla.stack.A(stack.cbi), compute=T),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
 control.fixed = list(
  prec = list(prior = "loggamma", param = c(0.5, 0.1))
 ),
 control.family = list(
  hyper = list(prec = list(prior = "loggamma", param = c(0.5, 0.1)))
 ),
 control.inla = list(strategy = "adaptive", int.strategy = "grid")
)
summary(ml.cbi.re.sp)
mean(ml.cbi.re.sp$cpo$cpo, na.rm = TRUE)



#===========POSTERIOR EFFECTS===========#

#########################################
# Plot all of the posterior fixed effects
# Extract fixed effect marginals
cbi_marginals <- ml.cbi.re.sp$marginals.fixed
# Tidy marginals for all fixed effects
tidy.effects.cbi <- tibble::tibble(
 parameter = names(cbi_marginals),
 data = purrr::map(cbi_marginals, ~ as.data.frame(.x))
) %>%
 unnest(data) %>%
 mutate(
  # Extract species names from parameters
  fortyp = case_when(
   str_detect(parameter, "fortypnm_gp") ~ str_extract(parameter, "(?<=fortypnm_gp)[a-zA-Z_ñ]+"),
   TRUE ~ NA_character_  # For non-species effects
  )
 ) %>%
 drop_na() %>%
 # tidy the parameter
 mutate(
  parameter = case_when(
   str_detect(parameter, "vpd:") ~ "VPD-mediated",
   TRUE ~ "Quaking aspen live BA"
  ),
  # handle forest type names
  fortyp = recode(
   fortyp,
   "mixed_conifer" = "Mixed-conifer",
   "lodgepole" = "Lodgepole",
   "ponderosa" = "Ponderosa",
   "spruce_fir" = "Spruce-fir",
   "piñon_juniper" = "Piñon-juniper",
   "quaking_aspen" = "Quaking aspen"
  )
 )

# sort the dataframe by mean effect (non VPD-mediated)
order.df <- tidy.effects.cbi %>%
 filter(parameter == "Quaking aspen live BA") %>%
 group_by(fortyp) %>%
 summarise(effect_mn = mean(x, na.rm = TRUE), .groups = "drop") %>%
 arrange(desc(effect_mn))  # Sort by mean effect
# reorder the factors
tidy.effects.cbi <- tidy.effects.cbi %>%
 left_join(order.df, by = "fortyp") %>%  # Merge mean effect values
 mutate(fortyp = factor(fortyp, levels = order.df$fortyp))  # Apply ordered factor levels
rm
# Verify ordering
levels(tidy.effects.cbi$fortyp)
head(tidy.effects.cbi)

# make the ridge plot
cbi.p1 <- ggplot(tidy.effects.cbi, 
                 aes(x = x, y = fortyp, height = y, 
                     fill = factor(parameter, levels = c("Quaking aspen live BA", "VPD-mediated")),
                     alpha = factor(parameter, levels = c("Quaking aspen live BA", "VPD-mediated")))) +
 geom_density_ridges(
  stat = "identity", scale = 1.5, show.legend=F) +
 geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
 labs(
  x = "Effect Size",
  y = "Predominant Forest Type",
 ) +
 scale_fill_manual(
  values = c(
   "Quaking aspen live BA" = "#1b9e77",
   "VPD-mediated" = "gray80"
  )
 ) +
 scale_alpha_manual(
  values = c(
   "Quaking aspen live BA" = 0.9,
   "VPD-mediated" = 0.7
  ),
  guide = "none"  # Hide alpha from legend
 ) +
 coord_cartesian(xlim=c(-0.26,0.51)) +
 theme_minimal() +
 theme(
  axis.text.y = element_text(angle = 0, hjust = 1, size=9),
  axis.text.x = element_text(angle = 0, hjust = 0, size=9),
  axis.title.y = element_text(size = 10, margin = margin(r = 12)),
  axis.title.x = element_text(size = 10, margin = margin(t = 12))
 )
cbi.p1


# Make a combined figure (FRP and CBI)
# Load the two ridge plots (assuming they are `frp_plot` and `cbi_plot`)
frp.cbi.p1 <- frp.p1 / cbi.p1 # "/" stacks them, "|" would place them side-by-side
frp.cbi.p1

# Save the plot
out_png <- paste0(maindir, 'figures/INLA_ASPEN_VPD_FRP-CBI_fortyp.png')
ggsave(out_png, plot = frp.cbi.p1,
       width = 6, height = 5, bg = 'white')


