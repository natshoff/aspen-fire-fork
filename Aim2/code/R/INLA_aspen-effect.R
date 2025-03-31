
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

#=========Load the Prepped gridcell data=========#

# gridcell_prep.R
fp <- paste0(maindir,"data/tabular/mod/model_data_cleaned.csv")
grid_tm <- read_csv(fp) %>%
 # filter to aspen fires
 filter(fire_aspen == 1) %>%
 # prep some of the attributes for modeling
 mutate(
  first_obs_date = as.factor(first_obs_date),
  fire_aspen = as.factor(fire_aspen),
  grid_aspen = as.factor(grid_aspen),
  log_fire_size = log(fire_acres) # log-scale fire size
 ) %>%
 # flag re-burn gridcells
 group_by(grid_index) %>%
 mutate(
  # Get earliest fire date for this gridcell
  first_fire_dt = min(fire_ig_dt, na.rm = TRUE),
  # Flag reburn: if current fire date is after earliest
  reburn = as.factor(fire_ig_dt > first_fire_dt)
 ) %>%
 ungroup() %>%
 # filter out re-burns
 filter(
  reburn == "FALSE"
 ) %>%
 # select the model attributes ...
 select(
  grid_idx, grid_index, Fire_ID, fortypnm_gp, first_obs_date,
  tpp_live_total, tpp_dead_total, tpp_live_pr,
  ba_live_total, ba_dead_total, ba_live_pr,  
  lf_canopy, lf_height, forest_pct, fortyp_pct,
  erc, erc_dv, vpd, vpd_dv, vs, elev, slope, northness, tpi,
  day_prop, overlap, dist_to_perim, log_fire_size, 
  log_frp_csum, CBIbc_p90, CBIbc_p95, CBIbc_p99, 
  x, y, fire_aspen, grid_aspen,
  aspen_ba_pr, aspen_tpp_pr, aspen_ba_live,
  aspen_tpp_live, aspen_ht_live, aspen_dia_live
 ) %>%
 # keep just one row per grid cell by predominant type
 distinct(grid_idx, fortypnm_gp, .keep_all = TRUE)
glimpse(grid_tm)

###################################################
# check how many grids are majority forested, etc #
print(length(unique(grid_tm$grid_idx))) # unique gridcells
# percent of gridcells majority forested
# Check how many grids are > 50% forested
print(paste0(
 "Percent forested gridcells (>50% gridcell area): ",
 round(dim(grid_tm %>% filter(forest_pct > 0.50) %>% distinct(grid_idx))[1]/
        dim(grid_tm %>% distinct(grid_idx))[1], 3) * 100, "%"
))
# # retain predominantly forested gridcells ...
# grid_tm <- grid_tm %>%
#  filter(forest_pct >= 0.50)

##########################################
# filter fires with not enough gridcells #
# should have at least 10 (?)
# check on the grid cell counts
gridcell_counts <- grid_tm %>%
 distinct(Fire_ID, grid_idx) %>% # keep only distinct rows
 group_by(Fire_ID) %>%
 summarise(n_gridcells = n())
# check the distribution
summary(gridcell_counts$n_gridcells)
# calculate the quantile distribution
(qt <- tibble::enframe(
 round(quantile(gridcell_counts$n_gridcells, probs = seq(.1, .9, by = .1))),
 name = "qt", value = "val"
))
(qt10 = qt[1,]$val) # 10%
# # filter fires below the 10th percentile
# grid_tm <- grid_tm %>%
#  group_by(Fire_ID) %>%
#  filter(n() >= qt10) %>%
#  ungroup()
# tidy up!
rm(gridcell_counts,qt)

##################################
# check on the aspen proportions #
summary(grid_tm$aspen_ba_pr)
summary(grid_tm$aspen_tpp_pr)

################################
# check how many grids and fires
length(unique(grid_tm$grid_idx))
length(unique(grid_tm$Fire_ID))

# relevel the grid_aspen factor
grid_tm$grid_aspen = fct_relevel(grid_tm$grid_aspen, c("1","0"))
levels(grid_tm$grid_aspen)


#===========MODEL SETUP==============#

# prep the model data frame
# center and scale fixed effects
da <- grid_tm %>%
 mutate(
  # force aspen to be the baseline factor
  fortypnm_gp = fct_relevel(fortypnm_gp, "quaking_aspen"),
  # center/scale metrics / fixed effects
  across(
   c(aspen_ba_live, aspen_tpp_live, aspen_ht_live, # gridcell aspen metrics
     aspen_dia_live, aspen_ba_pr, aspen_tpp_pr, # gridcell aspen metrics
     forest_pct, fortyp_pct, lf_canopy, lf_height, # gridcell forest characteristics
     ba_dead_total, ba_live_total, tpp_live_total, tpp_dead_total, # gridcell forest characteristics
     erc, erc_dv, vpd, vs, #climate
     elev, slope, northness, tpi, # topography
     day_prop, overlap, # VIIRS detection characteristics
     dist_to_perim 
   ), ~ as.numeric(scale(., center=T, scale=T))
  )
 ) %>%
 arrange(grid_idx)
rm(grid_tm)
gc()

# check the factor levels
# make sure aspen is first
levels(da$fortypnm_gp)



#===========MODEL FITTING==============#
set.seed(456)

####################
# define some priors
# fixed effects
pr.fixed <- list(
 prec.intercept = 1e-6, 
 prec = 0.001  # shrinkage prior
)
# for the Gaussian precision
pr.family <- list(
 hyper = list(
  prec = list(
   prior = "pc.prec", param = c(1, 0.5))
 )
)
# pc.prec for random effects, etc.
pr.pc_prec <- list(
 prec = list(
  prior = "pc.prec", param = c(1, 0.5))
)

########################
# FIRE RADIATIVE POWER #

#################################################
# 1. Baseline model (no random or latent effects)

# setup the model formula
mf.frp <- log_frp_csum ~ 1 +
 fortypnm_gp:aspen_ba_pr + 
 vpd:fortypnm_gp:aspen_ba_pr +
 fortyp_pct + # forest type percent cover
 lf_canopy + # gridcell forest canopy percent
 tpp_dead_total + # proportion live/dead basal area
 tpp_live_total + # total trees/pixel
 erc_dv + vs + # day-of-burn climate/weather
 elev + slope + northness + tpi + # topography 
 overlap + # gridcell VIIRS overlap (cumulative)
 day_prop + # gridcell proportion daytime detections 
 dist_to_perim + # gridcell distance to perimeter
 grid_aspen + # gridcell aspen presence
 # random effects
 f(Fire_ID, model = "iid", hyper = pr.pc_prec)

# fit the model
ml.frp <- inla(
 mf.frp, data = da,
 family = "gaussian",
 control.predictor = list(compute=T),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
 control.fixed = pr.fixed, # regularization on fixed effects
 control.family = pr.family # on the Gaussian observation
)
summary(ml.frp) # model summary table
# check on predictive power (CPO)
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

rm(coords) # tidy up

# Build the SPDE model
spde.ml <- inla.spde2.pcmatern(
 # Mesh and smoothness parameter
 mesh = mesh, 
 alpha = 2,
 # P(practic.range < 0.3) = 0.5
 prior.range = c(0.3, 0.5), # 50% certainty that range is below ~5km
 # P(sigma > 1) = 0.01
 prior.sigma = c(1, 0.01) # variance
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
 select(Fire_ID, first_obs_date, grid_index,
        aspen_ba_live, aspen_tpp_live, aspen_ht_live, 
        aspen_dia_live, aspen_ba_pr, aspen_tpp_pr, 
        fortypnm_gp, fortyp_pct, ba_live_pr,    
        lf_canopy, ba_dead_total, ba_live_total,
        tpp_dead_total, tpp_live_total,
        erc_dv, vpd, vs, elev, slope, northness, tpi,
        grid_aspen, day_prop, overlap, dist_to_perim)

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
 control.fixed = pr.fixed, # regularization on fixed effects
 control.family = pr.family # on the Gaussian observation
)
summary(ml.frp.re.sp)
mean(ml.frp.re.sp$cpo$cpo, na.rm = TRUE)

rm(stack.frp)
gc()


#=================MODEL STATEMENTS=================#

# compute the exponentiated effects
exp.frp <- ml.frp.re.sp$summary.fixed %>%
 rownames_to_column(var = "parameter") %>%
 mutate(
  exp_mean = exp(mean) - 1,  # Convert log(FRP) effect to % difference
  lower_ci = exp(`0.025quant`) - 1,  # 2.5% CI bound
  upper_ci = exp(`0.975quant`) - 1   # 97.5% CI bound
 )
# save this table
write_csv(exp.frp, paste0(maindir,"data/tabular/mod/results/INLA_exp_FRP_aspenEffect.csv"))
# check results
exp.frp%>%select(parameter,exp_mean,lower_ci,upper_ci)


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
 filter(!parameter %in% c("(Intercept)","day_prop","overlap")) %>%  
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
   str_detect(parameter, ":vpd") ~ "VPD-mediated",  
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

# Get 95% CI bounds for each parameter
ci_summary <- map_dfr(ml.frp.re.sp$marginals.fixed, function(m) {
 tibble(
  mean = inla.emarginal(identity, m),
  lower = inla.qmarginal(0.025, m),
  upper = inla.qmarginal(0.975, m)
 )
}, .id = "parameter") %>%
 left_join(
  tidy.effects.frp %>% 
   select(parameter, fill_species, effect) %>% 
   distinct(), 
  by = "parameter"
 ) %>%
 drop_na()


###########
# color map
color_map = c(
 "Aspen proportion" = "#E69F00",  # Keep orange for direct effect
 "VPD-mediated" = "gray80"  # Keep gray for mediation effect
)


#######################
# make the ridge plot #
# sort the factors by effect
aspen_order <- tidy.effects.frp %>%
 filter(effect == "Aspen proportion") %>%
 group_by(fill_species) %>%
 summarize(mean_effect = mean(x, na.rm = TRUE)) %>%
 arrange(-mean_effect) %>%  # for most negative at top
 pull(fill_species)

# plot it
(frp.p1 <- tidy.effects.frp %>%
  mutate(fill_species = factor(fill_species, levels = aspen_order)) %>%
  ggplot(., aes(
   x = x, y = fill_species, height = y, 
   fill = factor(effect, levels = c("VPD-mediated","Aspen proportion")),
   alpha = factor(effect, levels = c("VPD-mediated","Aspen proportion")))
  ) + 
  
  geom_density_ridges(stat = "identity", scale = 1.5, show.legend=T) +
  
  # geom_pointrange(
  #  data = ci_summary,
  #  aes(x = mean, xmin = lower, xmax = upper, y = fill_species),
  #  inherit.aes = FALSE,
  #  shape = 21, size = 0.3,
  #  stroke = 0.3,
  #  fill = "black",
  #  color = "black"
  # ) +
  
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
  # add a subplot label
  annotate("text", x = -0.48, y = 8.8,
           label = expression(bold("(A)") ~ "FRPc"),
           size = 4, hjust = 0.2) +
  coord_cartesian(xlim=c(-0.52,0.66)) +
  theme_classic() +
  theme(
   axis.text.y = element_text(angle = 0, hjust = 1, size=9),
   axis.text.x = element_text(angle = 0, hjust = 0, size=9),
   axis.title.y = element_text(size = 11, margin = margin(r = 12)),
   axis.title.x = element_blank(),
   legend.position = c(0.85, 0.88),
   legend.background = element_rect(
    fill = scales::alpha("white", 0.6), 
    color = NA, linewidth = 1),
   legend.title = element_text(size = 11),
   legend.text = element_text(size = 10)
  ))

#############
# version two

###########
# color map
color_map = c(
 "Aspen proportion" = "#E69F00",  # Keep orange for direct effect
 "VPD-mediated" = "#E69F00"  # Keep gray for mediation effect
)

(frp.p1.1 <- tidy.effects.frp %>%
  mutate(fill_species = factor(fill_species, levels = aspen_order)) %>%
  arrange(match(effect, c("VPD-mediated", "Aspen proportion"))) %>%
  ggplot(., aes(x = x, y = fill_species, height = y,
                fill = effect,
                alpha = effect,
                color = effect
  )) +
  geom_density_ridges(stat = "identity", scale = 1.5, aes(alpha = effect)) +
  # geom_density_ridges(stat = "identity", scale = 1.5, color = NA, alpha = 0) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(
   x = "Effect size",
   y = "Predominant Forest Type",
   fill = "Effect"
  ) +
  scale_fill_manual(
   values = c(
    "Aspen proportion" = "#E69F00",
    "VPD-mediated" = scales::alpha("#E69F00", 0.3)
   ),
   labels = c(
    "Aspen proportion" = "Aspen dominance (basal area)",
    "VPD-mediated" = "VPD × aspen dominance"
   ),
   guide = guide_legend(
    override.aes = list(
     fill = c("#E69F00", scales::alpha("#E69F00", 0.3)),
     alpha = c(1, 0.3),  # Match transparency
     color = c(
      scales::alpha("black", 0.7),
      scales::alpha("#E69F00", 0.3))
    ),
    order = 1  # Ensure this legend appears first
   )
  ) +
  scale_alpha_manual(
   values = c(
    "VPD-mediated" = 0.4,
    "Aspen proportion" = 0.98
   ),
   guide = "none"
  ) +
  scale_color_manual(
   values = c(
    "Aspen proportion" = scales::alpha("black", 0.7),
    "VPD-mediated" = scales::alpha("#E69F00", 0.5)
   ),
   guide = "none"
  ) +
  # add a subplot label
  annotate(
   "text", x = -0.48, y = 8.8,
   label = expression(bold("(A)") ~ "FRPc"),
   size = 4, hjust = 0.2
  ) +
  theme_classic() +
  coord_cartesian(xlim=c(-0.52,0.66)) +
  theme(
   axis.text.y = element_text(angle = 0, hjust = 1, size=9),
   axis.text.x = element_text(angle = 0, hjust = 0, size=9),
   axis.title.y = element_text(size = 10, margin = margin(r = 12)),
   axis.title.x = element_blank(),
   legend.position = c(0.80, 0.88),
   legend.background = element_rect(
    fill = scales::alpha("white", 0.6), 
    color = NA, linewidth = 1),
   legend.title = element_text(size = 10),
   legend.text = element_text(size = 9)
  ))
# save the plot
out_png <- paste0(maindir, 'figures/INLA_AspenEffect_FRP_fortyp.png')
ggsave(out_png, plot = frp.p1.1, dpi = 500, width = 7, height = 5, bg = 'white')


########################################
# make a version without VPD mediation #
(frp.p1.2 <- tidy.effects.frp %>%
  mutate(fill_species = factor(fill_species, levels = aspen_order)) %>%
  filter(effect == "Aspen proportion") %>%
  ggplot(., aes(x = x, y = fill_species, height = y,
                fill = effect,
                alpha = effect,
                color = effect
  )) +
  geom_density_ridges(stat = "identity", scale = 1.5, aes(alpha = effect)) +
  # geom_density_ridges(stat = "identity", scale = 1.5, color = NA, alpha = 0) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(
   x = "Effect size",
   y = "Predominant Forest Type",
   fill = "Effect"
  ) +
  scale_fill_manual(
   values = c("Aspen proportion" = "#E69F00"),
   labels = c("Aspen proportion" = "Aspen dominance (basal area)"),
   guide = guide_legend(
    override.aes = list(
     fill = c("#E69F00"),
     alpha = c(1),  # Match transparency
     color = c(scales::alpha("black", 0.7)),
     order = 1  # Ensure this legend appears first
    )
   )) +
  scale_alpha_manual(
   values = c("Aspen proportion" = 0.98),
   guide = "none"
  ) +
  scale_color_manual(
   values = c(
    "Aspen proportion" = scales::alpha("black", 0.7)
   ),
   guide = "none"
  ) +
  # add a subplot label
  annotate(
   "text", x = -0.48, y = 8.8,
   label = expression(bold("(A)") ~ "FRPc"),
   size = 4, hjust = 0.2
  ) +
  theme_classic() +
  coord_cartesian(xlim=c(-0.52,0.66)) +
  theme(
   axis.text.y = element_text(angle = 0, hjust = 1, size=9),
   axis.text.x = element_text(angle = 0, hjust = 0, size=9),
   axis.title.y = element_text(size = 10, margin = margin(r = 12)),
   axis.title.x = element_blank(),
   legend.position = c(0.80, 0.88),
   legend.background = element_rect(
    fill = scales::alpha("white", 0.6), 
    color = NA, linewidth = 1),
   legend.title = element_text(size = 10),
   legend.text = element_text(size = 9)
  ))
# save the plot
out_png <- paste0(maindir, 'figures/INLA_AspenEffect_FRP_fortyp_noVPD.png')
ggsave(out_png, plot = frp.p1.2, dpi = 500, width = 7, height = 5, bg = 'white')

# #######################################
# # Version 2: no Gambel oak or White fir
# (frp.p1.1 <- tidy.effects.frp %>%
#  filter(!fill_species %in% c("Gambel oak", "White fir")) %>%
#  ggplot(., aes(x = x, y = fill_species, height = y,
#                fill = factor(effect, levels = c("VPD-mediated","Aspen proportion")),
#                alpha = factor(effect, levels = c("VPD-mediated","Aspen proportion")))) +
#  geom_density_ridges(
#   stat = "identity", scale = 1.5, show.legend=T) +
#  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
#  labs(
#   x = "",
#   y = "Predominant Forest Type",
#   fill = "Model effect"
#  ) +
#  scale_fill_manual(
#   values = color_map,
#   guide = guide_legend(title = "Model effect")
#  ) +
#  scale_alpha_manual(
#   values = c(
#    "Aspen proportion" = 0.95,
#    "VPD-mediated" = 0.6
#   ),
#   guide = "none"  # Hide alpha from legend
#  ) +
#  # coord_cartesian(xlim=c(-0.268,0.51)) +
#  theme_classic() +
#  theme(
#   axis.text.y = element_text(angle = 0, hjust = 1, size=9),
#   axis.text.x = element_text(angle = 0, hjust = 0, size=9),
#   axis.title.y = element_text(size = 10, margin = margin(r = 12)),
#   axis.title.x = element_blank(),
#   legend.position = c(0.85, 0.88),
#   legend.background = element_rect(
#    fill = scales::alpha("white", 0.6), 
#    color = NA, linewidth = 1),
#   legend.title = element_text(size = 10),
#   legend.text = element_text(size = 9)
#  ))


#==========#

########################
# COMPOSITE BURN INDEX #

da <- da %>%
 # either add a constant or remove zeros
 # gamma requires strictly positive
 # CBI = 0 is likely "unburned" but had some FRP detection (?)
 # filter zeros
 filter(CBIbc_p90 > 0)
# # add a very small constant
# mutate(CBIbc_p90 = CBIbc_p90 + 1e-6)


#################################################
# 1. Baseline model (no random or latent effects)

# setup the model formula
mf.cbi <- CBIbc_p90 ~ 1 +
 fortypnm_gp:aspen_ba_pr + 
 vpd:fortypnm_gp:aspen_ba_pr +
 fortyp_pct + # forest type percent cover
 lf_canopy + # gridcell forest canopy percent
 tpp_dead_total + # proportion live/dead basal area
 tpp_live_total + # total trees/pixel
 erc_dv + vs + # day-of-burn climate/weather
 elev + slope + northness + tpi + # topography 
 overlap + # gridcell VIIRS overlap (cumulative)
 day_prop + # gridcell proportion daytime detections 
 dist_to_perim + # gridcell distance to perimeter
 grid_aspen + # gridcell aspen presence
 # random effects
 f(Fire_ID, model = "iid", hyper = pr.pc_prec)

# fit the model
ml.cbi <- inla(
 mf.cbi, data = da,
 family = "gamma",
 control.predictor = list(compute=T),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
 control.fixed = pr.fixed, # regularization on fixed effects
 control.family = pr.family # on the Gaussian observation
)
summary(ml.cbi) # model summary table
# check on predictive power (CPO)
mean(ml.cbi$cpo$cpo, na.rm = TRUE)


#######################
# 4. Spatial SPDE model

######################
# Spatial SPDE model #

# extract gridcell coordinates
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
# # Plot mesh to check
# plot(mesh, main = "SPDE Mesh for FRP Model")
# points(coords, col = "red", pch = 20)
rm(coords)

######################
# Build the SPDE model
spde.ml <- inla.spde2.pcmatern(
 # Mesh and smoothness parameter
 mesh = mesh, 
 alpha = 2,
 # P(practic.range < 0.3) = 0.5
 prior.range = c(0.3, 0.5), # 50% certainty that range is below ~5km
 # P(sigma > 1) = 0.01
 prior.sigma = c(1, 0.05) # variance
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
 select(Fire_ID, first_obs_date, grid_index,
        aspen_ba_live, aspen_tpp_live, aspen_ht_live, 
        aspen_dia_live, aspen_ba_pr, aspen_tpp_pr, 
        fortypnm_gp, fortyp_pct, lf_canopy, ba_dead_total,
        erc_dv, vpd, vs, elev, slope, northness, tpi,
        tpp_live_total, tpp_dead_total,
        grid_aspen, day_prop, overlap, dist_to_perim)

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
rm(grid_sf,X,coords_mat,field.idx,mesh)
gc()


##########################
# update the model formula
mf.cbi.sp <- update(
 mf.cbi, . ~ 1 + . + 
  f(mesh.idx, model = spde.ml) # spatial process model
)
# fit the model
ml.cbi.re.sp <- inla(
 mf.cbi.sp, 
 data = inla.stack.data(stack.cbi), ,
 family = "gaussian",
 control.predictor = list(A = inla.stack.A(stack.cbi), compute=T),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
 control.fixed = pr.fixed, # regularization on fixed effects
 control.family = pr.family # on the Gaussian observation
)
summary(ml.cbi.re.sp)
mean(ml.cbi.re.sp$cpo$cpo, na.rm = TRUE)

rm(A,spde.ml,stack.cbi) # tidy up


#=================MODEL STATEMENTS=================#

# compute the exponentiated effects
exp.cbi <- ml.cbi.re.sp$summary.fixed %>%
 rownames_to_column(var = "parameter") %>%
 mutate(
  exp_mean = exp(mean) - 1,  # Convert log(FRP) effect to % difference
  lower_ci = exp(`0.025quant`) - 1,  # 2.5% CI bound
  upper_ci = exp(`0.975quant`) - 1   # 97.5% CI bound
 )
# save this table
write_csv(exp.cbi, paste0(maindir,"data/tabular/mod/results/INLA_exp_CBI_aspenEffect.csv"))
# check results
exp.cbi%>%select(parameter,exp_mean,lower_ci,upper_ci)


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
 filter(!parameter %in% c("(Intercept)","day_prop","overlap")) %>%  
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
   str_detect(parameter, ":vpd") ~ "VPD-mediated",  
   str_detect(parameter, ":aspen_ba_pr") & 
    !str_detect(parameter, ":vpd")  ~ "Aspen proportion",
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
glimpse(tidy.effects.cbi)

# reorder the factors
tidy.effects.cbi <- tidy.effects.cbi %>%
 mutate(
  effect_order = case_when(
   !is.na(fortyp) ~ 2,  # Species-specific effects
   TRUE ~ 2              # Global effects
  ),
  effect = factor(effect, levels = tidy.effects.cbi %>%
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

###########
# color map
color_map = c(
 "Aspen proportion" = "#E69F00",  # Keep orange for direct effect
 "VPD-mediated" = "gray80"  # Keep gray for mediation effect
)

#######################
# make the ridge plot #

# plot it
(cbi.p1 <- tidy.effects.cbi %>%
  mutate(fill_species = factor(fill_species, levels = aspen_order)) %>%
  arrange(match(effect, c("VPD-mediated", "Aspen proportion"))) %>%
  ggplot(., aes(
   x = x, y = fill_species, height = y, 
   fill = factor(effect, levels = c("Aspen proportion", "VPD-mediated")),
   alpha = factor(effect, levels = c("Aspen proportion", "VPD-mediated")))
  ) + 
  geom_density_ridges(stat = "identity", scale = 1.5, show.legend=T) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(
   x = "",
   y = "Predominant Forest Type"
  ) +
  scale_fill_manual(
   values = color_map,
   guide = guide_legend(title = "Effect type"),
   labels = c(
    "Aspen proportion" = "Aspen dominance (basal area)",
    "VPD-mediated" = "VPD × aspen dominance"
   )
  ) +
  scale_alpha_manual(
   values = c(
    "VPD-mediated" = 0.6,
    "Aspen proportion" = 0.95
   ),
   guide = "none"  # Hide alpha from legend
  ) +
  # add a subplot label
  annotate("text", x = -0.48, y = 8.8,
           label = expression(bold("(B)") ~ "CBIbc"),
           size = 4, hjust = 0.2) +
  coord_cartesian(xlim=c(-0.52,0.66)) +
  theme_classic() +
  theme(
   axis.text.y = element_text(angle = 0, hjust = 1, size=9),
   axis.text.x = element_text(angle = 0, hjust = 0, size=9),
   axis.title.y = element_text(size = 10, margin = margin(r = 12)),
   axis.title.x = element_blank(),
   legend.position = c(0.85, 0.88),
   legend.background = element_rect(
    fill = scales::alpha("white", 0.6), 
    color = NA, linewidth = 1),
   legend.title = element_text(size = 10),
   legend.text = element_text(size = 9)
  ))

#############
# version two

###########
# color map
color_map = c(
 "Aspen proportion" = "#E69F00",  # Keep orange for direct effect
 "VPD-mediated" = "#E69F00"  # Keep gray for mediation effect
)

(cbi.p1.1 <- tidy.effects.cbi %>%
 mutate(fill_species = factor(fill_species, levels = aspen_order)) %>%
 arrange(match(effect, c("VPD-mediated", "Aspen proportion"))) %>%
 ggplot(., aes(x = x, y = fill_species, height = y,
               fill = effect,
               alpha = effect,
               color = effect
           )) +
 geom_density_ridges(stat = "identity", scale = 1.5, aes(alpha = effect)) +
 # geom_density_ridges(stat = "identity", scale = 1.5, color = NA, alpha = 0) +
 geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
 labs(
  x = "Effect size",
  y = "Predominant Forest Type",
  fill = "Effect"
 ) +
  scale_fill_manual(
   values = c(
    "Aspen proportion" = "#E69F00",
    "VPD-mediated" = scales::alpha("#E69F00", 0.3)
   ),
   labels = c(
    "Aspen proportion" = "Aspen dominance (basal area)",
    "VPD-mediated" = "VPD × aspen dominance"
   ),
   guide = guide_legend(
    override.aes = list(
     fill = c("#E69F00", scales::alpha("#E69F00", 0.3)),
     alpha = c(1, 0.3),  # Match transparency
     color = c(
      scales::alpha("black", 0.7),
      scales::alpha("#E69F00", 0.3))
    ),
    order = 1  # Ensure this legend appears first
   )
 ) +
 scale_alpha_manual(
  values = c(
   "VPD-mediated" = 0.4,
   "Aspen proportion" = 0.98
  ),
  guide = "none"
 ) +
 scale_color_manual(
  values = c(
   "Aspen proportion" = scales::alpha("black", 0.7),
   "VPD-mediated" = scales::alpha("#E69F00", 0.5)
  ),
  guide = "none"
 ) +
  # add a subplot label
  annotate(
   "text", x = -0.48, y = 8.8,
   label = expression(bold("(B)") ~ "CBIbc"),
   size = 4, hjust = 0.2
  ) +
  theme_classic() +
  coord_cartesian(xlim=c(-0.52,0.66)) +
  theme(
   axis.text.y = element_text(angle = 0, hjust = 1, size=9),
   axis.text.x = element_text(angle = 0, hjust = 0, size=9),
   axis.title.y = element_text(size = 10, margin = margin(r = 12)),
   axis.title.x = element_blank(),
   legend.position = c(0.80, 0.88),
   legend.background = element_rect(
    fill = scales::alpha("white", 0.6), 
    color = NA, linewidth = 1),
   legend.title = element_text(size = 10),
   legend.text = element_text(size = 9)
  ))
# save the plot
out_png <- paste0(maindir, 'figures/INLA_AspenEffect_CBI_fortyp.png')
ggsave(out_png, plot = cbi.p1.1, dpi = 500, width = 7, height = 5, bg = 'white')


########################################
# make a version without VPD mediation #
(cbi.p1.2 <- tidy.effects.cbi %>%
  mutate(fill_species = factor(fill_species, levels = aspen_order)) %>%
  filter(effect == "Aspen proportion") %>%
  ggplot(., aes(x = x, y = fill_species, height = y,
                fill = effect,
                alpha = effect,
                color = effect
  )) +
  geom_density_ridges(stat = "identity", scale = 1.5, aes(alpha = effect)) +
  # geom_density_ridges(stat = "identity", scale = 1.5, color = NA, alpha = 0) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(
   x = "Effect size",
   y = "Predominant Forest Type",
   fill = "Effect"
  ) +
  scale_fill_manual(
   values = c("Aspen proportion" = "#E69F00"),
   labels = c("Aspen proportion" = "Aspen dominance (basal area)"),
   guide = guide_legend(
    override.aes = list(
     fill = c("#E69F00"),
     alpha = c(1),  # Match transparency
     color = c(scales::alpha("black", 0.7)),
    order = 1  # Ensure this legend appears first
   )
  )) +
  scale_alpha_manual(
   values = c("Aspen proportion" = 0.98),
   guide = "none"
  ) +
  scale_color_manual(
   values = c(
    "Aspen proportion" = scales::alpha("black", 0.7)
   ),
   guide = "none"
  ) +
  # add a subplot label
  annotate(
   "text", x = -0.48, y = 8.8,
   label = expression(bold("(B)") ~ "CBIbc"),
   size = 4, hjust = 0.2
  ) +
  theme_classic() +
  coord_cartesian(xlim=c(-0.52,0.66)) +
  theme(
   axis.text.y = element_text(angle = 0, hjust = 1, size=9),
   axis.text.x = element_text(angle = 0, hjust = 0, size=9),
   axis.title.y = element_text(size = 10, margin = margin(r = 12)),
   axis.title.x = element_blank(),
   legend.position = c(0.80, 0.88),
   legend.background = element_rect(
    fill = scales::alpha("white", 0.6), 
    color = NA, linewidth = 1),
   legend.title = element_text(size = 10),
   legend.text = element_text(size = 9)
  ))
# save the plot
out_png <- paste0(maindir, 'figures/INLA_AspenEffect_CBI_fortyp_noVPD.png')
ggsave(out_png, plot = cbi.p1.2, dpi = 500, width = 7, height = 5, bg = 'white')

# #######################################
# # Version 2: no Gambel oak or White fir
# (cbi.p1.1 <- tidy.effects.cbi %>%
#   filter(!fill_species %in% c("Gambel oak", "White fir")) %>%
#   ggplot(., aes(x = x, y = fill_species, height = y,
#                 fill = factor(effect, levels = c("VPD-mediated","Aspen proportion")),
#                 alpha = factor(effect, levels = c("VPD-mediated","Aspen proportion")))) +
#   geom_density_ridges(
#    stat = "identity", scale = 1.5, show.legend=T) +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
#   labs(
#    x = "",
#    y = "Predominant Forest Type",
#    fill = "Model effect"
#   ) +
#   scale_fill_manual(
#    values = color_map,
#    guide = guide_legend(title = "Model effect")
#   ) +
#   scale_alpha_manual(
#    values = c(
#     "Aspen proportion" = 0.95,
#     "VPD-mediated" = 0.6
#    ),
#    guide = "none"  # Hide alpha from legend
#   ) +
#   # coord_cartesian(xlim=c(-0.268,0.51)) +
#   theme_classic() +
#   theme(
#    axis.text.y = element_text(angle = 0, hjust = 1, size=9),
#    axis.text.x = element_text(angle = 0, hjust = 0, size=9),
#    axis.title.y = element_text(size = 10, margin = margin(r = 12)),
#    axis.title.x = element_blank(),
#    legend.position = c(0.85, 0.88),
#    legend.background = element_rect(
#     fill = scales::alpha("white", 0.6), 
#     color = NA, size = 1),
#    legend.title = element_text(size = 10),
#    legend.text = element_text(size = 9)
#   ))


########################################
# Make a combined figure (FRP and CBI) #
# remove the legend from one of them
cbi.p1.1 <- cbi.p1.1 +
 theme(legend.position="none")
# merge them
frp.cbi.p1 <- frp.p1.1 / cbi.p1.1 # "/" stacks, "|" side-by-side
frp.cbi.p1
# Save the plot
out_png <- paste0(maindir, 'figures/INLA_AspenEffect_FRP-CBI_fortyp-VPD_panel.png')
ggsave(out_png, plot = frp.cbi.p1, dpi = 500, width = 7, height = 5, bg = 'white')

# without VPD
cbi.p1.2 <- cbi.p1.2 +
 theme(legend.position="none")
# merge them
frp.cbi.p2 <- frp.p1.2 / cbi.p1.2 # "/" stacks, "|" side-by-side
frp.cbi.p2
# Save the plot
out_png <- paste0(maindir, 'figures/INLA_AspenEffect_FRP-CBI_fortyp-VPD_panel_noVPD.png')
ggsave(out_png, plot = frp.cbi.p2, dpi = 500, width = 7, height = 5, bg = 'white')



###################################
# OPTION TWO: COMBINED RIDGE PLOT #

# Create the combined forest type relative to aspen
tidy.effects.frp$response <- "FRPc"
tidy.effects.cbi$response <- "CBIbc"
fortyp_effects <- bind_rows(tidy.effects.frp, tidy.effects.cbi) %>%
 filter(str_detect(parameter, "fortypnm_gp"),
        !str_detect(parameter, ":fortyp_pct")) %>%
 mutate(
  group = case_when(
   str_detect(parameter, "vpd:") ~ "VPD-mediated",
   TRUE ~ "Aspen live basal area"  # Default for all other fixed effects
  )) %>%
 group_by(fill_species) %>%
 mutate(mean_effect = mean(x, na.rm = TRUE)) %>%
 ungroup() %>%
 mutate(fill_species = fct_reorder(fill_species, -mean_effect))%>%
 # reorder the effects so main effect is drawn on top
 mutate(
  group = factor(group, levels = c("VPD-mediated", "Aspen live basal area")),
  fill_cat = factor(
   paste0(group, " (", response, ")"),
   levels = c(
    "VPD-mediated (CBIbc)", "VPD-mediated (FRPc)",
    "Aspen live basal area (CBIbc)", "Aspen live basal area (FRPc)"
   ))
 )
glimpse(fortyp_effects)

# order the forest types from low->high on the main effect
forest_order <- fortyp_effects %>%
 filter(fill_cat == "Aspen live basal area (FRPc)") %>%
 group_by(fill_species) %>%
 summarize(mean_effect = mean(x, na.rm = TRUE)) %>%
 arrange(mean_effect)
# apply the ordering
fortyp_effects$fill_species <- factor(
 fortyp_effects$fill_species, 
 levels = rev(forest_order$fill_species)
)

##########################
# Plot with VPD-mediation
(p3 <- ggplot(fortyp_effects, 
                aes(x = x, y = fill_species, height = y,
                    fill = fill_cat,
                    alpha = fill_cat,
                    color = fill_cat)) +
 geom_density_ridges(stat = "identity", scale = 1.5) +
 geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
 labs(
  x = "Effect size",
  y = "Predominant Forest Type",
  fill = "Effect Type",
 ) +
 scale_fill_manual(
  values = c(
   "Aspen live basal area (FRPc)" = "#FEB24C",
   "VPD-mediated (FRPc)" = scales::alpha("#FEB24C", 0.3),
   "Aspen live basal area (CBIbc)" = "#800026",
   "VPD-mediated (CBIbc)" = scales::alpha("#800026", 0.3)
  ),
  labels = c(
   "Aspen live basal area (FRPc)", "VPD-mediated (FRPc)",
   "Aspen live basal area (CBIbc)", "VPD-mediated (CBIbc)"
  ),
  guide = guide_legend(
   override.aes = list(
    fill = c("#FEB24C", scales::alpha("#FEB24C", 0.3),
             "#800026", scales::alpha("#800026", 0.3)),
    alpha = c(1, 0.4, 1, 0.4),  # Match transparency
    color = c(
     scales::alpha("black", 0.7),
     scales::alpha("#FEB24C", 0.4),
     scales::alpha("black", 0.7),
     scales::alpha("#800026", 0.4))
   ),
   order = 1  # Ensure this legend appears first
  )
 ) +
 scale_alpha_manual(
  values = c(
   "Aspen live basal area (FRPc)" = 1,
   "VPD-mediated (FRPc)" = 0.4,
   "Aspen live basal area (CBIbc)" = 1,
   "VPD-mediated (CBIbc)" = 0.4
  ),
  guide = "none"
 ) +
 scale_color_manual(
  values = c(
   "Aspen live basal area (FRPc)" = scales::alpha("black", 0.7),
   "VPD-mediated (FRPc)" = scales::alpha("#FEB24C", 0.4),
   "Aspen live basal area (CBIbc)" = scales::alpha("black", 0.7),
   "VPD-mediated (CBIbc)" = scales::alpha("#800026", 0.4)
  ),
  guide = "none"
 ) +
 # coord_cartesian(xlim=c(-0.38,0.44)) +
 theme_classic() +
 theme(axis.text.y = element_text(angle = 0, hjust = 1, size=9),
       axis.text.x = element_text(angle = 0, hjust = 0, size=9),
       axis.title.y = element_text(size = 10, margin = margin(r = 12)),
       axis.title.x = element_text(size = 10, margin = margin(t = 12)),
       legend.position = c(0.82, 0.84),
       legend.background = element_rect(
        fill = scales::alpha("white", 0.4),
        color = NA, size = 0.8),
       legend.title = element_text(size = 9),
       legend.text = element_text(size = 8)))
# Save the plot
out_png <- paste0(maindir, 'figures/INLA_AspenEffect_FRP-CBI_fortyp-VPD.png.png')
ggsave(out_png, plot = p3, dpi = 500, width = 7, height = 4, bg = 'white')
