
# ##################################################################
# # 4. Baseline + Temporal + Fire-level + Species-pair latent effect
# # requires building a precision matrix using species pairs
# # weighted by the observed co-occurrence frequency
# 
# # gather unique species pairs
# spp_pairs <- unique(da$spp_pair)
# n_pairs <- length(spp_pairs)
# # build the precision matrix
# prec_mat <- Matrix(0, nrow = n_pairs, ncol = n_pairs, sparse = TRUE)
# rownames(prec_mat) <- spp_pairs
# colnames(prec_mat) <- spp_pairs
# # populate diagonal with weights
# # appropriate for balancing rare/common pairings
# for (pair in spp_pairs) {
#  wt_norm <- mean(da$wt_norm[da$spp_pair == pair])  # Use mean weight for each spp_pair
#  prec_mat[pair, pair] <- wt_norm
# }
# # check positive-definiteness
# if (!all(eigen(prec_mat)$values > 0)) {
#  diag(prec_mat) <- diag(prec_mat) + 1e-4  # Add a small constant to diagonals
# }
# # ensure symmetry and positive-definiteness
# prec_mat <- Matrix::nearPD(prec_mat, corr = FALSE)$mat
# # match spp_pair levels in 'da' to precision matrix
# da$spp_pair <- factor(da$spp_pair, levels = rownames(prec_mat))
# # Verify alignment
# stopifnot(all(levels(da$spp_pair) %in% rownames(prec_mat)))
# stopifnot(all(rownames(prec_mat) %in% levels(da$spp_pair)))
# # Verify positive semi-definiteness
# isSymmetric(prec_mat)  # Should return TRUE
# all(eigen(prec_mat)$values > 0)  # Should return TRUE
# # Verify factor dimensions
# nrow(prec_mat) == length(unique(da$spp_pair))  # Should return TRUE
# 
# 
# ##########################
# # update the model formula 
# mf.frp.re2.spp <- update(
#  mf.frp.re2, 
#  . ~ 1 + . + 
#  f(spp_pair, model = "generic0", Cmatrix = prec_mat,
#    values = levels(da$spp_pair),
#    hyper = list(
#     prec = list(prior = "pc.prec", param = c(1, 0.01))
#    ))
# )
# 
# # fit the model
# ml.frp.re.spp <- inla(
#  mf.frp.re2.spp, data = da,
#  family = "gaussian",
#  control.predictor = list(compute=T),
#  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
#  control.fixed = list(
#   prec = list(prior = "pc.prec", param = c(1, 0.5))
#  ) 
# )
# summary(ml.frp.re.spp)
# 
# # check on predictive power of the random effects models
# mean(ml.frp.re2$cpo$cpo, na.rm = TRUE)
# mean(ml.frp.re.spp$cpo$cpo, na.rm = TRUE)
# 
# rm(prec_mat)
# gc()

##################################################
# Examine the semivariogram for spatial dependence
# Function to compute semivariogram for each fire
compute_variogram <- function(fire_df) {
 fire_id <- as.character(unique(fire_df$fire_id))  # Ensure fire_id is treated as character
 coordinates(fire_df) <- ~x + y  # Set spatial coordinates
 vgm_model <- variogram(frp_day ~ 1, data = fire_df)  # Compute semivariogram
 
 vfit <- tryCatch(
  fit.variogram(vgm_model, model = vgm("Sph")),  # Fit spherical model
  error = function(e) {
   return(data.frame(fire_id = fire_id, range = NA, status = "failed"))  # Mark failure
  }
 )
 
 if (!is.null(vfit) && "range" %in% names(vfit)) {
  range_val <- vfit$range[2]
  
  # Filter out invalid ranges (negative or excessively large)
  if (!is.na(range_val) && range_val > 0 & range_val < 5) {  
   return(data.frame(fire_id = fire_id, range = range_val, status = "valid"))
  } else {
   return(data.frame(fire_id = fire_id, range = range_val, status = "extreme"))  # Mark extreme values
  }
 } 
 return(data.frame(fire_id = fire_id, range = NA, status = "failed"))  # If something goes wrong
}

# Compute variograms for each fire and classify results
fire_results <- da %>%
 group_split(fire_id) %>%
 map_dfr(compute_variogram) %>%
 mutate(range_km = range * 111.32)  # Convert range to km

# Filter failed and extreme-range fires
failed_fires <- fire_results %>% filter(status == "failed") %>% pull(fire_id) %>% unique()
extreme_fires <- fire_results %>% filter(status == "extreme") %>% pull(fire_id) %>% unique()

# Print diagnostic information
cat("Fires that failed semivariogram fitting:", length(failed_fires), "\n")
print(failed_fires)

cat("\nFires with extreme range estimates (>5 degrees):", length(extreme_fires), "\n")
print(extreme_fires)

# Save results for review
write.csv(fire_results, "fire_results.csv", row.names = FALSE)

# Summary of within-fire spatial dependence (only valid values)
fire_ranges <- fire_results %>% filter(status == "valid")
summary(fire_ranges$range_km)
# Boxplot for visual inspection
ggplot(fire_ranges, aes(y = range_km)) +
 geom_boxplot(fill = "orange", alpha = 0.6) +
 labs(y = "Semivariogram Range (km)", title = "Distribution of Within-Fire Spatial Dependence") +
 theme_minimal()


# compute the average distance for different KNN
k_values <- c(3, 5, 7, 9)
avg_dist <- sapply(k_values, function(k) {
 nbs <- knearneigh(coords, k = k, longlat = TRUE)
 nb <- knn2nb(nbs, row.names = grid_sf$grid_index, sym = TRUE)
 mean(unlist(nbdists(nb, coords)))  # Average neighbor distance
})
data.frame(k = k_values, avg_distance = avg_dist)


#################################################
# Extract relevant metrics from the spatial model
sp_metrics <- data.frame(
 Model = "Baseline + Spatial",
 DIC = model_bl.frp.sp$dic$dic,
 WAIC = model_bl.frp.sp$waic$waic,
 Marginal_Log_Likelihood = model_bl.frp.sp$mlik[1,1],
 Num_Effective_Params = model_bl.frp.sp$waic$p.eff
)
# Add to the existing comparison table
model_comparison <- bind_rows(model_comparison, sp_metrics) %>% arrange(WAIC)
print(model_comparison)

# Tidy up !
rm(model_bl.frp, model_bl.frp.re, model_bl.frp.re2)
gc()

#################################################
# Check out the fire-level variability in effects
# Extract fire-level random effects
fire_effects <- model_bl.frp.sp$summary.random$fire_id %>%
 dplyr::select(ID, mean, sd, `0.025quant`, `0.5quant`, `0.975quant`) %>%
 rename(fire_id = ID)

# Merge with original fire metadata for plotting
fire_summary <- da %>%
 select(fire_id, fire_year, season) %>%
 distinct() %>%
 left_join(fire_effects, by = "fire_id")

# Plot fire-level random effects
ggplot(fire_summary, aes(x = reorder(fire_id, mean), y = mean)) +
 geom_point(aes(color = as.factor(fire_year))) +
 geom_errorbar(aes(ymin = `0.025quant`, ymax = `0.975quant`), width = 0.2) +
 labs(title = "Fire-level Random Effects",
      x = "Fire ID",
      y = "Estimated Fire Effect",
      color = "Fire Year") +
 theme_minimal() +
 coord_flip()  # To make the plot more readable

# Merge fire-level random effects with forest type data
da_forest_fire <- da %>%
 left_join(fire_effects, by = "fire_id")

# Boxplot to see distribution across fires
ggplot(da_forest_fire, aes(x = fortypnm_gp, y = frp, fill = fortypnm_gp)) +
 geom_boxplot(outlier.shape = NA, alpha = 0.6) +
 labs(title = "FRP Distribution by Forest Type Across Fires",
      x = "Forest Type",
      y = "Fire Radiative Power (FRP)",
      fill = "Forest Type") +
 theme_minimal() +
 coord_flip()

# Combine fixed effect estimates from models with and without fire effects
fixed_effects <- bind_rows(
 model_bl.frp$summary.fixed %>% mutate(Model = "No Fire Effect"),
 model_bl.frp.sp$summary.fixed %>% mutate(Model = "With Fire Effect")
) %>%
 filter(str_detect(row.names(.), "fortypnm_gp")) %>%
 mutate(ForestType = str_remove(row.names(.), "fortypnm_gp"))

# Plot comparison of forest type fixed effects
ggplot(fixed_effects, aes(x = ForestType, y = mean, fill = Model)) +
 geom_bar(stat = "identity", position = position_dodge()) +
 geom_errorbar(aes(ymin = `0.025quant`, ymax = `0.975quant`), position = position_dodge(0.9), width = 0.3) +
 labs(title = "Forest Type Effects Before and After Adding Fire-Level Random Effect",
      x = "Forest Type",
      y = "Effect Estimate",
      fill = "Model") +
 theme_minimal() +
 coord_flip()





# ################################
# # Co-occurrence precision matrix

# ####################################
# # Extract fixed effects for spp_pair
# fixed_effects <- as.data.frame(model_tm.re.sp$summary.fixed)
# fixed_effects$spp_pair <- rownames(fixed_effects)
# 
# # Filter spp_pair effects
# spp_pair_effects <- fixed_effects %>%
#  filter(str_detect(spp_pair, "spp_pair")) %>%
#  mutate(
#   spp_pair = str_remove(spp_pair, "spp_pair"),
#   lower_ci = mean - 1.96 * sd,
#   upper_ci = mean + 1.96 * sd,
#   significant = (lower_ci > 0 | upper_ci < 0)  # Significant if CI does not overlap 0
#  ) %>%
#  arrange(mean)
# 
# # Filter only significant effects
# significant_spp_pairs <- spp_pair_effects %>%
#  filter(significant)
# 
# # Create the ridge plot with bar heights representing the effect size
# ggplot(significant_spp_pairs, aes(x = mean, y = reorder(spp_pair, mean), fill = mean > 0)) +
#  geom_col(
#   aes(x = mean),  # Use the actual mean value
#   position = "identity",
#   width = 0.8,    # Bar width for the ridge
#   alpha = 0.7
#  ) +
#  geom_point(aes(x = mean), color = "black", size = 2) +
#  geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci), height = 0.2, color = "black") +
#  scale_fill_manual(values = c("TRUE" = "skyblue", "FALSE" = "coral"), guide = FALSE) +
#  labs(
#   x = "Effect Size (Mean ± 95% CI)",
#   y = "Species Pair",
#   title = "Significant Effects of Species Co-occurrence (spp_pair)"
#  ) +
#  theme_minimal() +
#  theme(
#   legend.position = "none",
#   axis.text.y = element_text(size = 8),
#   axis.text.x = element_text(size = 10),
#   axis.title = element_text(size = 12)
#  )


# spp_prec_mat <- sparseMatrix(
#  i = 1:nrow(spp_coo_gp),
#  j = 1:nrow(spp_coo_gp),
#  x = spp_coo_gp$wt
# )
# # assign col and row names
# rownames(spp_prec_mat) <- colnames(spp_prec_mat) <- spp_coo_gp$spp_pair
# rownames(spp_prec_mat) <- make.unique(rownames(spp_prec_mat))
# colnames(spp_prec_mat) <- make.unique(colnames(spp_prec_mat))
# # match factor levels to the precision matrix
# da$spp_pair <- factor(da$spp_pair, levels = rownames(spp_prec_mat))
# 
# # verify
# all(levels(da$spp_pair) %in% rownames(spp_prec_mat)) # should return TRUE
# is.factor(da$spp_pair)  # Should return TRUE
# levels(da$spp_pair) == rownames(spp_prec_mat)  # All should be TRUE
# setdiff(levels(da$spp_pair), rownames(spp_prec_mat))  # Should return an empty vector
# setdiff(rownames(spp_prec_mat), levels(da$spp_pair))  # Should also return an empty vector
# isSymmetric(spp_prec_mat)  # Should return TRUE
# all(eigen(spp_prec_mat)$values > 0)  # Should return TRUE
# 
# # subset
# used_levels <- levels(da$spp_pair)
# spp_prec_mat <- spp_prec_mat[used_levels, used_levels]


# # Ensure species pair columns are character strings
# spp_coo$spp_pair <- as.character(spp_coo$spp_pair)
# da$spp_pair <- as.character(da$spp_pair)
# 
# spp_coo <- spp_coo %>%
#  group_by(spp_pair) %>%
#  summarize(prob_cooccur = mean(prob_cooccur, na.rm = TRUE), .groups = "drop")
# 
# # Check for unmatched pairs between `da` and `spp_coo`
# unmatched <- setdiff(da$spp_pair, spp_coo$spp_pair)
# if (length(unmatched) > 0) {
#  cat("Unmatched pairs:\n")
#  print(unmatched)
# } else {
#  cat("All pairs match!\n")
# }
# 
# # Extract unique species pairs from spp_coo
# unique_pairs <- unique(spp_coo$spp_pair)
# # Initialize a similarity matrix (identity matrix for simplicity)
# pairwise_sim <- diag(1, nrow = length(unique_pairs))
# rownames(pairwise_sim) <- colnames(pairwise_sim) <- unique_pairs
# 
# # Populate the diagonal entries
# for (i in seq_along(unique_pairs)) {
#  pair <- unique_pairs[i]
#  prob <- spp_coo$prob_cooccur[spp_coo$spp_pair == pair]
#  if (length(prob) == 1 && !is.na(prob)) {
#   pairwise_sim[i, i] <- prob
#  } else {
#   pairwise_sim[i, i] <- 0.01  # Default for missing pairs
#  }
# }
# 
# # Regularize and invert to create the precision matrix
# Q <- solve(pairwise_sim + Diagonal(nrow(pairwise_sim)) * 1e-2)
# rownames(Q) <- colnames(Q) <- unique_pairs
# 
# # reorder to a factor
# da$spp_pair <- factor(da$spp_pair, levels = rownames(Q))
# 
# # check dimensions
# dim(Q) == c(length(levels(da$spp_pair)), length(levels(da$spp_pair)))  # Should return TRUE
# 
# # Check for mismatches between da$spp_pair levels and Q
# if (!all(levels(da$spp_pair) %in% rownames(Q))) {
#  stop("Mismatch between spp_pair levels in `da` and row names of Q")
# }
# 
# print("Levels in da$spp_pair:")
# print(levels(da$spp_pair))
# print("Row names in Q:")
# print(rownames(Q))
# 
# # Check for any mismatches
# mismatched <- setdiff(levels(da$spp_pair), rownames(Q))
# if (length(mismatched) > 0) {
#  print("Mismatched pairs:")
#  print(mismatched)
# } else {
#  print("All levels match!")
# }
# 
# is.factor(da$spp_pair)  # Should return TRUE
# levels(da$spp_pair) == rownames(Q)  # Should all return TRUE
# 
# isSymmetric(Q)  # Should return TRUE
# all(eigen(Q)$values > 0)  # Should return TRUE (positive-definite)
# 
# table(da$spp_pair, exclude = NULL)  # Check frequency of each level
# any(is.na(da$spp_pair))             # Ensure no missing values
# 
# eigenvalues <- eigen(Q)$values
# condition_number <- max(eigenvalues) / min(eigenvalues)
# print(condition_number)
# 
# unused_levels <- setdiff(levels(da$spp_pair), da$spp_pair)
# print(unused_levels)

# #####
# # CBI
# 
# # setup the model formula
# mf.cbi <- CBIbc_p90 ~
#  species_gp_n * (ba_live + qmd_live) + # Species and their structure/diversity
#  erc + vpd + elev + slope + tpi + chili + # Climate and topography predictors
#  f(Fire_ID, model = "iid") + # Fire-level random effect
#  f(grid_index, model = "iid") # Grid-level random effect
# 
# # fit the model
# model_bl_tm.cbi <- inla(
#  mf.cbi, data = grid_tm,
#  family = "gaussian",
#  control.predictor = list(compute=T),
#  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
# )
# summary(model_bl_tm.cbi)



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


#########################################################
# Extract posterior summaries for spp_pair latent effects

# create a lookup table for the species pair names
pair_lookup <- data.frame(
 spp_pair_id = seq_along(levels(da$spp_pair)),  # Numeric indices
 spp_pair_name = levels(da$spp_pair)  # Original species pair names
)

# Extract the posterior distribution of effects
spp_pair_effects <- as.data.frame(model_tm.re.sp$summary.random$spp_pair)
spp_pair_effects$spp_pair_id <- as.numeric(rownames(spp_pair_effects))

# Add confidence intervals and significance
spp_pair_effects <- spp_pair_effects %>%
 mutate(
  lower_ci = `0.025quant`,
  upper_ci = `0.975quant`,
  significant = (lower_ci > 0 | upper_ci < 0),  # Significant if CI does not overlap 0
  spp_pair = rownames(spp_pair_effects)  # Extract spp_pair names
 ) %>%
 arrange(mean) %>%
 filter(significant) %>% # keep significant pairs
 left_join(pair_lookup, by="spp_pair_id") %>%
 mutate(spp_pair_name = factor(spp_pair_name, levels = spp_pair_name[order(mean)]))

# Create the plot
ggplot(spp_pair_effects, aes(x = mean, y = spp_pair_name, fill = mean > 0)) +
 geom_col(
  aes(x = mean), 
  position = "identity",
  width = 0.8,
  alpha = 0.7
 ) +
 geom_point(aes(x = mean), color = "black", size = 2) +
 geom_errorbarh(aes(xmin = `0.025quant`, xmax = `0.975quant`), height = 0.2, color = "black") +
 scale_fill_manual(values = c("TRUE" = "coral", "FALSE" = "skyblue"), guide = "none") +
 labs(
  x = "Effect Size (Mean ± 95% CI)",
  y = "Species Pair",
 ) +
 theme_minimal() +
 theme(
  legend.position = "none",
  axis.text.y = element_text(size = 8),
  axis.text.x = element_text(size = 10),
  axis.title = element_text(size = 12)
 )

# Save the plot
out_png <- paste0(maindir, 'figures/INLA_Spp_COO_LatentEffects.png')
ggsave(out_png, dpi = 500, bg = 'white')

rm(spp_pair_effects, pair_lookup)
gc()
