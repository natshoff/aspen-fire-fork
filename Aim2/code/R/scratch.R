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