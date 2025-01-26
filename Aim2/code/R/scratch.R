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