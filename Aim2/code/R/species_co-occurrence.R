##################
library(tidyverse)
library(sf)
library(cooccur)

# Environment variables
maindir <- '/Users/max/Library/CloudStorage/OneDrive-Personal/mcook/aspen-fire/Aim2/'


#=========Prep the grid data=========#

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
  # interaction term between forest type and species co-occurring
  spp_int = interaction(fortypnm_gp, species_gp_n)
 ) %>%
 # be sure there are no duplicate rows for grid/fire/species
 distinct(grid_index, Fire_ID, species_gp_n, .keep_all = TRUE) # remove duplicates


#============PREP THE SPECIES DATA=============#

# calculate the species presence / absence data
# Define threshold for presence
dt = 0.10 # 5% of live BA or TPP

# Calculate presence/absence for each species group
# Assuming columns like species group contributions: balive_species, tpa_species, etc.
pres_abs <- grid_tm %>%
 select(grid_index, species_gp_n, ba_live_pr, tpp_live_pr) %>%
 mutate(
  grid_index = as.character(grid_index),
  dominance = ifelse(is.na(ba_live_pr) | ba_live_pr < dt, 0, 1),
  abundance = ifelse(is.na(tpp_live_pr) | tpp_live_pr < dt, 0, 1)
 ) %>%
 distinct(grid_index, species_gp_n, abundance)
# check on the table
head(pres_abs)

# Pivot the data to wide format for presence based on abundance
pmat.df <- pres_abs %>%
 select(grid_index, species_gp_n, abundance) %>%
 pivot_wider(
  names_from = species_gp_n, # Species become columns
  values_from = abundance, # Presence values (1/0)
  values_fill = 0
 ) %>%
 column_to_rownames("grid_index")  # Set grid_id as rownames
# View the matrix
print(head(pmat.df))


#=========MODEL FIT=========#

pmat <- as.matrix(pmat.df)

# run the co-occurrence model
coo.model <- cooccur(t(pmat), type = "spp_site", thresh = TRUE, spp_names = TRUE)
summary(coo.model)  # Overview of results

# Built-in plots:
plot(coo.model)
pair.profile(coo.model)


#=========PLOTTING OBSERVED PROBABILITIES=========#

# Extract probabilities of co-occurrence
probs <- as.data.frame(prob.table(coo.model))
head(probs)

# map numeric code to species name
spps <- c("aspen", "lodgepole", "mixed_conifer", "piñon_juniper", "ponderosa", "spruce_fir")

# Create a lookup table
species_lookup <- data.frame(
 id = 1:length(spps),  # Numeric IDs
 name = spps           # Species names
)

# Plot a heatmap of species co-occurrence probabilities
ggplot(probs, aes(x = sp1, y = sp2, fill = prob_cooccur)) +
 geom_tile(color = "white") +
 geom_text(aes(label = ifelse(prob_cooccur > 0.01, sprintf("%.2f", prob_cooccur), "")),
           color = "black", size = 3) +
 scale_fill_viridis_c(name = "Observed\nProbability", option = "viridis") +
 # Map species names to numeric IDs for axis labels
 scale_x_continuous(
  breaks = unique(probs$sp1),
  labels = unique(probs$sp1_name)
 ) +
 scale_y_continuous(
  breaks = unique(probs$sp2),
  labels = unique(probs$sp2_name)
 ) +
 labs(
  x = "Species A",
  y = "Species B"
 ) +
 theme_minimal() +
 theme(
  axis.text.x = element_text(angle = 45, hjust = 1),
  axis.text = element_text(size = 10),
  panel.grid = element_blank()
 )


#=========PLOTTING RELATIONSHIPS=========#

# Classify co-occurrence relationships
probs_cl <- probs %>%
 mutate(
  relationship = case_when(
   p_lt < 0.05 ~ "negative",
   p_gt < 0.05 ~ "positive",
   TRUE ~ "random"
  )
 )

# Reorder the species
spp_order <- unique(c(probs$sp1_name, probs_cl$sp2_name))
probs_cl <- probs_cl %>%
 mutate(
  sp1_name = factor(sp1_name, levels = spp_order),
  sp2_name = factor(sp2_name, levels = spp_order)
 )

# Plot using ggplot
ggplot(probs_cl, aes(x = sp1_name, y = sp2_name, fill = relationship)) +
 geom_tile(color = "white") +
 scale_fill_manual(values = c("positive" = "skyblue", 
                              "negative" = "orange", 
                              "random" = "grey80"),
                   name = "Relationship") +
 labs(
  x = "Species A",
  y = "Species B",
  title = "Species Co-occurrence Matrix"
 ) +
 theme_minimal() +
 theme(
  axis.text.x = element_text(angle = 45, hjust = 1),
  panel.grid = element_blank()
 )


#################################################
# Calculate observed probabilities for pure grids
# Add a column to mark pure stands in the grid data
pure_pa <- pmat.df %>%
 rowwise() %>%
 mutate(
  pure = ifelse(sum(c_across(where(is.numeric))) == 1, 1, 0)  # Only one species present
 ) %>%
 ungroup()
head(pure_pa)

n_pure_grids <- sum(pure_pa$pure == 1)
total_grids <- nrow(pure_pa)
pure_fraction <- n_pure_grids / total_grids

cat("Number of pure grids:", n_pure_grids, "\n")
cat("Fraction of pure grids:", pure_fraction, "\n")

# Identify species in pure stands
pure_species <- colnames(pmat)[colSums(pmat == 1 & pure_pa$pure == 1) > 0]

# Initialize a results table for pure stands
pure_probs <- tibble(
 species = pure_species,
 obs_pure = numeric(length(pure_species)),
 exp_pure = numeric(length(pure_species)),
 prob_pure = numeric(length(pure_species))
)
# Observed probability: Fraction of grids where the species is pure
pure_probs <- pure_probs %>%
 rowwise() %>%
 mutate(
  n_grids = sum(pmat.df[[species]] == 1 & pure_pa$pure == 1),  # Count grids where species is pure
  obs_pure = n_grids / nrow(pmat.df),                         # Fraction of grids where species is pure
  exp_pure = (sum(pmat.df[[species]]) / nrow(pmat.df))^2,     # Marginal probability squared
  prob_pure = ifelse(exp_pure > 0, obs_pure / exp_pure, NA)   # Avoid division by zero
 ) %>%
 ungroup() %>%
 mutate(
  species = paste0(species, ".pure")  # Label species with ".pure"
 )
head(pure_probs)


#############################################################################
# Identify the top two most dominant species by grid using abundance
# Create the species pair factor for each grid or label species:pure
df <- grid_tm %>%
 select(grid_index, fortypnm_gp, species_gp_n,
        ba_live_pr, tpp_live_pr) %>%
 filter(
  tpp_live_pr >= 0.10
 )

# Identify top two species by live basal area
top_spps <- df %>%
 mutate(species_gp_n = as.character(species_gp_n)) %>%
 group_by(grid_index) %>%
 arrange(desc(tpp_live_pr)) %>%  # Sort by live basal area
 slice_head(n = 2) %>%         # Select the top 2 species
 summarise(                    # Summarize into a single row per grid
  spp_1 = first(species_gp_n),         # Most dominant species
  spp_2 = nth(species_gp_n, 2, default = first(species_gp_n))  # Second most dominant species
 ) %>%
 mutate(
  spp_2 = ifelse(is.na(spp_2), spp_1, spp_2),  # Handle pure stands by filling spp_2 with spp_1
  spp_pair = if_else(
   spp_1 == spp_2,
   paste0(spp_1, ".pure"),  # Label pure stands
   paste(spp_1, spp_2, sep = ".")  # Label mixed stands
  )
 ) %>%
 # join to the probability of co-occurrence
 left_join(
  probs %>% 
   mutate(spp_pair = interaction(sp1_name, sp2_name)) %>%
   select(spp_pair, obs_cooccur, prob_cooccur),
  by = "spp_pair"
 ) %>%
 left_join(
  pure_probs %>%
   rename(spp_pair = species) %>%
   select(spp_pair, n_grids, prob_pure),  # Include n_grids and prob_pure for pure stands
  by = "spp_pair"
 ) %>%
 mutate(
  # Use prob_pure and n_grids for pure stands; retain prob_cooccur for mixed stands
  prob_cooccur = if_else(!is.na(prob_pure), prob_pure, prob_cooccur),
  obs_cooccur = if_else(!is.na(n_grids), n_grids, obs_cooccur)
 ) %>%
 select(-c(prob_pure, n_grids))
head(top_spps)

# write this file out
out_csv <- paste0(maindir,'data/tabular/mod/spp_cooccur_top_spp.csv')
write.csv(top_spps, out_csv)

