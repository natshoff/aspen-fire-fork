##################
library(tidyverse)
library(sf)
library(cooccur)

# Environment variables
maindir <- '/Users/max/Library/CloudStorage/OneDrive-Personal/mcook/aspen-fire/Aim2/'


#=========Prep the grid data=========#

fp <- paste0(maindir,'data/tabular/mod/gridstats_fortypnm_gp_tm_ct_frp-cbi.csv')
grid <- read_csv(fp)
glimpse(grid)

# calculate the species presence / absence data
# Define threshold for presence
dt = 0.01 # 1% of live basal area

# Calculate presence/absence for each species group
# Assuming columns like species group contributions: balive_species, tpa_species, etc.
pres_abs <- grid %>%
 select(grid_index, species_gp_n, ba_live_pr, tpp_live_pr) %>%
 mutate(
  grid_index = as.character(grid_index),
  dominance = ifelse(is.na(ba_live_pr) | ba_live_pr < dt, 0, 1),
  abundance = ifelse(is.na(tpp_live_pr) | tpp_live_pr < dt, 0, 1)
 ) %>%
 distinct(grid_index, species_gp_n, dominance, abundance)
# check on the table
head(pres_abs)

# Pivot the data to wide format for presence based on abundance
pmat <- pres_abs %>%
 select(grid_index, species_gp_n, abundance) %>%
 pivot_wider(
  names_from = species_gp_n, # Species become columns
  values_from = abundance, # Presence values (1/0)
  values_fill = 0
 ) %>%
 column_to_rownames("grid_index")  # Set grid_id as rownames
# View the matrix
print(head(pmat))


#=========MODEL FIT=========#

pmat <- as.matrix(pmat)

# run the co-occurrence model
coo.model <- cooccur(t(pmat), type = "spp_site", thresh = TRUE, spp_names = TRUE)
summary(coo.model)  # Overview of results

# Built-in plots:
plot(coo.model)
pair.profile(coo.model)

#=========PLOTTING PROBABILITIES=========#

# Extract probabilities of co-occurrence
probs <- as.data.frame(prob.table(coo.model))
head(probs)

# Tidy the data frame for plotting
# Ensure species names are factors in the same order
probs <- probs %>%
 mutate(
  sp1_name = factor(sp1_name, levels = sort(unique(c(sp1_name, sp2_name)))),
  sp2_name = factor(sp2_name, levels = sort(unique(c(sp1_name, sp2_name))))
 ) %>%
 # Filter to include only the lower triangle of the matrix
 filter(as.numeric(sp1_name) <= as.numeric(sp2_name))

# map numeric code to species name
spps <- c("Aspen", "Lodgepole", "Mixed-conifer", "Piñon-juniper", "Ponderosa", "Spruce-fir")

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

##############
# Network plot

# Create a data frame of the nodes in the network. 
nodes <- data.frame(id = 1:nrow(pmat),
                    label = rownames(pmat),
                    color = "#606482",
                    shadow = TRUE)