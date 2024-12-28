##################
library(tidyverse)
library(sf)
library(cooccur)

# Environment variables
maindir <- '/Users/max/Library/CloudStorage/OneDrive-Personal/mcook/aspen-fire/Aim2/'


#=========Prep the grid data=========#

fp <- paste0(maindir,'data/tabular/mod/viirs_gridstats_treemap.csv')
grid <- read_csv(fp)
glimpse(grid)

# calculate the species presence / absence data
# Define threshold for presence
dt <- 0.01  # Presence if species contributes ≥ 1% of BALIVE or TPA

# Calculate presence/absence for each species group
# Assuming columns like species group contributions: balive_species, tpa_species, etc.
pres_abs <- grid %>%
 select(grid_index, species_gp_n, sp_dominance_ld, sp_abundance_ld) %>%
 mutate(
  grid_index = as.character(grid_index),
  presence_rd = ifelse(is.na(sp_dominance_ld) | sp_dominance_ld < dt, 0, 1),
  presence_ab = ifelse(is.na(sp_abundance_ld) | sp_abundance_ld < dt, 0, 1)
 ) %>%
 distinct(grid_index, species_gp_n, presence_rd)

# check on the table
head(pres_abs)

# Pivot the data to wide format for presence based on dominance
pmat <- pres_abs %>%
 filter(presence_rd == 1) %>%  # Use presence based on dominance threshold
 select(grid_index, species_gp_n, presence_rd) %>%
 pivot_wider(
  names_from = species_gp_n, # Species become columns
  values_from = presence_rd, # Presence values (1/0)
  values_fill = list(presence_rd = 0)
 ) %>%
 column_to_rownames("grid_index")  # Set grid_id as rownames
# View the matrix
print(head(pmat))


#=========Run co-occurrence analysis=========#

pmat <- as.matrix(pmat)
pmat_t <- t(pmat)
rm(pmat)
gc()

# Initialize an empty matrix for probabilities
species <- rownames(pmat_t)
cooccur_prob <- matrix(0, nrow = length(species), ncol = length(species),
                       dimnames = list(species, species))

# Calculate observed and expected co-occurrences
for (i in 1:(nrow(pmat_t) - 1)) {
 for (j in (i + 1):nrow(pmat_t)) {
  # Observed co-occurrence
  observed <- sum(pmat_t[i, ] & pmat_t[j, ])
  
  # Marginal totals
  total_sites <- ncol(pmat_t)
  sp1_count <- sum(pmat_t[i, ])
  sp2_count <- sum(pmat_t[j, ])
  
  # Expected co-occurrence under independence
  expected <- (sp1_count * sp2_count) / total_sites
  
  # Hypergeometric test
  p_value <- phyper(observed - 1, sp1_count, total_sites - sp1_count, sp2_count, lower.tail = FALSE)
  
  # Assign probabilities
  cooccur_prob[i, j] <- p_value
  cooccur_prob[j, i] <- p_value  # Symmetry
 }
}

# Convert to a data frame for visualization (optional)
cooccur_prob_df <- as.data.frame(as.table(cooccur_prob))

# View results
print(head(cooccur_prob_df))

# Exclude self-co-occurrence (optional)
cooccur_prob_df <- cooccur_prob_df %>%
 filter(Var1 != Var2)

# Create a symmetric matrix (optional, if not already symmetric)
cooccur_prob_df <- cooccur_prob_df %>%
 mutate(Freq = as.numeric(Freq)) %>%  # Ensure Freq is numeric
 bind_rows(
  cooccur_prob_df %>%
   rename(Var1 = Var2, Var2 = Var1)  # Add reverse pair
 ) %>%
 distinct()  # Remove duplicates

# Ensure factor levels match species order
cooccur_prob_df <- cooccur_prob_df %>%
 mutate(Var1 = factor(Var1, levels = unique(Var1)),
        Var2 = factor(Var2, levels = unique(Var2)))

library(ggplot2)

# Heatmap plot
ggplot(cooccur_prob_df, aes(x = Var1, y = Var2, fill = Freq)) +
 geom_tile(color = "white") +
 scale_fill_viridis_c(option = "plasma", na.value = "grey50") +
 labs(
  title = "Species Co-Occurrence Heatmap",
  x = "Species A",
  y = "Species B",
  fill = "Co-Occurrence (p-value)"
 ) +
 theme_minimal() +
 theme(
  axis.text.x = element_text(angle = 45, hjust = 1),
  axis.text.y = element_text(size = 10)
 )
