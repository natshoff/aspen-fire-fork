# Self Organizing Maps (SOM) implementation
# for management priority landscapes

library(tidyverse)
library(sf)
library(kohonen)
library(reshape2)

projdir <- "/Users/max/Library/CloudStorage/OneDrive-Personal/mcook/aspen-fire/Aim3/"

# load and prep the firesheds data
firesheds <- st_read(paste0(projdir,'data/spatial/mod/srm_firesheds_model_data.gpkg'))
glimpse(firesheds)

# keep a version with just EVT information
firesheds_lf <- firesheds %>% 
 select(sfs_id, dom_evt1, dom_evt2, dom_evt3, dom_evt4)

# select numeric variables for SOM
X <- firesheds %>%
 st_drop_geometry() %>%
 # keep only the subset we want to use for SOM
 select(-c(sfs_id, fs_id, burned_area_c, ssp245, ssp585, 
           historic, aspen10_pixn, sfs_area, sfs_area_ha)) %>%  
 select(where(is.numeric)) %>%  # Keep only numeric columns
 # Fill NA in population data and scale everything
 mutate(
  across(starts_with("pop_"), ~replace_na(.x, 0)),  # Replace NA with 0 in population data
  across(everything(), ~as.numeric(scale(.)))  # Standardize (mean = 0, SD = 1)
 )
# Check structure
str(X)
summary(X)

#==============VARIABLE WEIGHTING================#

##########################
# WUI codes for reference:
# 1: 'Forest/Shrubland/Wetland-dominated Intermix WU'
# 2: 'Forest/Shrubland/Wetland-dominated Interface WUI'
# 3: 'Grassland-dominated Intermix WUI'
# 4: 'Grassland -dominated Interface WUI'


#==============SOM SETUP================#

# Define SOM grid size (adjust xdim & ydim as needed)
som.grid <- kohonen::somgrid(xdim = 10, ydim = 10, topo = "hexagonal")

# Train the SOM
set.seed(123)  # For reproducibility
som.model <- som(as.matrix(X), 
                 grid = som.grid, 
                 rlen = 1001,   # Number of training iterations
                 alpha = c(0.1, 0.02), # Learning rate (start, end)
                 keep.data = TRUE)
# View summary of the trained model
summary(som.model)
# Plot the nodes
plot(som.model, type = "counts")


#==============CLUSTER VIZ================#

# U-Matrix (shows distances between SOM nodes)
plot(som.model, type = "dist.neighbours", main = "SOM U-Matrix")
# Codebook vectors (cluster prototypes)
plot(som.model, type = "codes", main = "Cluster Prototypes")
# Component planes (influence of each variable)
par(mfrow = c(4, 4))  # Adjust layout for more variables
for (i in 1:ncol(X)) {
 plot(som.model, type = "property", property = som.model$codes[[1]][, i], 
      main = colnames(X)[i])
}


#==============MAP GRIDS================#

# Perform hierarchical clustering on SOM node weights
som.cluster <- kmeans(som.model$codes[[1]], centers = 6)$cluster  # Adjust 'k' as needed
# Assign clusters to each fireshed
firesheds$som.cluster <- som.cluster[som.model$unit.classif]
# Check distribution of clusters
table(firesheds$som.cluster)

# make the spatial map
ggplot(firesheds) +
 geom_sf(aes(fill = as.factor(som.cluster)), color = NA) +
 scale_fill_viridis_d(name = "SOM Cluster") +
 theme_minimal()



#==============VARIABLE CONTRIBUTION================#

codes_df <- as.data.frame(som.model$codes[[1]])
colnames(codes_df) <- colnames(X)
codes_df$cluster <- factor(kmeans(som.model$codes[[1]], centers = 5)$cluster)

# Reshape for visualization
codes_m <- melt(codes_df, id.vars = "cluster")

# Plot variable contributions by cluster
ggplot(codes_m, aes(x = cluster, y = value, fill = cluster)) +
 geom_boxplot() +
 facet_wrap(~ variable, scales = "free") +
 theme_minimal() +
 scale_fill_viridis_d(name = "SOM Cluster")


