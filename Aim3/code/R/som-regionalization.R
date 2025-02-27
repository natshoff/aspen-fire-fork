# Self Organizing Maps (SOM) implementation
# for management priority landscapes

library(tidyverse)
library(sf)
library(kohonen)

projdir <- "/Users/max/Library/CloudStorage/OneDrive-Personal/mcook/aspen-fire/Aim3/"

# load and prep the grid data
grid <- st_read(paste0(projdir,'data/spatial/mod/srm_model_data.gpkg'))
glimpse(grid)

# keep a version with just forest type information
grid_fortyp <- grid %>% select(grid_id, dom_spp1, dom_spp2, dom_spp3)

# select numeric variables for SOM
X <- grid %>%
 st_drop_geometry() %>%
 select(-c(grid_id, p_area, p_count)) %>%  # Remove ID and very small numeric values
 select(where(is.numeric)) %>%  # Keep only numeric columns
 # Fill NA in population data and scale everything
 mutate(
  across(starts_with("pop_"), ~replace_na(.x, 0)),  # Replace NA with 0 in population data
  across(everything(), ~as.numeric(scale(.)))  # Standardize (mean = 0, SD = 1)
 )
# Check structure
str(X)

# tidy up what we can
rm(grid)
gc()


#==============VARIABLE WEIGHTING================#

#####
# WUI codes:
# 1: 'Forest/Shrubland/Wetland-dominated Intermix WU'
# 2: 'Forest/Shrubland/Wetland-dominated Interface WUI'
# 3: 'Grassland-dominated Intermix WUI'
# 4: 'Grassland -dominated Interface WUI'
#####

#####################################################
# manually weight variables based on domain knowledge
# Define manual weights based on domain knowledge
wt.manual <- c(
 # future fire
 trend_area = 5, trend_count = 1, 
 # future aspen
 historic = 2, ssp245 = 1, ssp585 = 1, delta245 = 4, delta585 = 10,
 # current aspen
 aspen10_pct = 10, aspen10_pixn = 1,
 # wui designation
 wui1 = 5, wui2 = 5, wui3 = 5, wui4 = 5,  
 # combustible mass of the built environment/ building counts
 combust_sum = 3, msbf_count_sum = 3,
 # population density and count
 pop_density_mean = 4, pop_count_sum = 4,  # Strongly upweight population
 # fire history (burned area)
 burned_area = 1, burned_pct = 1,
 # fire risk (conditional flame length and burn probability)
 iFLP6_p90 = 3, iFLP5_p90 = 2, iFLP4_p90 = 1, iBP_p90 = 3,
 iFLP1_p90 = 1, iFLP3_p90 = 1, iFLP2_p90 = 1,
 # mean canopy percent and total live basal area
 canopypct_mean = 1, balive_sum = 1
)

# Apply manual weights
X.wt.m <- X * wt.manual[colnames(X)]


##############################################
# Data-driven weighting scheme (PCA, variance)
# Run PCA on the scaled data
pca.results <- prcomp(X, scale = TRUE)
# Extract absolute loadings from the first principal component (PC1)
pca.imp <- abs(pca.results$rotation[, 1])

# Compute weights: Invert importance (rare variables get higher weight)
# wt.quant <- 1 / (pca.imp + 0.01)  # Small constant to avoid division by zero
wt.quant <- pca.imp

# Normalize weights between 1 and 10 (or another range)
wt.quant <- scales::rescale(wt.quant, to = c(1, 10))
# Convert to named vector
names(wt.quant) <- colnames(X)
# Print weights for review
print(sort(wt.quant, decreasing = TRUE))

# apply the weights
X.wt.q <- X * wt.quant[colnames(X)]


#==============SOM SETUP================#

# Define SOM grid size (adjust xdim & ydim as needed)
som.grid <- kohonen::somgrid(xdim = 8, ydim = 8, topo = "hexagonal")

# Train the SOM
set.seed(123)  # For reproducibility
som.model <- som(as.matrix(X.wt.q), 
                 grid = som.grid, 
                 rlen = 101,   # Number of training iterations
                 alpha = c(0.05, 0.01), # Learning rate (start, end)
                 keep.data = TRUE)
# View summary of the trained model
summary(som.model)


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
som.cluster <- cutree(hclust(dist(som.model$codes[[1]])), k = 7)  # Adjust 'k' as needed
# Assign clusters to each grid cell
grid$som.cluster <- som.cluster[som.model$unit.classif]
# Check distribution of clusters
table(grid$som.cluster)

# make the spatial map
ggplot(grid) +
 geom_sf(aes(fill = as.factor(som.cluster)), color = NA) +
 scale_fill_viridis_d(name = "SOM Cluster") +
 theme_minimal()





