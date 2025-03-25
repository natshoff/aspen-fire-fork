##########################################################################
# Fit semivariograms to explore the spatial dependence in FRPc and CBIbc #
# Create the spatial mesh grid for R-INLA based on variogram results.    #
##########################################################################

library(tidyverse)
library(sf)
library(gstat)

maindir <- '/Users/max/Library/CloudStorage/OneDrive-Personal/mcook/aspen-fire/Aim2/'

#=========Load the Prepped gridcell data=========#

# gridcell_prep.R
fp <- paste0(maindir,"data/tabular/mod/model_data_cleaned.csv")
grid_sf <- read_csv(fp) %>%
 distinct(grid_index, .keep_all=T) %>%
 select(Fire_ID, grid_index, x, y,
        log_frp_csum, CBIbc_p90, CBIbc_p95) %>%
 arrange(grid_index) %>%
 st_as_sf(., coords = c("x", "y"), crs = 4326)
glimpse(grid_sf)

# extract a coordinate matrix
# for creating the mesh grid 
coords.mat <- grid_sf %>% 
 st_coordinates(.) %>%
 as.matrix()


#=========Fit the semivariograms=========#

###################################################
# Function to compute semivariogram for each fire #
fire_variogram <- function(Fire_ID, data, N) {
 fire_sp <- data %>% filter(Fire_ID == !!Fire_ID)
 # Only compute if there are enough points
 if (nrow(fire_sp) >= N) { 
  # compute the variogram
  vario <- fit.variogram(variogram(log_frp_csum ~ 1, data = fire_sp), model = vgm("Exp"))
  print(paste("Computed variogram for fire:", Fire_ID, "with", nrow(fire_sp), "points"))
  return(vario)
 } else {
  print(paste("Skipping fire:", Fire_ID, "- Not enough points"))
  return(NULL)
 }
}

# run for each fire separately
fire_ids <- unique(grid_sf$Fire_ID)
variograms <- lapply(fire_ids, fire_variogram, data = grid_sf, N = 10)
variograms <- do.call(rbind, variograms)  # Combine results

# Plot variograms grouped by Fire ID
(p1 <- ggplot(variograms, 
              aes(x = dist, y = gamma, 
                  group = Fire_ID, 
                  color = Fire_ID)) +
 geom_line(alpha = 0.3, size=0.9) +  # Light transparency to see overlap
 labs(x = "Distance (degrees)", y = "Semivariance",
      title = "Semivariograms for Individual Fires") +
 theme_minimal() +
 theme(legend.position = "none"))

# save the plot.
out_png <- paste0(maindir,'figures/INLA_Fire_SemiVariograms_FRP.png')
ggsave(out_png, plot=p1, dpi=500, bg = 'white')

# Find the approximate range (where semivariance plateaus) for each fire
range_est <- variograms %>%
 drop_na() %>%
 group_by(Fire_ID) %>%
 summarize(
  range_m = ifelse(any(gamma < max(gamma) * 0.9, na.rm = TRUE),
                   max(dist[gamma < max(gamma) * 0.9], na.rm = TRUE),
                   NA_real_)  # Assign NA when no valid values exist
 )
qt <- quantile(range_est$range_m, probs = seq(.1, .9, by = .1), na.rm=TRUE)
qt

# Compute and print the mean spatial range for FRP
mean_range_frp <- mean(range_est$range_m, na.rm = TRUE)
print(paste("Mean spatial range for FRP:", round(mean_range_frp, 2), "meters"))

# Histogram of estimated spatial ranges
ggplot(range_est, aes(x = range_m)) +
 geom_histogram(bins = 20, fill = "blue", alpha = 0.6) +
 labs(x = "Estimated Range (meters)", y = "Count",
      title = "Distribution of Within-Fire Spatial Correlation Ranges") +
 theme_minimal()

# save the plot.
out_png <- paste0(maindir,'figures/INLA_Fire_EstimatedRange_FRP.png')
ggsave(out_png, dpi=500, bg = 'white')

rm(sp, range_est, variograms)
gc()


#######################################
# Create the R-INLA SPDE Spatial Mesh #
#######################################

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

# ####### Save the mesh as spatial objects
# # Convert INLA mesh vertices to an sf object
# mesh_v <- data.frame(mesh$loc) # Extract coordinates
# colnames(mesh_v) <- c("x", "y") # set column names
# mesh_sf <- st_as_sf(mesh_v, coords = c("x", "y"), crs = 4326)
# # Export vertices to a gpkg
# st_write(mesh_sf, paste0(maindir,"data/spatial/mod/INLA_mesh_vertices.gpkg"),
#          layer = "mesh_vertices", driver = "GPKG", append=FALSE)
# # Function to convert triangles to polygons
# mesh_t <- mesh$graph$tv  # Triangle indices
# tri_coords <- mesh$loc  # Mesh coordinates
# # Create list of polygons
# poly_list <- lapply(1:nrow(mesh_t), function(i) {
#  coords <- tri_coords[mesh_t[i, ], ]
#  coords <- rbind(coords, coords[1, ])  # Close the polygon
#  st_polygon(list(coords))
# })
# # Convert to sf object
# mesh_poly_sf <- st_sfc(poly_list, crs = 4326) %>% st_sf(geometry = .)
# # Export mesh triangles to GPKG
# st_write(mesh_poly_sf, paste0(maindir,"data/spatial/mod/INLA_mesh_triangles.gpkg"),
#          layer = "mesh_triangles", driver = "GPKG", append=FALSE)
# rm(mesh_v, mesh_sf, mesh_t, tri_coords, poly_list, mesh_poly_sf)
# gc()