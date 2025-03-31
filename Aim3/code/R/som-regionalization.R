###########################################
# Self Organizing Maps (SOM) implementation
# for management priority landscapes

library(tidyverse)
library(sf)
library(kohonen)
library(reshape2)
library(scales)
library(ggspatial)
library(gridExtra)
library(ggradar)

projdir <- "/Users/max/Library/CloudStorage/OneDrive-Personal/mcook/aspen-fire/Aim3/"

# load and prep the firesheds data
firesheds <- st_read(paste0(projdir,'data/spatial/mod/srm_firesheds_model_data.gpkg'))
glimpse(firesheds)


####################################
# select numeric variables for SOM #
X <- firesheds %>%
 # drop the geometry
 st_drop_geometry() %>%
 rename(
  # tidy some of the column names
  n_patches = number_of_patches,
  patch_den = patch_density,
  big_patch = largest_patch_index,
  pop_den = pop_density_max,
  pop_n = pop_count_sum,
  wui_dist = wui_dist_mean,
  lf_canopy = forest_cc_mean
 ) %>%
 mutate(
  # calculate contemporary aspen hectares
  aspen10_ha = aspen10_pixn * 0.01,
  forest_ha = forest_pixels * 0.09,
  # percent of forested area that is aspen?
  aspen_forest = aspen10_ha / forest_ha,
  # combine the WUI classes (???)
  wui = wui1 + wui2 + wui3 + wui4
 ) %>%
 # keep only the subset we want to use for SOM
 select(
  trend_count, delta585, # future future / future aspen
  aspen_forest, aspen10_ha, n_patches, patch_den, big_patch, # contemporary aspen presence/patches
  combust_sum, pop_den, pop_n, wui, wui_dist, # built environment metrics
  burned_pct_c, # cumulative burned area percent of fireshed since 1984
  whp_p90, hui_p90, cfl_p90, # wildfire risk to communities
  forest_pct, lf_canopy, # forested area and mean canopy cover
  aspen, douglas_fir, lodgepole, gambel_oak, sagebrush, # forest type percents
  pinon_juniper, ponderosa, spruce_fir, white_fir # forest type percents
 ) %>%
 # Fill NAs with zero and scale variables
 mutate(
  across(everything(), ~ replace_na(.x, 0))
 ) %>%
 mutate(
  across(everything(), ~ as.numeric(scale(., center=T)))
 )
# Check structure
str(X)

#==============VARIABLE WEIGHTING================#

##########################
# WUI codes for reference:
# 1: 'Forest/Shrubland/Wetland-dominated Intermix WU'
# 2: 'Forest/Shrubland/Wetland-dominated Interface WUI'
# 3: 'Grassland-dominated Intermix WUI'
# 4: 'Grassland -dominated Interface WUI'


#==============SOM SETUP================#
set.seed(123)

# Define SOM grid size (adjust xdim & ydim as needed)
som.grid <- kohonen::somgrid(xdim = 20, ydim = 20, topo = "hexagonal")

# Train the SOM
som.model <- som(
 as.matrix(X),
 grid = som.grid,
 rlen = 5001,   # Number of training iterations
 alpha = c(0.1, 0.02), # Learning rate (start, end)
 keep.data = TRUE)

# View summary of the trained model
summary(som.model)
# Plot the nodes
plot(som.model)


#==============MAP GRIDS================#

# Compute distance matrix on SOM codebook vectors
dist_mat <- dist(som.model$codes[[1]]) 

# Perform hierarchical clustering on SOM node weights
# som.cluster <- kmeans(som.model$codes[[1]], centers = 6)$cluster  # Adjust 'k' as needed
hc <- hclust(dist_mat, method = "ward.D2")
som.cluster <- cutree(hc, k = 6)  

# Assign clusters to each fireshed
firesheds$som.cluster <- som.cluster[som.model$unit.classif]
# Check distribution of clusters
table(firesheds$som.cluster)

# make the spatial map
(cl.map <- ggplot(firesheds) +
 geom_sf(aes(fill = as.factor(som.cluster)), color = NA) +
 scale_fill_brewer(palette = "Accent", name = "SOM Cluster") +
 theme_void() +
 coord_sf(expand = FALSE) +
 ggspatial::annotation_scale(location = "br", width_hint = 0.1,
                             pad_x = unit(0.15,"in"), pad_y = unit(0.05,"in"),
                             line_width = 0.5, text_pad = unit(0.15,"cm"),
                             height = unit(0.15,"cm")) +
 ggspatial::annotation_north_arrow(location = "br", which_north = "true",
                                   pad_x = unit(0.25,"in"), pad_y= unit(0.2,"in"),
                                   width = unit(0.8,"cm"), height = unit(0.8,"cm"),
                                   style = north_arrow_fancy_orienteering) +
 theme(
  legend.title = element_text(angle = 90, size = 9, vjust = 1.2, hjust = 0.5),
  legend.text = element_text(size = 8),
  legend.key.size = unit(0.4, "cm"),  
  legend.position = c(0.98, 0.6)
 ) +
 guides(fill = guide_legend(title.position = "left", title.vjust = 1)))
# save the plot.
out_png <- paste0(projdir,'figures/SRM_Firesheds_SOMclusters.png')
ggsave(out_png, plot = cl.map, dpi = 500, width = 8, height = 6, bg = 'white')


#==============RADAR PLOTS================#

# Scale the SOM input variables (already done)
# Add SOM cluster assignments
# Get codebook vectors
som_codes <- som.model$codes[[1]]
# Run clustering on SOM units
dist_mat <- dist(som_codes)
hc <- hclust(dist_mat, method = "ward.D2")
som.cluster <- cutree(hc, k = 6)  
# Add cluster labels to each observation (in X)
X_cl.sc <- X %>%
 mutate(cluster = as.factor(som.cluster[som.model$unit.classif]))

# Compute cluster means of scaled variables
cl_means <- X_cl.sc %>%
 group_by(cluster) %>%
 summarise(
  # z-score means per cluster
  across(where(is.numeric), \(x) mean(x, na.rm = TRUE))
 ) 

##############
# Radar plot #
# Tidy the data:
cl_means_r <- cl_means %>%
 mutate(across(-cluster, ~ rescale(.x))) %>%
 mutate(group = paste0("Cluster ", cluster)) %>%
 select(-cluster) %>%
 relocate(group)

# Assign cluster colors same as map
cluster_colors <- RColorBrewer::brewer.pal(6, "Accent")
names(cluster_colors) <- paste0("Cluster ", 1:6)

# plot all clusters in one
ggradar(cl_means_r,
        values.radar = c("0", "0.5", "1"),
        grid.min = 0, grid.mid = 0.5, grid.max = 1,
        group.line.width = 0.8,
        group.point.size = 1,
        axis.label.size = 2.5, 
        grid.label.size = 3.5, 
        legend.text.size = 8,
        group.colours = cluster_colors,
        fill = TRUE, fill.alpha = 0.2,
        background.circle.colour = "grey90",
        gridline.mid.colour = "grey80",
        font.radar = "sans") +
 theme(
  plot.title = element_text(size = 12, hjust = 0.5)  
 )


# Function to generate ggplot-based radar chart per cluster
make_ggradar <- function(df_row, cluster_id, fill_color) {
 # Plot with ggradar
 ggradar(df_row,
         values.radar = c("0", "0.5", "1"),
         grid.min = 0, grid.mid = 0.5, grid.max = 1,
         group.line.width = 1.2,
         group.point.size = 2,
         group.colours = fill_color,
         fill = TRUE, fill.alpha = 0.3,
         axis.label.size = 1.8,
         grid.label.size = 3,
         background.circle.colour = "grey90",
         gridline.mid.colour = "grey80",
         font.radar = "sans",
         legend.position = "none")
}

# List of radar plots
(radar_plots <- lapply(1:6, function(i) {
 df_i <- cl_means_r %>% filter(group == paste0("Cluster ", i))
 make_ggradar(df_i, paste0("Cluster ", i), cluster_colors[[i]])
}))

# panel the radar plots
radar_grid <- patchwork::wrap_plots(radar_plots, ncol = 2)
radar_grid

# combine with the map
map_radar <- cl.map + radar_grid + 
 patchwork::plot_layout(
  ncol = 2, widths = c(1.2, 2)) 
map_radar
# save the plot.
out_png <- paste0(projdir,'figures/SRM_Firesheds_SOMclusters_wRadar.png')
ggsave(out_png, plot = map_radar, dpi = 500, 
       width = 8, height = 6, bg = 'white')



#==============ALTERNATE PLOTS================#

###########
# Heatmap #
# Long format for ggplot heatmap
cl_melt <- melt(cl_means, id.vars = "cluster")
# plot it.
ggplot(cl_melt, aes(x = variable, y = cluster, fill = value)) +
 geom_tile(color = "white") +
 scale_fill_viridis_c(option = "C", name = "Z-Score") +
 theme_minimal(base_size = 11) +
 theme(axis.text.x = element_text(angle = 45, hjust = 1),
       panel.grid = element_blank()) +
 labs(title = "SOM Cluster Characterization (Standardized)",
      x = "Variable", y = "Cluster")

###########
# Boxplot #
codes_df <- as.data.frame(som.model$codes[[1]])
colnames(codes_df) <- colnames(X)
codes_df$cluster <- factor(paste0("Cluster ", som.cluster))
codes_df$cluster_num <- factor(som.cluster) # keep a numeric version

# Reshape for visualization
codes_m <- melt(codes_df, id.vars = c("cluster", "cluster_num"))
# Set factor levels explicitly (in case of ordering)
codes_m$cluster <- factor(codes_m$cluster, levels = paste0("Cluster ", 1:6))
# Plot variable contributions by cluster with custom colors
(p.box <- ggplot(codes_m, aes(x = cluster_num, y = value, fill = cluster)) +
  geom_boxplot(outlier.size = 0.4, lwd = 0.2) +
  facet_wrap(~ variable, scales = "free") +
  theme_classic(base_size = 10) +
  labs(x = "Cluster", y = "Value") +
  scale_fill_manual(values = cluster_colors, name = "SOM Cluster") +
  theme(
   axis.text.x = element_text(angle = 0, hjust = 0.5),
   legend.position = c(0.78, 0.08),  # manual position
   legend.direction = "vertical",
   legend.title = element_text(angle = 0, vjust = 0.5, size = 11),
   legend.text = element_text(size = 10),
   legend.key.size = unit(0.4, "cm")
  ) +
  guides(fill = guide_legend(ncol = 2, byrow = TRUE)))
# save the plot.
out_png <- paste0(projdir,'figures/SRM_Firesheds_SOMclusters_VarMetrics.png')
ggsave(out_png, plot = p.box, dpi = 500, 
       width = 8, height = 6, bg = 'white')


# #==============CLUSTER VIZ================#
# 
# # Extract grid info
# som_codes <- som.model$codes[[1]]
# unit_coords <- som.model$grid$pts
# unit_classif <- som.model$unit.classif
# # Get SOM unit positions and codebook vectors
# xdim <- som.model$grid$xdim
# ydim <- som.model$grid$ydim
# 
# # Function to get neighbor indices
# get_neighbors <- function(i, coords) {
#  this_coord <- coords[i, ]
#  diffs <- sweep(coords, 2, this_coord)
#  dists <- sqrt(rowSums(diffs^2))
#  neighbors <- which(dists > 0 & dists <= 1.1)  # typical threshold for adjacent hexes
#  return(neighbors)
# }
# 
# # Calculate average distance to neighbors for each SOM unit
# u_matrix_vals <- map_dbl(1:nrow(som_codes), function(i) {
#  neighbors <- get_neighbors(i, unit_coords)
#  if (length(neighbors) == 0) return(NA)
#  dists <- apply(som_codes[neighbors, , drop = FALSE], 1, function(n) {
#   sqrt(sum((som_codes[i, ] - n)^2))
#  })
#  mean(dists)
# })
# 
# # Hierarchical clustering
# dist_mat <- dist(som_codes)
# hc <- hclust(dist_mat, method = "ward.D2")
# n_clusters <- 6
# som_cluster <- cutree(hc, k = n_clusters)
# 
# # Make a SOM dataframe
# som_df <- as_tibble(unit_coords) %>%
#  rename(x = x, y = y) %>%
#  mutate(node = row_number(),
#         uval = u_matrix_vals,
#         cluster = as.factor(som_cluster))
# 
# # Set hex size
# hex_size <- 0.5
# # Create hexagon coordinates
# hex_coords <- som_df %>%
#  rowwise() %>%
#  mutate(hex = list(tibble(
#   hx = x + hex_size * cos(seq(0, 2 * pi, length.out = 7)),
#   hy = y + hex_size * sin(seq(0, 2 * pi, length.out = 7))
#  ))) %>%
#  unnest(hex)
# 
# # make the plot:
# uplot <- ggplot(hex_coords, aes(x = hx, y = hy, group = node)) +
#  geom_polygon(aes(fill = uval), color = "white", size = 0.1) +
#  scale_fill_viridis_c(name = "U-Matrix Distance", option = "C", na.value = "grey90") +
#  geom_path(aes(group = cluster), color = "black", linewidth = 0.4) +
#  coord_equal() +
#  theme_void() +
#  theme(
#   legend.position = "right",
#   plot.title = element_text(hjust = 0.5)
#  ) +
#  labs(title = "SOM U-Matrix with Cluster Outlines")
# uplot

